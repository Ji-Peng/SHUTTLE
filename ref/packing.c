/*
 * packing.c - Serialization of keys and signatures for SHUTTLE.
 *
 * Public key format:
 *   rho (SEEDBYTES) || polypk_pack(b[0..M-1])
 *
 * Secret key format:
 *   rho || tr || key || polyeta_pack(s[0..L-1]) || polyeta_pack(e[0..M-1])
 *
 * Signature format (NGCC-Signature Alg 2 compressed form,
 *                   two-stream rANS layout per SHUTTLE_rANS.tex §3.3):
 *
 *   seedC || irs_signs
 *     || uint16 zhi_rans_len  || rANS(z-hi)      + pad to ZHI_RESERVED
 *     || polyz0_lo_pack(lo(z^(0)))
 *     || L * polyz1_lo_pack(lo(z^(1..lenS)))
 *     || uint16 hint_rans_len || rANS(hint)      + pad to HINT_RESERVED
 *
 * The z-hi rANS block carries n*(lenS+1) coefficients in concatenation
 * order [ HighBits(z^(0)), HighBits(z^(1)), ..., HighBits(z^(lenS)) ].
 * All share scale r/alpha_r so a single frequency table fits them all
 * (mode-128 also shares the table with the hint context — see
 * shuttle_rans.{c,h}).
 *
 * Each rANS block is a fixed-size slot (length prefix + encoded stream +
 * zero padding). The fixed slot is what lets the signature have a single
 * on-the-wire length. Reservation budgets live in params.h
 * (SHUTTLE_ZHI_RESERVED_BYTES / SHUTTLE_HINT_RESERVED_BYTES); the
 * derivation is in SHUTTLE_rANS.tex Tab 8 (p_rans^* = 2^-20).
 *
 * Two rANS streams (saves 6 B/sig of fixed overhead vs. the historical
 * three-stream design that packed z^(0) on its own stream):
 *   - stream 1 (z-hi): z^(0) hi  ⨁  z^(1..lenS) hi (single shared table)
 *   - stream 2 (hint): MakeHint output
 *
 * IRS sign bits: TAU bits packed into IRS_SIGNBYTES bytes. Bit i is 1 if
 * irs_signs[i] == +1, 0 if irs_signs[i] == -1. Combined with the challenge
 * c (from c_tilde), these bits define c_eff = sum_i irs_sign_i * x^{j_i}.
 *
 * Encoder-side overflow returns rc = -2 from pack_sig; OOV is impossible
 * by construction (theoretical vocabulary covers the 11*sigma tight bound).
 * sign.c::crypto_sign_signature treats rc != 0 as a signing-round
 * rejection. With p_rans^* = 2^-20 the throughput hit is invisible.
 */

#include <string.h>

#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "shuttle_rans.h"

/* ============================================================
 * Public / secret key packing (unchanged across scheme versions).
 * ============================================================ */

void pack_pk(uint8_t pk[SHUTTLE_PUBLICKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const polyveck *b)
{
  unsigned int i;

  memcpy(pk, rho, SHUTTLE_SEEDBYTES);
  pk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_M; ++i) {
    polypk_pack(pk, &b->vec[i]);
    pk += SHUTTLE_POLYPK_PACKEDBYTES;
  }
}

void unpack_pk(uint8_t rho[SHUTTLE_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[SHUTTLE_PUBLICKEYBYTES])
{
  unsigned int i;

  memcpy(rho, pk, SHUTTLE_SEEDBYTES);
  pk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_M; ++i) {
    polypk_unpack(&b->vec[i], pk);
    pk += SHUTTLE_POLYPK_PACKEDBYTES;
  }
}

void pack_sk(uint8_t sk[SHUTTLE_SECRETKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const uint8_t tr[SHUTTLE_TRBYTES],
             const uint8_t key[SHUTTLE_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e)
{
  unsigned int i;

  memcpy(sk, rho, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  memcpy(sk, tr, SHUTTLE_TRBYTES);
  sk += SHUTTLE_TRBYTES;

  memcpy(sk, key, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_L; ++i) {
    polyeta_pack(sk, &s->vec[i]);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }

  for(i = 0; i < SHUTTLE_M; ++i) {
    polyeta_pack(sk, &e->vec[i]);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }
}

void unpack_sk(uint8_t rho[SHUTTLE_SEEDBYTES],
               uint8_t tr[SHUTTLE_TRBYTES],
               uint8_t key[SHUTTLE_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[SHUTTLE_SECRETKEYBYTES])
{
  unsigned int i;

  memcpy(rho, sk, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  memcpy(tr, sk, SHUTTLE_TRBYTES);
  sk += SHUTTLE_TRBYTES;

  memcpy(key, sk, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_L; ++i) {
    polyeta_unpack(&s->vec[i], sk);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }

  for(i = 0; i < SHUTTLE_M; ++i) {
    polyeta_unpack(&e->vec[i], sk);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }
}

/* ============================================================
 * Signature packing (compressed form: two rANS streams + z bit-packs).
 * ============================================================ */

#define OFF_C_TILDE     0
#define OFF_IRS_SIGNS   (OFF_C_TILDE + SHUTTLE_CTILDEBYTES)
#define OFF_ZHI_LEN     (OFF_IRS_SIGNS + SHUTTLE_IRS_SIGNBYTES)
#define OFF_ZHI_DATA    (OFF_ZHI_LEN + 2)
#define OFF_Z0_LO       (OFF_ZHI_DATA + SHUTTLE_ZHI_RESERVED_BYTES)
#define OFF_Z1_LO       (OFF_Z0_LO + SHUTTLE_POLYZ0_LO_PACKEDBYTES)
#define OFF_HINT_LEN    (OFF_Z1_LO + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES)
#define OFF_HINT_DATA   (OFF_HINT_LEN + 2)

int pack_sig(uint8_t sig[SHUTTLE_BYTES],
             const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
             const int8_t irs_signs[SHUTTLE_TAU],
             const poly *z_1,
             const polyveck *h)
{
  unsigned int i, j;
  int rc;

  /* 1. seedC */
  memcpy(&sig[OFF_C_TILDE], c_tilde, SHUTTLE_CTILDEBYTES);

  /* 2. irs_signs bitmap: bit i = 1 iff irs_signs[i] > 0. */
  memset(&sig[OFF_IRS_SIGNS], 0, SHUTTLE_IRS_SIGNBYTES);
  for(i = 0; i < SHUTTLE_TAU; ++i) {
    if(irs_signs[i] > 0)
      sig[OFF_IRS_SIGNS + (i >> 3)] |= (uint8_t)(1u << (i & 7));
  }

  /* 3. Split z^(0) and z^(1..lenS) into hi (rANS) + lo (bit-pack).
   *    All hi arrays are concatenated into one buffer that feeds the
   *    unified z-hi rANS stream. Layout:
   *      z_hi_flat[0 .. n-1]                = HighBits_{alpha_0'}(z^(0))
   *      z_hi_flat[k*n .. (k+1)*n - 1]      = HighBits_{alpha_r}(z^(k))
   *                                           for k in {1..lenS}. */
  int32_t z_hi_flat[(SHUTTLE_L + 1) * SHUTTLE_N];
  int32_t lo_scratch[SHUTTLE_N];

  polyz0_split(&z_hi_flat[0], lo_scratch, &z_1[0]);
  polyz0_lo_pack(&sig[OFF_Z0_LO], lo_scratch);

  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz1_split(&z_hi_flat[(i + 1) * SHUTTLE_N], lo_scratch, &z_1[1 + i]);
    polyz1_lo_pack(&sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES],
                   lo_scratch);
  }

  /* 4. rANS-encode the concatenated z-hi stream (n*(lenS+1) coefs). */
  size_t zhi_rans_len = 0;
  rc = shuttle_rans_encode_zhi(&sig[OFF_ZHI_DATA], &zhi_rans_len,
                               SHUTTLE_ZHI_RESERVED_BYTES,
                               z_hi_flat, (SHUTTLE_L + 1) * SHUTTLE_N);
  if(rc != 0)
    return rc;

  sig[OFF_ZHI_LEN + 0] = (uint8_t)(zhi_rans_len & 0xFF);
  sig[OFF_ZHI_LEN + 1] = (uint8_t)((zhi_rans_len >> 8) & 0xFF);
  if(zhi_rans_len < SHUTTLE_ZHI_RESERVED_BYTES)
    memset(&sig[OFF_ZHI_DATA + zhi_rans_len], 0,
           SHUTTLE_ZHI_RESERVED_BYTES - zhi_rans_len);

  /* 5. rANS-encode hint h. */
  int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h_flat[i * SHUTTLE_N + j] = h->vec[i].coeffs[j];

  size_t hint_rans_len = 0;
  rc = shuttle_rans_encode_hint(&sig[OFF_HINT_DATA], &hint_rans_len,
                                SHUTTLE_HINT_RESERVED_BYTES,
                                h_flat, SHUTTLE_M * SHUTTLE_N);
  if(rc != 0)
    return rc;

  sig[OFF_HINT_LEN + 0] = (uint8_t)(hint_rans_len & 0xFF);
  sig[OFF_HINT_LEN + 1] = (uint8_t)((hint_rans_len >> 8) & 0xFF);
  if(hint_rans_len < SHUTTLE_HINT_RESERVED_BYTES)
    memset(&sig[OFF_HINT_DATA + hint_rans_len], 0,
           SHUTTLE_HINT_RESERVED_BYTES - hint_rans_len);

  return 0;
}

int unpack_sig(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
               int8_t irs_signs[SHUTTLE_TAU],
               poly *z_1,
               polyveck *h,
               const uint8_t sig[SHUTTLE_BYTES])
{
  unsigned int i, j;
  int rc;

  /* 1. seedC */
  memcpy(c_tilde, &sig[OFF_C_TILDE], SHUTTLE_CTILDEBYTES);

  /* 2. irs_signs: bit 1 -> +1, bit 0 -> -1. */
  for(i = 0; i < SHUTTLE_TAU; ++i) {
    if((sig[OFF_IRS_SIGNS + (i >> 3)] >> (i & 7)) & 1)
      irs_signs[i] = (int8_t)1;
    else
      irs_signs[i] = (int8_t)-1;
  }

  /* 3. z-hi rANS decode. */
  size_t zhi_rans_len = (size_t)sig[OFF_ZHI_LEN + 0]
                      | ((size_t)sig[OFF_ZHI_LEN + 1] << 8);
  if(zhi_rans_len == 0 || zhi_rans_len > SHUTTLE_ZHI_RESERVED_BYTES)
    return -1;

  int32_t z_hi_flat[(SHUTTLE_L + 1) * SHUTTLE_N];
  rc = shuttle_rans_decode_zhi(z_hi_flat, (SHUTTLE_L + 1) * SHUTTLE_N,
                               &sig[OFF_ZHI_DATA], zhi_rans_len);
  if(rc != 0)
    return rc;

  /* 4. z^(0): combine hi (from z-hi stream) + lo (bit-packed). */
  int32_t lo_scratch[SHUTTLE_N];
  polyz0_lo_unpack(lo_scratch, &sig[OFF_Z0_LO]);
  polyz0_combine(&z_1[0], &z_hi_flat[0], lo_scratch);

  /* 5. z^(1..lenS): combine hi + lo. */
  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz1_lo_unpack(lo_scratch,
                     &sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES]);
    polyz1_combine(&z_1[1 + i],
                   &z_hi_flat[(i + 1) * SHUTTLE_N],
                   lo_scratch);
  }

  /* 6. hint rANS. */
  size_t hint_rans_len = (size_t)sig[OFF_HINT_LEN + 0]
                       | ((size_t)sig[OFF_HINT_LEN + 1] << 8);
  if(hint_rans_len == 0 || hint_rans_len > SHUTTLE_HINT_RESERVED_BYTES)
    return -1;

  int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
  rc = shuttle_rans_decode_hint(h_flat, SHUTTLE_M * SHUTTLE_N,
                                &sig[OFF_HINT_DATA], hint_rans_len);
  if(rc != 0)
    return rc;

  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h->vec[i].coeffs[j] = h_flat[i * SHUTTLE_N + j];

  return 0;
}
