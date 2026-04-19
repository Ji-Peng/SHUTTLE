/*
 * packing.c - Serialization of keys and signatures for SHUTTLE.
 *
 * Public key format:
 *   rho (SEEDBYTES) || polypk_pack(b[0..M-1])
 *
 * Secret key format:
 *   rho || tr || key || polyeta_pack(s[0..L-1]) || polyeta_pack(e[0..M-1])
 *
 * Signature format (NGCC-Signature Alg 2 compressed form with full rANS
 * entropy coding for Z_0, HighBits(z[1..L]) and hint h; see packing.h
 * for the byte-slot diagram and per-mode reservation budgets):
 *   seedC || irs_signs
 *     || uint16 z0_rans_len || rANS(Z_0) + pad
 *     || L * polyz1_lo_pack(lo(z[1..L]))
 *     || uint16 z1_rans_len || rANS(hi(z[1..L])) + pad
 *     || uint16 hint_rans_len || rANS(h) + pad
 *
 * IRS sign bits: TAU bits packed into IRS_SIGNBYTES bytes. Bit i is 1 if
 * irs_signs[i] == +1, 0 if irs_signs[i] == -1. Combined with the challenge
 * c (from c_tilde), these bits define c_eff = sum_i irs_sign_i * x^{j_i}.
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
 * Signature packing (compressed form, rANS for Z_0 + z1 highs + hint h).
 * ============================================================ */

#define OFF_C_TILDE     0
#define OFF_IRS_SIGNS   (OFF_C_TILDE + SHUTTLE_CTILDEBYTES)
#define OFF_Z0_LEN      (OFF_IRS_SIGNS + SHUTTLE_IRS_SIGNBYTES)
#define OFF_Z0_DATA     (OFF_Z0_LEN + 2)
#define OFF_Z1_LO       (OFF_Z0_DATA + SHUTTLE_Z0_RANS_RESERVED_BYTES)
#define OFF_Z1_HI_LEN   (OFF_Z1_LO + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES)
#define OFF_Z1_HI_DATA  (OFF_Z1_HI_LEN + 2)
#define OFF_HINT_LEN    (OFF_Z1_HI_DATA + SHUTTLE_Z1_RANS_RESERVED_BYTES)
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

  /* 3. Z_0 rANS. */
  int32_t z0_flat[SHUTTLE_N];
  for(i = 0; i < SHUTTLE_N; ++i)
    z0_flat[i] = z_1[0].coeffs[i];

  size_t z0_rans_len = 0;
  rc = shuttle_rans_encode_z0(&sig[OFF_Z0_DATA], &z0_rans_len,
                              SHUTTLE_Z0_RANS_RESERVED_BYTES,
                              z0_flat, SHUTTLE_N);
  if(rc != 0)
    return rc;

  sig[OFF_Z0_LEN + 0] = (uint8_t)(z0_rans_len & 0xFF);
  sig[OFF_Z0_LEN + 1] = (uint8_t)((z0_rans_len >> 8) & 0xFF);
  if(z0_rans_len < SHUTTLE_Z0_RANS_RESERVED_BYTES)
    memset(&sig[OFF_Z0_DATA + z0_rans_len], 0,
           SHUTTLE_Z0_RANS_RESERVED_BYTES - z0_rans_len);

  /* 4. z_1[1..L]: split each poly, write lo parts, accumulate hi parts. */
  int32_t z1_hi_flat[SHUTTLE_L * SHUTTLE_N];
  int32_t lo_scratch[SHUTTLE_N];
  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz1_split(&z1_hi_flat[i * SHUTTLE_N], lo_scratch, &z_1[1 + i]);
    polyz1_lo_pack(&sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES],
                   lo_scratch);
  }

  /* 5. rANS-encode the L*N high values into the z1 reservation. */
  size_t z1_rans_len = 0;
  rc = shuttle_rans_encode_z1(&sig[OFF_Z1_HI_DATA], &z1_rans_len,
                              SHUTTLE_Z1_RANS_RESERVED_BYTES,
                              z1_hi_flat, SHUTTLE_L * SHUTTLE_N);
  if(rc != 0)
    return rc;

  sig[OFF_Z1_HI_LEN + 0] = (uint8_t)(z1_rans_len & 0xFF);
  sig[OFF_Z1_HI_LEN + 1] = (uint8_t)((z1_rans_len >> 8) & 0xFF);
  if(z1_rans_len < SHUTTLE_Z1_RANS_RESERVED_BYTES)
    memset(&sig[OFF_Z1_HI_DATA + z1_rans_len], 0,
           SHUTTLE_Z1_RANS_RESERVED_BYTES - z1_rans_len);

  /* 6. rANS-encode hint h. */
  int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h_flat[i * SHUTTLE_N + j] = h->vec[i].coeffs[j];

  size_t hint_rans_len = 0;
  rc = shuttle_rans_encode(&sig[OFF_HINT_DATA], &hint_rans_len,
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

  /* 3. Z_0 rANS decode. */
  size_t z0_rans_len = (size_t)sig[OFF_Z0_LEN + 0]
                     | ((size_t)sig[OFF_Z0_LEN + 1] << 8);
  if(z0_rans_len == 0 || z0_rans_len > SHUTTLE_Z0_RANS_RESERVED_BYTES)
    return -1;

  int32_t z0_flat[SHUTTLE_N];
  rc = shuttle_rans_decode_z0(z0_flat, SHUTTLE_N,
                              &sig[OFF_Z0_DATA], z0_rans_len);
  if(rc != 0)
    return rc;
  for(i = 0; i < SHUTTLE_N; ++i)
    z_1[0].coeffs[i] = z0_flat[i];

  /* 4. z1 rANS high stream. */
  size_t z1_rans_len = (size_t)sig[OFF_Z1_HI_LEN + 0]
                     | ((size_t)sig[OFF_Z1_HI_LEN + 1] << 8);
  if(z1_rans_len == 0 || z1_rans_len > SHUTTLE_Z1_RANS_RESERVED_BYTES)
    return -1;

  int32_t z1_hi_flat[SHUTTLE_L * SHUTTLE_N];
  rc = shuttle_rans_decode_z1(z1_hi_flat, SHUTTLE_L * SHUTTLE_N,
                              &sig[OFF_Z1_HI_DATA], z1_rans_len);
  if(rc != 0)
    return rc;

  /* 5. z1 low parts: unpack and combine with the decoded high parts. */
  int32_t lo_scratch[SHUTTLE_N];
  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz1_lo_unpack(lo_scratch,
                     &sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES]);
    polyz1_combine(&z_1[1 + i], &z1_hi_flat[i * SHUTTLE_N], lo_scratch);
  }

  /* 6. hint rANS. */
  size_t hint_rans_len = (size_t)sig[OFF_HINT_LEN + 0]
                       | ((size_t)sig[OFF_HINT_LEN + 1] << 8);
  if(hint_rans_len == 0 || hint_rans_len > SHUTTLE_HINT_RESERVED_BYTES)
    return -1;

  int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
  rc = shuttle_rans_decode(h_flat, SHUTTLE_M * SHUTTLE_N,
                           &sig[OFF_HINT_DATA], hint_rans_len);
  if(rc != 0)
    return rc;

  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h->vec[i].coeffs[j] = h_flat[i * SHUTTLE_N + j];

  return 0;
}
