/*
 * packing_v2.c - NGCC-Signature Alg 2 compressed signature packing.
 *
 * See packing_v2.h for layout and error semantics.
 */

#include <string.h>

#include "packing_v2.h"
#include "shuttle_rans.h"
#include "poly.h"
#include "polyvec.h"

/* ============================================================
 * Byte-slot offsets within the fixed-size signature.
 * ============================================================ */
#define OFF_C_TILDE     0
#define OFF_IRS_SIGNS   (OFF_C_TILDE + SHUTTLE_CTILDEBYTES)
#define OFF_Z0          (OFF_IRS_SIGNS + SHUTTLE_IRS_SIGNBYTES)
#define OFF_Z_MID       (OFF_Z0 + SHUTTLE_POLYZ0_PACKEDBYTES)
#define OFF_HINT_LEN    (OFF_Z_MID + SHUTTLE_L * SHUTTLE_POLYZ_PACKEDBYTES)
#define OFF_HINT_DATA   (OFF_HINT_LEN + 2)
/* End offset = OFF_HINT_DATA + HINT_RESERVED_BYTES = SHUTTLE_BYTES_V2. */

/* ============================================================
 * pack_sig_v2
 * ============================================================ */
int pack_sig_v2(uint8_t sig[SHUTTLE_BYTES_V2],
                const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                const int8_t irs_signs[SHUTTLE_TAU],
                const poly *z_1,
                const polyveck *h)
{
  unsigned int i, j;

  /* 1. seedC */
  memcpy(&sig[OFF_C_TILDE], c_tilde, SHUTTLE_CTILDEBYTES);

  /* 2. irs_signs bitmap: bit i = 1 iff irs_signs[i] > 0. */
  memset(&sig[OFF_IRS_SIGNS], 0, SHUTTLE_IRS_SIGNBYTES);
  for(i = 0; i < SHUTTLE_TAU; ++i) {
    if(irs_signs[i] > 0)
      sig[OFF_IRS_SIGNS + (i >> 3)] |= (uint8_t)(1u << (i & 7));
  }

  /* 3. Z_0 at 11 bits/coef via polyz0_pack. */
  polyz0_pack(&sig[OFF_Z0], &z_1[0]);

  /* 4. z_1[1..L] at 14 bits/coef via polyz_pack. */
  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz_pack(&sig[OFF_Z_MID + i * SHUTTLE_POLYZ_PACKEDBYTES], &z_1[1 + i]);
  }

  /* 5. rANS-encode hint h. Flatten h->vec[i][j] into a 1D array. */
  int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h_flat[i * SHUTTLE_N + j] = h->vec[i].coeffs[j];

  size_t rans_len = 0;
  int rc = shuttle_rans_encode(&sig[OFF_HINT_DATA], &rans_len,
                               SHUTTLE_HINT_RESERVED_BYTES,
                               h_flat, SHUTTLE_M * SHUTTLE_N);
  if(rc != 0)
    return rc;  /* -1 OOV, -2 overflow; caller rejects signing round */

  /* Length prefix, little-endian uint16. Sanity: rans_len <= 65535. */
  sig[OFF_HINT_LEN + 0] = (uint8_t)(rans_len & 0xFF);
  sig[OFF_HINT_LEN + 1] = (uint8_t)((rans_len >> 8) & 0xFF);

  /* Zero-pad the unused tail. */
  if(rans_len < SHUTTLE_HINT_RESERVED_BYTES)
    memset(&sig[OFF_HINT_DATA + rans_len], 0,
           SHUTTLE_HINT_RESERVED_BYTES - rans_len);

  return 0;
}

/* ============================================================
 * unpack_sig_v2
 * ============================================================ */
int unpack_sig_v2(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                  int8_t irs_signs[SHUTTLE_TAU],
                  poly *z_1,
                  polyveck *h,
                  const uint8_t sig[SHUTTLE_BYTES_V2])
{
  unsigned int i, j;

  /* 1. seedC */
  memcpy(c_tilde, &sig[OFF_C_TILDE], SHUTTLE_CTILDEBYTES);

  /* 2. irs_signs: bit 1 -> +1, bit 0 -> -1. */
  for(i = 0; i < SHUTTLE_TAU; ++i) {
    if((sig[OFF_IRS_SIGNS + (i >> 3)] >> (i & 7)) & 1)
      irs_signs[i] = (int8_t)1;
    else
      irs_signs[i] = (int8_t)-1;
  }

  /* 3. Z_0 */
  polyz0_unpack(&z_1[0], &sig[OFF_Z0]);

  /* 4. z[1..L] */
  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz_unpack(&z_1[1 + i], &sig[OFF_Z_MID + i * SHUTTLE_POLYZ_PACKEDBYTES]);
  }

  /* 5. rANS-decode hint. */
  size_t rans_len = (size_t)sig[OFF_HINT_LEN + 0]
                  | ((size_t)sig[OFF_HINT_LEN + 1] << 8);
  if(rans_len == 0 || rans_len > SHUTTLE_HINT_RESERVED_BYTES)
    return -1;  /* malformed */

  int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
  int rc = shuttle_rans_decode(h_flat, SHUTTLE_M * SHUTTLE_N,
                               &sig[OFF_HINT_DATA], rans_len);
  if(rc != 0)
    return rc;

  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h->vec[i].coeffs[j] = h_flat[i * SHUTTLE_N + j];

  return 0;
}
