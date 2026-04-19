/*
 * packing_v3.c - Phase 6c compressed signature packing.
 *
 * See packing_v3.h for layout, reserved-size budgets, and error semantics.
 */

#include <string.h>

#include "packing_v3.h"
#include "shuttle_rans.h"
#include "poly.h"
#include "polyvec.h"

/* ============================================================
 * Byte-slot offsets within SHUTTLE_BYTES_V3.
 * ============================================================ */
#define OFF_C_TILDE     0
#define OFF_IRS_SIGNS   (OFF_C_TILDE + SHUTTLE_CTILDEBYTES)
#define OFF_Z0          (OFF_IRS_SIGNS + SHUTTLE_IRS_SIGNBYTES)
#define OFF_Z1_LO       (OFF_Z0 + SHUTTLE_POLYZ0_PACKEDBYTES)
#define OFF_Z1_HI_LEN   (OFF_Z1_LO + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES)
#define OFF_Z1_HI_DATA  (OFF_Z1_HI_LEN + 2)
#define OFF_HINT_LEN    (OFF_Z1_HI_DATA + SHUTTLE_Z1_RANS_RESERVED_BYTES)
#define OFF_HINT_DATA   (OFF_HINT_LEN + 2)
/* End offset = OFF_HINT_DATA + HINT_RESERVED_BYTES = SHUTTLE_BYTES_V3. */

/* ============================================================
 * pack_sig_v3
 * ============================================================ */
int pack_sig_v3(uint8_t sig[SHUTTLE_BYTES_V3],
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

  /* 3. Z_0 at 11 bits/coef via polyz0_pack. */
  polyz0_pack(&sig[OFF_Z0], &z_1[0]);

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
    return rc;  /* -1 OOV, -2 overflow; caller rejects the signing round. */

  sig[OFF_Z1_HI_LEN + 0] = (uint8_t)(z1_rans_len & 0xFF);
  sig[OFF_Z1_HI_LEN + 1] = (uint8_t)((z1_rans_len >> 8) & 0xFF);
  if(z1_rans_len < SHUTTLE_Z1_RANS_RESERVED_BYTES)
    memset(&sig[OFF_Z1_HI_DATA + z1_rans_len], 0,
           SHUTTLE_Z1_RANS_RESERVED_BYTES - z1_rans_len);

  /* 6. rANS-encode hint h into the hint reservation (unchanged from v2). */
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

/* ============================================================
 * unpack_sig_v3
 * ============================================================ */
int unpack_sig_v3(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                  int8_t irs_signs[SHUTTLE_TAU],
                  poly *z_1,
                  polyveck *h,
                  const uint8_t sig[SHUTTLE_BYTES_V3])
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

  /* 3. Z_0 */
  polyz0_unpack(&z_1[0], &sig[OFF_Z0]);

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
