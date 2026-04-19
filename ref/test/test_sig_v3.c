/*
 * test_sig_v3.c - round-trip tests for pack_sig_v3 / unpack_sig_v3 (Phase 6c-3).
 *
 * Constructs synthetic inputs matching the post-CompressY / post-MakeHint
 * ranges and validates byte-for-byte round-trip. The z[1..L] coefficients
 * are drawn from a triangular-sum distribution that keeps HighBits inside
 * the z1 rANS vocabulary; rare OOV throws are rejected before round-trip.
 *
 * Tests:
 *   1. All-zero minimal round-trip.
 *   2. 10 random rounds (z_1 from narrow-triangular distribution, h from
 *      the hint-table distribution).
 *   3. OOV rejection on the z1 path (a z_1[1][0] coefficient so large its
 *      HighBit falls outside the z1 vocabulary).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../params.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../packing_v3.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_MODE == 128
#  define HINT_NUM_SYMS   SHUTTLE128_RANS_NUM_SYMS
#  define hint_syms       shuttle128_rans_syms
#  define hint_freqs      shuttle128_rans_freqs
#  define HINT_PROB_BITS  SHUTTLE128_RANS_PROB_BITS
#  define Z1_NUM_SYMS     SHUTTLE128_RANS_Z1_NUM_SYMS
#  define z1_syms         shuttle128_rans_z1_syms
#  define z1_freqs        shuttle128_rans_z1_freqs
#  define Z1_PROB_BITS    SHUTTLE128_RANS_Z1_PROB_BITS
#  define Z1_SYM_MAX_ABS  SHUTTLE128_RANS_Z1_SYM_MAX
#else
#  define HINT_NUM_SYMS   SHUTTLE256_RANS_NUM_SYMS
#  define hint_syms       shuttle256_rans_syms
#  define hint_freqs      shuttle256_rans_freqs
#  define HINT_PROB_BITS  SHUTTLE256_RANS_PROB_BITS
#  define Z1_NUM_SYMS     SHUTTLE256_RANS_Z1_NUM_SYMS
#  define z1_syms         shuttle256_rans_z1_syms
#  define z1_freqs        shuttle256_rans_z1_freqs
#  define Z1_PROB_BITS    SHUTTLE256_RANS_Z1_PROB_BITS
#  define Z1_SYM_MAX_ABS  SHUTTLE256_RANS_Z1_SYM_MAX
#endif

#define HINT_PROB_TOTAL   (1u << HINT_PROB_BITS)
#define Z1_PROB_TOTAL     (1u << Z1_PROB_BITS)

static int32_t sample_hint_sym(void) {
  uint32_t r = (uint32_t)rand() & (HINT_PROB_TOTAL - 1);
  uint32_t cum = 0;
  for (unsigned i = 0; i < HINT_NUM_SYMS; ++i) {
    cum += hint_freqs[i];
    if (r < cum) return hint_syms[i];
  }
  return hint_syms[HINT_NUM_SYMS - 1];
}

/* Draw a z1 HighBit from the rANS table distribution. */
static int32_t sample_z1_hi(void) {
  uint32_t r = (uint32_t)rand() & (Z1_PROB_TOTAL - 1);
  uint32_t cum = 0;
  for (unsigned i = 0; i < Z1_NUM_SYMS; ++i) {
    cum += z1_freqs[i];
    if (r < cum) return z1_syms[i];
  }
  return z1_syms[Z1_NUM_SYMS - 1];
}

/* Fill z_1[1..L] matching the actual z1 distribution used by the signer:
 *   z[i] = hi * alpha_h + lo,  hi ~ z1 table,  lo ~ uniform [-alpha_h/2, alpha_h/2).
 * This guarantees pack_sig_v3 never hits OOV or overflow. */
static void fill_z1(poly *z_1) {
  /* Z_0 coefficients: fit polyz0_pack (11 bits signed). */
  for (unsigned j = 0; j < SHUTTLE_N; ++j)
    z_1[0].coeffs[j] = (rand() % 341) - 170;

  for (unsigned i = 1; i <= SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j) {
      int32_t hi = sample_z1_hi();
      int32_t lo = (rand() % SHUTTLE_ALPHA_H) - SHUTTLE_HALF_ALPHA_H;
      z_1[i].coeffs[j] = (hi << SHUTTLE_ALPHA_H_BITS) + lo;
    }
}

/* ------------------------------------------------------------ */
static int test_minimal(void) {
  printf("Test 1: minimal hand-crafted round-trip\n");

  uint8_t sig[SHUTTLE_BYTES_V3];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
  int8_t  irs[SHUTTLE_TAU];
  poly    z_1[1 + SHUTTLE_L];
  polyveck h, h_rec;

  for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)(i * 17 + 3);
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = (i & 1) ? (int8_t)1 : (int8_t)-1;
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  memset(&h, 0, sizeof h);

  int rc = pack_sig_v3(sig, c_tilde, irs, z_1, &h);
  if (rc != 0) FAIL("pack_sig_v3 rc=%d", rc);

  uint8_t c_rec[SHUTTLE_CTILDEBYTES];
  int8_t  irs_rec[SHUTTLE_TAU];
  poly    z_rec[1 + SHUTTLE_L];
  rc = unpack_sig_v3(c_rec, irs_rec, z_rec, &h_rec, sig);
  if (rc != 0) FAIL("unpack_sig_v3 rc=%d", rc);

  if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0)
    FAIL("c_tilde mismatch%s", "");
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
    if (irs[i] != irs_rec[i])
      FAIL("irs_signs[%u]: %d vs %d", i, irs[i], irs_rec[i]);
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j)
      if (z_1[i].coeffs[j] != z_rec[i].coeffs[j])
        FAIL("z_1[%u][%u]: %d vs %d",
             i, j, z_1[i].coeffs[j], z_rec[i].coeffs[j]);
  for (unsigned i = 0; i < SHUTTLE_M; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j)
      if (h.vec[i].coeffs[j] != h_rec.vec[i].coeffs[j])
        FAIL("h[%u][%u]: %d vs %d",
             i, j, h.vec[i].coeffs[j], h_rec.vec[i].coeffs[j]);

  printf("    SHUTTLE_BYTES_V3 = %d\n", SHUTTLE_BYTES_V3);
  PASS("minimal round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_random_rounds(void) {
  printf("Test 2: 10 random rounds\n");

  size_t min_z1_len = (size_t)-1, max_z1_len = 0, sum_z1_len = 0;
  size_t min_h_len  = (size_t)-1, max_h_len  = 0, sum_h_len  = 0;

  for (unsigned round = 0; round < 10; ++round) {
    uint8_t sig[SHUTTLE_BYTES_V3];
    uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
    int8_t  irs[SHUTTLE_TAU];
    poly    z_1[1 + SHUTTLE_L];
    polyveck h, h_rec;

    for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)rand();
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = (rand() & 1) ? (int8_t)1 : (int8_t)-1;

    fill_z1(z_1);

    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        h.vec[i].coeffs[j] = sample_hint_sym();

    int rc = pack_sig_v3(sig, c_tilde, irs, z_1, &h);
    if (rc != 0) FAIL("round %u: pack rc=%d", round, rc);

    uint8_t c_rec[SHUTTLE_CTILDEBYTES];
    int8_t  irs_rec[SHUTTLE_TAU];
    poly    z_rec[1 + SHUTTLE_L];
    rc = unpack_sig_v3(c_rec, irs_rec, z_rec, &h_rec, sig);
    if (rc != 0) FAIL("round %u: unpack rc=%d", round, rc);

    if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0) FAIL("round %u: c_tilde", round);
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
      if (irs[i] != irs_rec[i])
        FAIL("round %u: irs[%u]", round, i);
    for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        if (z_1[i].coeffs[j] != z_rec[i].coeffs[j])
          FAIL("round %u: z_1[%u][%u] %d vs %d",
               round, i, j, z_1[i].coeffs[j], z_rec[i].coeffs[j]);
    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        if (h.vec[i].coeffs[j] != h_rec.vec[i].coeffs[j])
          FAIL("round %u: h[%u][%u]", round, i, j);

    /* Peek at z1 rANS + hint rANS length prefixes. */
    size_t off_z1_len = SHUTTLE_CTILDEBYTES + SHUTTLE_IRS_SIGNBYTES
                      + SHUTTLE_POLYZ0_PACKEDBYTES
                      + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES;
    size_t z1_len = (size_t)sig[off_z1_len] | ((size_t)sig[off_z1_len + 1] << 8);
    size_t off_h_len  = off_z1_len + 2 + SHUTTLE_Z1_RANS_RESERVED_BYTES;
    size_t h_len  = (size_t)sig[off_h_len]  | ((size_t)sig[off_h_len  + 1] << 8);

    if (z1_len < min_z1_len) min_z1_len = z1_len;
    if (z1_len > max_z1_len) max_z1_len = z1_len;
    sum_z1_len += z1_len;
    if (h_len  < min_h_len)  min_h_len  = h_len;
    if (h_len  > max_h_len)  max_h_len  = h_len;
    sum_h_len  += h_len;
  }

  printf("    z1 rANS len: min=%zu max=%zu avg=%zu (reserved %d)\n",
         min_z1_len, max_z1_len, sum_z1_len / 10,
         SHUTTLE_Z1_RANS_RESERVED_BYTES);
  printf("    hint rANS len: min=%zu max=%zu avg=%zu (reserved %d)\n",
         min_h_len, max_h_len, sum_h_len / 10,
         SHUTTLE_HINT_RESERVED_BYTES);

  PASS("10 round round-trips OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_oov_z1(void) {
  printf("Test 3: pack_sig_v3 rejects out-of-vocabulary z1 HighBit\n");

  uint8_t sig[SHUTTLE_BYTES_V3];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES] = {0};
  int8_t  irs[SHUTTLE_TAU];
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = 1;
  poly    z_1[1 + SHUTTLE_L];
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  polyveck h = {0};

  /* Force a HighBit well outside the z1 vocabulary. */
  z_1[1].coeffs[0] = (int32_t)(Z1_SYM_MAX_ABS + 10) * SHUTTLE_ALPHA_H;

  int rc = pack_sig_v3(sig, c_tilde, irs, z_1, &h);
  if (rc == 0) FAIL("expected pack to fail, got rc=0%s", "");
  printf("    pack_sig_v3 returned rc=%d as expected\n", rc);
  PASS("OOV z1 rejection");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xC0DECAFE ^ SHUTTLE_MODE);

  printf("=== test_sig_v3 (MODE=%d, SHUTTLE_BYTES_V3=%d,\n"
         "                 Z1_RANS_RESERVED=%d, HINT_RESERVED=%d) ===\n",
         SHUTTLE_MODE, SHUTTLE_BYTES_V3,
         SHUTTLE_Z1_RANS_RESERVED_BYTES, SHUTTLE_HINT_RESERVED_BYTES);

  ret |= test_minimal();
  ret |= test_random_rounds();
  ret |= test_oov_z1();

  if (ret == 0) printf("\n=== All sig_v3 tests PASSED ===\n");
  else          printf("\n=== Some sig_v3 tests FAILED ===\n");
  return ret;
}
