/*
 * test_sig_v2.c - round-trip tests for pack_sig_v2 / unpack_sig_v2 (Phase 6b-5b).
 *
 * This test does NOT exercise sign.c (Phase 6a) or any end-to-end signing
 * pipeline. It constructs synthetic inputs matching the ranges expected
 * from Alg 2 (after CompressY and MakeHint) and validates that pack +
 * unpack returns the original inputs byte-for-byte in the exposed fields.
 *
 * Tests:
 *   1. Pack/unpack a hand-crafted minimal input.
 *   2. Pack/unpack with hints drawn from the rANS frequency-table
 *      distribution (~4K random samples per round, 10 rounds).
 *   3. Pack/unpack with irs_signs alternating +/-1.
 *   4. Verify that OOV hint coefficients are rejected by pack_sig_v2.
 *
 * Phase 6b-5c will then wire the actual Alg 2 signing / verifying flow to
 * use these functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../params.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../packing_v2.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_MODE == 128
#  define T_NUM_SYMS   SHUTTLE128_RANS_NUM_SYMS
#  define T_SYM_MIN    SHUTTLE128_RANS_SYM_MIN
#  define T_SYM_MAX    SHUTTLE128_RANS_SYM_MAX
#  define t_syms       shuttle128_rans_syms
#  define t_freqs      shuttle128_rans_freqs
#  define T_PROB_BITS  SHUTTLE128_RANS_PROB_BITS
#else
#  define T_NUM_SYMS   SHUTTLE256_RANS_NUM_SYMS
#  define T_SYM_MIN    SHUTTLE256_RANS_SYM_MIN
#  define T_SYM_MAX    SHUTTLE256_RANS_SYM_MAX
#  define t_syms       shuttle256_rans_syms
#  define t_freqs      shuttle256_rans_freqs
#  define T_PROB_BITS  SHUTTLE256_RANS_PROB_BITS
#endif

#define T_PROB_TOTAL   (1u << T_PROB_BITS)

/* Sample a symbol from the rANS table distribution. */
static int32_t sample_table_sym(void) {
  uint32_t r = (uint32_t)rand() & (T_PROB_TOTAL - 1);
  uint32_t cum = 0;
  for (unsigned i = 0; i < T_NUM_SYMS; ++i) {
    cum += t_freqs[i];
    if (r < cum) return t_syms[i];
  }
  return t_syms[T_NUM_SYMS - 1];
}

/* Fill z_1 (length 1 + L) with inputs matching the expected post-CompressY /
 * post-IRS ranges. z_1[0] has |Z_0| <= ~170 (11-bit signed), z_1[i >= 1]
 * has |z| <= 11*sigma + tau*eta, within 14-bit polyz_pack range. */
static void fill_z1(poly *z_1) {
  /* Z_0 coefficients: symmetric small integers fitting polyz0_pack (11 bits
   * signed). Range ~ [-170, 170]. */
  for (unsigned j = 0; j < SHUTTLE_N; ++j)
    z_1[0].coeffs[j] = (rand() % 341) - 170;

  /* z[1..L]: range ~ [-SHUTTLE_Z_BOUND, SHUTTLE_Z_BOUND] fits polyz_pack. */
  for (unsigned i = 1; i <= SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j)
      z_1[i].coeffs[j] = (rand() % (2 * SHUTTLE_Z_BOUND + 1)) - SHUTTLE_Z_BOUND;
}

/* ------------------------------------------------------------ */
static int test_minimal(void) {
  printf("Test 1: minimal hand-crafted round-trip\n");

  uint8_t sig[SHUTTLE_BYTES_V2];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
  int8_t  irs[SHUTTLE_TAU];
  poly    z_1[1 + SHUTTLE_L];
  polyveck h, h_rec;

  /* c_tilde: deterministic pattern. */
  for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)(i * 17 + 3);
  /* irs_signs: alternating +/-. */
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = (i & 1) ? (int8_t)1 : (int8_t)-1;
  /* z_1: all zeros (simplest case). */
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  /* h: all zeros (table always contains symbol 0). */
  memset(&h, 0, sizeof h);

  int rc = pack_sig_v2(sig, c_tilde, irs, z_1, &h);
  if (rc != 0) FAIL("pack_sig_v2 rc=%d", rc);

  /* Reverse. */
  uint8_t c_rec[SHUTTLE_CTILDEBYTES];
  int8_t  irs_rec[SHUTTLE_TAU];
  poly    z_rec[1 + SHUTTLE_L];

  rc = unpack_sig_v2(c_rec, irs_rec, z_rec, &h_rec, sig);
  if (rc != 0) FAIL("unpack_sig_v2 rc=%d", rc);

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

  printf("    SHUTTLE_BYTES_V2 = %d\n", SHUTTLE_BYTES_V2);
  PASS("minimal round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_random_rounds(void) {
  printf("Test 2: 10 random rounds (z_1 full-range, h from rANS table)\n");

  for (unsigned round = 0; round < 10; ++round) {
    uint8_t sig[SHUTTLE_BYTES_V2];
    uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
    int8_t  irs[SHUTTLE_TAU];
    poly    z_1[1 + SHUTTLE_L];
    polyveck h, h_rec;

    for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)rand();
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = (rand() & 1) ? (int8_t)1 : (int8_t)-1;

    fill_z1(z_1);

    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        h.vec[i].coeffs[j] = sample_table_sym();

    int rc = pack_sig_v2(sig, c_tilde, irs, z_1, &h);
    if (rc != 0) FAIL("round %u: pack rc=%d", round, rc);

    uint8_t c_rec[SHUTTLE_CTILDEBYTES];
    int8_t  irs_rec[SHUTTLE_TAU];
    poly    z_rec[1 + SHUTTLE_L];
    rc = unpack_sig_v2(c_rec, irs_rec, z_rec, &h_rec, sig);
    if (rc != 0) FAIL("round %u: unpack rc=%d", round, rc);

    if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0) FAIL("round %u: c_tilde", round);
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
      if (irs[i] != irs_rec[i])
        FAIL("round %u: irs[%u]", round, i);
    for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        if (z_1[i].coeffs[j] != z_rec[i].coeffs[j])
          FAIL("round %u: z_1[%u][%u]", round, i, j);
    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        if (h.vec[i].coeffs[j] != h_rec.vec[i].coeffs[j])
          FAIL("round %u: h[%u][%u]", round, i, j);

    /* Peek at the rANS length prefix to report average compression. */
    size_t rans_len = (size_t)sig[SHUTTLE_CTILDEBYTES
                                   + SHUTTLE_IRS_SIGNBYTES
                                   + SHUTTLE_POLYZ0_PACKEDBYTES
                                   + SHUTTLE_L * SHUTTLE_POLYZ_PACKEDBYTES]
                    | ((size_t)sig[SHUTTLE_CTILDEBYTES
                                   + SHUTTLE_IRS_SIGNBYTES
                                   + SHUTTLE_POLYZ0_PACKEDBYTES
                                   + SHUTTLE_L * SHUTTLE_POLYZ_PACKEDBYTES + 1] << 8);
    if (round == 0)
      printf("    round %u rANS hint bytes: %zu (reserved %d)\n",
             round, rans_len, SHUTTLE_HINT_RESERVED_BYTES);
  }

  PASS("10 round round-trips OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_oov_rejected(void) {
  printf("Test 3: pack_sig_v2 rejects out-of-vocabulary hint symbol\n");

  uint8_t sig[SHUTTLE_BYTES_V2];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES] = {0};
  int8_t  irs[SHUTTLE_TAU];
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = 1;
  poly    z_1[1 + SHUTTLE_L];
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  polyveck h = {0};
  h.vec[0].coeffs[0] = T_SYM_MAX + 100;   /* deliberate OOV */

  int rc = pack_sig_v2(sig, c_tilde, irs, z_1, &h);
  if (rc == 0) FAIL("expected pack to fail, got rc=0%s", "");
  printf("    pack_sig_v2 returned rc=%d as expected\n", rc);
  PASS("OOV hint rejection");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xDECAF0 ^ SHUTTLE_MODE);

  printf("=== test_sig_v2 (MODE=%d, SHUTTLE_BYTES_V2=%d, HINT_RESERVED=%d) ===\n",
         SHUTTLE_MODE, SHUTTLE_BYTES_V2, SHUTTLE_HINT_RESERVED_BYTES);

  ret |= test_minimal();
  ret |= test_random_rounds();
  ret |= test_oov_rejected();

  if (ret == 0)
    printf("\n=== All sig_v2 tests PASSED ===\n");
  else
    printf("\n=== Some sig_v2 tests FAILED ===\n");

  return ret;
}
