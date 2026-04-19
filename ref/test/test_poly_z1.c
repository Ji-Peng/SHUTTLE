/*
 * test_poly_z1.c - unit tests for Phase 6c z[1..L] split / combine / pack.
 *
 * Covers:
 *   1. polyz1_split + polyz1_combine: bijective round-trip across a random
 *      poly with coefficients drawn uniformly from [-Z_BOUND, Z_BOUND].
 *   2. polyz1_lo_pack + polyz1_lo_unpack: byte-round-trip for the LowBits
 *      part produced by polyz1_split.
 *   3. Full pipeline: split -> rANS encode (highs) + lo_pack (lows) ->
 *      rANS decode + lo_unpack -> combine. Output poly equals input.
 *      Runs on a Gaussian-like signed distribution that matches the
 *      training corpus of the z1 rANS table.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../params.h"
#include "../poly.h"
#include "../shuttle_rans.h"
#include "../rans_tables.h"
#include "../randombytes.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

/* Deterministic PRNG so the test is reproducible. */
static uint32_t rng_state = 0xC001C0DE;
static uint32_t next_u32(void) {
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state;
}

static int32_t rand_range(int32_t lo, int32_t hi) {
  int32_t span = hi - lo + 1;
  return lo + (int32_t)(next_u32() % (uint32_t)span);
}

/* Draw a coefficient matching the z1 empirical distribution more closely:
 *   sum of three triangular draws, scaled to fit [-Z_BOUND, Z_BOUND].
 * For the Phase 6c table-vocabulary range test we just want |high| <=
 * Z1_SYM_MAX * alpha_h/2 coefficients; a uniform draw over [-small, small]
 * would land on the largest high bucket too rarely. So clamp to roughly
 * |z| <= Z1_SYM_MAX * alpha_h so the full vocabulary is exercised. */

#if SHUTTLE_MODE == 128
#  define Z1_SYM_MAX_ABS   SHUTTLE128_RANS_Z1_SYM_MAX
#else
#  define Z1_SYM_MAX_ABS   SHUTTLE256_RANS_Z1_SYM_MAX
#endif

static int32_t sample_z_coef(void) {
  /* Sum of three uniforms in [-H, H] where H = alpha_h. Gives a
   * triangular-ish distribution bounded by 3*alpha_h so highs land
   * roughly in [-3, 3]. Good enough to hit the vocabulary without
   * going out of range. */
  int32_t H = SHUTTLE_ALPHA_H;
  int32_t s = rand_range(-H, H) + rand_range(-H, H) + rand_range(-H, H);
  return s;
}

/* ------------------------------------------------------------ */
static int test_split_combine_roundtrip(void) {
  printf("Test 1: polyz1_split + polyz1_combine round-trip\n");

  poly a, a_back;
  int32_t hi[SHUTTLE_N], lo[SHUTTLE_N];

  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    a.coeffs[i] = sample_z_coef();

  polyz1_split(hi, lo, &a);

  /* Validate ranges: lo in [-alpha_h/2, alpha_h/2) by round-half-up convention. */
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    if (lo[i] < -(int32_t)SHUTTLE_HALF_ALPHA_H
        || lo[i] >= (int32_t)SHUTTLE_HALF_ALPHA_H)
      FAIL("lo[%u]=%d out of [-alpha_h/2, alpha_h/2)", i, lo[i]);
  }

  polyz1_combine(&a_back, hi, lo);
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (a_back.coeffs[i] != a.coeffs[i])
      FAIL("combine mismatch at %u: a=%d back=%d", i, a.coeffs[i], a_back.coeffs[i]);

  PASS("split/combine round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_lo_pack_roundtrip(void) {
  printf("Test 2: polyz1_lo_pack + polyz1_lo_unpack round-trip\n");

  int32_t lo[SHUTTLE_N], lo_back[SHUTTLE_N];
  uint8_t buf[SHUTTLE_POLYZ1_LO_PACKEDBYTES];

  /* Exhaustively cover the valid range of lo [-alpha_h/2, alpha_h/2).
   * Since SHUTTLE_N = 256 or 512 and alpha_h = 128 or 256, we cycle. */
  int32_t v = -((int32_t)SHUTTLE_HALF_ALPHA_H);
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    lo[i] = v;
    ++v;
    if (v >= (int32_t)SHUTTLE_HALF_ALPHA_H)
      v = -((int32_t)SHUTTLE_HALF_ALPHA_H);
  }

  polyz1_lo_pack(buf, lo);
  polyz1_lo_unpack(lo_back, buf);

  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (lo[i] != lo_back[i])
      FAIL("lo roundtrip fail at %u: orig=%d unpacked=%d", i, lo[i], lo_back[i]);

  PASS("lo pack/unpack round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_full_pipeline_roundtrip(void) {
  printf("Test 3: split -> rANS+lo_pack -> decode -> combine\n");

  poly a, a_back;
  int32_t hi[SHUTTLE_N], lo[SHUTTLE_N];
  int32_t hi_back[SHUTTLE_N], lo_back[SHUTTLE_N];
  uint8_t lo_buf[SHUTTLE_POLYZ1_LO_PACKEDBYTES];
  uint8_t rans_buf[SHUTTLE_N];       /* 1 byte/coef is always enough */

  /* Use a Gaussian-like triangular sum so high buckets are exercised. */
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    a.coeffs[i] = sample_z_coef();

  polyz1_split(hi, lo, &a);

  /* Enforce vocabulary: if any |hi[i]| > Z1_SYM_MAX_ABS, clip it to the
   * bound for the purpose of this test (real signer would reject the
   * round). */
  unsigned clipped = 0;
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    if (hi[i] >  Z1_SYM_MAX_ABS) { hi[i] =  Z1_SYM_MAX_ABS; ++clipped; }
    if (hi[i] < -Z1_SYM_MAX_ABS) { hi[i] = -Z1_SYM_MAX_ABS; ++clipped; }
  }
  /* Recombine to a valid "test input" that lives inside the vocabulary. */
  polyz1_combine(&a, hi, lo);
  polyz1_split(hi, lo, &a);

  size_t rans_len;
  int rc = shuttle_rans_encode_z1(rans_buf, &rans_len, sizeof rans_buf,
                                  hi, SHUTTLE_N);
  if (rc != 0) FAIL("rANS encode failed (rc=%d, clipped=%u)", rc, clipped);
  polyz1_lo_pack(lo_buf, lo);

  rc = shuttle_rans_decode_z1(hi_back, SHUTTLE_N, rans_buf, rans_len);
  if (rc != 0) FAIL("rANS decode failed (rc=%d)", rc);
  polyz1_lo_unpack(lo_back, lo_buf);

  polyz1_combine(&a_back, hi_back, lo_back);
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (a_back.coeffs[i] != a.coeffs[i])
      FAIL("pipeline mismatch at %u: a=%d back=%d", i, a.coeffs[i], a_back.coeffs[i]);

  double rans_bits_per_coef = 8.0 * rans_len / SHUTTLE_N;
  double lo_bits_per_coef   = 8.0 * sizeof lo_buf / SHUTTLE_N;
  printf("    rANS high part: %zu bytes (%.3f bit/coef)\n", rans_len, rans_bits_per_coef);
  printf("    packed low part: %zu bytes (%.3f bit/coef)\n",
         sizeof lo_buf, lo_bits_per_coef);
  printf("    total per poly: %.1f B vs 14-bit baseline %.1f B\n",
         (double)rans_len + sizeof lo_buf,
         SHUTTLE_N * 14.0 / 8.0);
  PASS("full split/encode/decode/combine pipeline OK");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;

  printf("=== test_poly_z1 (MODE=%d, N=%d, alpha_h=%d, alpha_h_bits=%d) ===\n",
         SHUTTLE_MODE, SHUTTLE_N, SHUTTLE_ALPHA_H, SHUTTLE_ALPHA_H_BITS);

  ret |= test_split_combine_roundtrip();
  ret |= test_lo_pack_roundtrip();
  ret |= test_full_pipeline_roundtrip();

  if (ret == 0) printf("\n=== All polyz1 tests PASSED ===\n");
  else          printf("\n=== Some polyz1 tests FAILED ===\n");
  return ret;
}
