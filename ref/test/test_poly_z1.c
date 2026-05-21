/*
 * test_poly_z1.c - unit tests for polyz1_split (alpha_r) and polyz0_split
 *                  (alpha_0' = alpha_r/alpha_1).
 *
 * Covers:
 *   1. polyz1_split + polyz1_combine: bijective round-trip.
 *   2. polyz1_lo_pack + polyz1_lo_unpack: byte-round-trip across the
 *      full [-alpha_r/2, alpha_r/2) range.
 *   3. polyz0_split + polyz0_combine: bijective round-trip.
 *   4. polyz0_lo_pack + polyz0_lo_unpack: byte-round-trip.
 *   5. Full z-hi pipeline (z^(0) + z^(1..lenS)): split -> rANS encode +
 *      bit-pack -> rANS decode + bit-unpack -> combine recovers all
 *      polynomials. Exercises the unified z-hi rANS stream.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../params.h"
#include "../poly.h"
#include "../shuttle_rans.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_MODE == 128
#  define ZHI_SYM_MAX_ABS  SHUTTLE128_RANS_ZHI_SYM_MAX
#elif SHUTTLE_MODE == 256
#  define ZHI_SYM_MAX_ABS  SHUTTLE256_RANS_ZHI_SYM_MAX
#else
#  define ZHI_SYM_MAX_ABS  SHUTTLE512_RANS_ZHI_SYM_MAX
#endif

/* Deterministic PRNG so the test is reproducible. */
static uint32_t rng_state = 0xC001C0DEu;
static uint32_t next_u32(void) {
  rng_state = rng_state * 1664525u + 1013904223u;
  return rng_state;
}
static int32_t rand_range(int32_t lo, int32_t hi) {
  int32_t span = hi - lo + 1;
  return lo + (int32_t)(next_u32() % (uint32_t)span);
}

/* Draw a z^(i) coefficient: sum of three uniforms in [-alpha_r, alpha_r]
 * yields a triangular-ish bounded distribution that exercises the
 * z-hi vocabulary without straying past it (|hi| stays inside +/-3). */
static int32_t sample_zi_coef(void) {
  int32_t H = SHUTTLE_ALPHA_R;
  return rand_range(-H, H) + rand_range(-H, H) + rand_range(-H, H);
}

/* Draw a z^(0) coefficient: smaller magnitude after CompressY's alpha_1
 * compression. Sum-of-three with H = alpha_0' gives |hi| in low range. */
static int32_t sample_z0_coef(void) {
  int32_t H = SHUTTLE_ALPHA_0P;
  return rand_range(-H, H) + rand_range(-H, H) + rand_range(-H, H);
}

/* ------------------------------------------------------------ */
static int test_z1_split_combine(void) {
  printf("Test 1: polyz1_split + polyz1_combine round-trip\n");

  poly a, a_back;
  int32_t hi[SHUTTLE_N], lo[SHUTTLE_N];

  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    a.coeffs[i] = sample_zi_coef();

  polyz1_split(hi, lo, &a);

  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    if (lo[i] < -(int32_t)SHUTTLE_HALF_ALPHA_R
        || lo[i] >= (int32_t)SHUTTLE_HALF_ALPHA_R)
      FAIL("lo[%u]=%d out of [-alpha_r/2, alpha_r/2)", i, lo[i]);
  }

  polyz1_combine(&a_back, hi, lo);
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (a_back.coeffs[i] != a.coeffs[i])
      FAIL("combine mismatch at %u: a=%d back=%d", i, a.coeffs[i], a_back.coeffs[i]);

  PASS("z1 split/combine round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_z1_lo_pack(void) {
  printf("Test 2: polyz1_lo_pack + polyz1_lo_unpack round-trip\n");
  int32_t lo[SHUTTLE_N], lo_back[SHUTTLE_N];
  uint8_t buf[SHUTTLE_POLYZ1_LO_PACKEDBYTES];

  int32_t v = -((int32_t)SHUTTLE_HALF_ALPHA_R);
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    lo[i] = v;
    ++v;
    if (v >= (int32_t)SHUTTLE_HALF_ALPHA_R)
      v = -((int32_t)SHUTTLE_HALF_ALPHA_R);
  }

  polyz1_lo_pack(buf, lo);
  polyz1_lo_unpack(lo_back, buf);
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (lo[i] != lo_back[i])
      FAIL("z1 lo roundtrip fail at %u: orig=%d unpacked=%d",
           i, lo[i], lo_back[i]);

  PASS("z1 lo pack/unpack round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_z0_split_combine(void) {
  printf("Test 3: polyz0_split + polyz0_combine round-trip\n");

  poly a, a_back;
  int32_t hi[SHUTTLE_N], lo[SHUTTLE_N];

  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    a.coeffs[i] = sample_z0_coef();

  polyz0_split(hi, lo, &a);

  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    if (lo[i] < -(int32_t)SHUTTLE_HALF_ALPHA_0P
        || lo[i] >= (int32_t)SHUTTLE_HALF_ALPHA_0P)
      FAIL("z0 lo[%u]=%d out of [-alpha_0p/2, alpha_0p/2)", i, lo[i]);
  }

  polyz0_combine(&a_back, hi, lo);
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (a_back.coeffs[i] != a.coeffs[i])
      FAIL("z0 combine mismatch at %u: a=%d back=%d",
           i, a.coeffs[i], a_back.coeffs[i]);

  PASS("z0 split/combine round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_z0_lo_pack(void) {
  printf("Test 4: polyz0_lo_pack + polyz0_lo_unpack round-trip\n");
  int32_t lo[SHUTTLE_N], lo_back[SHUTTLE_N];
  uint8_t buf[SHUTTLE_POLYZ0_LO_PACKEDBYTES];

  int32_t v = -((int32_t)SHUTTLE_HALF_ALPHA_0P);
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    lo[i] = v;
    ++v;
    if (v >= (int32_t)SHUTTLE_HALF_ALPHA_0P)
      v = -((int32_t)SHUTTLE_HALF_ALPHA_0P);
  }

  polyz0_lo_pack(buf, lo);
  polyz0_lo_unpack(lo_back, buf);
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    if (lo[i] != lo_back[i])
      FAIL("z0 lo roundtrip fail at %u: orig=%d unpacked=%d",
           i, lo[i], lo_back[i]);

  PASS("z0 lo pack/unpack round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_zhi_pipeline(void) {
  printf("Test 5: full z-hi pipeline (z^(0) + z^(1..lenS) unified stream)\n");

  poly z[SHUTTLE_L + 1], z_back[SHUTTLE_L + 1];
  int32_t hi[(SHUTTLE_L + 1) * SHUTTLE_N];
  int32_t hi_back[(SHUTTLE_L + 1) * SHUTTLE_N];
  int32_t lo_scratch[SHUTTLE_N], lo_back[SHUTTLE_N];
  uint8_t z0_lo_buf[SHUTTLE_POLYZ0_LO_PACKEDBYTES];
  uint8_t z1_lo_bufs[SHUTTLE_L][SHUTTLE_POLYZ1_LO_PACKEDBYTES];
  uint8_t rans_buf[SHUTTLE_ZHI_RESERVED_BYTES];

  /* Populate z^(0) and z^(1..lenS) with controlled samples, then split. */
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    z[0].coeffs[i] = sample_z0_coef();
  for (unsigned k = 1; k <= SHUTTLE_L; ++k)
    for (unsigned i = 0; i < SHUTTLE_N; ++i)
      z[k].coeffs[i] = sample_zi_coef();

  /* Encoder side: split each poly, accumulate hi into one flat array,
   * bit-pack lo per poly. */
  polyz0_split(&hi[0], lo_scratch, &z[0]);
  polyz0_lo_pack(z0_lo_buf, lo_scratch);
  for (unsigned k = 0; k < SHUTTLE_L; ++k) {
    polyz1_split(&hi[(k + 1) * SHUTTLE_N], lo_scratch, &z[1 + k]);
    polyz1_lo_pack(z1_lo_bufs[k], lo_scratch);
  }

  /* Sanity: every hi value must already lie in the vocabulary by the
   * 11-sigma tight bound. Our triangular sampler may sometimes overshoot,
   * so we clamp + recombine + re-split as the real signer would not. */
  unsigned clipped = 0;
  for (unsigned i = 0; i < (SHUTTLE_L + 1) * SHUTTLE_N; ++i) {
    if (hi[i] >  ZHI_SYM_MAX_ABS) { hi[i] =  ZHI_SYM_MAX_ABS; ++clipped; }
    if (hi[i] < -ZHI_SYM_MAX_ABS) { hi[i] = -ZHI_SYM_MAX_ABS; ++clipped; }
  }
  if (clipped > 0) {
    /* Recombine the test input from clipped hi + original lo so we end up
     * comparing against the right reference (the signer rejection path
     * is irrelevant to this unit test). */
    polyz0_lo_unpack(lo_scratch, z0_lo_buf);
    polyz0_combine(&z[0], &hi[0], lo_scratch);
    for (unsigned k = 0; k < SHUTTLE_L; ++k) {
      polyz1_lo_unpack(lo_scratch, z1_lo_bufs[k]);
      polyz1_combine(&z[1 + k], &hi[(k + 1) * SHUTTLE_N], lo_scratch);
    }
  }

  size_t rans_len;
  int rc = shuttle_rans_encode_zhi(rans_buf, &rans_len, sizeof rans_buf,
                                   hi, (SHUTTLE_L + 1) * SHUTTLE_N);
  if (rc != 0)
    FAIL("rANS encode failed (rc=%d, clipped=%u)", rc, clipped);

  rc = shuttle_rans_decode_zhi(hi_back, (SHUTTLE_L + 1) * SHUTTLE_N,
                               rans_buf, rans_len);
  if (rc != 0) FAIL("rANS decode failed (rc=%d)", rc);

  /* Decoder side: unpack lo, recombine, compare. */
  polyz0_lo_unpack(lo_back, z0_lo_buf);
  polyz0_combine(&z_back[0], &hi_back[0], lo_back);
  for (unsigned k = 0; k < SHUTTLE_L; ++k) {
    polyz1_lo_unpack(lo_back, z1_lo_bufs[k]);
    polyz1_combine(&z_back[1 + k], &hi_back[(k + 1) * SHUTTLE_N], lo_back);
  }

  for (unsigned k = 0; k <= SHUTTLE_L; ++k)
    for (unsigned i = 0; i < SHUTTLE_N; ++i)
      if (z_back[k].coeffs[i] != z[k].coeffs[i])
        FAIL("pipeline mismatch at z[%u].coeffs[%u]: orig=%d back=%d",
             k, i, z[k].coeffs[i], z_back[k].coeffs[i]);

  printf("    z-hi rANS: %zu bytes for %u coefs (%.3f bit/coef)\n",
         rans_len, (SHUTTLE_L + 1) * SHUTTLE_N,
         8.0 * rans_len / ((SHUTTLE_L + 1) * SHUTTLE_N));
  printf("    z0 lo pack: %zu bytes, z1 lo pack: %zu bytes/poly\n",
         (size_t)SHUTTLE_POLYZ0_LO_PACKEDBYTES,
         (size_t)SHUTTLE_POLYZ1_LO_PACKEDBYTES);
  PASS("unified z-hi pipeline OK");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;

  printf("=== test_poly_z1 (MODE=%d, N=%d, alpha_r=%d, alpha_0'=%d) ===\n",
         SHUTTLE_MODE, SHUTTLE_N, SHUTTLE_ALPHA_R, SHUTTLE_ALPHA_0P);

  ret |= test_z1_split_combine();
  ret |= test_z1_lo_pack();
  ret |= test_z0_split_combine();
  ret |= test_z0_lo_pack();
  ret |= test_zhi_pipeline();

  if (ret == 0) printf("\n=== All polyz1 tests PASSED ===\n");
  else          printf("\n=== Some polyz1 tests FAILED ===\n");
  return ret;
}
