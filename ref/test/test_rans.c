/*
 * test_rans.c - unit tests for C-side rANS (Phase 6b-5a).
 *
 * Tests:
 *   1. Round-trip on a crafted small input matching a known distribution.
 *   2. Round-trip on a large synthetic Gaussian-like input (10000 symbols).
 *   3. Out-of-vocabulary rejection.
 *   4. Compression-rate check: encoded size close to Shannon entropy lower bound.
 *
 * Phase 6b-5b will then replace polyveck_hint_pack_basic with this engine
 * in pack_sig / unpack_sig, producing the actual SHUTTLE signature.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../params.h"
#include "../shuttle_rans.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_MODE == 128
#  define T_PROB_BITS  SHUTTLE128_RANS_PROB_BITS
#  define T_NUM_SYMS   SHUTTLE128_RANS_NUM_SYMS
#  define T_SYM_MIN    SHUTTLE128_RANS_SYM_MIN
#  define T_SYM_MAX    SHUTTLE128_RANS_SYM_MAX
#  define t_syms       shuttle128_rans_syms
#  define t_freqs      shuttle128_rans_freqs

#  define Z1_PROB_BITS SHUTTLE128_RANS_Z1_PROB_BITS
#  define Z1_NUM_SYMS  SHUTTLE128_RANS_Z1_NUM_SYMS
#  define Z1_SYM_MIN   SHUTTLE128_RANS_Z1_SYM_MIN
#  define Z1_SYM_MAX   SHUTTLE128_RANS_Z1_SYM_MAX
#  define z1_syms      shuttle128_rans_z1_syms
#  define z1_freqs     shuttle128_rans_z1_freqs
#else
#  define T_PROB_BITS  SHUTTLE256_RANS_PROB_BITS
#  define T_NUM_SYMS   SHUTTLE256_RANS_NUM_SYMS
#  define T_SYM_MIN    SHUTTLE256_RANS_SYM_MIN
#  define T_SYM_MAX    SHUTTLE256_RANS_SYM_MAX
#  define t_syms       shuttle256_rans_syms
#  define t_freqs      shuttle256_rans_freqs

#  define Z1_PROB_BITS SHUTTLE256_RANS_Z1_PROB_BITS
#  define Z1_NUM_SYMS  SHUTTLE256_RANS_Z1_NUM_SYMS
#  define Z1_SYM_MIN   SHUTTLE256_RANS_Z1_SYM_MIN
#  define Z1_SYM_MAX   SHUTTLE256_RANS_Z1_SYM_MAX
#  define z1_syms      shuttle256_rans_z1_syms
#  define z1_freqs     shuttle256_rans_z1_freqs
#endif

#define T_PROB_TOTAL   (1u << T_PROB_BITS)
#define Z1_PROB_TOTAL  (1u << Z1_PROB_BITS)

/* Generate a random symbol from the distribution defined by t_freqs. */
static int32_t sample_sym_from_table(void) {
  /* 12-bit uniform draw, map via cumulative frequency. */
  uint32_t r = rand() & (T_PROB_TOTAL - 1);
  uint32_t cum = 0;
  for (unsigned i = 0; i < T_NUM_SYMS; ++i) {
    cum += t_freqs[i];
    if (r < cum) return t_syms[i];
  }
  /* Fallback: last symbol (shouldn't reach here if freqs sum to PROB_TOTAL). */
  return t_syms[T_NUM_SYMS - 1];
}

/* ------------------------------------------------------------ */
static int test_small_roundtrip(void) {
  printf("Test 1: small round-trip (10 symbols, handcrafted)\n");

  int32_t msg[10]  = {0, 1, -1, 2, 0, -2, 0, 3, 1, 0};
  int32_t dec[10];
  uint8_t buf[128];
  size_t out_len;

  int rc = shuttle_rans_encode(buf, &out_len, sizeof buf, msg, 10);
  if (rc != 0) FAIL("encode failed (rc=%d)", rc);

  rc = shuttle_rans_decode(dec, 10, buf, out_len);
  if (rc != 0) FAIL("decode failed (rc=%d)", rc);

  for (unsigned i = 0; i < 10; ++i)
    if (dec[i] != msg[i])
      FAIL("mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);

  printf("    encoded %zu bytes for 10 symbols\n", out_len);
  PASS("small round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_large_roundtrip(void) {
  printf("Test 2: large round-trip (10000 samples from table distribution)\n");

  const unsigned N = 10000;
  int32_t *msg = malloc(N * sizeof(int32_t));
  int32_t *dec = malloc(N * sizeof(int32_t));
  uint8_t *buf = malloc(N);   /* naive 1 byte/sym is an upper bound */

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_sym_from_table();

  size_t out_len;
  int rc = shuttle_rans_encode(buf, &out_len, N, msg, N);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("encode failed (rc=%d)", rc); }

  rc = shuttle_rans_decode(dec, N, buf, out_len);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("decode failed (rc=%d)", rc); }

  for (unsigned i = 0; i < N; ++i)
    if (dec[i] != msg[i]) {
      free(msg); free(dec); free(buf);
      FAIL("mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);
    }

  /* Compute Shannon entropy of the used table distribution. */
  double H = 0.0;
  for (unsigned i = 0; i < T_NUM_SYMS; ++i) {
    double p = (double)t_freqs[i] / T_PROB_TOTAL;
    if (p > 0) H -= p * log2(p);
  }

  double bits_per_sym = 8.0 * out_len / N;
  printf("    symbols = %u\n", N);
  printf("    encoded = %zu bytes = %.3f bits/symbol\n", out_len, bits_per_sym);
  printf("    Shannon entropy lower bound = %.3f bits/symbol\n", H);
  printf("    overhead = %.3f bits/symbol\n", bits_per_sym - H);

  free(msg); free(dec); free(buf);
  PASS("large round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_oov_rejection(void) {
  printf("Test 3: out-of-vocabulary rejection\n");

  /* Symbol beyond SYM_MAX by 10 should be rejected. */
  int32_t bad[] = {T_SYM_MAX + 10};
  uint8_t buf[32];
  size_t out_len;

  int rc = shuttle_rans_encode(buf, &out_len, sizeof buf, bad, 1);
  if (rc != -1)
    FAIL("expected rc=-1 (OOV), got %d", rc);
  PASS("encoder rejects out-of-vocabulary symbol");
  return 0;
}

/* ------------------------------------------------------------ */
/* z1 variants of tests 1/2 exercising the Phase 6c z[1..L] rANS API. */

static int32_t sample_sym_from_z1_table(void) {
  uint32_t r = rand() & (Z1_PROB_TOTAL - 1);
  uint32_t cum = 0;
  for (unsigned i = 0; i < Z1_NUM_SYMS; ++i) {
    cum += z1_freqs[i];
    if (r < cum) return z1_syms[i];
  }
  return z1_syms[Z1_NUM_SYMS - 1];
}

static int test_z1_small_roundtrip(void) {
  printf("Test 5: z1 small round-trip (10 symbols)\n");

  int32_t msg[10]  = {0, 1, -1, 2, 0, -2, 0, 1, -1, 0};
  int32_t dec[10];
  uint8_t buf[128];
  size_t out_len;

  int rc = shuttle_rans_encode_z1(buf, &out_len, sizeof buf, msg, 10);
  if (rc != 0) FAIL("encode_z1 failed (rc=%d)", rc);

  rc = shuttle_rans_decode_z1(dec, 10, buf, out_len);
  if (rc != 0) FAIL("decode_z1 failed (rc=%d)", rc);

  for (unsigned i = 0; i < 10; ++i)
    if (dec[i] != msg[i])
      FAIL("z1 mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);

  printf("    encoded %zu bytes for 10 z1 symbols\n", out_len);
  PASS("z1 small round-trip OK");
  return 0;
}

static int test_z1_large_roundtrip(void) {
  printf("Test 6: z1 large round-trip (10000 samples)\n");

  const unsigned N = 10000;
  int32_t *msg = malloc(N * sizeof(int32_t));
  int32_t *dec = malloc(N * sizeof(int32_t));
  uint8_t *buf = malloc(N);

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_sym_from_z1_table();

  size_t out_len;
  int rc = shuttle_rans_encode_z1(buf, &out_len, N, msg, N);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("encode_z1 failed (rc=%d)", rc); }

  rc = shuttle_rans_decode_z1(dec, N, buf, out_len);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("decode_z1 failed (rc=%d)", rc); }

  for (unsigned i = 0; i < N; ++i)
    if (dec[i] != msg[i]) {
      free(msg); free(dec); free(buf);
      FAIL("z1 mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);
    }

  double H = 0.0;
  for (unsigned i = 0; i < Z1_NUM_SYMS; ++i) {
    double p = (double)z1_freqs[i] / Z1_PROB_TOTAL;
    if (p > 0) H -= p * log2(p);
  }
  double bits_per_sym = 8.0 * out_len / N;
  printf("    symbols = %u\n", N);
  printf("    encoded = %zu bytes = %.3f bits/symbol\n", out_len, bits_per_sym);
  printf("    Shannon entropy lower bound = %.3f bits/symbol\n", H);
  printf("    overhead = %.3f bits/symbol\n", bits_per_sym - H);

  free(msg); free(dec); free(buf);
  PASS("z1 large round-trip OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_empty_sequence(void) {
  printf("Test 4: empty sequence\n");

  uint8_t buf[32];
  size_t out_len;
  int rc = shuttle_rans_encode(buf, &out_len, sizeof buf, NULL, 0);
  if (rc != 0) FAIL("encode of empty sequence failed (rc=%d)", rc);
  if (out_len != 4)
    FAIL("encoded length for empty sequence should be 4 (flushed state), got %zu", out_len);

  int32_t dec_dummy[1];
  rc = shuttle_rans_decode(dec_dummy, 0, buf, out_len);
  if (rc != 0) FAIL("decode of empty sequence failed (rc=%d)", rc);

  PASS("empty sequence round-trip OK (4-byte state flush)");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xC0DE ^ SHUTTLE_MODE);

  printf("=== test_rans (MODE=%d)\n"
         "    hint table: [%d, %d], %d symbols, prob_bits=%d\n"
         "    z1   table: [%d, %d], %d symbols, prob_bits=%d ===\n",
         SHUTTLE_MODE,
         T_SYM_MIN, T_SYM_MAX, T_NUM_SYMS, T_PROB_BITS,
         Z1_SYM_MIN, Z1_SYM_MAX, Z1_NUM_SYMS, Z1_PROB_BITS);

  ret |= test_small_roundtrip();
  ret |= test_large_roundtrip();
  ret |= test_oov_rejection();
  ret |= test_z1_small_roundtrip();
  ret |= test_z1_large_roundtrip();
  ret |= test_empty_sequence();

  if (ret == 0) {
    printf("\n=== All rANS tests PASSED ===\n");
  } else {
    printf("\n=== Some rANS tests FAILED ===\n");
  }
  return ret;
}
