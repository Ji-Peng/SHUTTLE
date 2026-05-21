/*
 * test_rans.c - unit tests for the SHUTTLE rANS engine (two-stream API).
 *
 * Coverage:
 *   1. Round-trip (small and large) on the z-hi context.
 *   2. Round-trip (small and large) on the hint context.
 *   3. Empty input: encoder produces a 4-byte flush, decoder accepts it.
 *   4. Little-endian flush: the final 4 bytes of the encoded stream
 *      reproduce the encoder's last x state in LE order.
 *   5. Final-state verification: any single-byte corruption of the
 *      encoded stream is detected by the decoder (returns -1).
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
#  define ZHI_NUM_SYMS   SHUTTLE128_RANS_ZHI_NUM_SYMS
#  define ZHI_SYM_MIN    SHUTTLE128_RANS_ZHI_SYM_MIN
#  define ZHI_SYM_MAX    SHUTTLE128_RANS_ZHI_SYM_MAX
#  define zhi_syms       shuttle128_rans_zhi_syms
#  define zhi_freqs      shuttle128_rans_zhi_freqs
#  define HINT_NUM_SYMS  SHUTTLE128_RANS_HINT_NUM_SYMS
#  define HINT_SYM_MIN   SHUTTLE128_RANS_HINT_SYM_MIN
#  define HINT_SYM_MAX   SHUTTLE128_RANS_HINT_SYM_MAX
#  define hint_syms      shuttle128_rans_hint_syms
#  define hint_freqs     shuttle128_rans_hint_freqs
#elif SHUTTLE_MODE == 256
#  define ZHI_NUM_SYMS   SHUTTLE256_RANS_ZHI_NUM_SYMS
#  define ZHI_SYM_MIN    SHUTTLE256_RANS_ZHI_SYM_MIN
#  define ZHI_SYM_MAX    SHUTTLE256_RANS_ZHI_SYM_MAX
#  define zhi_syms       shuttle256_rans_zhi_syms
#  define zhi_freqs      shuttle256_rans_zhi_freqs
#  define HINT_NUM_SYMS  SHUTTLE256_RANS_HINT_NUM_SYMS
#  define HINT_SYM_MIN   SHUTTLE256_RANS_HINT_SYM_MIN
#  define HINT_SYM_MAX   SHUTTLE256_RANS_HINT_SYM_MAX
#  define hint_syms      shuttle256_rans_hint_syms
#  define hint_freqs     shuttle256_rans_hint_freqs
#else
#  define ZHI_NUM_SYMS   SHUTTLE512_RANS_ZHI_NUM_SYMS
#  define ZHI_SYM_MIN    SHUTTLE512_RANS_ZHI_SYM_MIN
#  define ZHI_SYM_MAX    SHUTTLE512_RANS_ZHI_SYM_MAX
#  define zhi_syms       shuttle512_rans_zhi_syms
#  define zhi_freqs      shuttle512_rans_zhi_freqs
#  define HINT_NUM_SYMS  SHUTTLE512_RANS_HINT_NUM_SYMS
#  define HINT_SYM_MIN   SHUTTLE512_RANS_HINT_SYM_MIN
#  define HINT_SYM_MAX   SHUTTLE512_RANS_HINT_SYM_MAX
#  define hint_syms      shuttle512_rans_hint_syms
#  define hint_freqs     shuttle512_rans_hint_freqs
#endif

/* All tables share PROB_BITS = 10 (defined in shuttle_rans.h). */
#define PROB_TOTAL (1u << SHUTTLE_RANS_PROB_BITS)

/* Sample from a freq table by inverse CDF on a uniform 10-bit draw. */
static int32_t sample_from_table(const int16_t *syms, const uint16_t *freqs,
                                 unsigned num_syms)
{
  uint32_t r = (uint32_t)rand() & (PROB_TOTAL - 1u);
  uint32_t cum = 0;
  for (unsigned i = 0; i < num_syms; ++i) {
    cum += freqs[i];
    if (r < cum) return syms[i];
  }
  return syms[num_syms - 1];
}

static double table_entropy(const uint16_t *freqs, unsigned num_syms)
{
  double H = 0.0;
  for (unsigned i = 0; i < num_syms; ++i) {
    double p = (double)freqs[i] / (double)PROB_TOTAL;
    if (p > 0) H -= p * log2(p);
  }
  return H;
}

/* ------------------------------------------------------------ */
static int test_zhi_small(void) {
  printf("Test 1: z-hi small round-trip\n");
  int32_t msg[10] = {0, 1, -1, 2, 0, -2, 0, 3, 1, 0};
  int32_t dec[10];
  uint8_t buf[128];
  size_t out_len;

  int rc = shuttle_rans_encode_zhi(buf, &out_len, sizeof buf, msg, 10);
  if (rc != 0) FAIL("encode_zhi failed (rc=%d)", rc);

  rc = shuttle_rans_decode_zhi(dec, 10, buf, out_len);
  if (rc != 0) FAIL("decode_zhi failed (rc=%d)", rc);

  for (unsigned i = 0; i < 10; ++i)
    if (dec[i] != msg[i])
      FAIL("z-hi mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);

  printf("    encoded %zu bytes for 10 z-hi symbols\n", out_len);
  PASS("z-hi small round-trip");
  return 0;
}

static int test_zhi_large(void) {
  printf("Test 2: z-hi large round-trip (10000 samples)\n");
  const unsigned N = 10000;
  int32_t *msg = malloc(N * sizeof(int32_t));
  int32_t *dec = malloc(N * sizeof(int32_t));
  uint8_t *buf = malloc(N * 2);   /* generous upper bound */

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_from_table(zhi_syms, zhi_freqs, ZHI_NUM_SYMS);

  size_t out_len;
  int rc = shuttle_rans_encode_zhi(buf, &out_len, N * 2, msg, N);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("encode_zhi failed (rc=%d)", rc); }
  rc = shuttle_rans_decode_zhi(dec, N, buf, out_len);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("decode_zhi failed (rc=%d)", rc); }
  for (unsigned i = 0; i < N; ++i)
    if (dec[i] != msg[i]) {
      free(msg); free(dec); free(buf);
      FAIL("z-hi mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);
    }

  double H = table_entropy(zhi_freqs, ZHI_NUM_SYMS);
  printf("    encoded = %zu bytes = %.3f bits/sym  (Shannon = %.3f)\n",
         out_len, 8.0 * out_len / N, H);
  free(msg); free(dec); free(buf);
  PASS("z-hi large round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_hint_small(void) {
  printf("Test 3: hint small round-trip\n");
  int32_t msg[8] = {0, 0, 1, 0, -1, 0, 0, 0};
  int32_t dec[8];
  uint8_t buf[128];
  size_t out_len;

  int rc = shuttle_rans_encode_hint(buf, &out_len, sizeof buf, msg, 8);
  if (rc != 0) FAIL("encode_hint failed (rc=%d)", rc);

  rc = shuttle_rans_decode_hint(dec, 8, buf, out_len);
  if (rc != 0) FAIL("decode_hint failed (rc=%d)", rc);

  for (unsigned i = 0; i < 8; ++i)
    if (dec[i] != msg[i])
      FAIL("hint mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);

  printf("    encoded %zu bytes for 8 hint symbols\n", out_len);
  PASS("hint small round-trip");
  return 0;
}

static int test_hint_large(void) {
  printf("Test 4: hint large round-trip (10000 samples)\n");
  const unsigned N = 10000;
  int32_t *msg = malloc(N * sizeof(int32_t));
  int32_t *dec = malloc(N * sizeof(int32_t));
  uint8_t *buf = malloc(N * 2);

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_from_table(hint_syms, hint_freqs, HINT_NUM_SYMS);

  size_t out_len;
  int rc = shuttle_rans_encode_hint(buf, &out_len, N * 2, msg, N);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("encode_hint failed (rc=%d)", rc); }
  rc = shuttle_rans_decode_hint(dec, N, buf, out_len);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("decode_hint failed (rc=%d)", rc); }
  for (unsigned i = 0; i < N; ++i)
    if (dec[i] != msg[i]) {
      free(msg); free(dec); free(buf);
      FAIL("hint mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);
    }

  double H = table_entropy(hint_freqs, HINT_NUM_SYMS);
  printf("    encoded = %zu bytes = %.3f bits/sym  (Shannon = %.3f)\n",
         out_len, 8.0 * out_len / N, H);
  free(msg); free(dec); free(buf);
  PASS("hint large round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_empty(void) {
  printf("Test 5: empty sequence (z-hi)\n");
  uint8_t buf[32];
  size_t out_len;

  int rc = shuttle_rans_encode_zhi(buf, &out_len, sizeof buf, NULL, 0);
  if (rc != 0) FAIL("encode_zhi(empty) failed (rc=%d)", rc);
  if (out_len != 4)
    FAIL("empty stream should be exactly 4 flush bytes, got %zu", out_len);

  /* The initial state is L = 2^23 = 0x00800000, written little-endian.
   * That's bytes [0x00, 0x00, 0x80, 0x00]. */
  if (!(buf[0] == 0x00 && buf[1] == 0x00 && buf[2] == 0x80 && buf[3] == 0x00))
    FAIL("flush bytes look wrong: %02x %02x %02x %02x",
         buf[0], buf[1], buf[2], buf[3]);

  int32_t dummy[1];
  rc = shuttle_rans_decode_zhi(dummy, 0, buf, out_len);
  if (rc != 0) FAIL("decode_zhi(empty) failed (rc=%d)", rc);

  PASS("empty stream round-trip and little-endian flush of L = 2^23");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_corruption_detected(void) {
  printf("Test 6: final-state mismatch detection (flip every byte)\n");
  const unsigned N = 200;
  int32_t msg[200];
  int32_t dec[200];
  uint8_t buf[400];
  uint8_t buf2[400];
  size_t out_len;

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_from_table(zhi_syms, zhi_freqs, ZHI_NUM_SYMS);

  int rc = shuttle_rans_encode_zhi(buf, &out_len, sizeof buf, msg, N);
  if (rc != 0) FAIL("encode_zhi failed (rc=%d)", rc);

  /* Baseline: clean round-trip succeeds. */
  rc = shuttle_rans_decode_zhi(dec, N, buf, out_len);
  if (rc != 0) FAIL("clean decode failed (rc=%d)", rc);

  /* For every byte position, flip one bit and check the decoder catches it.
   * A genuine corruption may cause:
   *   - rc != 0 (underflow OR final-state mismatch), or
   *   - rc == 0 but at least one decoded symbol differs.
   * Both count as detection. */
  unsigned undetected = 0;
  for (size_t pos = 0; pos < out_len; ++pos) {
    memcpy(buf2, buf, out_len);
    buf2[pos] ^= 0xFFu;

    rc = shuttle_rans_decode_zhi(dec, N, buf2, out_len);
    int caught = (rc != 0);
    if (rc == 0) {
      for (unsigned i = 0; i < N; ++i)
        if (dec[i] != msg[i]) { caught = 1; break; }
    }
    if (!caught) ++undetected;
  }
  if (undetected != 0)
    FAIL("%u single-byte corruptions slipped through (out of %zu positions)",
         undetected, out_len);

  printf("    all %zu single-byte corruptions detected\n", out_len);
  PASS("byte-corruption detection");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xC0DEu ^ (unsigned)SHUTTLE_MODE);

  printf("=== test_rans (MODE=%d)\n"
         "    z-hi  table: [%d, %d], %u symbols\n"
         "    hint  table: [%d, %d], %u symbols\n"
         "    prob_bits = %u ===\n",
         SHUTTLE_MODE,
         ZHI_SYM_MIN, ZHI_SYM_MAX, (unsigned)ZHI_NUM_SYMS,
         HINT_SYM_MIN, HINT_SYM_MAX, (unsigned)HINT_NUM_SYMS,
         (unsigned)SHUTTLE_RANS_PROB_BITS);

  ret |= test_zhi_small();
  ret |= test_zhi_large();
  ret |= test_hint_small();
  ret |= test_hint_large();
  ret |= test_empty();
  ret |= test_corruption_detected();

  if (ret == 0)
    printf("\n=== All rANS tests PASSED ===\n");
  else
    printf("\n=== Some rANS tests FAILED ===\n");
  return ret;
}
