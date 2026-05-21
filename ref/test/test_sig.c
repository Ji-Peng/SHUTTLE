/*
 * test_sig.c - round-trip tests for pack_sig / unpack_sig under the
 *              two-stream rANS layout (z-hi + hint).
 *
 * Tests:
 *   1. Minimal hand-crafted round-trip (all-zero z and h).
 *   2. 10 random rounds: z_1 drawn so HighBits live in the z-hi
 *      vocabulary, h drawn from the hint table distribution; verifies
 *      byte-for-byte round-trip and reports per-stream rANS lengths.
 *
 * Note: there is no OOV test here. The theoretical vocabulary covers
 * every legal coefficient (SHUTTLE_rANS.tex §3 tight bound), and the
 * encoder asserts on out-of-range symbols in debug builds rather than
 * returning a soft error.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../params.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../packing.h"
#include "../shuttle_rans.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_MODE == 128
#  define ZHI_NUM_SYMS   SHUTTLE128_RANS_ZHI_NUM_SYMS
#  define zhi_syms       shuttle128_rans_zhi_syms
#  define zhi_freqs      shuttle128_rans_zhi_freqs
#  define HINT_NUM_SYMS  SHUTTLE128_RANS_HINT_NUM_SYMS
#  define hint_syms      shuttle128_rans_hint_syms
#  define hint_freqs     shuttle128_rans_hint_freqs
#elif SHUTTLE_MODE == 256
#  define ZHI_NUM_SYMS   SHUTTLE256_RANS_ZHI_NUM_SYMS
#  define zhi_syms       shuttle256_rans_zhi_syms
#  define zhi_freqs      shuttle256_rans_zhi_freqs
#  define HINT_NUM_SYMS  SHUTTLE256_RANS_HINT_NUM_SYMS
#  define hint_syms      shuttle256_rans_hint_syms
#  define hint_freqs     shuttle256_rans_hint_freqs
#else
#  define ZHI_NUM_SYMS   SHUTTLE512_RANS_ZHI_NUM_SYMS
#  define zhi_syms       shuttle512_rans_zhi_syms
#  define zhi_freqs      shuttle512_rans_zhi_freqs
#  define HINT_NUM_SYMS  SHUTTLE512_RANS_HINT_NUM_SYMS
#  define hint_syms      shuttle512_rans_hint_syms
#  define hint_freqs     shuttle512_rans_hint_freqs
#endif

#define PROB_TOTAL (1u << SHUTTLE_RANS_PROB_BITS)

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

/* Pick the per-context vocab radii. */
#if SHUTTLE_MODE == 128
#  define ZHI_M_VOC   SHUTTLE128_RANS_ZHI_SYM_MAX
#  define HINT_M_VOC  SHUTTLE128_RANS_HINT_SYM_MAX
#elif SHUTTLE_MODE == 256
#  define ZHI_M_VOC   SHUTTLE256_RANS_ZHI_SYM_MAX
#  define HINT_M_VOC  SHUTTLE256_RANS_HINT_SYM_MAX
#else
#  define ZHI_M_VOC   SHUTTLE512_RANS_ZHI_SYM_MAX
#  define HINT_M_VOC  SHUTTLE512_RANS_HINT_SYM_MAX
#endif

/* Discrete-Gaussian sampler over [-M_voc, M_voc] via inverse CDF on the
 * theoretical PMF p(k) ∝ exp(-k^2 / (2 sigma^2)).
 *
 * The real signer samples y ~ D_{Z, r}; after the alpha_r split the
 * HighBits follow this discrete distribution exactly (up to integer
 * truncation). A continuous Gaussian rounded to integer is NOT
 * equivalent — for narrow sigma (e.g. the hint context in mode-256/512)
 * the rounded distribution overshoots p(±1) by an order of magnitude,
 * which would blow the hint reservation in this test even though the
 * real signer never gets near it. */
static int32_t sample_discrete_gaussian(double sigma, int M_voc) {
  /* Cache CDFs keyed by (sigma, M_voc). 4 slots is enough for the
   * mode's z-hi + hint pair. */
  enum { CACHE_SLOTS = 4 };
  static struct {
    double sigma;
    int    M_voc;
    int    valid;
    double cdf[256];   /* 2*M_voc+1 entries; M_voc <= 50 across all modes. */
  } cache[CACHE_SLOTS];
  static int next_slot;

  int slot = -1;
  for (int s = 0; s < CACHE_SLOTS; ++s)
    if (cache[s].valid && cache[s].sigma == sigma && cache[s].M_voc == M_voc) {
      slot = s; break;
    }
  if (slot < 0) {
    slot = next_slot;
    next_slot = (next_slot + 1) % CACHE_SLOTS;
    cache[slot].sigma = sigma;
    cache[slot].M_voc = M_voc;
    double Z = 0.0;
    for (int k = -M_voc; k <= M_voc; ++k)
      Z += exp(-(double)k * k / (2.0 * sigma * sigma));
    double running = 0.0;
    for (int k = -M_voc; k <= M_voc; ++k) {
      running += exp(-(double)k * k / (2.0 * sigma * sigma)) / Z;
      cache[slot].cdf[k + M_voc] = running;
    }
    cache[slot].valid = 1;
  }

  double u = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
  int len = 2 * M_voc + 1;
  for (int i = 0; i < len; ++i)
    if (u <= cache[slot].cdf[i]) return (int32_t)(i - M_voc);
  return (int32_t)M_voc;
}

static inline int32_t sample_gaussian_zhi(int M_voc) {
  return sample_discrete_gaussian(
      (double)SHUTTLE_SIGMA / (double)SHUTTLE_ALPHA_R, M_voc);
}

static inline int32_t sample_gaussian_hint(int M_voc) {
  return sample_discrete_gaussian(
      2.0 * (double)SHUTTLE_SIGMA / (double)SHUTTLE_ALPHA_H, M_voc);
}

/* Fill z_1[0..L] with coefficients whose split (z0 at alpha_0', z1 at
 * alpha_r) yields HighBits drawn from the theoretical z-hi PMF. */
static void fill_z1(poly *z_1) {
  for (unsigned j = 0; j < SHUTTLE_N; ++j) {
    int32_t hi = sample_gaussian_zhi(ZHI_M_VOC);
    int32_t lo = (rand() % SHUTTLE_ALPHA_0P) - SHUTTLE_HALF_ALPHA_0P;
    z_1[0].coeffs[j] = (hi << SHUTTLE_ALPHA_0P_BITS) + lo;
  }
  for (unsigned i = 1; i <= SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j) {
      int32_t hi = sample_gaussian_zhi(ZHI_M_VOC);
      int32_t lo = (rand() % SHUTTLE_ALPHA_R) - SHUTTLE_HALF_ALPHA_R;
      z_1[i].coeffs[j] = (hi << SHUTTLE_ALPHA_R_BITS) + lo;
    }
}

/* ------------------------------------------------------------ */
static int test_minimal(void) {
  printf("Test 1: minimal hand-crafted round-trip\n");

  uint8_t sig[SHUTTLE_BYTES];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
  int8_t  irs[SHUTTLE_TAU];
  poly    z_1[1 + SHUTTLE_L];
  polyveck h, h_rec;

  for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)(i * 17 + 3);
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = (i & 1) ? (int8_t)1 : (int8_t)-1;
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  memset(&h, 0, sizeof h);

  int rc = pack_sig(sig, c_tilde, irs, z_1, &h);
  if (rc != 0) FAIL("pack_sig rc=%d", rc);

  uint8_t c_rec[SHUTTLE_CTILDEBYTES];
  int8_t  irs_rec[SHUTTLE_TAU];
  poly    z_rec[1 + SHUTTLE_L];
  rc = unpack_sig(c_rec, irs_rec, z_rec, &h_rec, sig);
  if (rc != 0) FAIL("unpack_sig rc=%d", rc);

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

  printf("    SHUTTLE_BYTES = %d\n", SHUTTLE_BYTES);
  PASS("minimal round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_random_rounds(void) {
  printf("Test 2: 10 random rounds\n");

  size_t min_zhi  = (size_t)-1, max_zhi  = 0, sum_zhi  = 0;
  size_t min_hint = (size_t)-1, max_hint = 0, sum_hint = 0;

  /* On-the-wire offsets, mirroring packing.c::OFF_* constants. */
  const size_t off_zhi_len  = SHUTTLE_CTILDEBYTES + SHUTTLE_IRS_SIGNBYTES;
  const size_t off_hint_len = off_zhi_len + 2 + SHUTTLE_ZHI_RESERVED_BYTES
                            + SHUTTLE_POLYZ0_LO_PACKEDBYTES
                            + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES;

  for (unsigned round = 0; round < 10; ++round) {
    uint8_t sig[SHUTTLE_BYTES];
    uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
    int8_t  irs[SHUTTLE_TAU];
    poly    z_1[1 + SHUTTLE_L];
    polyveck h, h_rec;

    for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)rand();
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
      irs[i] = (rand() & 1) ? (int8_t)1 : (int8_t)-1;

    fill_z1(z_1);
    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        h.vec[i].coeffs[j] = sample_gaussian_hint(HINT_M_VOC);
    (void)sample_from_table;  /* keep available for future tests */

    int rc = pack_sig(sig, c_tilde, irs, z_1, &h);
    if (rc != 0) FAIL("round %u: pack rc=%d", round, rc);

    uint8_t c_rec[SHUTTLE_CTILDEBYTES];
    int8_t  irs_rec[SHUTTLE_TAU];
    poly    z_rec[1 + SHUTTLE_L];
    rc = unpack_sig(c_rec, irs_rec, z_rec, &h_rec, sig);
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

    size_t zhi_len  = (size_t)sig[off_zhi_len]  | ((size_t)sig[off_zhi_len  + 1] << 8);
    size_t hint_len = (size_t)sig[off_hint_len] | ((size_t)sig[off_hint_len + 1] << 8);

    if (zhi_len  < min_zhi)  min_zhi  = zhi_len;
    if (zhi_len  > max_zhi)  max_zhi  = zhi_len;
    sum_zhi  += zhi_len;
    if (hint_len < min_hint) min_hint = hint_len;
    if (hint_len > max_hint) max_hint = hint_len;
    sum_hint += hint_len;
  }

  printf("    z-hi  rANS len: min=%zu max=%zu avg=%zu (reserved %d)\n",
         min_zhi, max_zhi, sum_zhi / 10, SHUTTLE_ZHI_RESERVED_BYTES);
  printf("    hint  rANS len: min=%zu max=%zu avg=%zu (reserved %d)\n",
         min_hint, max_hint, sum_hint / 10, SHUTTLE_HINT_RESERVED_BYTES);
  PASS("10 random rounds OK");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xC0DECAFEu ^ (unsigned)SHUTTLE_MODE);

  printf("=== test_sig (MODE=%d, SHUTTLE_BYTES=%d,\n"
         "                 ZHI_RESERVED=%d, HINT_RESERVED=%d) ===\n",
         SHUTTLE_MODE, SHUTTLE_BYTES,
         SHUTTLE_ZHI_RESERVED_BYTES, SHUTTLE_HINT_RESERVED_BYTES);

  ret |= test_minimal();
  ret |= test_random_rounds();

  if (ret == 0) printf("\n=== All sig tests PASSED ===\n");
  else          printf("\n=== Some sig tests FAILED ===\n");
  return ret;
}
