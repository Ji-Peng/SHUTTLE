/*
 * sampler.h - Discrete Gaussian sampler for SHUTTLE.
 *
 * Algorithm (BLISS-style convolution), parameterized by SHUTTLE_SIGMA:
 *   sigma = small_sigma * k, with k = 64 (SHUTTLE_K_BITS = 6).
 *   1) CDT-sample x ~ D_{small_sigma} using 93-bit RCDT (RCDT_ENTRIES limbs,
 *      3x31-bit each). RCDT_ENTRIES and the table values depend on
 *      SHUTTLE_SIGMA (see params.h):
 *        sigma=128 -> 22 entries, small_sigma = 2.0000
 *        sigma=101 -> 17 entries, small_sigma = 1.5781   (SHUTTLE-128)
 *        sigma=149 -> 26 entries, small_sigma = 2.3281   (SHUTTLE-256)
 *   2) Sample y uniformly from {0, ..., 2^SHUTTLE_Y_BITS - 1}.
 *      y is generated in batch via sampler_y: 16 six-bit values from
 *      exactly Y_RAND_BYTES = 12 stream bytes (= 16*6/8), zero waste.
 *   3) Candidate r = (x << SHUTTLE_K_BITS) | y.
 *   4) Accept with probability exp(-y*(y + 2*k*x) / (2*sigma^2)).
 *   5) Assign random sign.
 *
 * Output: signed int16_t (max |r| = 11*sigma; fits in int16_t for all modes).
 *
 * Unified buffer design:
 *   One SHAKE-256 stream per sample_gauss_N call produces ALL randomness.
 *   Byte order per stream:
 *       signs (len/8 bytes, one-off)
 *    then repeating mini-batches of:
 *       CDT randomness      = SIGMA2_RAND_BYTES (= 192)
 *       y  randomness       = Y_RAND_BYTES      (= 12)
 *       per-attempt tails   = GAUSS_BATCH * GAUSS_RAND_BYTES (= 16 * 8 = 128)
 *   Total per mini-batch    = MINIBATCH_RAND_BYTES (= 332)
 *
 *   This schedule is identical between ref and AVX2, so both implementations
 *   are KAT-compatible bit-for-bit.
 *
 * Renyi divergence from ideal (1025-th order):
 *   sigma=128 : 1 + 2^{-93.85}
 *   sigma=101 : 1 + 2^{-93.40}
 *   sigma=149 : 1 + 2^{-93.75}
 */

#ifndef SHUTTLE_SAMPLER_H
#define SHUTTLE_SAMPLER_H

#include <stdint.h>
#include <stddef.h>
#include "symmetric.h"
#include "params.h"

/* RCDT_ENTRIES is defined in params.h (mode-dependent). */

/* ============================================================
 * Mini-batch size (number of attempts per batch).
 * Must be a multiple of 8 for AVX2 register alignment.
 * ============================================================ */
#ifndef GAUSS_BATCH
#define GAUSS_BATCH 16
#endif

/* CDT random bytes per mini-batch: 12 bytes/sample in AVX2 3x31 layout.
 *   Layout: groups of 8 samples, each group = 3 x 32 bytes (low, mid, high).
 *   BATCH=16: 2 groups x 96 bytes = 192 bytes. */
#define SIGMA2_RAND_BYTES   (GAUSS_BATCH * 12)

/* y random bytes per mini-batch: 16 six-bit fields packed into 96 bits = 12
 * bytes, consumed by sampler_y. No wasted bits (vs. the previous 1 byte per y
 * which discarded the top 2 bits). */
#define Y_RAND_BYTES        ((GAUSS_BATCH * SHUTTLE_Y_BITS + 7) / 8)

/* Per-attempt random bytes AFTER y is pulled out:
 *   8 bytes carrying 1 sign-for-zero bit + 63 rejection bits. */
#define GAUSS_RAND_BYTES    8

/* Total random bytes per mini-batch: CDT + y + attempt tails. */
#define MINIBATCH_RAND_BYTES \
    (SIGMA2_RAND_BYTES + Y_RAND_BYTES + GAUSS_BATCH * GAUSS_RAND_BYTES)

/* SHAKE-256 blocks per refill (enough for one mini-batch). */
#define SQUEEZE_NBLOCKS \
    ((MINIBATCH_RAND_BYTES + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES)

/* Internal buffer size: residual + one squeeze. */
#define GAUSS_BUF_SIZE \
    (MINIBATCH_RAND_BYTES + SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES)

/* ============================================================
 * Portable little-endian load helpers
 * ============================================================ */
static inline uint64_t load_le64(const uint8_t *p) {
    return (uint64_t)p[0]
         | ((uint64_t)p[1] << 8)
         | ((uint64_t)p[2] << 16)
         | ((uint64_t)p[3] << 24)
         | ((uint64_t)p[4] << 32)
         | ((uint64_t)p[5] << 40)
         | ((uint64_t)p[6] << 48)
         | ((uint64_t)p[7] << 56);
}

static inline uint32_t load_le32(const uint8_t *p) {
    return (uint32_t)p[0]
         | ((uint32_t)p[1] << 8)
         | ((uint32_t)p[2] << 16)
         | ((uint32_t)p[3] << 24);
}

/* Branchless unsigned less-than: returns 1 if a < b */
static inline uint32_t ct_lt_u32(uint32_t a, uint32_t b) {
    return (a - b) >> 31;
}

/* ============================================================
 * Batched small-sigma CDT sampler: sampler_sigma2 (pure function).
 *
 * Despite the historical name, the CDT targets the mode-specific
 * small_sigma = sigma / k (not literally sigma_2 = 2). The table and
 * bound RCDT_ENTRIES are selected by SHUTTLE_SIGMA in params.h.
 *
 * Input:  rand - SIGMA2_RAND_BYTES random bytes in AVX2-friendly layout.
 * Output: z_out - GAUSS_BATCH int16_t values in [0, RCDT_ENTRIES].
 *         Internal computation in 32-bit for precision.
 * Returns: GAUSS_BATCH (always).
 * ============================================================ */
int sampler_sigma2(int16_t *z_out, const uint8_t *rand);

/* ============================================================
 * Batch Gaussian sampling from one unified SHAKE-256 stream.
 *
 * Output: r[len] - signed Gaussian samples as int16_t.
 *         Signs applied directly. |r[i]| <= 11*sigma = 1408.
 *
 * Stream layout: signs (len/8 bytes) | [CDT | y+rej] x N mini-batches
 * One nonce per call. KAT-compatible between ref and AVX2.
 * ============================================================ */
void sample_gauss_N(int16_t *r,
                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                    uint64_t nonce, size_t len);

/* ============================================================
 * 4-way parallel Gaussian sampling.
 *
 * Each lane runs an independent sample_gauss_N with its own nonce.
 * Ref: wrapper calling sample_gauss_N 4 times.
 * AVX2: exploits keccak4x for 4-way parallel SHAKE-256 squeeze.
 *
 * Unused lanes: set len=0.
 * ============================================================ */
void sample_gauss_N_4x(int16_t *r0, int16_t *r1,
                        int16_t *r2, int16_t *r3,
                        const uint8_t seed[SHUTTLE_SEEDBYTES],
                        uint64_t nonce0, uint64_t nonce1,
                        uint64_t nonce2, uint64_t nonce3,
                        size_t len0, size_t len1,
                        size_t len2, size_t len3);

#endif /* SHUTTLE_SAMPLER_H */
