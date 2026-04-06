/*
 * sampler.h - Discrete Gaussian sampler for NGCC_SIGN.
 *
 * Algorithm (BLISS-style convolution):
 *   sigma = 128 = 2^7, decomposed as sigma = k * sigma_2 where k = 2^6, sigma_2 = 2.
 *   1) CDT-sample x ~ D_{sigma_2} using 93-bit RCDT table (22 entries)
 *   2) Sample y uniformly from {0, ..., 2^6 - 1}
 *   3) Candidate r = x * 2^6 + y
 *   4) Accept with probability exp(-y(y + 2*k*x) / (2*sigma^2))
 *   5) Assign random sign
 *
 * Output: signed int16_t (max |r| = 11*128 = 1408, fits in int16_t).
 *
 * Unified buffer design:
 *   One SHAKE-256 stream per sample_gauss_N call produces ALL randomness
 *   (signs, CDT, y, rejection). Processing in mini-batches of NGCC_GAUSS_BATCH
 *   attempts for low memory footprint. This ensures KAT compatibility between
 *   ref and AVX2 implementations.
 *
 * Renyi divergence from ideal: 1 + 2^{-93.85}
 */

#ifndef NGCC_SIGN_SAMPLER_H
#define NGCC_SIGN_SAMPLER_H

#include <stdint.h>
#include <stddef.h>
#include "symmetric.h"
#include "params.h"

/* ============================================================
 * RCDT table parameters
 * ============================================================ */
#define RCDT_ENTRIES 22

/* ============================================================
 * Mini-batch size (number of attempts per batch).
 * Must be a multiple of 8 for AVX2 register alignment.
 * ============================================================ */
#ifndef NGCC_GAUSS_BATCH
#define NGCC_GAUSS_BATCH 16
#endif

/* CDT random bytes per mini-batch: 12 bytes/sample in AVX2 3x31 layout.
 *   Layout: groups of 8 samples, each group = 3 x 32 bytes (low, mid, high).
 *   BATCH=16: 2 groups x 96 bytes = 192 bytes. */
#define SIGMA2_RAND_BYTES   (NGCC_GAUSS_BATCH * 12)

/* Per-attempt random bytes (excluding CDT):
 *   1 byte for y (6 bits) + 8 bytes for sign_r0+rejection = 9 bytes. */
#define GAUSS_RAND_BYTES    9

/* Total random bytes per mini-batch: CDT + attempts. */
#define MINIBATCH_RAND_BYTES (SIGMA2_RAND_BYTES + NGCC_GAUSS_BATCH * GAUSS_RAND_BYTES)

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
 * Batched CDT sampler: sampler_sigma2 (pure function)
 *
 * Input:  rand - SIGMA2_RAND_BYTES random bytes in AVX2-friendly layout.
 * Output: z_out - NGCC_GAUSS_BATCH int16_t values in [0, 22].
 *         Internal computation in 32-bit for precision.
 * Returns: NGCC_GAUSS_BATCH (always).
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
                    const uint8_t seed[NGCC_SIGN_SEEDBYTES],
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
                        const uint8_t seed[NGCC_SIGN_SEEDBYTES],
                        uint64_t nonce0, uint64_t nonce1,
                        uint64_t nonce2, uint64_t nonce3,
                        size_t len0, size_t len1,
                        size_t len2, size_t len3);

#endif /* NGCC_SIGN_SAMPLER_H */
