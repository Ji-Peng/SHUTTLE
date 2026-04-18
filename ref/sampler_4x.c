/*
 * sampler_4x.c - 4-way parallel Gaussian sampling (ref version).
 *
 * Simple wrapper calling sample_gauss_N 4 times sequentially.
 * KAT-compatible with AVX2's keccak4x implementation because each
 * lane uses an independent SHAKE-256 stream.
 */

#include "sampler.h"

void sample_gauss_N_4x(int16_t *r0, int16_t *r1,
                        int16_t *r2, int16_t *r3,
                        const uint8_t seed[SHUTTLE_SEEDBYTES],
                        uint64_t nonce0, uint64_t nonce1,
                        uint64_t nonce2, uint64_t nonce3,
                        size_t len0, size_t len1,
                        size_t len2, size_t len3) {
    if (len0 > 0) sample_gauss_N(r0, seed, nonce0, len0);
    if (len1 > 0) sample_gauss_N(r1, seed, nonce1, len1);
    if (len2 > 0) sample_gauss_N(r2, seed, nonce2, len2);
    if (len3 > 0) sample_gauss_N(r3, seed, nonce3, len3);
}
