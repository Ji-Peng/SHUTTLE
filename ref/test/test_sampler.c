/*
 * test_sampler.c - Test the NGCC_SIGN discrete Gaussian sampler.
 *
 * Tests:
 *   1) Basic statistical properties (mean ~0, variance ~sigma^2)
 *   2) Bound check: |r| <= 11*sigma = 1408
 *   3) KAT reproducibility: same seed+nonce -> same output
 *   4) 4x KAT: sample_gauss_N_4x lane 0 matches sample_gauss_N
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sampler.h"
#include "../randombytes.h"

#define N 256

int main(void) {
    uint8_t seed[NGCC_SIGN_SEEDBYTES];
    int16_t r[N];
    int ok = 1;

    /* Fixed seed for reproducibility */
    memset(seed, 0x42, NGCC_SIGN_SEEDBYTES);

    /* Test 1: Basic sampling and statistics */
    printf("=== NGCC_SIGN Gaussian Sampler Test (sigma=%d, N=%d) ===\n\n",
           NGCC_SIGN_SIGMA, N);

    sample_gauss_N(r, seed, 0, N);

    double mean = 0, var = 0;
    int max_abs = 0;
    for (int i = 0; i < N; i++) {
        mean += r[i];
        var += (double)r[i] * r[i];
        int a = abs(r[i]);
        if (a > max_abs) max_abs = a;
    }
    mean /= N;
    var = var / N - mean * mean;

    printf("Statistics:\n");
    printf("  Mean:     %.2f (expected ~0)\n", mean);
    printf("  Variance: %.2f (expected ~%.2f)\n", var, (double)NGCC_SIGN_SIGMA * NGCC_SIGN_SIGMA);
    printf("  Max|r|:   %d (bound: %d = 11*sigma)\n", max_abs, NGCC_SIGN_BOUND);

    /* Test 2: Bound check */
    printf("\nBound check: ");
    for (int i = 0; i < N; i++) {
        if (r[i] > NGCC_SIGN_BOUND || r[i] < -NGCC_SIGN_BOUND) {
            printf("FAIL (r[%d] = %d out of bounds!)\n", i, r[i]);
            ok = 0;
        }
    }
    if (ok) printf("PASS\n");

    /* Test 3: KAT reproducibility */
    int16_t r2[N];
    sample_gauss_N(r2, seed, 0, N);
    printf("KAT test:   ");
    if (memcmp(r, r2, N * sizeof(int16_t)) != 0) {
        printf("FAIL (different output for same seed+nonce)\n");
        ok = 0;
    } else {
        printf("PASS\n");
    }

    /* Test 4: 4x KAT - lane 0 must match single call */
    int16_t r4x[4][N];
    sample_gauss_N_4x(r4x[0], r4x[1], r4x[2], r4x[3],
                       seed, 0, 1, 2, 3, N, N, N, N);

    printf("4x KAT L0:  ");
    if (memcmp(r4x[0], r, N * sizeof(int16_t)) != 0) {
        printf("FAIL (4x lane 0 != single)\n");
        ok = 0;
    } else {
        printf("PASS\n");
    }

    /* Check other lanes match individual calls */
    for (int lane = 1; lane < 4; lane++) {
        int16_t rs[N];
        sample_gauss_N(rs, seed, (uint64_t)lane, N);
        printf("4x KAT L%d:  ", lane);
        if (memcmp(r4x[lane], rs, N * sizeof(int16_t)) != 0) {
            printf("FAIL\n");
            ok = 0;
        } else {
            printf("PASS\n");
        }
    }

    /* Test 5: Larger statistical test with random seed */
    printf("\nLarge statistical test (10 x N=%d):\n", N);
    randombytes(seed, NGCC_SIGN_SEEDBYTES);
    double total_mean = 0, total_var = 0;
    int total_max = 0;
    int total_samples = 10 * N;

    for (int round = 0; round < 10; round++) {
        sample_gauss_N(r, seed, (uint64_t)round, N);
        for (int i = 0; i < N; i++) {
            total_mean += r[i];
            total_var += (double)r[i] * r[i];
            int a = abs(r[i]);
            if (a > total_max) total_max = a;
            if (r[i] > NGCC_SIGN_BOUND || r[i] < -NGCC_SIGN_BOUND) {
                printf("  ERROR: round %d, r[%d] = %d out of bounds!\n",
                       round, i, r[i]);
                ok = 0;
            }
        }
    }
    total_mean /= total_samples;
    total_var = total_var / total_samples - total_mean * total_mean;

    printf("  Mean:     %.2f\n", total_mean);
    printf("  Variance: %.2f (expected ~%.2f)\n", total_var,
           (double)NGCC_SIGN_SIGMA * NGCC_SIGN_SIGMA);
    printf("  Max|r|:   %d\n", total_max);

    printf("\n=== Overall: %s ===\n", ok ? "ALL PASSED" : "SOME FAILED");
    return ok ? 0 : 1;
}
