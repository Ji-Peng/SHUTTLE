/*
 * speed_sampler.c - Performance benchmarks for NGCC_SIGN components.
 *
 * Benchmarks:
 *   1) sampler_sigma2 (batched CDT)
 *   2) approx_exp (V2: degree 9, table 2)
 *   3) stream256 (SHAKE-256 per-byte throughput)
 *   4) sample_gauss_N (full sampler, N=256)
 *   5) sample_gauss_N_4x (4-way parallel, 4 x N=256)
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "../sampler.h"
#include "../approx_exp.h"
#include "../randombytes.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 10000
#define N 256

int main(void) {
    uint8_t seed[NGCC_SIGN_SEEDBYTES];
    int16_t r[N];
    uint64_t t[NTESTS];

    randombytes(seed, NGCC_SIGN_SEEDBYTES);

    printf("=== NGCC_SIGN Component Benchmarks (sigma=%d) ===\n\n",
           NGCC_SIGN_SIGMA);

    /* ---- 1. sampler_sigma2 (batched CDT, BATCH=16) ---- */
    {
        uint8_t cdt_rand[SIGMA2_RAND_BYTES];
        int16_t z[NGCC_GAUSS_BATCH];
        randombytes(cdt_rand, sizeof(cdt_rand));

        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            sampler_sigma2(z, cdt_rand);
        }
        print_results("sampler_sigma2 (BATCH=16):", t, NTESTS);
    }

    /* ---- 2. approx_exp (V2: degree 9, table 2) ---- */
    {
        /* Test with representative exponent values from the sampler range.
         * max exponent ~5.535, in Q60: 5.535 * 2^60 ~ 6.38e18 */
        uint64_t test_inputs[8] = {
            0,                                      /* exp(-0) = 1 */
            UINT64_C(576460752303423488),            /* ~0.5 in Q60 */
            UINT64_C(1152921504606846976),           /* ~1.0 in Q60 */
            UINT64_C(2305843009213693952),           /* ~2.0 in Q60 */
            UINT64_C(3458764513820540928),           /* ~3.0 in Q60 */
            UINT64_C(4611686018427387904),           /* ~4.0 in Q60 */
            UINT64_C(5764607523034234880),           /* ~5.0 in Q60 */
            UINT64_C(6381655765338685440),           /* ~5.535 in Q60 */
        };
        volatile uint64_t sink = 0;

        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            sink += approx_exp(test_inputs[i & 7]);
        }
        print_results("approx_exp (V2, single call):", t, NTESTS);
        (void)sink;
    }

    /* ---- 3. stream256 (SHAKE-256 per-byte throughput) ---- */
    {
        /* Measure cost of squeezing N_BLOCKS blocks */
        #define STREAM_NBLOCKS 16
        #define STREAM_BYTES (STREAM_NBLOCKS * STREAM256_BLOCKBYTES)
        uint8_t stream_buf[STREAM_BYTES];
        stream256_state state;

        for (int i = 0; i < NTESTS; i++) {
            stream256_init(&state, seed, (uint64_t)i);
            t[i] = cpucycles();
            stream256_squeezeblocks(stream_buf, STREAM_NBLOCKS, &state);
        }
        print_results("stream256_squeezeblocks (16 blocks = 2176 bytes):", t, NTESTS);

        /* Compute and print per-byte cost using median cycles/sample */
        uint64_t med;
        stream256_init(&state, seed, 0);
        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            stream256_squeezeblocks(stream_buf, STREAM_NBLOCKS, &state);
        }
        med = cycles_median(t, NTESTS);
        double avg_per_byte = (double)med / (double)STREAM_BYTES;
        printf("stream256 avg per-byte: %.2f cycles/byte (%d bytes/squeeze)\n\n",
               avg_per_byte, STREAM_BYTES);
        #undef STREAM_NBLOCKS
        #undef STREAM_BYTES
    }

    /* ---- 4. sample_gauss_N (full sampler, N=256) ---- */
    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        sample_gauss_N(r, seed, (uint64_t)i, N);
    }
    print_results("sample_gauss_N (N=256):", t, NTESTS);

    /* ---- 5. sample_gauss_N_4x (4-way parallel) ---- */
    {
        int16_t r0[N], r1[N], r2[N], r3[N];
        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            sample_gauss_N_4x(r0, r1, r2, r3,
                               seed,
                               (uint64_t)(4*i), (uint64_t)(4*i+1),
                               (uint64_t)(4*i+2), (uint64_t)(4*i+3),
                               N, N, N, N);
        }
        print_results("sample_gauss_N_4x (4 x N=256):", t, NTESTS);
    }

    return 0;
}
