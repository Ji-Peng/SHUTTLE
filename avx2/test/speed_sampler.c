/*
 * speed_sampler.c - Performance benchmarks for NGCC_SIGN components
 * (AVX2).
 *
 * Benchmarks:
 *   1) sampler_sigma2 (AVX2 8-wide CDT)
 *   2) approx_exp (V2: degree 9, table 2)
 *   3) shake256x4_squeezeblocks (4-way parallel SHAKE-256, per-byte
 * throughput) 4) sample_gauss_N (full sampler, N=256) 5) sample_gauss_N_4x
 * (4-way parallel via keccak4x, 4 x N=256)
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "../approx_exp.h"
#include "../randombytes.h"
#include "../sampler.h"
#include "../symmetric.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 10000
#define N 256

int main(void)
{
    uint8_t seed[NGCC_SIGN_SEEDBYTES];
    int16_t r[N];
    uint64_t t[NTESTS];

    randombytes(seed, NGCC_SIGN_SEEDBYTES);

    printf("=== NGCC_SIGN Component Benchmarks - AVX2 (sigma=%d) ===\n\n",
           NGCC_SIGN_SIGMA);

    /* ---- 1. sampler_sigma2 (AVX2 batched CDT, BATCH=16) ---- */
    {
        uint8_t cdt_rand[SIGMA2_RAND_BYTES];
        int16_t z[NGCC_GAUSS_BATCH];
        randombytes(cdt_rand, sizeof(cdt_rand));

        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            sampler_sigma2(z, cdt_rand);
        }
        print_results("sampler_sigma2 (AVX2, BATCH=16):", t, NTESTS);
    }

    /* ---- 2. approx_exp (V2: degree 9, table 2) ---- */
    {
        uint64_t test_inputs[8] = {
            0,
            UINT64_C(576460752303423488),  /* ~0.5 in Q60 */
            UINT64_C(1152921504606846976), /* ~1.0 in Q60 */
            UINT64_C(2305843009213693952), /* ~2.0 in Q60 */
            UINT64_C(3458764513820540928), /* ~3.0 in Q60 */
            UINT64_C(4611686018427387904), /* ~4.0 in Q60 */
            UINT64_C(5764607523034234880), /* ~5.0 in Q60 */
            UINT64_C(6381655765338685440), /* ~5.535 in Q60 */
        };
        volatile uint64_t sink = 0;

        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            sink += approx_exp(test_inputs[i & 7]);
        }
        print_results("approx_exp (V2, single call):", t, NTESTS);
        (void)sink;
    }

    /* ---- 3. shake256x4_squeezeblocks (4-way parallel SHAKE-256) ---- */
    {
#define STREAM_NBLOCKS 128
#define STREAM_BYTES (STREAM_NBLOCKS * STREAM256_BLOCKBYTES)
#define TOTAL_BYTES_4X (4 * STREAM_BYTES)

        uint8_t buf0[STREAM_BYTES], buf1[STREAM_BYTES];
        uint8_t buf2[STREAM_BYTES], buf3[STREAM_BYTES];
        uint8_t inbuf[4][NGCC_SIGN_SEEDBYTES + 8];
        keccakx4_state state4x;

        /* Prepare 4 input buffers */
        for (int j = 0; j < 4; j++) {
            memcpy(inbuf[j], seed, NGCC_SIGN_SEEDBYTES);
            uint64_t nonce = (uint64_t)j;
            for (int k = 0; k < 8; k++)
                inbuf[j][NGCC_SIGN_SEEDBYTES + k] =
                    (uint8_t)(nonce >> (8 * k));
        }

        shake256x4_absorb_once(&state4x, inbuf[0], inbuf[1], inbuf[2],
                               inbuf[3], NGCC_SIGN_SEEDBYTES + 8);

        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            shake256x4_squeezeblocks(buf0, buf1, buf2, buf3,
                                     STREAM_NBLOCKS, &state4x);
        }
        print_results(
            "shake256x4_squeezeblocks (16 blk x 4 lanes = 8704 bytes):", t,
            NTESTS);

        /* Per-byte cost */
        uint64_t med;
        shake256x4_absorb_once(&state4x, inbuf[0], inbuf[1], inbuf[2],
                               inbuf[3], NGCC_SIGN_SEEDBYTES + 8);
        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            shake256x4_squeezeblocks(buf0, buf1, buf2, buf3,
                                     STREAM_NBLOCKS, &state4x);
        }
        med = cycles_median(t, NTESTS);
        double avg_per_byte = (double)med / (double)TOTAL_BYTES_4X;
        printf(
            "shake256x4 avg per-byte: %.2f cycles/byte (%d total "
            "bytes/squeeze)\n\n",
            avg_per_byte, TOTAL_BYTES_4X);

        /* Also benchmark single-lane stream256 for comparison */
        {
            uint8_t sbuf[STREAM_BYTES];
            stream256_state state;

            for (int i = 0; i < NTESTS; i++) {
                stream256_init(&state, seed, (uint64_t)i);
                t[i] = cpucycles();
                stream256_squeezeblocks(sbuf, STREAM_NBLOCKS, &state);
            }
            print_results(
                "stream256_squeezeblocks (16 blocks = 2176 bytes, "
                "single-lane ref):",
                t, NTESTS);

            stream256_init(&state, seed, 0);
            for (int i = 0; i < NTESTS; i++) {
                t[i] = cpucycles();
                stream256_squeezeblocks(sbuf, STREAM_NBLOCKS, &state);
            }
            med = cycles_median(t, NTESTS);
            avg_per_byte = (double)med / (double)STREAM_BYTES;
            printf(
                "stream256 (single-lane) avg per-byte: %.2f cycles/byte "
                "(%d bytes/squeeze)\n\n",
                avg_per_byte, STREAM_BYTES);
        }

#undef STREAM_NBLOCKS
#undef STREAM_BYTES
#undef TOTAL_BYTES_4X
    }

    /* ---- 4. sample_gauss_N (full sampler, N=256) ---- */
    for (int i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        sample_gauss_N(r, seed, (uint64_t)i, N);
    }
    print_results("sample_gauss_N (N=256):", t, NTESTS);

    /* ---- 5. sample_gauss_N_4x (4-way parallel via keccak4x) ---- */
    {
        int16_t r0[N], r1[N], r2[N], r3[N];
        for (int i = 0; i < NTESTS; i++) {
            t[i] = cpucycles();
            sample_gauss_N_4x(r0, r1, r2, r3, seed, (uint64_t)(4 * i),
                              (uint64_t)(4 * i + 1), (uint64_t)(4 * i + 2),
                              (uint64_t)(4 * i + 3), N, N, N, N);
        }
        print_results("sample_gauss_N_4x (4 x N=256):", t, NTESTS);
    }

    return 0;
}
