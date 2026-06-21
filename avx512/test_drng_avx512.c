/*
 * Correctness test for the 16-way AVX512 SM3 Hash-DRBG vs the scalar
 * reference (SHUTTLE/ref/drng.c), driven per-lane through instantiate + a
 * sequence of generate calls; both output and internal state must match
 * byte-for-byte.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "drng.h"        /* reference */
#include "drng_avx512.h" /* 16-way */

#define NWAY DRNG_WAY_AVX512

static uint64_t rng_state = 0xfeedface12345678ULL;
static unsigned char rb(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned char)(rng_state >> 24);
}

static int state_matches(const DRNG_ctx *ref, const DRNG_ctx_avx512 *avx,
                         int k)
{
    return memcmp(ref->V, avx->V[k], SEEDLEN) == 0 &&
           memcmp(ref->C, avx->C[k], SEEDLEN) == 0 &&
           memcmp(ref->reseed_counter, avx->reseed_counter[k], SEEDLEN) ==
               0;
}

static int run_case(size_t seed_len, const unsigned long long *gens,
                    int ngen)
{
    unsigned char *seed[NWAY];
    const unsigned char *cseed[NWAY];
    DRNG_ctx ref[NWAY];
    DRNG_ctx_avx512 avx;
    int ok = 1;

    for (int k = 0; k < NWAY; k++) {
        seed[k] = malloc(seed_len ? seed_len : 1);
        for (size_t i = 0; i < seed_len; i++)
            seed[k][i] = rb();
        cseed[k] = seed[k];
    }

    for (int k = 0; k < NWAY; k++)
        init_random_number(&ref[k], seed[k], seed_len);
    init_random_number_avx512(&avx, cseed, seed_len);

    for (int k = 0; k < NWAY && ok; k++)
        if (!state_matches(&ref[k], &avx, k)) {
            printf("  [FAIL] init state lane=%d seed_len=%zu\n", k,
                   seed_len);
            ok = 0;
        }

    for (int g = 0; g < ngen && ok; g++) {
        unsigned long long bits = gens[g];
        size_t obytes = (size_t)((bits + 7) / 8);
        unsigned char *oref[NWAY], *oavx[NWAY];
        unsigned char *coavx[NWAY];
        for (int k = 0; k < NWAY; k++) {
            oref[k] = malloc(obytes ? obytes : 1);
            oavx[k] = malloc(obytes ? obytes : 1);
            memset(oref[k], 0xAA, obytes ? obytes : 1);
            memset(oavx[k], 0x55, obytes ? obytes : 1);
            coavx[k] = oavx[k];
        }
        for (int k = 0; k < NWAY; k++)
            get_random_number(&ref[k], oref[k], bits);
        get_random_number_avx512(&avx, coavx, bits);

        for (int k = 0; k < NWAY; k++) {
            if (obytes && memcmp(oref[k], oavx[k], obytes) != 0) {
                printf(
                    "  [FAIL] gen output lane=%d seed_len=%zu bits=%llu\n",
                    k, seed_len, bits);
                ok = 0;
            }
            if (!state_matches(&ref[k], &avx, k)) {
                printf(
                    "  [FAIL] gen state lane=%d seed_len=%zu bits=%llu\n",
                    k, seed_len, bits);
                ok = 0;
            }
        }
        for (int k = 0; k < NWAY; k++) {
            free(oref[k]);
            free(oavx[k]);
        }
    }

    for (int k = 0; k < NWAY; k++)
        free(seed[k]);
    return ok;
}

int main(void)
{
    int pass = 0, total = 0;
    size_t seed_lens[] = {0, 1, 16, 32, 48, 55, 56, 64, 100, 256};
    unsigned long long gens[] = {1,   8,    255,  256,  257, 440,
                                 512, 1000, 2047, 2048, 4096};
    int nseed = sizeof(seed_lens) / sizeof(seed_lens[0]);
    int ngen = sizeof(gens) / sizeof(gens[0]);

    printf("== correctness: AVX512 16-way DRNG vs reference ==\n");
    for (int s = 0; s < nseed; s++) {
        total++;
        pass += run_case(seed_lens[s], gens, ngen);
    }
    for (int s = 0; s < nseed; s++)
        for (int g = 0; g < ngen; g++) {
            total++;
            pass += run_case(seed_lens[s], &gens[g], 1);
        }

    printf("DRNG correctness: %d/%d scenarios passed\n", pass, total);
    if (pass != total) {
        printf("\nRESULT: FAIL\n");
        return 1;
    }
    printf("\nRESULT: PASS\n");
    return 0;
}
