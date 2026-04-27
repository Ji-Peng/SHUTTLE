/*
 * rng.c - SHAKE-256 based deterministic DRBG for test-vector generation.
 *
 * Single 48-byte seed; all subsequent randombytes() calls squeeze from
 * the same stream. The byte schedule depends only on the seed.
 */
#include <string.h>

#include "rng.h"
#include "fips202.h"

static keccak_state g_state;
static int          g_seeded = 0;

void randombytes_seed(const uint8_t seed[48])
{
    shake256_init(&g_state);
    shake256_absorb(&g_state, seed, 48);
    shake256_finalize(&g_state);
    g_seeded = 1;
}

void randombytes(uint8_t *out, size_t outlen)
{
    if (!g_seeded) {
        uint8_t zeros[48];
        memset(zeros, 0, sizeof zeros);
        randombytes_seed(zeros);
    }
    shake256_squeeze(out, outlen, &g_state);
}
