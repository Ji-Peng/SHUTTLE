/*
 * symmetric.h - Symmetric primitives for NGCC_SIGN (AVX2).
 *
 * Single-lane: fips202.h (AVX2 Keccak permutation)
 * Four-lane:   fips202x4.h (4-way parallel SIMD256)
 */

#ifndef NGCC_SIGN_SYMMETRIC_H
#define NGCC_SIGN_SYMMETRIC_H

#include <stdint.h>
#include <stddef.h>
#include "fips202.h"
#include "fips202x4.h"
#include "params.h"

/* Stream256: SHAKE-256 based XOF (single-lane) */
typedef keccak_state stream256_state;

#define STREAM256_BLOCKBYTES SHAKE256_RATE  /* 136 */

/*
 * Initialize SHAKE-256 stream from seed (32 bytes) + nonce (8 bytes LE).
 */
void ngcc_sign_shake256_stream_init(stream256_state *state,
                                    const uint8_t seed[NGCC_SIGN_SEEDBYTES],
                                    uint64_t nonce);

/* Convenience macros */
#define stream256_init(STATE, SEED, NONCE) \
    ngcc_sign_shake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, NBLOCKS, STATE) \
    shake256_squeezeblocks(OUT, NBLOCKS, STATE)

#endif /* NGCC_SIGN_SYMMETRIC_H */
