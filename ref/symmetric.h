/*
 * symmetric.h - Symmetric primitives for SHUTTLE.
 *
 * Wraps Dilithium's FIPS 202 (SHAKE-256) as the stream cipher.
 * Uses keccak_state (stack-allocated, no heap) from Dilithium.
 */

#ifndef SHUTTLE_SYMMETRIC_H
#define SHUTTLE_SYMMETRIC_H

#include <stdint.h>
#include <stddef.h>
#include "fips202.h"
#include "params.h"

/* Stream256: SHAKE-256 based XOF */
typedef keccak_state stream256_state;

#define STREAM256_BLOCKBYTES SHAKE256_RATE  /* 136 */

/*
 * Initialize SHAKE-256 stream from seed (32 bytes) + nonce (8 bytes LE).
 */
void shuttle_shake256_stream_init(stream256_state *state,
                                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                                    uint64_t nonce);

/* Convenience macros */
#define stream256_init(STATE, SEED, NONCE) \
    shuttle_shake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, NBLOCKS, STATE) \
    shake256_squeezeblocks(OUT, NBLOCKS, STATE)

#endif /* SHUTTLE_SYMMETRIC_H */
