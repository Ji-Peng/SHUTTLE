/*
 * symmetric-shake.c - SHAKE-256 stream initialization for SHUTTLE.
 */

#include "symmetric.h"
#include <string.h>

void shuttle_shake256_stream_init(stream256_state *state,
                                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                                    uint64_t nonce) {
    uint8_t buf[SHUTTLE_SEEDBYTES + 8];
    memcpy(buf, seed, SHUTTLE_SEEDBYTES);
    buf[SHUTTLE_SEEDBYTES + 0] = (uint8_t)(nonce);
    buf[SHUTTLE_SEEDBYTES + 1] = (uint8_t)(nonce >> 8);
    buf[SHUTTLE_SEEDBYTES + 2] = (uint8_t)(nonce >> 16);
    buf[SHUTTLE_SEEDBYTES + 3] = (uint8_t)(nonce >> 24);
    buf[SHUTTLE_SEEDBYTES + 4] = (uint8_t)(nonce >> 32);
    buf[SHUTTLE_SEEDBYTES + 5] = (uint8_t)(nonce >> 40);
    buf[SHUTTLE_SEEDBYTES + 6] = (uint8_t)(nonce >> 48);
    buf[SHUTTLE_SEEDBYTES + 7] = (uint8_t)(nonce >> 56);
    shake256_absorb_once(state, buf, sizeof(buf));
}
