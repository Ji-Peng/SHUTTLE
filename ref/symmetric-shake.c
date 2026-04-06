/*
 * symmetric-shake.c - SHAKE-256 stream initialization for NGCC_SIGN.
 */

#include "symmetric.h"
#include <string.h>

void ngcc_sign_shake256_stream_init(stream256_state *state,
                                    const uint8_t seed[NGCC_SIGN_SEEDBYTES],
                                    uint64_t nonce) {
    uint8_t buf[NGCC_SIGN_SEEDBYTES + 8];
    memcpy(buf, seed, NGCC_SIGN_SEEDBYTES);
    buf[NGCC_SIGN_SEEDBYTES + 0] = (uint8_t)(nonce);
    buf[NGCC_SIGN_SEEDBYTES + 1] = (uint8_t)(nonce >> 8);
    buf[NGCC_SIGN_SEEDBYTES + 2] = (uint8_t)(nonce >> 16);
    buf[NGCC_SIGN_SEEDBYTES + 3] = (uint8_t)(nonce >> 24);
    buf[NGCC_SIGN_SEEDBYTES + 4] = (uint8_t)(nonce >> 32);
    buf[NGCC_SIGN_SEEDBYTES + 5] = (uint8_t)(nonce >> 40);
    buf[NGCC_SIGN_SEEDBYTES + 6] = (uint8_t)(nonce >> 48);
    buf[NGCC_SIGN_SEEDBYTES + 7] = (uint8_t)(nonce >> 56);
    shake256_absorb_once(state, buf, sizeof(buf));
}
