/*
 * api.h - NIST-style API header for NGCC_SIGN signature scheme.
 */

#ifndef NGCC_SIGN_API_H
#define NGCC_SIGN_API_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"

#define CRYPTO_PUBLICKEYBYTES NGCC_SIGN_PUBLICKEYBYTES
#define CRYPTO_SECRETKEYBYTES NGCC_SIGN_SECRETKEYBYTES
#define CRYPTO_BYTES          NGCC_SIGN_BYTES
#define CRYPTO_ALGNAME        "NGCC_SIGN"

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

#endif /* NGCC_SIGN_API_H */
