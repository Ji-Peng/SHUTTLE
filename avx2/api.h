/*
 * api.h - NIST-style API header for SHUTTLE signature scheme.
 */

#ifndef SHUTTLE_API_H
#define SHUTTLE_API_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"

#define CRYPTO_PUBLICKEYBYTES SHUTTLE_PUBLICKEYBYTES
#define CRYPTO_SECRETKEYBYTES SHUTTLE_SECRETKEYBYTES
#define CRYPTO_BYTES          SHUTTLE_BYTES
#define CRYPTO_ALGNAME        "SHUTTLE"

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

#endif /* SHUTTLE_API_H */
