/*
 * sign.h - Signature API for SHUTTLE.
 *
 * Provides key generation, signing, and verification.
 */

#ifndef SHUTTLE_SIGN_H
#define SHUTTLE_SIGN_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"

#define crypto_sign_keypair SHUTTLE_NAMESPACE(keypair)
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature SHUTTLE_NAMESPACE(signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_verify SHUTTLE_NAMESPACE(verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

#define crypto_sign SHUTTLE_NAMESPACE(sign)
int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

#define crypto_sign_open SHUTTLE_NAMESPACE(open)
int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

#endif /* SHUTTLE_SIGN_H */
