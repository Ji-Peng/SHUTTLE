/*
 * sign_v3.h - NGCC-Signature Alg 2 compressed signature, v3 (Phase 6c).
 *
 * v3 shares the exact same KeyGen and Sign/Verify math as v2; the only
 * difference is the on-wire signature layout, which additionally entropy-
 * codes HighBits(z[1..L]) via the Phase 6c z1 rANS table. Keys generated
 * by v3 are interoperable with v2 (same pk format); signatures are NOT.
 *
 * Sizes:
 *   SHUTTLE_PUBLICKEYBYTES   (unchanged, shared with v2/v1)
 *   SHUTTLE_SECRETKEYBYTES   (unchanged)
 *   SHUTTLE_BYTES_V3         (defined in packing_v3.h)
 */

#ifndef SHUTTLE_SIGN_V3_H
#define SHUTTLE_SIGN_V3_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "packing_v3.h"

#define crypto_sign_keypair_v3   SHUTTLE_NAMESPACE(keypair_v3)
int crypto_sign_keypair_v3(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature_v3 SHUTTLE_NAMESPACE(signature_v3)
int crypto_sign_signature_v3(uint8_t *sig, size_t *siglen,
                             const uint8_t *m, size_t mlen,
                             const uint8_t *sk);

#define crypto_sign_verify_v3    SHUTTLE_NAMESPACE(verify_v3)
int crypto_sign_verify_v3(const uint8_t *sig, size_t siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *pk);

#define crypto_sign_v3           SHUTTLE_NAMESPACE(sign_v3)
int crypto_sign_v3(uint8_t *sm, size_t *smlen,
                   const uint8_t *m, size_t mlen,
                   const uint8_t *sk);

#define crypto_sign_open_v3      SHUTTLE_NAMESPACE(open_v3)
int crypto_sign_open_v3(uint8_t *m, size_t *mlen,
                        const uint8_t *sm, size_t smlen,
                        const uint8_t *pk);

#endif /* SHUTTLE_SIGN_V3_H */
