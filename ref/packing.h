/*
 * packing.h - Serialization of keys and signatures for NGCC_SIGN.
 */

#ifndef NGCC_SIGN_PACKING_H
#define NGCC_SIGN_PACKING_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define pack_pk NGCC_SIGN_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES],
             const uint8_t rho[NGCC_SIGN_SEEDBYTES],
             const polyveck *b);

#define unpack_pk NGCC_SIGN_NAMESPACE(unpack_pk)
void unpack_pk(uint8_t rho[NGCC_SIGN_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES]);

#define pack_sk NGCC_SIGN_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[NGCC_SIGN_SECRETKEYBYTES],
             const uint8_t rho[NGCC_SIGN_SEEDBYTES],
             const uint8_t tr[NGCC_SIGN_TRBYTES],
             const uint8_t key[NGCC_SIGN_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e);

#define unpack_sk NGCC_SIGN_NAMESPACE(unpack_sk)
void unpack_sk(uint8_t rho[NGCC_SIGN_SEEDBYTES],
               uint8_t tr[NGCC_SIGN_TRBYTES],
               uint8_t key[NGCC_SIGN_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[NGCC_SIGN_SECRETKEYBYTES]);

#define pack_sig NGCC_SIGN_NAMESPACE(pack_sig)
void pack_sig(uint8_t sig[NGCC_SIGN_BYTES],
              const uint8_t c_tilde[NGCC_SIGN_CTILDEBYTES],
              const int8_t irs_signs[NGCC_SIGN_TAU],
              const polyvec *z);

#define unpack_sig NGCC_SIGN_NAMESPACE(unpack_sig)
int unpack_sig(uint8_t c_tilde[NGCC_SIGN_CTILDEBYTES],
               int8_t irs_signs[NGCC_SIGN_TAU],
               polyvec *z,
               const uint8_t sig[NGCC_SIGN_BYTES]);

#endif /* NGCC_SIGN_PACKING_H */
