/*
 * packing.h - Serialization of keys and signatures for SHUTTLE.
 */

#ifndef SHUTTLE_PACKING_H
#define SHUTTLE_PACKING_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define pack_pk SHUTTLE_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[SHUTTLE_PUBLICKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const polyveck *b);

#define unpack_pk SHUTTLE_NAMESPACE(unpack_pk)
void unpack_pk(uint8_t rho[SHUTTLE_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[SHUTTLE_PUBLICKEYBYTES]);

#define pack_sk SHUTTLE_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[SHUTTLE_SECRETKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const uint8_t tr[SHUTTLE_TRBYTES],
             const uint8_t key[SHUTTLE_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e);

#define unpack_sk SHUTTLE_NAMESPACE(unpack_sk)
void unpack_sk(uint8_t rho[SHUTTLE_SEEDBYTES],
               uint8_t tr[SHUTTLE_TRBYTES],
               uint8_t key[SHUTTLE_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[SHUTTLE_SECRETKEYBYTES]);

#define pack_sig SHUTTLE_NAMESPACE(pack_sig)
void pack_sig(uint8_t sig[SHUTTLE_BYTES],
              const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
              const int8_t irs_signs[SHUTTLE_TAU],
              const polyvec *z);

#define unpack_sig SHUTTLE_NAMESPACE(unpack_sig)
int unpack_sig(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
               int8_t irs_signs[SHUTTLE_TAU],
               polyvec *z,
               const uint8_t sig[SHUTTLE_BYTES]);

#endif /* SHUTTLE_PACKING_H */
