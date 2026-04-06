/*
 * packing.c - Serialization of keys and signatures for NGCC_SIGN.
 *
 * Public key format (928 bytes):
 *   rho (32B) || polypk_pack(b[0]) (448B) || polypk_pack(b[1]) (448B)
 *
 * Secret key format (448 bytes):
 *   rho (32B) || tr (64B) || key (32B) ||
 *   polyeta_pack(s[0..2]) (3*64B) || polyeta_pack(e[0..1]) (2*64B)
 *
 * Signature format (2724 bytes):
 *   c_tilde (32B) || irs_signs (4B) || polyz_pack(z[0..5]) (6*448B)
 *
 * IRS sign bits: TAU=30 bits packed into ceil(30/8)=4 bytes.
 * Bit i is 1 if the effective coefficient for monomial i is positive,
 * 0 if negative. Combined with the challenge c (from c_tilde), these
 * bits define c_eff = sum_i irs_sign_i * x^{j_i}.
 */

#include <string.h>
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"

/*************************************************
* Name:        pack_pk
**************************************************/
void pack_pk(uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES],
             const uint8_t rho[NGCC_SIGN_SEEDBYTES],
             const polyveck *b)
{
  unsigned int i;

  memcpy(pk, rho, NGCC_SIGN_SEEDBYTES);
  pk += NGCC_SIGN_SEEDBYTES;

  for(i = 0; i < NGCC_SIGN_M; ++i) {
    polypk_pack(pk, &b->vec[i]);
    pk += NGCC_SIGN_POLYPK_PACKEDBYTES;
  }
}

/*************************************************
* Name:        unpack_pk
**************************************************/
void unpack_pk(uint8_t rho[NGCC_SIGN_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES])
{
  unsigned int i;

  memcpy(rho, pk, NGCC_SIGN_SEEDBYTES);
  pk += NGCC_SIGN_SEEDBYTES;

  for(i = 0; i < NGCC_SIGN_M; ++i) {
    polypk_unpack(&b->vec[i], pk);
    pk += NGCC_SIGN_POLYPK_PACKEDBYTES;
  }
}

/*************************************************
* Name:        pack_sk
**************************************************/
void pack_sk(uint8_t sk[NGCC_SIGN_SECRETKEYBYTES],
             const uint8_t rho[NGCC_SIGN_SEEDBYTES],
             const uint8_t tr[NGCC_SIGN_TRBYTES],
             const uint8_t key[NGCC_SIGN_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e)
{
  unsigned int i;

  memcpy(sk, rho, NGCC_SIGN_SEEDBYTES);
  sk += NGCC_SIGN_SEEDBYTES;

  memcpy(sk, tr, NGCC_SIGN_TRBYTES);
  sk += NGCC_SIGN_TRBYTES;

  memcpy(sk, key, NGCC_SIGN_SEEDBYTES);
  sk += NGCC_SIGN_SEEDBYTES;

  for(i = 0; i < NGCC_SIGN_L; ++i) {
    polyeta_pack(sk, &s->vec[i]);
    sk += NGCC_SIGN_POLYETA_PACKEDBYTES;
  }

  for(i = 0; i < NGCC_SIGN_M; ++i) {
    polyeta_pack(sk, &e->vec[i]);
    sk += NGCC_SIGN_POLYETA_PACKEDBYTES;
  }
}

/*************************************************
* Name:        unpack_sk
**************************************************/
void unpack_sk(uint8_t rho[NGCC_SIGN_SEEDBYTES],
               uint8_t tr[NGCC_SIGN_TRBYTES],
               uint8_t key[NGCC_SIGN_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[NGCC_SIGN_SECRETKEYBYTES])
{
  unsigned int i;

  memcpy(rho, sk, NGCC_SIGN_SEEDBYTES);
  sk += NGCC_SIGN_SEEDBYTES;

  memcpy(tr, sk, NGCC_SIGN_TRBYTES);
  sk += NGCC_SIGN_TRBYTES;

  memcpy(key, sk, NGCC_SIGN_SEEDBYTES);
  sk += NGCC_SIGN_SEEDBYTES;

  for(i = 0; i < NGCC_SIGN_L; ++i) {
    polyeta_unpack(&s->vec[i], sk);
    sk += NGCC_SIGN_POLYETA_PACKEDBYTES;
  }

  for(i = 0; i < NGCC_SIGN_M; ++i) {
    polyeta_unpack(&e->vec[i], sk);
    sk += NGCC_SIGN_POLYETA_PACKEDBYTES;
  }
}

/*************************************************
* Name:        pack_sig
*
* Description: Serialize signature as byte array.
*              sig = c_tilde || irs_signs || polyz_pack(z[0..5])
*
*              IRS signs: TAU bits packed into IRS_SIGNBYTES bytes.
*              Bit i is 1 if irs_signs[i] == +1, 0 if irs_signs[i] == -1.
*
* Arguments:   - uint8_t sig[]: output byte array (BYTES)
*              - const uint8_t c_tilde[]: challenge hash (CTILDEBYTES)
*              - const int8_t irs_signs[]: per-monomial sign choices from IRS
*              - const polyvec *z: full response vector (VECLEN polys)
**************************************************/
void pack_sig(uint8_t sig[NGCC_SIGN_BYTES],
              const uint8_t c_tilde[NGCC_SIGN_CTILDEBYTES],
              const int8_t irs_signs[NGCC_SIGN_TAU],
              const polyvec *z)
{
  unsigned int i;

  memcpy(sig, c_tilde, NGCC_SIGN_CTILDEBYTES);
  sig += NGCC_SIGN_CTILDEBYTES;

  /* Pack IRS sign bits: bit i is 1 if irs_signs[i] > 0, 0 if < 0 */
  memset(sig, 0, NGCC_SIGN_IRS_SIGNBYTES);
  for(i = 0; i < NGCC_SIGN_TAU; ++i) {
    if(irs_signs[i] > 0)
      sig[i / 8] |= (uint8_t)(1u << (i % 8));
  }
  sig += NGCC_SIGN_IRS_SIGNBYTES;

  for(i = 0; i < NGCC_SIGN_VECLEN; ++i) {
    polyz_pack(sig, &z->vec[i]);
    sig += NGCC_SIGN_POLYZ_PACKEDBYTES;
  }
}

/*************************************************
* Name:        unpack_sig
*
* Description: Deserialize signature from byte array.
*
* Arguments:   - uint8_t c_tilde[]: output challenge hash (CTILDEBYTES)
*              - int8_t irs_signs[]: output per-monomial sign choices
*              - polyvec *z: output full response vector (VECLEN polys)
*              - const uint8_t sig[]: input byte array (BYTES)
*
* Returns 0 on success.
**************************************************/
int unpack_sig(uint8_t c_tilde[NGCC_SIGN_CTILDEBYTES],
               int8_t irs_signs[NGCC_SIGN_TAU],
               polyvec *z,
               const uint8_t sig[NGCC_SIGN_BYTES])
{
  unsigned int i;

  memcpy(c_tilde, sig, NGCC_SIGN_CTILDEBYTES);
  sig += NGCC_SIGN_CTILDEBYTES;

  /* Unpack IRS sign bits */
  for(i = 0; i < NGCC_SIGN_TAU; ++i) {
    if((sig[i / 8] >> (i % 8)) & 1)
      irs_signs[i] = 1;
    else
      irs_signs[i] = -1;
  }
  sig += NGCC_SIGN_IRS_SIGNBYTES;

  for(i = 0; i < NGCC_SIGN_VECLEN; ++i) {
    polyz_unpack(&z->vec[i], sig);
    sig += NGCC_SIGN_POLYZ_PACKEDBYTES;
  }

  return 0;
}
