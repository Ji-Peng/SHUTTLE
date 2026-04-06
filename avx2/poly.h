#ifndef NGCC_SIGN_POLY_H
#define NGCC_SIGN_POLY_H

#include <stdint.h>
#include "params.h"

typedef struct {
  int32_t coeffs[NGCC_SIGN_N];
} poly;

#define poly_reduce NGCC_SIGN_NAMESPACE(poly_reduce)
void poly_reduce(poly *a);
#define poly_caddq NGCC_SIGN_NAMESPACE(poly_caddq)
void poly_caddq(poly *a);

#define poly_add NGCC_SIGN_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub NGCC_SIGN_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

#define poly_ntt NGCC_SIGN_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont NGCC_SIGN_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);
#define poly_pointwise_montgomery NGCC_SIGN_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_chknorm NGCC_SIGN_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);

#define poly_sq_norm NGCC_SIGN_NAMESPACE(poly_sq_norm)
int64_t poly_sq_norm(const poly *a);

#define poly_uniform NGCC_SIGN_NAMESPACE(poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[NGCC_SIGN_SEEDBYTES],
                  uint16_t nonce);

#define poly_uniform_eta NGCC_SIGN_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[NGCC_SIGN_CRHBYTES],
                      uint16_t nonce);

#define poly_challenge NGCC_SIGN_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[NGCC_SIGN_CTILDEBYTES]);

#define poly_decompose NGCC_SIGN_NAMESPACE(poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);
#define poly_make_hint NGCC_SIGN_NAMESPACE(poly_make_hint)
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
#define poly_use_hint NGCC_SIGN_NAMESPACE(poly_use_hint)
void poly_use_hint(poly *b, const poly *a, const poly *h);

#define polyeta_pack NGCC_SIGN_NAMESPACE(polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_unpack NGCC_SIGN_NAMESPACE(polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);

#define polypk_pack NGCC_SIGN_NAMESPACE(polypk_pack)
void polypk_pack(uint8_t *r, const poly *a);
#define polypk_unpack NGCC_SIGN_NAMESPACE(polypk_unpack)
void polypk_unpack(poly *r, const uint8_t *a);

#define polyz_pack NGCC_SIGN_NAMESPACE(polyz_pack)
void polyz_pack(uint8_t *r, const poly *a);
#define polyz_unpack NGCC_SIGN_NAMESPACE(polyz_unpack)
void polyz_unpack(poly *r, const uint8_t *a);

#define polyw1_pack NGCC_SIGN_NAMESPACE(polyw1_pack)
void polyw1_pack(uint8_t *r, const poly *a);

#endif
