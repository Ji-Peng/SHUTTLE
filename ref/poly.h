#ifndef SHUTTLE_POLY_H
#define SHUTTLE_POLY_H

#include <stdint.h>
#include "params.h"

typedef struct {
  int32_t coeffs[SHUTTLE_N];
} poly;

#define poly_reduce SHUTTLE_NAMESPACE(poly_reduce)
void poly_reduce(poly *a);
#define poly_caddq SHUTTLE_NAMESPACE(poly_caddq)
void poly_caddq(poly *a);

#define poly_add SHUTTLE_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub SHUTTLE_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

#define poly_ntt SHUTTLE_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont SHUTTLE_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);
#define poly_pointwise_montgomery SHUTTLE_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_chknorm SHUTTLE_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);

#define poly_sq_norm SHUTTLE_NAMESPACE(poly_sq_norm)
int64_t poly_sq_norm(const poly *a);

#define poly_uniform SHUTTLE_NAMESPACE(poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SHUTTLE_SEEDBYTES],
                  uint16_t nonce);

#define poly_uniform_eta SHUTTLE_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[SHUTTLE_CRHBYTES],
                      uint16_t nonce);

#define poly_challenge SHUTTLE_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[SHUTTLE_CTILDEBYTES]);

#define poly_decompose SHUTTLE_NAMESPACE(poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);
#define poly_make_hint SHUTTLE_NAMESPACE(poly_make_hint)
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
#define poly_use_hint SHUTTLE_NAMESPACE(poly_use_hint)
void poly_use_hint(poly *b, const poly *a, const poly *h);

#define polyeta_pack SHUTTLE_NAMESPACE(polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_unpack SHUTTLE_NAMESPACE(polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);

#define polypk_pack SHUTTLE_NAMESPACE(polypk_pack)
void polypk_pack(uint8_t *r, const poly *a);
#define polypk_unpack SHUTTLE_NAMESPACE(polypk_unpack)
void polypk_unpack(poly *r, const uint8_t *a);

#define polyz_pack SHUTTLE_NAMESPACE(polyz_pack)
void polyz_pack(uint8_t *r, const poly *a);
#define polyz_unpack SHUTTLE_NAMESPACE(polyz_unpack)
void polyz_unpack(poly *r, const uint8_t *a);

#define polyz0_pack SHUTTLE_NAMESPACE(polyz0_pack)
void polyz0_pack(uint8_t *r, const poly *a);
#define polyz0_unpack SHUTTLE_NAMESPACE(polyz0_unpack)
void polyz0_unpack(poly *r, const uint8_t *a);

#define polyw1_pack SHUTTLE_NAMESPACE(polyw1_pack)
void polyw1_pack(uint8_t *r, const poly *a);

#endif
