#ifndef SHUTTLE_POLYVEC_H
#define SHUTTLE_POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

/**************************************************************/
/*** Polynomial vector types                                ***/
/**************************************************************/

/* Vector of L polynomials (secret s-vector) */
typedef struct {
  poly vec[SHUTTLE_L];
} polyvecl;

/* Vector of M polynomials (secret e-vector, public b) */
typedef struct {
  poly vec[SHUTTLE_M];
} polyveck;

/* Vector of VECLEN=1+L+M polynomials (full secret key [1,s,e]) */
typedef struct {
  poly vec[SHUTTLE_VECLEN];
} polyvec;

/**************************************************************/
/*** polyvecl operations (length L = 3)                     ***/
/**************************************************************/

#define polyvecl_uniform_eta SHUTTLE_NAMESPACE(polyvecl_uniform_eta)
void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[SHUTTLE_CRHBYTES],
                          uint16_t nonce);

#define polyvecl_reduce SHUTTLE_NAMESPACE(polyvecl_reduce)
void polyvecl_reduce(polyvecl *v);

#define polyvecl_caddq SHUTTLE_NAMESPACE(polyvecl_caddq)
void polyvecl_caddq(polyvecl *v);

#define polyvecl_add SHUTTLE_NAMESPACE(polyvecl_add)
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v);

#define polyvecl_sub SHUTTLE_NAMESPACE(polyvecl_sub)
void polyvecl_sub(polyvecl *w, const polyvecl *u, const polyvecl *v);

#define polyvecl_ntt SHUTTLE_NAMESPACE(polyvecl_ntt)
void polyvecl_ntt(polyvecl *v);

#define polyvecl_invntt_tomont SHUTTLE_NAMESPACE(polyvecl_invntt_tomont)
void polyvecl_invntt_tomont(polyvecl *v);

#define polyvecl_pointwise_poly_montgomery \
        SHUTTLE_NAMESPACE(polyvecl_pointwise_poly_montgomery)
void polyvecl_pointwise_poly_montgomery(polyvecl *r, const poly *a,
                                        const polyvecl *v);

#define polyvecl_pointwise_acc_montgomery \
        SHUTTLE_NAMESPACE(polyvecl_pointwise_acc_montgomery)
void polyvecl_pointwise_acc_montgomery(poly *w,
                                       const polyvecl *u,
                                       const polyvecl *v);

#define polyvecl_chknorm SHUTTLE_NAMESPACE(polyvecl_chknorm)
int polyvecl_chknorm(const polyvecl *v, int32_t B);

#define polyvecl_sq_norm SHUTTLE_NAMESPACE(polyvecl_sq_norm)
int64_t polyvecl_sq_norm(const polyvecl *v);

/**************************************************************/
/*** polyveck operations (length M = 2)                     ***/
/**************************************************************/

#define polyveck_uniform_eta SHUTTLE_NAMESPACE(polyveck_uniform_eta)
void polyveck_uniform_eta(polyveck *v, const uint8_t seed[SHUTTLE_CRHBYTES],
                          uint16_t nonce);

#define polyveck_reduce SHUTTLE_NAMESPACE(polyveck_reduce)
void polyveck_reduce(polyveck *v);

#define polyveck_caddq SHUTTLE_NAMESPACE(polyveck_caddq)
void polyveck_caddq(polyveck *v);

#define polyveck_add SHUTTLE_NAMESPACE(polyveck_add)
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v);

#define polyveck_sub SHUTTLE_NAMESPACE(polyveck_sub)
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v);

#define polyveck_ntt SHUTTLE_NAMESPACE(polyveck_ntt)
void polyveck_ntt(polyveck *v);

#define polyveck_invntt_tomont SHUTTLE_NAMESPACE(polyveck_invntt_tomont)
void polyveck_invntt_tomont(polyveck *v);

#define polyveck_pointwise_poly_montgomery \
        SHUTTLE_NAMESPACE(polyveck_pointwise_poly_montgomery)
void polyveck_pointwise_poly_montgomery(polyveck *r, const poly *a,
                                        const polyveck *v);

#define polyveck_chknorm SHUTTLE_NAMESPACE(polyveck_chknorm)
int polyveck_chknorm(const polyveck *v, int32_t B);

#define polyveck_sq_norm SHUTTLE_NAMESPACE(polyveck_sq_norm)
int64_t polyveck_sq_norm(const polyveck *v);

#define polyveck_decompose SHUTTLE_NAMESPACE(polyveck_decompose)
void polyveck_decompose(polyveck *v1, polyveck *v0, const polyveck *v);

#define polyveck_make_hint SHUTTLE_NAMESPACE(polyveck_make_hint)
unsigned int polyveck_make_hint(polyveck *h,
                                const polyveck *v0,
                                const polyveck *v1);

#define polyveck_use_hint SHUTTLE_NAMESPACE(polyveck_use_hint)
void polyveck_use_hint(polyveck *w, const polyveck *v, const polyveck *h);

#define polyveck_pack_w1 SHUTTLE_NAMESPACE(polyveck_pack_w1)
void polyveck_pack_w1(uint8_t r[SHUTTLE_M * SHUTTLE_POLYW1_PACKEDBYTES],
                      const polyveck *w1);

/**************************************************************/
/*** polyvec operations (full length VECLEN = 6)            ***/
/**************************************************************/

#define polyvec_ntt SHUTTLE_NAMESPACE(polyvec_ntt)
void polyvec_ntt(polyvec *v);

#define polyvec_invntt_tomont SHUTTLE_NAMESPACE(polyvec_invntt_tomont)
void polyvec_invntt_tomont(polyvec *v);

#define polyvec_reduce SHUTTLE_NAMESPACE(polyvec_reduce)
void polyvec_reduce(polyvec *v);

#define polyvec_caddq SHUTTLE_NAMESPACE(polyvec_caddq)
void polyvec_caddq(polyvec *v);

#define polyvec_add SHUTTLE_NAMESPACE(polyvec_add)
void polyvec_add(polyvec *w, const polyvec *u, const polyvec *v);

#define polyvec_sub SHUTTLE_NAMESPACE(polyvec_sub)
void polyvec_sub(polyvec *w, const polyvec *u, const polyvec *v);

#define polyvec_chknorm SHUTTLE_NAMESPACE(polyvec_chknorm)
int polyvec_chknorm(const polyvec *v, int32_t B);

#define polyvec_sq_norm SHUTTLE_NAMESPACE(polyvec_sq_norm)
int64_t polyvec_sq_norm(const polyvec *v);

/**************************************************************/
/*** Matrix operations                                      ***/
/**************************************************************/

/*
 * SHUTTLE public matrix structure:
 * The public matrix B = [a | A | I_m] is m x (1+l+m) = 2 x 6.
 * - a_gen: polyveck (m=2 polys), the first column
 * - A_gen: polyveck[L] (3 columns of m=2 polys each)
 * - The identity I_m is implicit.
 *
 * For matrix expansion, we generate a_gen and A_gen from seedA.
 */

#define polyvec_matrix_expand SHUTTLE_NAMESPACE(polyvec_matrix_expand)
void polyvec_matrix_expand(polyveck *a_gen,
                           polyveck A_gen[SHUTTLE_L],
                           const uint8_t rho[SHUTTLE_SEEDBYTES]);

/*
 * polyvec_matrix_pointwise_montgomery:
 * Compute w = B * v_ntt where v_ntt = [v0, v1..vL, vL+1..vL+M]
 * is the full secret vector in NTT domain.
 *
 * For the i-th component of w (i in [0,M)):
 *   w[i] = a_gen[i]*v[0] + sum_j A_gen[j][i]*v[1+j] + v[1+L+i]
 *
 * All inputs must be in NTT domain.
 */
#define polyvec_matrix_pointwise_montgomery \
        SHUTTLE_NAMESPACE(polyvec_matrix_pointwise_montgomery)
void polyvec_matrix_pointwise_montgomery(polyveck *w,
                                         const polyveck *a_gen,
                                         const polyveck A_gen[SHUTTLE_L],
                                         const polyvec *v);

#endif
