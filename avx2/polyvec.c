#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

/**************************************************************/
/*** Matrix expansion and multiplication                    ***/
/**************************************************************/

/*************************************************
* Name:        polyvec_matrix_expand
*
* Description: Generate the public matrix components a_gen (m polys)
*              and A_gen (L columns, each m polys) from seed rho.
*              Nonce encoding: (row << 8) + col, where col 0 = a_gen,
*              cols 1..L = A_gen[0..L-1].
*
* Arguments:   - polyveck *a_gen: output first column (m polys)
*              - polyveck A_gen[L]: output columns 1..L
*              - const uint8_t rho[]: seed of length SEEDBYTES
**************************************************/
void polyvec_matrix_expand(polyveck *a_gen,
                           polyveck A_gen[SHUTTLE_L],
                           const uint8_t rho[SHUTTLE_SEEDBYTES])
{
  unsigned int i, j;

  /* Generate a_gen: column 0 */
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_uniform(&a_gen->vec[i], rho, (uint16_t)((i << 8) + 0));

  /* Generate A_gen: columns 1..L */
  for(j = 0; j < SHUTTLE_L; ++j)
    for(i = 0; i < SHUTTLE_M; ++i)
      poly_uniform(&A_gen[j].vec[i], rho, (uint16_t)((i << 8) + 1 + j));
}

/*************************************************
* Name:        polyvec_matrix_pointwise_montgomery
*
* Description: Compute w = B * v in NTT domain, where B = [a | A | I_m].
*              For each row i (0 <= i < M):
*                w[i] = a_gen[i]*v[0] + sum_{j=0}^{L-1} A_gen[j][i]*v[1+j]
*                       + v[1+L+i]
*
* Arguments:   - polyveck *w: output (m polys)
*              - const polyveck *a_gen: first column
*              - const polyveck A_gen[L]: columns 1..L
*              - const polyvec *v: full vector (1+L+M polys) in NTT domain
**************************************************/
void polyvec_matrix_pointwise_montgomery(polyveck *w,
                                         const polyveck *a_gen,
                                         const polyveck A_gen[SHUTTLE_L],
                                         const polyvec *v)
{
  unsigned int i, j;
  poly t;

  for(i = 0; i < SHUTTLE_M; ++i) {
    /* w[i] = a_gen[i] * v[0] */
    poly_pointwise_montgomery(&w->vec[i], &a_gen->vec[i], &v->vec[0]);

    /* w[i] += A_gen[j][i] * v[1+j] for j = 0..L-1 */
    for(j = 0; j < SHUTTLE_L; ++j) {
      poly_pointwise_montgomery(&t, &A_gen[j].vec[i], &v->vec[1 + j]);
      poly_add(&w->vec[i], &w->vec[i], &t);
    }

    /* w[i] += v[1+L+i] (identity block) */
    poly_add(&w->vec[i], &w->vec[i], &v->vec[1 + SHUTTLE_L + i]);
  }
}

/**************************************************************/
/*** Vectors of polynomials of length L                     ***/
/**************************************************************/

void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[SHUTTLE_CRHBYTES],
                          uint16_t nonce)
{
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_uniform_eta(&v->vec[i], seed, nonce++);
}

void polyvecl_reduce(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_reduce(&v->vec[i]);
}

void polyvecl_caddq(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_caddq(&v->vec[i]);
}

/*************************************************
* Name:        polyvecl_add
*
* Description: Add vectors of polynomials of length L.
*              No modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - const polyvecl *u: pointer to first summand
*              - const polyvecl *v: pointer to second summand
**************************************************/
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_sub
*
* Description: Subtract vectors of polynomials of length L.
*              No modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/
void polyvecl_sub(polyvecl *w, const polyvecl *u, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_ntt
*
* Description: Forward NTT of all polynomials in vector of length L.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
void polyvecl_ntt(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_ntt(&v->vec[i]);
}

void polyvecl_invntt_tomont(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

void polyvecl_pointwise_poly_montgomery(polyvecl *r, const poly *a,
                                        const polyvecl *v)
{
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_pointwise_acc_montgomery
*
* Description: Pointwise multiply vectors of polynomials of length L, multiply
*              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/
void polyvecl_pointwise_acc_montgomery(poly *w,
                                       const polyvecl *u,
                                       const polyvecl *v)
{
  unsigned int i;
  poly t;

  poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
  for(i = 1; i < SHUTTLE_L; ++i) {
    poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
    poly_add(w, w, &t);
  }
}

/*************************************************
* Name:        polyvecl_chknorm
*
* Description: Check infinity norm of polynomials in vector of length L.
*
* Arguments:   - const polyvecl *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials is strictly smaller than B
* and 1 otherwise.
**************************************************/
int polyvecl_chknorm(const polyvecl *v, int32_t bound) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_L; ++i)
    if(poly_chknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

/*************************************************
* Name:        polyvecl_sq_norm
*
* Description: Compute sum of squared L2 norms of all polynomials
*              in vector of length L.
*
* Arguments:   - const polyvecl *v: pointer to vector
*
* Returns sum of squared norms as int64_t.
**************************************************/
int64_t polyvecl_sq_norm(const polyvecl *v) {
  unsigned int i;
  int64_t norm = 0;

  for(i = 0; i < SHUTTLE_L; ++i)
    norm += poly_sq_norm(&v->vec[i]);

  return norm;
}

/**************************************************************/
/*** Vectors of polynomials of length K (= M = 2)          ***/
/**************************************************************/

void polyveck_uniform_eta(polyveck *v, const uint8_t seed[SHUTTLE_CRHBYTES],
                          uint16_t nonce)
{
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_uniform_eta(&v->vec[i], seed, nonce++);
}

/*************************************************
* Name:        polyveck_reduce
*
* Description: Reduce coefficients of polynomials in vector of length M.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_reduce(polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_reduce(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_caddq
*
* Description: For all coefficients of polynomials in vector of length M
*              add Q if coefficient is negative.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_caddq(polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_caddq(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_add
*
* Description: Add vectors of polynomials of length M.
*              No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - const polyveck *u: pointer to first summand
*              - const polyveck *v: pointer to second summand
**************************************************/
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_sub
*
* Description: Subtract vectors of polynomials of length M.
*              No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - const polyveck *u: pointer to first input vector
*              - const polyveck *v: pointer to second input vector
**************************************************/
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_ntt
*
* Description: Forward NTT of all polynomials in vector of length M.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_ntt(polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_ntt(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_invntt_tomont
*
* Description: Inverse NTT and multiplication by 2^{32} of polynomials
*              in vector of length M.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_invntt_tomont(polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

void polyveck_pointwise_poly_montgomery(polyveck *r, const poly *a,
                                        const polyveck *v)
{
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
* Name:        polyveck_chknorm
*
* Description: Check infinity norm of polynomials in vector of length M.
*
* Arguments:   - const polyveck *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials are strictly smaller than B
* and 1 otherwise.
**************************************************/
int polyveck_chknorm(const polyveck *v, int32_t bound) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    if(poly_chknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

/*************************************************
* Name:        polyveck_sq_norm
*
* Description: Compute sum of squared L2 norms of all polynomials
*              in vector of length M.
*
* Arguments:   - const polyveck *v: pointer to vector
*
* Returns sum of squared norms as int64_t.
**************************************************/
int64_t polyveck_sq_norm(const polyveck *v) {
  unsigned int i;
  int64_t norm = 0;

  for(i = 0; i < SHUTTLE_M; ++i)
    norm += poly_sq_norm(&v->vec[i]);

  return norm;
}

/*************************************************
* Name:        polyveck_decompose
*
* Description: For all coefficients a of polynomials in vector of length M,
*              compute high and low bits a1, a0 such that
*              a mod^+ Q = a1 * 2*ALPHA_H + a0.
*
* Arguments:   - polyveck *v1: pointer to output vector with high bits
*              - polyveck *v0: pointer to output vector with low bits
*              - const polyveck *v: pointer to input vector
**************************************************/
void polyveck_decompose(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_decompose(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_make_hint
*
* Description: Compute hint vector.
*
* Arguments:   - polyveck *h: pointer to output vector
*              - const polyveck *v0: pointer to low part of input vector
*              - const polyveck *v1: pointer to high part of input vector
*
* Returns number of 1 bits.
**************************************************/
unsigned int polyveck_make_hint(polyveck *h,
                                const polyveck *v0,
                                const polyveck *v1)
{
  unsigned int i, s = 0;

  for(i = 0; i < SHUTTLE_M; ++i)
    s += poly_make_hint(&h->vec[i], &v0->vec[i], &v1->vec[i]);

  return s;
}

/*************************************************
* Name:        polyveck_use_hint
*
* Description: Use hint vector to correct the high bits of input vector.
*
* Arguments:   - polyveck *w: pointer to output vector with corrected high bits
*              - const polyveck *v: pointer to input vector
*              - const polyveck *h: pointer to input hint vector
**************************************************/
void polyveck_use_hint(polyveck *w, const polyveck *v, const polyveck *h) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    poly_use_hint(&w->vec[i], &v->vec[i], &h->vec[i]);
}

void polyveck_pack_w1(uint8_t r[SHUTTLE_M * SHUTTLE_POLYW1_PACKEDBYTES],
                      const polyveck *w1)
{
  unsigned int i;

  for(i = 0; i < SHUTTLE_M; ++i)
    polyw1_pack(&r[i * SHUTTLE_POLYW1_PACKEDBYTES], &w1->vec[i]);
}

/**************************************************************/
/*** Full vector operations (length VECLEN = 1+L+M = 6)    ***/
/**************************************************************/

void polyvec_ntt(polyvec *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    poly_ntt(&v->vec[i]);
}

void polyvec_invntt_tomont(polyvec *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

void polyvec_reduce(polyvec *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    poly_reduce(&v->vec[i]);
}

void polyvec_caddq(polyvec *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    poly_caddq(&v->vec[i]);
}

void polyvec_add(polyvec *w, const polyvec *u, const polyvec *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyvec_sub(polyvec *w, const polyvec *u, const polyvec *v) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

int polyvec_chknorm(const polyvec *v, int32_t B) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    if(poly_chknorm(&v->vec[i], B))
      return 1;

  return 0;
}

int64_t polyvec_sq_norm(const polyvec *v) {
  unsigned int i;
  int64_t norm = 0;

  for(i = 0; i < SHUTTLE_VECLEN; ++i)
    norm += poly_sq_norm(&v->vec[i]);

  return norm;
}
