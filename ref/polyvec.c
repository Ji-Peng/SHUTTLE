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

/**************************************************************/
/*********** mod 2q polyveck helpers (Phase 6b-1) ************/
/**************************************************************/

void polyveck_highbits_mod_2q(polyveck *out, const polyveck *w) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_highbits_mod_2q(&out->vec[i], &w->vec[i]);
}

void polyveck_lsb_extract(uint8_t bitmap[SHUTTLE_M * SHUTTLE_N / 8],
                          const polyveck *w)
{
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_lsb_extract(&bitmap[i * (SHUTTLE_N / 8)], &w->vec[i]);
}

void polyveck_lift_to_2q(polyveck *out,
                         const polyveck *u_mod_q,
                         const poly *Y0)
{
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_lift_to_2q(&out->vec[i], &u_mod_q->vec[i], Y0);
}

void polyveck_sub_2z2_mod2q(polyveck *w, const polyveck *z2) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_sub_2z2_mod2q(&w->vec[i], &z2->vec[i]);
}

/*************************************************
* Name:        compute_commitment_mod2q
*
* Description: comY = hat_A * v (mod 2q) via mod q NTT + lift.
*              See polyvec.h for argument semantics.
**************************************************/
void compute_commitment_mod2q(polyveck *comY,
                              const polyveck *a_gen_hat,
                              const polyveck A_gen_hat[SHUTTLE_L],
                              const polyveck *b_hat,
                              const polyvec *v)
{
  polyveck U;
  polyvec v_hat;
  poly tmp;
  unsigned int i, j;

  /* NTT-transform the input (private copy). */
  v_hat = *v;
  polyvec_ntt(&v_hat);

  /* Step 1: Build U_modq in NTT domain:
   *   U[i] = (a_gen[i] - b[i]) * v_hat[0] + sum_j A_gen_hat[j][i] * v_hat[1+j]
   *
   * Note: we DO NOT include the identity-block term v_hat[1+L+i] here because
   * poly_pointwise_montgomery introduces a 1/MONT factor while a bare addition
   * of NTT-form values does not. Mixing the two would force invntt_tomont to
   * apply an inconsistent Montgomery correction. We add the identity block in
   * the coefficient domain after invntt_tomont, where both operands are
   * plain integers mod q. */
  for(i = 0; i < SHUTTLE_M; ++i) {
    poly_pointwise_montgomery(&U.vec[i], &a_gen_hat->vec[i], &v_hat.vec[0]);

    poly_pointwise_montgomery(&tmp, &b_hat->vec[i], &v_hat.vec[0]);
    poly_sub(&U.vec[i], &U.vec[i], &tmp);

    for(j = 0; j < SHUTTLE_L; ++j) {
      poly_pointwise_montgomery(&tmp, &A_gen_hat[j].vec[i], &v_hat.vec[1 + j]);
      poly_add(&U.vec[i], &U.vec[i], &tmp);
    }
  }

  /* INTT and canonicalize to [0, q). */
  polyveck_invntt_tomont(&U);
  polyveck_reduce(&U);
  polyveck_caddq(&U);

  /* Step 2: Add the identity-block contribution v[1+L+i] in coefficient
   * domain. Uses the ORIGINAL (pre-NTT) v, not v_hat. */
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_add(&U.vec[i], &U.vec[i], &v->vec[1 + SHUTTLE_L + i]);
  polyveck_reduce(&U);
  polyveck_caddq(&U);

  /* Step 3: Lift to mod 2q using pre-NTT v->vec[0] (Y_0) parity bits. */
  polyveck_lift_to_2q(comY, &U, &v->vec[0]);
}

/*************************************************
* Name:        polyveck_make_hint_mod2q
**************************************************/
void polyveck_make_hint_mod2q(polyveck *h,
                              const polyveck *w,
                              const polyveck *z2)
{
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_make_hint_mod2q(&h->vec[i], &w->vec[i], &z2->vec[i]);
}

/*************************************************
* Name:        polyveck_use_hint_wh_mod2q
**************************************************/
void polyveck_use_hint_wh_mod2q(polyveck *w_h,
                                const polyveck *tilde_w,
                                const polyveck *h)
{
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_use_hint_wh_mod2q(&w_h->vec[i], &tilde_w->vec[i], &h->vec[i]);
}

/*************************************************
* Name:        polyveck_recover_z2_mod2q
**************************************************/
void polyveck_recover_z2_mod2q(polyveck *z2_out,
                               const polyveck *w_h,
                               const polyveck *w_0,
                               const polyveck *tilde_w)
{
  unsigned int i;
  for(i = 0; i < SHUTTLE_M; ++i)
    poly_recover_z2_mod2q(&z2_out->vec[i], &w_h->vec[i],
                          &w_0->vec[i], &tilde_w->vec[i]);
}

/*************************************************
* Name:        polyveck_lsb_from_bitmap
*
* Description: Inverse of polyveck_lsb_extract. Each bit becomes a 0/1
*              coefficient of the output polyveck.
**************************************************/
void polyveck_lsb_from_bitmap(polyveck *out,
                              const uint8_t bitmap[SHUTTLE_M * SHUTTLE_N / 8])
{
  unsigned int i, j;
  for(i = 0; i < SHUTTLE_M; ++i) {
    const uint8_t *row = &bitmap[i * (SHUTTLE_N / 8)];
    for(j = 0; j < SHUTTLE_N; ++j)
      out->vec[i].coeffs[j] = (int32_t)((row[j >> 3] >> (j & 7)) & 1);
  }
}

/*************************************************
* Name:        polyveck_hint_pack_basic
*
* Description: Pack hint polyveck using 1 signed byte per coefficient.
*              For current modes, |h|_inf <= 18 (mode-128) or 14 (mode-256),
*              well within int8 range. Total size M * N bytes.
**************************************************/
void polyveck_hint_pack_basic(uint8_t out[SHUTTLE_HINT_PACKEDBYTES_BASIC],
                              const polyveck *h)
{
  unsigned int i, j;
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      out[i * SHUTTLE_N + j] = (uint8_t)(int8_t)h->vec[i].coeffs[j];
}

/*************************************************
* Name:        polyveck_hint_unpack_basic
**************************************************/
void polyveck_hint_unpack_basic(polyveck *h,
                                const uint8_t in[SHUTTLE_HINT_PACKEDBYTES_BASIC])
{
  unsigned int i, j;
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h->vec[i].coeffs[j] = (int32_t)(int8_t)in[i * SHUTTLE_N + j];
}

/*************************************************
* Name:        compute_b
*
* Description: Clean keygen b = a_gen + A_gen*s + e (mod q). Identity-block
*              addition is done in the coefficient domain to avoid the
*              mixed-Montgomery scaling of Phase 6a's compute_commitment.
**************************************************/
void compute_b(polyveck *b,
                  const polyveck *a_gen,
                  const polyveck A_gen_hat[SHUTTLE_L],
                  const polyvecl *s,
                  const polyveck *e)
{
  unsigned int i, j;
  polyvecl s_hat;
  poly tmp, acc;

  /* NTT(s) — private copy. */
  s_hat = *s;
  polyvecl_ntt(&s_hat);

  for(i = 0; i < SHUTTLE_M; ++i) {
    /* acc = sum_j A_gen_hat[j][i] * s_hat[j]  (NTT domain, /MONT per product) */
    poly_pointwise_montgomery(&acc, &A_gen_hat[0].vec[i], &s_hat.vec[0]);
    for(j = 1; j < SHUTTLE_L; ++j) {
      poly_pointwise_montgomery(&tmp, &A_gen_hat[j].vec[i], &s_hat.vec[j]);
      poly_add(&acc, &acc, &tmp);
    }
    /* INTT + canonicalize to [0, q). */
    poly_invntt_tomont(&acc);
    poly_reduce(&acc);
    poly_caddq(&acc);

    /* b[i] = a_gen[i] + acc + e[i]  (coefficient domain). */
    poly_add(&b->vec[i], &a_gen->vec[i], &acc);
    poly_add(&b->vec[i], &b->vec[i], &e->vec[i]);
    poly_reduce(&b->vec[i]);
    poly_caddq(&b->vec[i]);
  }
}

/*************************************************
* Name:        compute_tilde_w_mod2q
*
* Description: Verifier's tilde_w = hat_A_1 * z_1 - q*c*j (mod 2q). See
*              polyvec.h for the full formula and argument semantics.
*
*              Implemented via:
*                U_ver = (a - b) * Z_0 + A_gen * z_1[1..L]  (mod q)
*                tilde_w = lift(U_ver, Z_0 - c)             (mod 2q)
*              The identity-block addition is absent (verifier has no z_2).
**************************************************/
void compute_tilde_w_mod2q(polyveck *tilde_w,
                           const polyveck *a_gen_hat,
                           const polyveck A_gen_hat[SHUTTLE_L],
                           const polyveck *b_hat,
                           const poly *z_1,
                           const poly *c_poly)
{
  unsigned int i, j;
  polyveck U_ver;
  poly tmp;

  /* NTT(z_1) on a private array of 1 + L polys. */
  poly z_hat[1 + SHUTTLE_L];
  for(i = 0; i < 1 + SHUTTLE_L; ++i) {
    z_hat[i] = z_1[i];
    poly_ntt(&z_hat[i]);
  }

  for(i = 0; i < SHUTTLE_M; ++i) {
    /* U_ver[i] = (a_gen[i] - b[i]) * Z_0  (NTT domain, /MONT per product). */
    poly_pointwise_montgomery(&U_ver.vec[i], &a_gen_hat->vec[i], &z_hat[0]);
    poly_pointwise_montgomery(&tmp, &b_hat->vec[i], &z_hat[0]);
    poly_sub(&U_ver.vec[i], &U_ver.vec[i], &tmp);

    /* + sum_j A_gen_hat[j][i] * z_hat[1+j]. */
    for(j = 0; j < SHUTTLE_L; ++j) {
      poly_pointwise_montgomery(&tmp, &A_gen_hat[j].vec[i], &z_hat[1 + j]);
      poly_add(&U_ver.vec[i], &U_ver.vec[i], &tmp);
    }
  }

  /* INTT + canonicalize to [0, q). */
  polyveck_invntt_tomont(&U_ver);
  polyveck_reduce(&U_ver);
  polyveck_caddq(&U_ver);

  /* Lift parity = (Z_0 - c) mod 2. Compute a throwaway poly for parity. */
  poly parity;
  for(j = 0; j < SHUTTLE_N; ++j)
    parity.coeffs[j] = z_1[0].coeffs[j] - c_poly->coeffs[j];

  polyveck_lift_to_2q(tilde_w, &U_ver, &parity);
}
