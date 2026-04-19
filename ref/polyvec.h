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

/* ============================================================
 * mod 2q polyveck helpers (Phase 6b-1). See docs/NGCC_Sign/SHUTTLE_draft.md
 * sections 6 and 13.
 * ============================================================ */

#define polyveck_highbits_mod_2q SHUTTLE_NAMESPACE(polyveck_highbits_mod_2q)
void polyveck_highbits_mod_2q(polyveck *out, const polyveck *w);

/* Extract LSB of each of M*N coefficients. Output layout: M * (N/8) bytes.
 * Bitmap for row i starts at offset i * N/8. */
#define polyveck_lsb_extract SHUTTLE_NAMESPACE(polyveck_lsb_extract)
void polyveck_lsb_extract(uint8_t bitmap[SHUTTLE_M * SHUTTLE_N / 8],
                          const polyveck *w);

/* Lift: out <- 2 * u_mod_q + q * (Y0 & 1) * j  (mod 2q).
 * u_mod_q is a polyveck; Y0 is a single polynomial whose LSB is broadcast
 * to all M rows (since q*j*Y0 contributes identically to each row). */
#define polyveck_lift_to_2q SHUTTLE_NAMESPACE(polyveck_lift_to_2q)
void polyveck_lift_to_2q(polyveck *out,
                         const polyveck *u_mod_q,
                         const poly *Y0);

/* In-place w <- (w - 2*z2) mod 2q, component-wise. */
#define polyveck_sub_2z2_mod2q SHUTTLE_NAMESPACE(polyveck_sub_2z2_mod2q)
void polyveck_sub_2z2_mod2q(polyveck *w, const polyveck *z2);

/* Compute the full mod 2q commitment:
 *   comY = hat_A * v (mod 2q)
 *         where hat_A = [2(a_gen - b) + q*j | 2*A_gen | 2*I_M].
 *
 * Implemented via the existing mod q NTT + lift trick:
 *   1. Compute U_mod_q = (a_gen - b)*v[0] + A_gen*v[1..L] + v[L+1..L+M]
 *      by standard NTT-domain multiplication, INTT, reduce, caddq.
 *   2. comY = 2 * U_mod_q + q * (v[0] & 1) * j  (mod 2q), per coefficient.
 *
 * Arguments:
 *   comY       : output polyveck, coefficients in [0, 2q).
 *   a_gen_hat  : NTT form of a_gen (polyveck).
 *   A_gen_hat  : NTT forms of A_gen columns (polyveck[L]).
 *   b_hat      : NTT form of public-key vector b (polyveck).
 *   v          : input polyvec (VECLEN) in coefficient form. v->vec[0] is Y_0
 *                (already CompressY'd), v->vec[1..L+M] are the remaining
 *                input slots. Caller retains ownership; function NTTs a
 *                private copy. */
#define compute_commitment_mod2q SHUTTLE_NAMESPACE(compute_commitment_mod2q)
void compute_commitment_mod2q(polyveck *comY,
                              const polyveck *a_gen_hat,
                              const polyveck A_gen_hat[SHUTTLE_L],
                              const polyveck *b_hat,
                              const polyvec *v);

/* ============================================================
 * Polyveck-level MakeHint / UseHint + basic hint packing (Phase 6b-3).
 * ============================================================ */

#define polyveck_make_hint_mod2q SHUTTLE_NAMESPACE(polyveck_make_hint_mod2q)
void polyveck_make_hint_mod2q(polyveck *h,
                              const polyveck *w,
                              const polyveck *z2);

#define polyveck_use_hint_wh_mod2q SHUTTLE_NAMESPACE(polyveck_use_hint_wh_mod2q)
void polyveck_use_hint_wh_mod2q(polyveck *w_h,
                                const polyveck *tilde_w,
                                const polyveck *h);

#define polyveck_recover_z2_mod2q SHUTTLE_NAMESPACE(polyveck_recover_z2_mod2q)
void polyveck_recover_z2_mod2q(polyveck *z2_out,
                               const polyveck *w_h,
                               const polyveck *w_0,
                               const polyveck *tilde_w);

/* LSB-polyveck helper: take a packed bitmap (M * N/8 bytes) and expand it
 * back into a polyveck where each coefficient is 0 or 1. Inverse of
 * polyveck_lsb_extract. Used at verify time to turn the LSB bitmap
 * comY_0 back into a polyveck for recover_z2. */
#define polyveck_lsb_from_bitmap SHUTTLE_NAMESPACE(polyveck_lsb_from_bitmap)
void polyveck_lsb_from_bitmap(polyveck *out,
                              const uint8_t bitmap[SHUTTLE_M * SHUTTLE_N / 8]);

/* ============================================================
 * Phase 6b-5c helpers: keygen b and verify-side tilde_w under Alg 2.
 * ============================================================ */

/* Compute b = a_gen + A_gen*s + e (mod q), with the identity-block
 * addition performed in the coefficient domain to avoid the mixed-
 * Montgomery scaling that Phase 6a's compute_commitment has.
 *
 * Inputs:
 *   a_gen      : polyveck (coefficient form).
 *   A_gen_hat  : polyveck[L] in NTT form.
 *   s          : polyvecl (coefficient form, CBD-sampled secret).
 *   e          : polyveck (coefficient form, CBD-sampled secret).
 * Output:
 *   b          : polyveck in [0, q). */
#define compute_b SHUTTLE_NAMESPACE(compute_b)
void compute_b(polyveck *b,
                  const polyveck *a_gen,
                  const polyveck A_gen_hat[SHUTTLE_L],
                  const polyvecl *s,
                  const polyveck *e);

/* Compute verifier's tilde_w = hat_A_1 * z_1 - q*c*j (mod 2q), where
 *   hat_A_1 = [2(a_gen - b) + q*j | 2 * A_gen_hat]
 *   z_1     = (Z_0, z[1..L])
 *
 * Algebraically:
 *   tilde_w = 2 * U_ver + q * (Z_0 - c) * j     (mod 2q)
 * with U_ver = (a_gen - b) * Z_0 + A_gen * z[1..L]  (mod q).
 *
 * Inputs (all in NTT form except z_1 which is coefficient form):
 *   a_gen_hat, b_hat : polyveck.
 *   A_gen_hat        : polyveck[L].
 *   z_1              : array of 1 + L polys (Z_0 in z_1[0], rest in z_1[1..L]).
 *   c_poly           : challenge polynomial (coefficient form).
 *
 * Output:
 *   tilde_w          : polyveck in [0, 2q). */
#define compute_tilde_w_mod2q SHUTTLE_NAMESPACE(compute_tilde_w_mod2q)
void compute_tilde_w_mod2q(polyveck *tilde_w,
                           const polyveck *a_gen_hat,
                           const polyveck A_gen_hat[SHUTTLE_L],
                           const polyveck *b_hat,
                           const poly *z_1,   /* length 1 + L */
                           const poly *c_poly);

/* Basic (non-rANS) hint encoding: one signed int8 per coefficient.
 * Size = M * N bytes. Safe for current modes where |h|_inf << 128.
 * Phase 6b-5 will replace this with rANS compression. */
#define SHUTTLE_HINT_PACKEDBYTES_BASIC (SHUTTLE_M * SHUTTLE_N)

#define polyveck_hint_pack_basic SHUTTLE_NAMESPACE(polyveck_hint_pack_basic)
void polyveck_hint_pack_basic(uint8_t out[SHUTTLE_HINT_PACKEDBYTES_BASIC],
                              const polyveck *h);

#define polyveck_hint_unpack_basic SHUTTLE_NAMESPACE(polyveck_hint_unpack_basic)
void polyveck_hint_unpack_basic(polyveck *h,
                                const uint8_t in[SHUTTLE_HINT_PACKEDBYTES_BASIC]);

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
