/*
 * ntt.c - Forward / inverse NTT for SHUTTLE.
 *
 * Two NTT flavours selected by SHUTTLE_BASE_DEG:
 *
 *   base_deg = 1  (SHUTTLE-128, full NTT)
 *     log2(N) Cooley-Tukey layers. After ntt() the polynomial is in N
 *     point-value form: a[i] = f(rho^{2*brv(i)+1}) for primitive 2N-th
 *     root rho. Pointwise multiplication is plain int * int.
 *
 *   base_deg = 2  (SHUTTLE-256, SHUTTLE-512, incomplete NTT)
 *     log2(N) - 1 Cooley-Tukey layers. After ntt() the polynomial is in
 *     N/2 basecase form: pairs (a[2k], a[2k+1]) live in
 *     Z_q[X]/(X^2 - zetas_basemul[k]) for a primitive N-th root rho.
 *     "Pointwise" multiplication becomes a basemul over each basecase.
 *
 * The zetas[] table (size SHUTTLE_NTT_ZETAS_SIZE) is shared between
 * butterflies and basemul (Kyber convention):
 *   - butterflies use zetas[1..ZETAS_SIZE-1]
 *   - basemul uses zetas[ZETAS_SIZE/2 .. ZETAS_SIZE-1] with +/- pairs.
 *
 * Constants are generated offline by gen_zetas.py.
 */

#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"
#include "ntt_constants.h"

/*************************************************
* Name:        ntt
*
* Description: Forward NTT, in-place. No modular reduction is performed
*              after additions or subtractions. Output ordering is the
*              standard Cooley-Tukey bit-reversed layout.
*
*              For base_deg = 1 we run log2(N) layers (full NTT).
*              For base_deg = 2 we run log2(N) - 1 layers (skipping the
*              last len=1 layer); the result is a sequence of N/2 deg-1
*              polynomials in basecase rings Z_q[X]/(X^2 - zeta).
*
* Arguments:   - int32_t a[SHUTTLE_N]: input/output coefficient array
**************************************************/
void ntt(int32_t a[SHUTTLE_N]) {
  unsigned int len, start, j, k;
  int32_t zeta, t;

  k = 0;
  for(len = SHUTTLE_N / 2; len >= (unsigned)SHUTTLE_BASE_DEG; len >>= 1) {
    for(start = 0; start < SHUTTLE_N; start = j + len) {
      zeta = zetas[++k];
      for(j = start; j < start + len; ++j) {
        t = montgomery_reduce((int64_t)zeta * a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}

/*************************************************
* Name:        invntt_tomont
*
* Description: Inverse NTT and multiplication by Montgomery factor 2^32.
*              In-place. Input coefficients need to be smaller than Q in
*              absolute value. Output coefficients are smaller than Q in
*              absolute value.
*
*              The number of inverse butterfly layers matches the forward
*              NTT (log2(N) for full, log2(N)-1 for incomplete). The final
*              scale factor INVNTT_F = MONT^2 / NTT_ZETAS_SIZE folds the
*              standard 1/N normalization with the Montgomery correction
*              and the residual N/2 factor (for incomplete NTT) into a
*              single multiply.
*
* Arguments:   - int32_t a[SHUTTLE_N]: input/output coefficient array
**************************************************/
void invntt_tomont(int32_t a[SHUTTLE_N]) {
  unsigned int start, len, j, k;
  int32_t t, zeta;
  const int32_t f = SHUTTLE_INVNTT_F;

  k = SHUTTLE_NTT_ZETAS_SIZE;
  for(len = SHUTTLE_BASE_DEG; len < SHUTTLE_N; len <<= 1) {
    for(start = 0; start < SHUTTLE_N; start = j + len) {
      zeta = -zetas[--k];
      for(j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
      }
    }
  }

  for(j = 0; j < SHUTTLE_N; ++j) {
    a[j] = montgomery_reduce((int64_t)f * a[j]);
  }
}

#if SHUTTLE_BASE_DEG == 2
/*************************************************
* Name:        basemul
*
* Description: One Kyber-style basemul over Z_q[X]/(X^2 - zeta):
*                r0 = a0*b0 + a1*b1*zeta   (Montgomery-reduced)
*                r1 = a0*b1 + a1*b0
*              All operands int32_t; output stays in Montgomery form.
*
* Arguments:   - int32_t r[2]: output basecase pair
*              - const int32_t a[2]: first input basecase pair
*              - const int32_t b[2]: second input basecase pair
*              - int32_t zeta: basecase modulus exponent in Montgomery form
**************************************************/
static inline void basemul(int32_t r[2],
                           const int32_t a[2],
                           const int32_t b[2],
                           int32_t zeta)
{
  int32_t t;
  t    = montgomery_reduce((int64_t)a[1] * b[1]);
  t    = montgomery_reduce((int64_t)t * zeta);
  t   += montgomery_reduce((int64_t)a[0] * b[0]);
  r[0] = t;
  t    = montgomery_reduce((int64_t)a[0] * b[1]);
  t   += montgomery_reduce((int64_t)a[1] * b[0]);
  r[1] = t;
}

/*************************************************
* Name:        poly_basemul_montgomery_native
*
* Description: Apply basemul to all N/2 basecase pairs of a NTT-domain
*              polynomial. Paired basecases share a |zeta| value with
*              opposite signs:
*                pair (4i, 4i+1) uses  zetas[ZETAS_SIZE/2 + i]
*                pair (4i+2, 4i+3) uses -zetas[ZETAS_SIZE/2 + i]
*
* Arguments:   - int32_t r[SHUTTLE_N]: output coefficient array
*              - const int32_t a[SHUTTLE_N]: first input (NTT domain)
*              - const int32_t b[SHUTTLE_N]: second input (NTT domain)
**************************************************/
void poly_basemul_montgomery_native(int32_t r[SHUTTLE_N],
                                    const int32_t a[SHUTTLE_N],
                                    const int32_t b[SHUTTLE_N])
{
  unsigned int i;
  const unsigned int half = SHUTTLE_NTT_ZETAS_SIZE / 2;
  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    int32_t zeta = zetas[half + i];
    basemul(&r[4*i],     &a[4*i],     &b[4*i],     zeta);
    basemul(&r[4*i + 2], &a[4*i + 2], &b[4*i + 2], -zeta);
  }
}
#endif /* SHUTTLE_BASE_DEG == 2 */
