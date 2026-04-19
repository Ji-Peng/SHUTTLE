#include <stdint.h>
#include "params.h"
#include "reduce.h"
#include "rounding.h"

/*
 * Decomposition for SHUTTLE.
 *
 * For a standard representative a in [0, Q-1]:
 *   a1 = (a + alpha_h) / (2*alpha_h)              (rounded)
 *   a0 = a - a1 * 2*alpha_h                        (in (-alpha_h, alpha_h])
 *
 * With q=13313 the natural range of a1 is [0, W1_MAX] and the edge case
 * "a1 overflows to W1_MAX+1" that existed for q=15361 never occurs because
 * (q-1) is exactly divisible by 2*alpha_h for both modes here (q-1 = 13312 =
 * 52*256 = 26*512), so the old ad-hoc wrap guard is unnecessary.
 *
 * 2*alpha_h is always a power of two (256 for mode-128, 512 for mode-256),
 * so the compiler lowers the division to a right shift.
 */

/*************************************************
* Name:        decompose
*
* Description: For finite field element a (standard representative in
*              [0,Q-1]), compute high and low bits a0, a1 such that
*              a mod^+ Q = a1 * 2*ALPHA_H + a0 with -ALPHA_H < a0 <= ALPHA_H.
*
* Arguments:   - int32_t *a0: pointer to output low part
*              - int32_t a:   input element (standard representative)
*
* Returns a1 (high part).
**************************************************/
int32_t decompose(int32_t *a0, int32_t a) {
  int32_t a1;

  a1 = (a + SHUTTLE_ALPHA_H) / (2 * SHUTTLE_ALPHA_H);
  *a0 = a - a1 * (2 * SHUTTLE_ALPHA_H);
  return a1;
}

/*************************************************
* Name:        make_hint
*
* Description: Compute hint bit indicating whether the low bits of the
*              input element overflow into the high bits.
*              Returns 1 if adding z2 to comY would change the high bits.
*
* Arguments:   - int32_t z2:   low bits perturbation
*              - int32_t comY: low part of commitment
*
* Returns 1 if overflow.
**************************************************/
unsigned int make_hint(int32_t z2, int32_t comY) {
  if(z2 > (int32_t)SHUTTLE_ALPHA_H
     || z2 < -(int32_t)SHUTTLE_ALPHA_H
     || (z2 == -(int32_t)SHUTTLE_ALPHA_H && comY != 0))
    return 1;

  return 0;
}

/*************************************************
* Name:        use_hint
*
* Description: Correct high bits according to hint. High bits live in
*              [0, W1_MAX]; the ±1 adjustment wraps modulo (W1_MAX + 1).
*
* Arguments:   - int32_t a: input element (standard representative)
*              - unsigned int hint: hint bit
*
* Returns corrected high bits.
**************************************************/
int32_t use_hint(int32_t a, unsigned int hint) {
  int32_t a0, a1;

  a1 = decompose(&a0, a);

  if(hint == 0)
    return a1;

  if(a0 > 0)
    return (a1 == SHUTTLE_W1_MAX) ? 0 : a1 + 1;
  else
    return (a1 == 0) ? SHUTTLE_W1_MAX : a1 - 1;
}

/*************************************************
* Name:        highbits_mod_2q
*
* Description: Round-half-up HighBits for mod 2q commitment. Mirrors the
*              spec's <x/alpha_h> and HAETAE's decompose_hint.
*
*              Edge case: round(w/alpha_h) can reach HINT_MAX (at the top
*              of [0, 2q)); wrap it back to 0 to stay in [0, HINT_MAX).
*              Implemented with bit-mask for constant-time.
*
* Arguments:   - int32_t w: input, expected in [0, 2q).
*
* Returns bucket index in [0, HINT_MAX).
**************************************************/
int32_t highbits_mod_2q(int32_t w) {
  int32_t hb;
  int32_t edgecase;

  /* (w + alpha_h/2) >> log2(alpha_h) == (w + alpha_h/2) / alpha_h. */
  hb = (w + SHUTTLE_HALF_ALPHA_H) >> SHUTTLE_ALPHA_H_BITS;

  /* edgecase = -1 if hb >= HINT_MAX else 0. Using sign extension of the
   * difference: if (HINT_MAX - 1 - hb) < 0 then hb >= HINT_MAX. */
  edgecase = (SHUTTLE_HINT_MAX - 1 - hb) >> 31;
  hb -= SHUTTLE_HINT_MAX & edgecase;
  return hb;
}

/*************************************************
* Name:        lift_to_2q
*
* Description: Combine the "even part" (2 * u, with u in [0, q)) and the
*              "q-parity part" (q * bit) into a single [0, 2q) residue.
*              This is the per-coefficient operation behind compute_commitment
*              when moving the output of standard mod q NTT into the mod 2q
*              domain demanded by the derived public matrix hat_A.
*
* Arguments:   - int32_t u_mod_q: input in [0, q).
*              - int32_t parity:  0 or 1.
*
* Returns (2*u + q*parity) mod 2q in [0, 2q). Constant-time.
**************************************************/
int32_t lift_to_2q(int32_t u_mod_q, int32_t parity) {
  int32_t r;
  int32_t over;

  r = (u_mod_q << 1) + ((parity & 1) ? SHUTTLE_Q : 0);

  /* r is in [0, 2q) + [0, q) = [0, 3q). If r >= 2q, subtract 2q.
   * Constant-time using sign bit of (2q - 1 - r). */
  over = (SHUTTLE_DQ - 1 - r) >> 31;   /* -1 if r >= 2q else 0 */
  r -= SHUTTLE_DQ & over;
  return r;
}

/*************************************************
* Name:        make_hint_mod2q
*
* Description: Compute MakeHint difference per coefficient (spec Alg 9).
*              See rounding.h for the exact formula.
*
*              Returns the centered representative in (-HINT_MAX/2, HINT_MAX/2]:
*              the raw difference w_h - tilde_w_h lives in (-HINT_MAX, HINT_MAX),
*              and near bucket-wrap boundaries it would otherwise drift to
*              magnitude ~HINT_MAX. Centering unifies the two equivalent
*              representations (e.g. -207 == +1 mod HINT_MAX for mode-128)
*              into a single small-magnitude distribution suitable for
*              rANS encoding.
**************************************************/
int32_t make_hint_mod2q(int32_t w, int32_t z2_coef) {
  int32_t w_h = highbits_mod_2q(w);
  int32_t tilde_w = reduce_mod_2q(w - 2 * z2_coef);
  int32_t tilde_w_h = highbits_mod_2q(tilde_w);
  int32_t h = w_h - tilde_w_h;

  /* Center to (-HINT_MAX/2, HINT_MAX/2], constant-time. */
  int32_t half = SHUTTLE_HINT_MAX >> 1;
  int32_t upper_mask = (half - h) >> 31;      /* -1 if h > half */
  h -= SHUTTLE_HINT_MAX & upper_mask;
  int32_t lower_mask = (h + half) >> 31;      /* -1 if h < -half */
  h += SHUTTLE_HINT_MAX & lower_mask;
  return h;
}

/*************************************************
* Name:        use_hint_wh_mod2q
*
* Description: Recover signer's w_h from tilde_w and hint h (spec Alg 10
*              first half). Returns w_h in [0, HINT_MAX). Constant-time.
**************************************************/
int32_t use_hint_wh_mod2q(int32_t tilde_w, int32_t h) {
  int32_t tilde_w_h = highbits_mod_2q(tilde_w);
  int32_t w_h = tilde_w_h + h;
  int32_t over_mask, under_mask;

  /* Constant-time modular normalization to [0, HINT_MAX):
   *   if w_h >= HINT_MAX  ->  w_h -= HINT_MAX
   *   if w_h < 0          ->  w_h += HINT_MAX
   * We use sign-bit masking of (HINT_MAX - 1 - w_h) and (w_h) respectively. */
  over_mask  = (SHUTTLE_HINT_MAX - 1 - w_h) >> 31;   /* -1 if w_h >= HINT_MAX */
  w_h -= SHUTTLE_HINT_MAX & over_mask;

  under_mask = w_h >> 31;                             /* -1 if w_h < 0 */
  w_h += SHUTTLE_HINT_MAX & under_mask;
  return w_h;
}

/*************************************************
* Name:        recover_z2_coef_mod2q
*
* Description: Recover z_2 coefficient (spec Alg 10 second half).
*              See rounding.h for the exact formula and the +/-alpha_h/2
*              rounding-error bound.
**************************************************/
int32_t recover_z2_coef_mod2q(int32_t w_h, int32_t w_0, int32_t tilde_w) {
  int32_t w_app = SHUTTLE_ALPHA_H * w_h + w_0;
  int32_t diff  = w_app - tilde_w;
  int32_t over_mask, under_mask;

  /* diff currently in (-2q, 2q). Canonicalize to (-q, q]:
   *   if diff >  q  ->  diff -= 2q
   *   if diff <= -q ->  diff += 2q */
  over_mask  = (SHUTTLE_Q - diff) >> 31;          /* -1 if diff > q */
  diff -= SHUTTLE_DQ & over_mask;

  under_mask = (diff + SHUTTLE_Q) >> 31;          /* -1 if diff < -q */
  diff += SHUTTLE_DQ & under_mask;

  /* Integer division by 2 (arithmetic shift for two's complement). */
  return diff >> 1;
}
