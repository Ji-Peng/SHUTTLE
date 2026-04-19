#include <stdint.h>
#include "params.h"
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
