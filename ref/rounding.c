#include <stdint.h>
#include "params.h"
#include "rounding.h"

/*
 * Decomposition for NGCC_SIGN with alpha_h = 128, so 2*alpha_h = 256 = 2^8.
 *
 * For a standard representative a in [0, Q-1]:
 *   a1 = floor((a + 128) / 256) = (a + 128) >> 8
 *   a0 = a - a1 * 256
 *
 * Maximum a1 value: (15360 + 128) >> 8 = 60.
 * When a1 == 60, we wrap: set a1 = 0 and a0 = a - Q (negative).
 * This ensures a = a1 * 256 + a0  (mod Q) always holds.
 */

/*************************************************
* Name:        decompose
*
* Description: For finite field element a (standard representative in [0,Q-1]),
*              compute high and low bits a0, a1 such that
*              a mod^+ Q = a1 * 2*ALPHA_H + a0
*              with -ALPHA_H < a0 <= ALPHA_H,
*              except when a1 would equal (Q-1)/(2*ALPHA_H)+1 = 60,
*              in which case we set a1 = 0 and a0 = a - Q < 0.
*
* Arguments:   - int32_t *a0: pointer to output low part
*              - int32_t a: input element (standard representative)
*
* Returns a1 (high part).
**************************************************/
int32_t decompose(int32_t *a0, int32_t a) {
  int32_t a1;

  /* a1 = (a + alpha_h) / (2*alpha_h) = (a + 128) >> 8 */
  a1 = (a + NGCC_SIGN_ALPHA_H) >> 8;

  /* Handle wrap-around: Q = 15361, (Q-1)/256 = 60.
   * When a1 reaches 60, the high bits would exceed
   * the valid range, so wrap to 0. */
  a1 ^= ((59 - a1) >> 31) & a1;  /* if a1 == 60, set a1 = 0 */

  *a0 = a - a1 * 2 * NGCC_SIGN_ALPHA_H;
  /* When a1 was set to 0 but a >= 60*256 = 15360,
   * we need a0 = a - Q to get negative remainder */
  *a0 -= (((NGCC_SIGN_Q - 1) / 2 - *a0) >> 31) & NGCC_SIGN_Q;

  return a1;
}

/*************************************************
* Name:        make_hint
*
* Description: Compute hint bit indicating whether the low bits of the
*              input element overflow into the high bits.
*              Returns 1 if adding z2 to comY changes the high bits.
*
* Arguments:   - int32_t z2: low bits perturbation
*              - int32_t comY: high bits of commitment
*
* Returns 1 if overflow.
**************************************************/
unsigned int make_hint(int32_t z2, int32_t comY) {
  if(z2 > (int32_t)NGCC_SIGN_ALPHA_H
     || z2 < -(int32_t)NGCC_SIGN_ALPHA_H
     || (z2 == -(int32_t)NGCC_SIGN_ALPHA_H && comY != 0))
    return 1;

  return 0;
}

/*************************************************
* Name:        use_hint
*
* Description: Correct high bits according to hint.
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

  /* Number of valid high-bit values: floor((Q-1) / (2*ALPHA_H)) = 60 */
  if(a0 > 0)
    return (a1 == 59) ? 0 : a1 + 1;
  else
    return (a1 == 0) ? 59 : a1 - 1;
}
