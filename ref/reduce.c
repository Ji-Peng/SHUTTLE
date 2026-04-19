#include <stdint.h>
#include "params.h"
#include "reduce.h"

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a with -2^{31}*Q <= a <= Q*2^{31},
*              compute r \equiv a*2^{-32} (mod Q) such that -Q < r < Q.
*
* Arguments:   - int64_t: finite field element a
*
* Returns r.
**************************************************/
int32_t montgomery_reduce(int64_t a) {
  int32_t t;

  t = (int32_t)((uint32_t)a * SHUTTLE_QINV);
  t = (a - (int64_t)t * SHUTTLE_Q) >> 32;
  return t;
}

/*************************************************
* Name:        reduce32
*
* Description: For finite field element a in [-2^31, 2^31), compute r
*              congruent to a mod Q with r in (-Q, Q).
*
*              Because Q = 13313 is not close to a power of two, the
*              (a + 2^13) >> 14 approximation used for q=15361 is too loose
*              and a Barrett estimator has off-by-one cases around r=Q.
*              The direct integer-% operator is exact and fast enough for
*              the reference implementation (the AVX2 port will carry a
*              dedicated Barrett routine).
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t reduce32(int32_t a) {
  return a % SHUTTLE_Q;
}

/*************************************************
* Name:        caddq
*
* Description: Add Q if input coefficient is negative.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t caddq(int32_t a) {
  a += (a >> 31) & SHUTTLE_Q;
  return a;
}

/*************************************************
* Name:        freeze
*
* Description: For finite field element a, compute standard representative
*              r = a mod^+ Q.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t freeze(int32_t a) {
  a = reduce32(a);
  a = caddq(a);
  return a;
}
