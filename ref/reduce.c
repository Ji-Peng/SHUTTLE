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
/* ============================================================
 * Constant-time Barrett reduction (division-free) for the scheme-domain
 * reductions, which run on SECRET data (NTT residues of s/y/z, commitment w).
 * `a % d` must not become a hardware/variable-time divide, and must not rely on
 * the compiler strength-reducing `% const`.  REC = floor(2^32 / d) is a
 * COMPILE-TIME constant (const/const, folded -- no runtime divide); SH=32 keeps
 * au*REC < 2^63 for any q < 2^16 with |a| <= 2^31, the shortfall < 0.5, so
 * (au*REC)>>32 + one branchless correction is the exact floor.  Bit-identical to
 * `a % d` (verified over the full int32 range for all SHUTTLE q).
 * ============================================================ */
#define SHUTTLE_BARRETT_SH 32
#define SHUTTLE_BARRETT_Q  (((uint64_t)1 << SHUTTLE_BARRETT_SH) / (uint64_t)SHUTTLE_Q)
#define SHUTTLE_BARRETT_2Q (((uint64_t)1 << SHUTTLE_BARRETT_SH) / (uint64_t)SHUTTLE_DQ)

static inline uint32_t barrett_mod_u(uint32_t au, uint32_t d, uint64_t REC) {
  uint32_t qh = (uint32_t)(((uint64_t)au * REC) >> SHUTTLE_BARRETT_SH);
  uint32_t r  = au - qh * d;                          /* in [0, 2d) */
  uint32_t ge = (uint32_t)0 - (uint32_t)(r >= d);     /* all-ones iff r >= d */
  return r - (ge & d);                                /* in [0, d) */
}

int32_t reduce32(int32_t a) {
  uint32_t sa = (uint32_t)(a >> 31);                  /* 0 or 0xFFFFFFFF */
  uint32_t au = ((uint32_t)a ^ sa) - sa;              /* |a| */
  uint32_t r  = barrett_mod_u(au, (uint32_t)SHUTTLE_Q, SHUTTLE_BARRETT_Q);
  return (int32_t)((r ^ sa) - sa);                    /* reapply sign (== a % Q) */
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

/*************************************************
* Name:        caddq2
*
* Description: Conditional add 2q. If input is negative, returns a + 2q;
*              otherwise returns a. Constant-time (sign-bit mask).
*
* Arguments:   - int32_t: input
*
* Returns a if a >= 0, else a + 2q.
**************************************************/
int32_t caddq2(int32_t a) {
  a += (a >> 31) & SHUTTLE_DQ;
  return a;
}

/*************************************************
* Name:        reduce_mod_2q
*
* Description: Reduce arbitrary int32_t to canonical residue in [0, 2q).
*              Uses integer modulus then positive-corrects.
*
*              This is the reference implementation; it relies on % which
*              the compiler lowers efficiently for 14-bit q. A Barrett
*              version is possible but unnecessary for ref code since this
*              helper is called only post-NTT (already small range).
*
* Arguments:   - int32_t: input
*
* Returns r in [0, 2q) with r congruent to a mod 2q.
**************************************************/
int32_t reduce_mod_2q(int32_t a) {
  uint32_t sa = (uint32_t)(a >> 31);
  uint32_t au = ((uint32_t)a ^ sa) - sa;              /* |a| */
  uint32_t ru = barrett_mod_u(au, (uint32_t)SHUTTLE_DQ, SHUTTLE_BARRETT_2Q);
  int32_t  r  = (int32_t)((ru ^ sa) - sa);            /* a % 2q (truncate) */
  r += (r >> 31) & SHUTTLE_DQ;                        /* make non-negative -> [0,2q) */
  return r;
}
