/*
 * approx_exp.c - Fixed-point exp(-a) for the SHUTTLE rejection sampler.
 *
 * This implements Algorithm 1 of agent/ApproxExp/ApproxExp.tex,
 * a six-step "exp(+)" pipeline with parameters chosen so the relative
 * precision matches the K = 53 security target across all SHUTTLE
 * parameter sets:
 *
 *   Input:  a_q60 -- Q60 fixed-point, a in [0, 11*ln2) ~ [0, 7.624)
 *   Output: round(exp(-a) * 2^63) -- Q63 unsigned (in [0, 2^63])
 *
 * The mathematical identity exploited is
 *
 *      exp(-a) = 2^{-m/B} * exp(r_l),
 *
 *   m   = ceil(B*a/ln2) + (B - 1),
 *   r_l = (m/B) * ln2 - a       in   [ (B-1)/B * ln2 ,  ln2 ),
 *
 * which has three useful structural properties:
 *
 *   (a) r_l is bounded away from zero -- the polynomial input never
 *       enters the "1 - exp(r_l) -> 0" angle of the old design, so the
 *       rejection error is no longer amplified by 1/r_l;
 *   (b) exp(r_l) lives in [exp(31/32 ln2), 2) ~ [1.957, 2), so on the
 *       fitting interval [31/32 ln2, ln2] every Sobolev-optimal degree-6
 *       coefficient is strictly positive -- the entire Horner chain runs
 *       unsigned (every Horner accumulator stays < 2^64 because
 *       p < exp(ln2) = 2 in Q63);
 *   (c) 2^{-m/B} factors exactly into  2^{-n} * 2^{-j/B}  (n = m / B,
 *       j = m mod B). The integer shift by n is loss-less, and the
 *       fractional factor 2^{-j/B} comes from the B-entry table T[].
 *
 * Step-by-step error budget (relative to exp(-a); see ApproxExp.tex
 * Section 4 for the full derivation):
 *
 *      E_aq60   <= 2^-59   (computed by the sampler -- see sampler.c)
 *      E_range  <= 2^-61   (m*LB and Q(64+beta) -> Q64 truncation)
 *      E_poly   <= 2^-61.31 (Sage Horner-aware Sobolev fit)
 *      E_table  <= 2^-63
 *      E_mul    <= 2^-62   (one mulh64 truncation in Step 5)
 *      E_shift  <= 2^-64 / exp(-a)
 *
 * The first five terms are a-independent and sum to ~2^-58. E_shift
 * grows with a; at the worst case a ~ 7.37 (SHUTTLE-512) it dominates
 * at ~2^-53.3, still leaving ~0.3 bit of margin over K = 53.
 *
 * Constant-time guarantees: there is no data-dependent control flow.
 *   - Step 1 has no division, only mulh64 + add + shift.
 *   - Step 4 scans the entire B-entry table with bit masks.
 *   - Step 5's final shift is a runtime variable, but only depends on
 *     n = m >> beta which is itself a non-secret function of a; the
 *     sampler treats acceptance as the only secret signal.
 */

#include "approx_exp.h"
#include "approx_exp_constants.h"

/* ============================================================
 * Sobolev-optimal degree-6 polynomial coefficients (Q63 unsigned).
 *
 *   P(r_l) = sum_{k=0..6} C_k * r_l^k  ~  exp(r_l)
 *
 * Region of validity: r_l in [31/32 * ln2, ln2].
 * Source: tools/LogPolyApprox/expx_31div32mulln2_ln2_53_64.txt
 *         (GALACTICS Sobolev fit + LLL Babai-rounding to Q63).
 *
 * Note that C_0 ~ 1.0000249 * 2^63 > 2^63, so it must live in a
 * uint64 (it does not fit in int64). Every coefficient is positive,
 * so Horner runs as a pure  mulh64 + add  chain in unsigned arithmetic.
 *
 * Measured precision (Sage, Q63, 4096 samples on the fit interval):
 *      max relative Horner error  = 2^-61.31  (worst x ~ 0.6798)
 *      max rejection Horner error = 2^-60.29
 *      Sobolev best-fit error     = 2^-60
 * ============================================================ */
#define APPROX_EXP_POLY_C0  UINT64_C(9223601430956470991)
#define APPROX_EXP_POLY_C1  UINT64_C(9221045155417009835)
#define APPROX_EXP_POLY_C2  UINT64_C(4621763433873803717)
#define APPROX_EXP_POLY_C3  UINT64_C(1513127296998826967)
#define APPROX_EXP_POLY_C4  UINT64_C( 418531209137423100)
#define APPROX_EXP_POLY_C5  UINT64_C(  48310490603772128)
#define APPROX_EXP_POLY_C6  UINT64_C(  25344359382375445)

/* ============================================================
 * Constant-time helpers
 * ============================================================ */

/* Constant-time B-entry table lookup.
 *
 * Iterates the entire table and OR-masks the matching slot. The
 * sampler treats j as private (it is a deterministic function of a,
 * which is in turn derived from the secret seed via the CDT), so a
 * data-dependent indexed load would leak via the cache. */
static inline uint64_t ct_lookup_table(const uint64_t *table, uint64_t idx,
                                       uint64_t size)
{
    uint64_t r = 0;
    for (uint64_t i = 0; i < size; ++i) {
        /* mask = -((i == idx) ? 1 : 0)  computed without branching. */
        uint64_t diff = i ^ idx;
        /* If diff == 0:  diff | -diff  has bit 63 = 0  -> mask = -1.
         * If diff != 0:  diff | -diff  has bit 63 = 1  -> mask =  0. */
        uint64_t mask = ((diff | (uint64_t)(-(int64_t)diff)) >> 63) - 1U;
        r |= table[i] & mask;
    }
    return r;
}

/* ============================================================
 * Main entry point
 * ============================================================ */
uint64_t approx_exp(uint64_t a_q60) {
    /* --------------------------------------------------------------
     * Step 1. Compute  m = ceil(B*a/ln2) + (B - 1).
     *
     * Using the Q57 reciprocal  RCP = floor((B/ln2) * 2^57):
     *   V = mulh64(a_q60, RCP)  =  floor((B*a/ln2) * 2^{60+57-64})
     *                          =  floor((B*a/ln2) * 2^53).
     * The next "+(2^53 - 1)) >> 53" converts floor to ceiling.
     *
     * Both summands (a_q60 < 2^63 and RCP < 2^63) are < 2^64, so the
     * 128-bit intermediate is fully captured by the standard mulh64
     * primitive.
     * -------------------------------------------------------------- */
    uint64_t V = mulh64(a_q60, APPROX_EXP_RCP);
    uint64_t m = ((V + ((UINT64_C(1) << 53) - 1U)) >> 53) + (APPROX_EXP_B - 1U);
    uint64_t n = m >> APPROX_EXP_BETA;
    uint64_t j = m & (APPROX_EXP_B - 1U);

    /* --------------------------------------------------------------
     * Step 2. Range-reduce to r_l in [(B-1)/B * ln2, ln2).
     *
     * Mathematically  r_l = (m/B) * ln2 - a.
     *
     * With  LB = round((ln2/B) * 2^(64+beta))  and
     * a_q60 in Q60, the inner 128-bit subtract assembles
     *
     *   R     = m * LB - (a_q60 << (beta + 4))
     *         = r_l * 2^(64 + beta)
     *
     * (the `<< (beta+4)` lifts a_q60 from Q60 to Q(64+beta) = Q69),
     * and a final  >> beta  shift drops it to Q64.
     *
     * Sizes (worst case across all SHUTTLE modes):
     *   m * LB         <  372 * 2^63.47  <  2^72.4
     *   a_q60 << 9     <  2^63 * 2^9     =  2^72
     *   R              <  2^68.5   (=> well within 128 bits)
     *   X = R >> beta  <  2^64     (Q64, fits uint64 exactly)
     * -------------------------------------------------------------- */
    {
        unsigned __int128 R = (unsigned __int128)m * APPROX_EXP_LB
                            - ((unsigned __int128)a_q60 << (APPROX_EXP_BETA + 4));
        uint64_t X = (uint64_t)(R >> APPROX_EXP_BETA);

        /* ----------------------------------------------------------
         * Step 3. Horner evaluation of P(r_l) ~ exp(r_l) in Q63.
         *
         * Each step: p_Q63 = mulh64(p_Q63, X_Q64) + C_i_Q63 (Q63).
         * Because X = r_l * 2^64 < 2^64 (X is Q64) and every
         * Horner accumulator is bounded by exp(r_l) < 2 (Q63 ~ 2^64),
         * the chain never overflows uint64.
         * ---------------------------------------------------------- */
        uint64_t p = APPROX_EXP_POLY_C6;
        p = mulh64(p, X) + APPROX_EXP_POLY_C5;
        p = mulh64(p, X) + APPROX_EXP_POLY_C4;
        p = mulh64(p, X) + APPROX_EXP_POLY_C3;
        p = mulh64(p, X) + APPROX_EXP_POLY_C2;
        p = mulh64(p, X) + APPROX_EXP_POLY_C1;
        p = mulh64(p, X) + APPROX_EXP_POLY_C0;

        /* ----------------------------------------------------------
         * Step 4. Look up s = 2^{-j/B} (Q63), constant-time.
         * ---------------------------------------------------------- */
        uint64_t s = ct_lookup_table(APPROX_EXP_TABLE_T, j, APPROX_EXP_B);

        /* ----------------------------------------------------------
         * Step 5. Combine and shift.
         *
         *   mulh64(p_Q63, s_Q63) gives p*s in Q62; shift left by one
         *   to land in Q63. The product 2^{-j/B} * exp(r_l) is exactly
         *   2^{r_l/ln2 - j/B} in [1, 2), so W lives in [2^63, 2^64).
         *
         *   The 2^{-n} factor is then applied by a round-to-nearest
         *   right shift by n: add (1 << n) >> 1 in a 128-bit container
         *   to avoid the (rare) carry past bit 63, then shift by n.
         *
         *   Result is exp(-a) * 2^63 to Q63 precision, with absolute
         *   error bounded by 2^-64 from this last quantisation step.
         * ---------------------------------------------------------- */
        uint64_t ps_q62 = mulh64(p, s);
        uint64_t W = ps_q62 << 1;
        uint64_t round_bit = (UINT64_C(1) << n) >> 1;   /* 0 when n == 0 */
        unsigned __int128 W_rounded = (unsigned __int128)W + round_bit;
        return (uint64_t)(W_rounded >> n);
    }
}
