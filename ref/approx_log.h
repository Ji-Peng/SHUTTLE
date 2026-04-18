/*
 * approx_log.h - ApproxNegLn module for SHUTTLE rejection sampling.
 *
 * Computes -ln(r_rand) in Q62 fixed-point using:
 *   1) CLZ range reduction: r_rand = 2^{-(e+1)} * m_frac, m_frac in [1, 2)
 *   2) 64-entry table subdivision: ln(1 + j/64) and 1/(1 + j/64)
 *   3) Degree-8 polynomial: Q(s) = ln(1+s)/s on [0, 1/64)
 *
 * Input:  r_rand in (0, 2^64), representing a uniform value in (0, 1)
 *         (r_rand = 0 is handled by clamping to a large output)
 * Output: ell = -ln(r_rand / 2^64) in Q62 fixed-point
 *         Range: [0, 64*ln(2)] ~ [0, 44.36], representable in Q62
 *
 * Precision: ~57 bits absolute (see docs/NGCC_Sign/ApproxLog.md)
 *
 * Method: CLZ range reduction + 64-entry table + degree-8 polynomial
 * All operations are constant-time (no data-dependent branches).
 */

#ifndef SHUTTLE_APPROX_LOG_H
#define SHUTTLE_APPROX_LOG_H

#include <stdint.h>

/*
 * approx_neg_ln - Compute -ln(r_rand / 2^64) in Q62 fixed-point.
 *
 * Input:  r_rand in (0, 2^64), representing a uniform value in (0, 1)
 *         r_rand = 0 is clamped (returns 64*ln(2) in Q62).
 * Output: ell = -ln(r_rand * 2^{-64}) in Q62 fixed-point.
 *         Guaranteed non-negative.
 */
int64_t approx_neg_ln(uint64_t r_rand);

#endif /* SHUTTLE_APPROX_LOG_H */
