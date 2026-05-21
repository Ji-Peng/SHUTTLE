/*
 * approx_exp.h - Fixed-point exp(-a) primitive for the SHUTTLE rejection
 *                sampler.
 *
 * Public interface to the implementation in approx_exp.c, which realises
 * Algorithm 1 of agent/ApproxExp/ApproxExp.tex. The function approx_exp()
 * computes round(exp(-a) * 2^63) given a in Q60 fixed-point.
 *
 * Caller responsibility: the input must satisfy 0 <= a < 11 * ln(2) ~
 * 7.624. Every SHUTTLE parameter set's worst-case a_max stays below this
 * bound (SHUTTLE-128: 6.914, SHUTTLE-256: 4.812, SHUTTLE-512: 7.369).
 *
 * Output precision: 53-bit relative (a-independent terms ~ 2^-58 plus
 * an output-quantisation term 2^-64 / exp(-a) that dominates near
 * a = a_max). See approx_exp.c and ApproxExp.tex Section 4 for the full
 * error budget.
 *
 * Side-channel: no data-dependent control flow; the B-entry table is
 * always scanned end-to-end with bit masks.
 */

#ifndef SHUTTLE_APPROX_EXP_H
#define SHUTTLE_APPROX_EXP_H

#include <stdint.h>

/* ============================================================
 * 64-bit multiply-high primitive (mulh64).
 *
 * Returns the upper 64 bits of the unsigned 128-bit product a*b.
 * Used pervasively inside approx_exp (Horner chain, range reduction,
 * Q57 reciprocal) and inside the sampler's a_q60 computation.
 *
 * On x86-64 GCC / Clang the unsigned __int128 form compiles to a
 * single `mulq` instruction (~3 cycles).
 * ============================================================ */
#if defined(__GNUC__) || defined(__clang__)

__extension__ typedef unsigned __int128 wide_uint128;
__extension__ typedef __int128 wide_int128;

static inline uint64_t mulh64(uint64_t a, uint64_t b) {
    return (uint64_t)(((wide_uint128)a * b) >> 64);
}

static inline int64_t smulh64(int64_t a, int64_t b) {
    return (int64_t)(((wide_int128)a * b) >> 64);
}

#elif defined(_MSC_VER)

#include <intrin.h>

static inline uint64_t mulh64(uint64_t a, uint64_t b) {
    return __umulh(a, b);
}

static inline int64_t smulh64(int64_t a, int64_t b) {
    return __mulh(a, b);
}

#else
#  error "Unsupported compiler: need 128-bit multiply or intrinsics for 64-bit mulh"
#endif

/* ============================================================
 * approx_exp -- core primitive.
 *
 *   Input:  a_q60 = round(a * 2^60), with a in [0, 11 * ln 2).
 *   Output: round(exp(-a) * 2^63) as uint64_t, in [0, 2^63].
 *
 * Constant-time with respect to a_q60 and the implied internal index
 * j = m mod B (table is scanned in full with masks). All operations are
 * 64- or 128-bit unsigned integer arithmetic; no division, no branches.
 * ============================================================ */
uint64_t approx_exp(uint64_t a_q60);

#endif /* SHUTTLE_APPROX_EXP_H */
