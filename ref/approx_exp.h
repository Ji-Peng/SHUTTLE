/*
 * approx_exp.h - ApproxExp module for SHUTTLE discrete Gaussian sampler.
 *
 * Computes exp(-x) in Q63 fixed-point using the exp(+) approach (V2 variant):
 *   1) Transform: u = N*ln2 - x, so exp(-x) = exp(u) * 2^{-N}
 *   2) Decompose: u = nh*ln2 + rh*(ln2/2) + rl
 *   3) Precomputed table: 2^{+rh/2} for rh in {0,1}
 *   4) Degree-9 polynomial: exp(+rl) for rl in [0, ln2/2)
 *   5) Combine: exp(-x) = 2^{nh-N} * 2^{rh/2} * exp(rl)
 *
 * V2 variant: [0, ln2/2) interval, degree 9, 2-entry table.
 * Precision: ~55.9 bits.
 *
 * Input range: x in [0, N_SHIFT * ln2) = [0, ~7.62).
 * Sufficient for all SHUTTLE sampler modes (max ~6.914 at sigma=101).
 */

#ifndef SHUTTLE_APPROX_EXP_H
#define SHUTTLE_APPROX_EXP_H

#include <stdint.h>

/* ============================================================
 * 64-bit multiply-high primitives (platform-specific)
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
#error "Unsupported compiler: need 128-bit multiply or intrinsics for 64-bit mulh"
#endif

/*
 * N_SHIFT: the integer offset for the exp(+) transformation.
 * exp(-x) = exp(N_SHIFT*ln2 - x) * 2^{-N_SHIFT}
 * Supports x in [0, N_SHIFT * ln2) = [0, ~7.62).
 *
 * Sized for the largest worst-case rejection exponent across all three
 * SHUTTLE modes:
 *   sigma=128: max 63*(63+128*22)/32768 ~ 5.535   (fits in N=9)
 *   sigma=149: max 63*(63+128*26)/44402 ~ 4.812   (fits in N=9)
 *   sigma=101: max 63*(63+128*17)/20402 ~ 6.914   (needs N>=10; we pick 11)
 *
 * Bumping N from 9 to 11 is mathematically invariant for any previously
 * supported input: nh = floor((N*ln2 - x)/ln2) shifts by exactly +2 when N
 * increases by +2, so shift = N - nh is unchanged and so is the returned
 * value. Existing KATs are therefore preserved.
 */
#define APPROX_EXP_N_SHIFT 11

/*
 * Compute exp(-x) in Q63 fixed-point (V2: degree 9, table 2).
 *
 * Input:  x_q60 = round(x * 2^60), where x in [0, N_SHIFT * ln2).
 * Output: round(exp(-x) * 2^63) as uint64_t.
 *         Returns 0 if exp(-x) underflows to 0.
 *
 * All operations are constant-time (no data-dependent branches/loops).
 */
uint64_t approx_exp(uint64_t x_q60);

#endif /* SHUTTLE_APPROX_EXP_H */
