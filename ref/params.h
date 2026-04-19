/*
 * params.h - Parameters for SHUTTLE signature scheme.
 *
 * Ring: R_q = Z_q[x]/(x^n + 1), n and q selected per SHUTTLE_MODE.
 * Module: secret key [alpha_1, s, e] with s in R^l, e in R^m.
 *
 * Three parameter sets are supported (see NGCC-Signature Table 2):
 *   SHUTTLE_MODE=128 : n=256, q=13313, sigma=101
 *   SHUTTLE_MODE=256 : n=512, q=13313, sigma=149
 *   SHUTTLE_MODE=512 : parameters unspecified; #error in this file
 *
 * SHUTTLE_MODE is selected via config.h.
 */

#ifndef SHUTTLE_PARAMS_H
#define SHUTTLE_PARAMS_H

#include <stdint.h>
#include "config.h"

/* ============================================================
 * Shared constants (hash / seed sizes, common across all modes)
 * ============================================================ */
#define SHUTTLE_Q            13313   /* Prime modulus, shared across all modes. q = 13 * 1024 + 1 */
#define SHUTTLE_QBITS        14      /* ceil(log2(q)) = 14 */

/* Seed / hash byte lengths. Currently fixed at the 32-byte "seedBytes"
 * baseline used by the NGCC-Signature draft: pending spec confirmation of
 * whether these should later scale with lambda.
 *
 * ExpandSeeds in KeyGen (Alg, line 115) takes xi || I2B(lenE, 1) and
 * produces 4 * seedBytes bytes, split as
 *   seedA      : seedBytes            (SHUTTLE_SEEDBYTES    = 32)
 *   seedsk     : 2 * seedBytes        (SHUTTLE_SKSEEDBYTES  = 64)
 *   masterSeed : seedBytes            (SHUTTLE_SEEDBYTES    = 32)
 */
#define SHUTTLE_SEEDBYTES    32
#define SHUTTLE_SKSEEDBYTES  64
#define SHUTTLE_CRHBYTES     64
#define SHUTTLE_TRBYTES      64
#define SHUTTLE_RNDBYTES     32
#define SHUTTLE_CTILDEBYTES  32

/* ============================================================
 * Mode-specific parameters (Table 2 of NGCC-Signature)
 * ============================================================ */
#if SHUTTLE_MODE == 128

#define SHUTTLE_N          256   /* Ring dimension */
#define SHUTTLE_L          3     /* # components in s */
#define SHUTTLE_M          2     /* # components in e */
#define SHUTTLE_ETA        1     /* CBD parameter for s and e */
#define SHUTTLE_TAU        30    /* Hamming weight of challenge c */
#define SHUTTLE_ALPHA_H    128   /* Hint compression parameter */
#define SHUTTLE_ALPHA_1    8     /* Compression parameter for first response slot */
#define SHUTTLE_W1_BITS    6     /* bits/coeff for high-part packing (floor((q-1)/(2*alpha_h)) = 52 fits in 6) */

/* Norm bounds. Spec: B_k = 25.79, B_s = 3893.66, B_v = 4663.0.
 * Stored as integer squares (floor) for cheap comparisons. */
#define SHUTTLE_BK         26            /* ceil(B_k) */
#define SHUTTLE_BK_LOW     0
#define SHUTTLE_BS_SQ      15160508UL    /* floor(3893.66^2) */
#define SHUTTLE_BV_SQ      21743569UL    /* floor(4663.0^2) */

#elif SHUTTLE_MODE == 256

#define SHUTTLE_N          512
#define SHUTTLE_L          3
#define SHUTTLE_M          2
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        58
#define SHUTTLE_ALPHA_H    256
#define SHUTTLE_ALPHA_1    16
#define SHUTTLE_W1_BITS    5   /* floor((q-1)/(2*alpha_h)) = 26 fits in 5 */

#define SHUTTLE_BK         37           /* ceil(36.26) */
#define SHUTTLE_BK_LOW     0
#define SHUTTLE_BS_SQ      82628100UL   /* 9090^2 */
#define SHUTTLE_BV_SQ      125484804UL  /* 11202^2 */

#elif SHUTTLE_MODE == 512

#error "SHUTTLE-512 parameters not yet specified; see NGCC-Signature Table 2"

#endif

/* Full sk vector length = [alpha_1, s, e] */
#define SHUTTLE_VECLEN     (1 + SHUTTLE_L + SHUTTLE_M)

/* Decompose helpers.
 *   2 * alpha_h is a power of two across all modes (256 or 512), so
 *   division by it compiles to a shift.
 *   W1_MAX = floor((q-1) / (2*alpha_h)) is the highest valid high-bit value;
 *   valid bucket count is W1_MAX + 1.
 *   For q=13313:  mode-128 -> W1_MAX=52 (53 buckets, 6 bits);
 *                 mode-256 -> W1_MAX=26 (27 buckets, 5 bits). */
#define SHUTTLE_W1_MAX     ((SHUTTLE_Q - 1) / (2 * SHUTTLE_ALPHA_H))

/* ============================================================
 * Derived packing sizes (all in bytes)
 * ============================================================ */
/* eta=1: coefficients in {-1,0,1}, encode as 2 bits/coeff */
#define SHUTTLE_POLYETA_PACKEDBYTES   (SHUTTLE_N * 2 / 8)

/* Public key b: 14 bits/coeff (unsigned range [0, q-1]) */
#define SHUTTLE_POLYPK_PACKEDBYTES    (SHUTTLE_N * SHUTTLE_QBITS / 8)

/* z[0] after CompressY: coefficients in [-q/(2*alpha_1), q/(2*alpha_1)],
 * packed as signed integers at Z0_BITS bits/coeff.
 * For alpha_1=8, q/(2*8) ~ 832, so 11 signed bits.
 * For alpha_1=16, q/(2*16) ~ 416, so 10 signed bits.
 * Fixed at 11 across modes to keep the pack routine simple. */
#define SHUTTLE_Z0_BITS               11
#define SHUTTLE_POLYZ0_PACKEDBYTES    (SHUTTLE_N * SHUTTLE_Z0_BITS / 8)

/* z[1..L]: full-range signed coefficients, 14 bits/coeff */
#define SHUTTLE_POLYZ_PACKEDBYTES     (SHUTTLE_N * SHUTTLE_QBITS / 8)

/* High-bits (w1) packing, 5 or 6 bits/coeff per mode */
#define SHUTTLE_POLYW1_PACKEDBYTES    (SHUTTLE_N * SHUTTLE_W1_BITS / 8)

/* Hint encoding (sparse-index list). M polys; up to ~OMEGA=alpha_h/2 hint bits
 * per poly plus M index-offset counters. Exact bound depends on MakeHint analysis.
 * Conservative upper bound reserved here; final size must be asserted at
 * compile time by packing.c once hint encoding is implemented. */
#define SHUTTLE_POLYVECH_PACKEDBYTES  (SHUTTLE_M * (SHUTTLE_N / 8 + 1))

/* IRS sign bits: ceil(TAU/8) */
#define SHUTTLE_IRS_SIGNBYTES         ((SHUTTLE_TAU + 7) / 8)

/* Public / secret key sizes */
#define SHUTTLE_PUBLICKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYPK_PACKEDBYTES)

#define SHUTTLE_SECRETKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_TRBYTES \
                                   + SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_L * SHUTTLE_POLYETA_PACKEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYETA_PACKEDBYTES)

/* Signature layout.
 *
 * CURRENT (uncompressed, what packing.c writes today):
 *   c_tilde (32B)
 *   irs_signs (ceil(TAU/8) B)
 *   polyz_pack(z[0..VECLEN-1])               -- VECLEN polys, 14 bits/coeff
 *
 * TARGET (compressed per NGCC-Signature Table 2, to be implemented in phase 6):
 *   c_tilde (32B)
 *   irs_signs (ceil(TAU/8) B)
 *   polyz0_pack(z[0])                        -- 11 bits/coeff, N coeffs
 *   polyz_pack(z[1..L])                      -- L polys, 14 bits/coeff
 *   polyvech_pack(h)                         -- hint over the M "z2" polys
 *   -> SHUTTLE-128: 1560 B, SHUTTLE-256: 2311 B (spec targets)
 *
 * SHUTTLE_BYTES currently reflects the uncompressed layout. Phase 6 will
 * switch both pack_sig and SHUTTLE_BYTES to the compressed form atomically. */
#define SHUTTLE_BYTES  (SHUTTLE_CTILDEBYTES \
                          + SHUTTLE_IRS_SIGNBYTES \
                          + SHUTTLE_VECLEN * SHUTTLE_POLYZ_PACKEDBYTES)

/* ============================================================
 * Gaussian sampler parameters (driven by SHUTTLE_SIGMA, set in config.h)
 *
 *   sigma | small sigma | RCDT | 2*sigma^2 |   Renyi (1025)  | parameter set
 *   ------+------------+------+-----------+-----------------+---------------
 *    128  |    2.0000  |  22  |  32768    | 1 + 2^-93.85    | legacy / reference (not tied to any mode)
 *    101  |    1.5781  |  17  |  20402    | 1 + 2^-93.40    | SHUTTLE-128
 *    149  |    2.3281  |  26  |  44402    | 1 + 2^-93.75    | SHUTTLE-256
 *
 * All three modes share k = 64 and truncation = 11*sigma.
 * ============================================================ */

#if SHUTTLE_SIGMA == 128
#  define RCDT_ENTRIES 22
#elif SHUTTLE_SIGMA == 101
#  define RCDT_ENTRIES 17
#elif SHUTTLE_SIGMA == 149
#  define RCDT_ENTRIES 26
#else
#  error "Unsupported SHUTTLE_SIGMA (expected 128, 101, or 149)"
#endif

#define SHUTTLE_GAUSS_K      64
#define SHUTTLE_K_BITS       6
#define SHUTTLE_Y_BITS       6
#define SHUTTLE_TWO_K_BITS   7
#define SHUTTLE_TRUNC        11
#define SHUTTLE_BOUND        (SHUTTLE_TRUNC * SHUTTLE_SIGMA)

/* Packing bound for z coefficients after the alpha_1 stretch.
 *   z[i] = y[i] + c_eff * sk_full[i]
 * with sk_full[0] = alpha_1 (constant poly), so |z[0]|_inf <= 11*sigma + alpha_1*tau.
 * For i >= 1, sk_full[i] has eta-bounded coefficients so |z[i]|_inf <= 11*sigma + tau.
 * Use the larger of the two as a uniform bound so polyz_pack/unpack stays
 * a single code path over all six slots.
 * Both modes satisfy 2 * SHUTTLE_Z_BOUND < 2^14 so 14-bit packing still fits. */
#define SHUTTLE_Z_BOUND      (SHUTTLE_BOUND + SHUTTLE_ALPHA_1 * SHUTTLE_TAU)

/* Reciprocal of 2*sigma^2 in Q64. See rejsample.c. */
#define SHUTTLE_INV_2SIGMA2_Q64 \
    ((uint64_t)(-1ULL / (2ULL * (uint64_t)SHUTTLE_SIGMA * (uint64_t)SHUTTLE_SIGMA)))

/* Reciprocal of sigma^2 in Q62: floor(2^62 / sigma^2), fits uint64_t. Used
 * by rejsample.c to convert an integer inner-product into the Q62 quantity
 *   u_q62 = (inner / sigma^2) * 2^62 ~= inner * SHUTTLE_INV_SIGMA2_Q62
 * without assuming sigma^2 is a power of two. */
#define SHUTTLE_INV_SIGMA2_Q62 \
    ((uint64_t)(((uint64_t)1 << 62) / ((uint64_t)SHUTTLE_SIGMA * (uint64_t)SHUTTLE_SIGMA)))

#endif /* SHUTTLE_PARAMS_H */
