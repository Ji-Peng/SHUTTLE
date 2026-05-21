/*
 * params.h - Parameters for SHUTTLE signature scheme.
 *
 * Three parameter sets, aligned to SHUTTLE-Spec/sections/Description.tex Table 2:
 *
 *   SHUTTLE_MODE=128 : n=256,  q=13313, sigma=101, full NTT (base_deg=1)
 *   SHUTTLE_MODE=256 : n=512,  q=32257, sigma=149, incomplete NTT (base_deg=2)
 *   SHUTTLE_MODE=512 : n=1024, q=64513, sigma=202, incomplete NTT (base_deg=2)
 *
 * For SHUTTLE-256 and SHUTTLE-512, q does not satisfy q == 1 (mod 2n), so
 * a "Kyber-style" incomplete NTT is used: log2(n)-1 butterfly layers leave
 * the polynomial in n/2 deg-2 basecase rings Z_q[X]/(X^2 - psi_k), and
 * pointwise multiplication is replaced by basemul (see ntt.c).
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
#define SHUTTLE_SEEDBYTES    32
#define SHUTTLE_SKSEEDBYTES  64
#define SHUTTLE_CRHBYTES     64
#define SHUTTLE_TRBYTES      64
#define SHUTTLE_RNDBYTES     32
#define SHUTTLE_CTILDEBYTES  32

/* ============================================================
 * Mode-specific parameters (Spec Table 2)
 * ============================================================ */
#if SHUTTLE_MODE == 128

#define SHUTTLE_N          256
#define SHUTTLE_Q          13313
#define SHUTTLE_QBITS      14            /* ceil(log2(13313)) = 14 */
#define SHUTTLE_BASE_DEG   1             /* full NTT */
#define SHUTTLE_L          3
#define SHUTTLE_M          2
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        30
#define SHUTTLE_ALPHA_H    128
#define SHUTTLE_ALPHA_1    8

/* Fixed-point bounds (Q2 = original * 100, then squared).
 *   B_k = 25.79  -> B_k * 100 = 2579,   (B_k*100)^2 = 6651241
 *   B_s = 3893.66 -> B_s*100 = 389366,  (B_s*100)^2 = 151605881956
 *   B_v = 4663.0 -> B_v*100 = 466300,   (B_v*100)^2 = 217435690000   */
#define SHUTTLE_BK_Q2          2579L
#define SHUTTLE_BS_Q2          389366L
#define SHUTTLE_BV_Q2          466300L
#define SHUTTLE_BK_LOW_Q2      0L

#elif SHUTTLE_MODE == 256

#define SHUTTLE_N          512
#define SHUTTLE_Q          32257
#define SHUTTLE_QBITS      15            /* ceil(log2(32257)) = 15 */
#define SHUTTLE_BASE_DEG   2             /* incomplete NTT, basecase deg 2 */
#define SHUTTLE_L          3
#define SHUTTLE_M          2
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        58
#define SHUTTLE_ALPHA_H    1024
#define SHUTTLE_ALPHA_1    16

/*   B_k = 36.26   -> 3626,   3626^2 = 13147876
 *   B_s = 8391    -> 839100, 839100^2 = 704088810000
 *   B_v = 16647   -> 1664700, 1664700^2 = 2771226090000 */
#define SHUTTLE_BK_Q2          3626L
#define SHUTTLE_BS_Q2          839100L
#define SHUTTLE_BV_Q2          1664700L
#define SHUTTLE_BK_LOW_Q2      0L

#elif SHUTTLE_MODE == 512

#define SHUTTLE_N          1024
#define SHUTTLE_Q          64513
#define SHUTTLE_QBITS      16            /* ceil(log2(64513)) = 16 */
#define SHUTTLE_BASE_DEG   2             /* incomplete NTT */
#define SHUTTLE_L          3
#define SHUTTLE_M          2
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        115
#define SHUTTLE_ALPHA_H    2048
#define SHUTTLE_ALPHA_1    16

/*   B_k = 51.08   -> 5108,    5108^2 = 26091664
 *   B_s = 31639   -> 3163900, 3163900^2 = 10010264410000
 *   B_v = 54900   -> 5490000, 5490000^2 = 30140100000000 */
#define SHUTTLE_BK_Q2          5108L
#define SHUTTLE_BS_Q2          3163900L
#define SHUTTLE_BV_Q2          5490000L
#define SHUTTLE_BK_LOW_Q2      0L

#else
#  error "Unsupported SHUTTLE_MODE (expected 128, 256, or 512)"
#endif

/* Pull in MONT / QINV / BARRETT_V / INVNTT_F / zetas table per mode. */
#include "ntt_constants.h"

/* ============================================================
 * Fixed-point bound scale.
 *
 * Bounds B_k, B_s, B_v in the spec are decimals (e.g. B_k = 25.79).
 * We store them as integers in Q2 (multiplied by 100) and compare a
 * regular int64 norm_sq against the squared Q2 value by scaling norm_sq
 * by SHUTTLE_BOUND_SCALE_SQ:
 *
 *   accept iff  SHUTTLE_BOUND_SCALE_SQ * norm_sq < SHUTTLE_*_Q2 * SHUTTLE_*_Q2
 *
 * For the worst-case mode (SHUTTLE-512):
 *   norm_sq <= B_v^2 ~ 3.0e9
 *   * 10000 = 3.0e13   (fits in int64, max ~9.2e18).
 * ============================================================ */
#define SHUTTLE_BOUND_SCALE    100L
#define SHUTTLE_BOUND_SCALE_SQ (SHUTTLE_BOUND_SCALE * SHUTTLE_BOUND_SCALE)

#define SHUTTLE_BK_SQ_FX   ((int64_t)SHUTTLE_BK_Q2 * (int64_t)SHUTTLE_BK_Q2)
#define SHUTTLE_BS_SQ_FX   ((int64_t)SHUTTLE_BS_Q2 * (int64_t)SHUTTLE_BS_Q2)
#define SHUTTLE_BV_SQ_FX   ((int64_t)SHUTTLE_BV_Q2 * (int64_t)SHUTTLE_BV_Q2)
#define SHUTTLE_BK_LOW_SQ_FX ((int64_t)SHUTTLE_BK_LOW_Q2 * (int64_t)SHUTTLE_BK_LOW_Q2)

/* Test a non-negative integer norm_sq against a fixed-point squared bound.
 * Returns 1 iff norm_sq * scale_sq < bound_sq_fx. */
#define SHUTTLE_NORM_LT_FX(norm_sq, bound_sq_fx) \
    ((int64_t)(norm_sq) * SHUTTLE_BOUND_SCALE_SQ < (bound_sq_fx))

#define SHUTTLE_NORM_GE_FX(norm_sq, bound_sq_fx) \
    ((int64_t)(norm_sq) * SHUTTLE_BOUND_SCALE_SQ >= (bound_sq_fx))

#define SHUTTLE_NORM_GT_FX(norm_sq, bound_sq_fx) \
    ((int64_t)(norm_sq) * SHUTTLE_BOUND_SCALE_SQ > (bound_sq_fx))

/* Full sk vector length = [alpha_1, s, e] */
#define SHUTTLE_VECLEN     (1 + SHUTTLE_L + SHUTTLE_M)

/* Decompose helpers.
 *   The number of valid HighBits buckets is W1_MAX+1; W1_BITS is
 *   ceil(log2(W1_MAX+1)). Since 2*alpha_h is a power of 2 the divisor
 *   becomes a shift. */
#define SHUTTLE_W1_MAX     ((SHUTTLE_Q - 1) / (2 * SHUTTLE_ALPHA_H))

#if SHUTTLE_MODE == 128
#  define SHUTTLE_W1_BITS  6        /* W1_MAX = 52, fits in 6 bits */
#elif SHUTTLE_MODE == 256
#  define SHUTTLE_W1_BITS  5        /* W1_MAX = 15 (with rounding edge 16), 5 bits */
#elif SHUTTLE_MODE == 512
#  define SHUTTLE_W1_BITS  5        /* W1_MAX = 15 (with rounding edge 16), 5 bits */
#endif

/* ============================================================
 * mod 2q infrastructure (Phase 6b)
 * ============================================================ */
#define SHUTTLE_DQ          (2 * SHUTTLE_Q)
#define SHUTTLE_HALF_ALPHA_H (SHUTTLE_ALPHA_H / 2)
#define SHUTTLE_HINT_MOD    (2 * (SHUTTLE_Q - 1))
#define SHUTTLE_HINT_MAX    (SHUTTLE_HINT_MOD / SHUTTLE_ALPHA_H)

#if SHUTTLE_ALPHA_H == 128
#  define SHUTTLE_ALPHA_H_BITS 7
#elif SHUTTLE_ALPHA_H == 256
#  define SHUTTLE_ALPHA_H_BITS 8
#elif SHUTTLE_ALPHA_H == 1024
#  define SHUTTLE_ALPHA_H_BITS 10
#elif SHUTTLE_ALPHA_H == 2048
#  define SHUTTLE_ALPHA_H_BITS 11
#else
#  error "Unsupported SHUTTLE_ALPHA_H"
#endif

#if SHUTTLE_ALPHA_1 == 8
#  define SHUTTLE_ALPHA_1_BITS 3
#elif SHUTTLE_ALPHA_1 == 16
#  define SHUTTLE_ALPHA_1_BITS 4
#else
#  error "Unsupported SHUTTLE_ALPHA_1 (expected 8 or 16)"
#endif

/* ============================================================
 * rANS hi/lo split parameters (implementation-defined; see
 *   agent/rANS/SHUTTLE_rANS.tex Tab 1).
 *
 *   alpha_r  : the rANS hi/lo split step for z^(1..lenS); a power of two,
 *              fixed across modes. The HighBits of z^(i) live in the
 *              z-hi rANS vocabulary at scale r/alpha_r; the LowBits are
 *              uniform in [-alpha_r/2, alpha_r/2) and bit-packed at
 *              log2(alpha_r) bits/coef.
 *   alpha_0' : the second-level split step for z^(0). Chosen as
 *              alpha_r / alpha_1 so that HighBits_{alpha_0'}(z^(0))
 *              shares the same effective scale r/alpha_r as the rest
 *              of the z-hi stream (see tex §2.5 and Cor 3). LowBits
 *              are bit-packed at log2(alpha_0') bits/coef.
 *
 *   alpha_h  : NOT used by rANS; it is the spec-level HighBits/MakeHint
 *              parameter, separate from alpha_r / alpha_0'.
 * ============================================================ */
#define SHUTTLE_ALPHA_R       64
#define SHUTTLE_ALPHA_R_BITS  6
#define SHUTTLE_HALF_ALPHA_R  (SHUTTLE_ALPHA_R / 2)

#define SHUTTLE_ALPHA_0P      (SHUTTLE_ALPHA_R / SHUTTLE_ALPHA_1)
#if   SHUTTLE_ALPHA_0P == 8
#  define SHUTTLE_ALPHA_0P_BITS 3
#elif SHUTTLE_ALPHA_0P == 4
#  define SHUTTLE_ALPHA_0P_BITS 2
#else
#  error "Unsupported SHUTTLE_ALPHA_0P (expected 4 or 8)"
#endif
#define SHUTTLE_HALF_ALPHA_0P (SHUTTLE_ALPHA_0P / 2)

/* ============================================================
 * Derived packing sizes (all in bytes)
 * ============================================================ */
/* eta=1: coefficients in {-1,0,1}, encode as 2 bits/coeff */
#define SHUTTLE_POLYETA_PACKEDBYTES   (SHUTTLE_N * 2 / 8)

/* Public key b: QBITS bits/coeff (unsigned range [0, q-1]) */
#define SHUTTLE_POLYPK_PACKEDBYTES    ((SHUTTLE_N * SHUTTLE_QBITS) / 8)

/* z[0] after CompressY: coefficients in [-q/(2*alpha_1), q/(2*alpha_1)].
 * For each mode:
 *   mode-128:  q/(2*8)  = 832,  11 signed bits
 *   mode-256:  q/(2*16) = 1008, 11 signed bits
 *   mode-512:  q/(2*16) = 2016, 12 signed bits */
#if SHUTTLE_MODE == 512
#  define SHUTTLE_Z0_BITS  12
#else
#  define SHUTTLE_Z0_BITS  11
#endif
#define SHUTTLE_POLYZ0_PACKEDBYTES    ((SHUTTLE_N * SHUTTLE_Z0_BITS + 7) / 8)

/* z[1..L]: full-range signed coefficients. Z_BOUND <= 11*sigma + alpha_1*tau.
 *   mode-128:  1351, 2*Z = 2702 -> 12 bits
 *   mode-256:  2567, 2*Z = 5134 -> 13 bits
 *   mode-512:  4062, 2*Z = 8124 -> 14 bits
 * Use 14 bits for all modes to keep one code path. */
#define SHUTTLE_POLYZ_BITS  14
#define SHUTTLE_POLYZ_PACKEDBYTES     ((SHUTTLE_N * SHUTTLE_POLYZ_BITS) / 8)

/* High-bits (w1) packing */
#define SHUTTLE_POLYW1_PACKEDBYTES    ((SHUTTLE_N * SHUTTLE_W1_BITS + 7) / 8)

/* Hint encoding (sparse-index list) -- placeholder; see polyveck_hint_pack_basic */
#define SHUTTLE_POLYVECH_PACKEDBYTES  (SHUTTLE_M * (SHUTTLE_N / 8 + 1))

/* IRS sign bits */
#define SHUTTLE_IRS_SIGNBYTES         ((SHUTTLE_TAU + 7) / 8)

/* ============================================================
 * rANS reservation budgets (consumed by packing.c).
 *
 * Analytically derived from agent/rANS/SHUTTLE_rANS.tex Tab 8
 * (tab:reserve-2stream): per-stream overflow probability p_k = 2^-21,
 * total p_rans^* = 2^-20 by union bound. z-hi stream uses Gauss CLT;
 * hint stream uses Gauss CLT for mode-128 (lambda >> 30) and the
 * Poisson model for mode-256/512 (narrow sigma).
 *
 * Regenerate with: python3 tools/SigSize.py
 *
 * OOV failure does not exist: the rANS vocabulary covers the full
 * |sym| <= M_voc tight bound derived from the 11*sigma truncation
 * of SampleY, so the only rejection cause is the rare overflow event.
 * ============================================================ */
#if SHUTTLE_MODE == 128
/* mode-128 plumbing: sigma_zhi = sigma_hint = r/alpha_r = 101/64 ~= 1.58,
 * so z-hi and hint share the same frequency table but stay in two streams
 * (kept independent for decoder simplicity). */
#  define SHUTTLE_ZHI_RESERVED_BYTES    375
#  define SHUTTLE_HINT_RESERVED_BYTES   194
#elif SHUTTLE_MODE == 256
#  define SHUTTLE_ZHI_RESERVED_BYTES    883
#  define SHUTTLE_HINT_RESERVED_BYTES    28
#elif SHUTTLE_MODE == 512
#  define SHUTTLE_ZHI_RESERVED_BYTES   1987
#  define SHUTTLE_HINT_RESERVED_BYTES     9
#endif

#define SHUTTLE_ZHI_BLOCK_BYTES     (2 + SHUTTLE_ZHI_RESERVED_BYTES)
#define SHUTTLE_HINT_BLOCK_BYTES    (2 + SHUTTLE_HINT_RESERVED_BYTES)

/* Packed size of the LowBits part of z^(0): ALPHA_0P_BITS per coef. */
#define SHUTTLE_POLYZ0_LO_PACKEDBYTES ((SHUTTLE_N * SHUTTLE_ALPHA_0P_BITS + 7) / 8)

/* Packed size of the LowBits part of z[1..lenS]: ALPHA_R_BITS per coef. */
#define SHUTTLE_POLYZ1_LO_PACKEDBYTES ((SHUTTLE_N * SHUTTLE_ALPHA_R_BITS + 7) / 8)

/* ============================================================
 * Public / secret key + signature sizes
 * ============================================================ */
#define SHUTTLE_PUBLICKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYPK_PACKEDBYTES)

#define SHUTTLE_SECRETKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_TRBYTES \
                                   + SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_L * SHUTTLE_POLYETA_PACKEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYETA_PACKEDBYTES)

/* Signature layout (two rANS streams, see SHUTTLE_rANS.tex §3.3):
 *   seedC || irs_signs
 *   || uint16 zhi_rans_len || rANS(z-hi) + pad to ZHI_RESERVED
 *   || polyz0_lo_pack(lo(z^(0)))
 *   || L * polyz1_lo_pack(lo(z^(1..lenS)))
 *   || uint16 hint_rans_len || rANS(hint) + pad to HINT_RESERVED
 *
 * The z-hi rANS block carries n*(lenS+1) coefficients: the HighBits of
 * z^(0) at step alpha_0' followed by the HighBits of z^(1..lenS) at
 * step alpha_r. All share the same effective scale r/alpha_r, hence the
 * same frequency table (tex §2.5). */
#define SHUTTLE_BYTES  ( SHUTTLE_CTILDEBYTES \
                       + SHUTTLE_IRS_SIGNBYTES \
                       + SHUTTLE_ZHI_BLOCK_BYTES \
                       + SHUTTLE_POLYZ0_LO_PACKEDBYTES \
                       + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES \
                       + SHUTTLE_HINT_BLOCK_BYTES )

/* ============================================================
 * Gaussian sampler parameters (driven by SHUTTLE_SIGMA, set in config.h)
 *
 *   sigma | small sigma | k=sigma/small | RCDT bits | RCDT entries
 *   ------+------------+---------------+-----------+--------------
 *    128  |   2.0000    |  64           |   93      |   22
 *    101  |   1.5781    |  64           |   93      |   17  (SHUTTLE-128)
 *    149  |   2.3281    |  64           |   93      |   26  (SHUTTLE-256)
 *    202  |   1.5781    | 128           |   96      |   18  (SHUTTLE-512)
 *
 * SHUTTLE-512 doubles k from 64 to 128 and uses 96-bit RCDT (3x32-bit
 * limbs) instead of 93-bit (3x31-bit limbs), per tools/BaseSampler.ipynb.
 * ============================================================ */

#if SHUTTLE_SIGMA == 128
#  define RCDT_ENTRIES 22
#  define SHUTTLE_RCDT_LIMB_BITS 31
#  define SHUTTLE_RCDT_BITS 93
#elif SHUTTLE_SIGMA == 101
#  define RCDT_ENTRIES 17
#  define SHUTTLE_RCDT_LIMB_BITS 31
#  define SHUTTLE_RCDT_BITS 93
#elif SHUTTLE_SIGMA == 149
#  define RCDT_ENTRIES 26
#  define SHUTTLE_RCDT_LIMB_BITS 31
#  define SHUTTLE_RCDT_BITS 93
#elif SHUTTLE_SIGMA == 202
#  define RCDT_ENTRIES 18
#  define SHUTTLE_RCDT_LIMB_BITS 32
#  define SHUTTLE_RCDT_BITS 96
#else
#  error "Unsupported SHUTTLE_SIGMA (expected 101, 149, 202, or legacy 128)"
#endif

#if SHUTTLE_SIGMA == 202
#  define SHUTTLE_GAUSS_K      128
#  define SHUTTLE_K_BITS       7
#  define SHUTTLE_Y_BITS       7
#  define SHUTTLE_TWO_K_BITS   8
#else
#  define SHUTTLE_GAUSS_K      64
#  define SHUTTLE_K_BITS       6
#  define SHUTTLE_Y_BITS       6
#  define SHUTTLE_TWO_K_BITS   7
#endif

#define SHUTTLE_TRUNC        11
#define SHUTTLE_BOUND        (SHUTTLE_TRUNC * SHUTTLE_SIGMA)

/* Packing bound for z coefficients (post alpha_1 stretch).
 *   |z[0]|_inf <= 11*sigma + alpha_1*tau
 *   |z[i]|_inf <= 11*sigma + tau    for i>=1
 * Use the larger of the two as a uniform bound. */
#define SHUTTLE_Z_BOUND      (SHUTTLE_BOUND + SHUTTLE_ALPHA_1 * SHUTTLE_TAU)

/* Q64 reciprocal of 2*sigma^2 was previously used by sampler.c to build
 * the approx_exp input a_q60 via mulh64(num<<44, INV_Q64). That path's
 * error chain bottomed out at 2^-46 and dominated the total budget; it
 * has been replaced by the Q80 reciprocal (APPROX_EXP_I_HI/_LO in
 * approx_exp_constants.h) coupled with a 128-bit multiply, which drops
 * the chain to 2^-59. The legacy macro is intentionally not kept here. */

/* Reciprocal of sigma^2 in Q62 (used by rejsample.c). */
#define SHUTTLE_INV_SIGMA2_Q62 \
    ((uint64_t)(((uint64_t)1 << 62) / ((uint64_t)SHUTTLE_SIGMA * (uint64_t)SHUTTLE_SIGMA)))

#endif /* SHUTTLE_PARAMS_H */
