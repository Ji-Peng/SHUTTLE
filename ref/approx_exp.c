/*
 * approx_exp.c - ApproxExp V2 implementation for SHUTTLE.
 *
 * V2 variant: [0, ln2/2) sub-interval, degree 9, 2-entry table.
 * 9 mulh64 (Horner) + 1 mulh64 (combine) + 2 ct_lookup entries.
 * Precision: ~55.9 bits.
 *
 * Algorithm (exp(+) approach):
 *   Given x >= 0 in Q60, compute exp(-x) * 2^63.
 *   1) Transform: u = N*ln2 - x > 0
 *   2) Decompose: u = nh*ln2 + rh*(ln2/2) + rl
 *   3) Polynomial: exp(+rl) via degree-9 unsigned Horner
 *   4) Table: 2^{rh/2} from 2-entry table (constant-time lookup)
 *   5) Combine: exp(-x) = 2^{nh-N} * 2^{rh/2} * exp(rl)
 */

#include "approx_exp.h"

/* ============================================================
 * Constants (Q60 range reduction)
 * ============================================================ */

#define LN2_Q60         UINT64_C(799144290325165979)
#define LN2_Q63         UINT64_C(6393154322601327830)
#define LN2_2_Q60       UINT64_C(399572145162582990)   /* ln2/2 * 2^60 */
#define N_LN2_Q60       UINT64_C(7192298612926493811)  /* 9 * LN2_Q60 */
#define RECIP_LN2_Q63   UINT64_C(13306513097844322492) /* round(1/ln2 * 2^63) */
#define RECIP_LN2_2_Q62 UINT64_C(13306513097844322492) /* round(2/ln2 * 2^62) */
#define N_SHIFT          APPROX_EXP_N_SHIFT  /* 9 */

/* ============================================================
 * Precomputed table: 2^{+i/2} in Q63, i in {0, 1}
 * ============================================================ */

static const uint64_t TABLE_2[2] = {
    UINT64_C( 9223372036854775808),  /* 2^{0/2} = 1.000000 */
    UINT64_C(13043817825332782212),  /* 2^{1/2} = sqrt(2) = 1.414214 */
};

/* ============================================================
 * V2 polynomial coefficients: exp(+x) on [0, ln2/2), degree 9, Q63
 * ============================================================ */

#define V2_C9   UINT64_C(30246879334790)
#define V2_C8   UINT64_C(225090733290635)
#define V2_C7   UINT64_C(1831549596915338)
#define V2_C6   ((uint64_t)(UINT64_C(6404931195738088) << 1))
#define V2_C5   ((uint64_t)(UINT64_C(4803843210556027) << 4))
#define V2_C4   ((uint64_t)(UINT64_C(6004799419128215) << 6))
#define V2_C3   ((uint64_t)(UINT64_C(6004799504283434) << 8))
#define V2_C2   ((uint64_t)(UINT64_C(9007199254725769) << 9))
#define V2_C1   ((uint64_t)(UINT64_C(4503599627370536) << 11))

/* ============================================================
 * Constant-time primitives
 * ============================================================ */

/* Constant-time unsigned greater-or-equal: returns 1 if a >= b */
static inline uint64_t ct_ge(uint64_t a, uint64_t b) {
    return 1 - ((a - b) >> 63);
}

/* Constant-time 2-entry table lookup (full scan for side-channel resistance) */
static inline uint64_t ct_lookup2(const uint64_t table[2], uint64_t idx) {
    uint64_t result = 0;
    for (uint64_t i = 0; i < 2; i++) {
        uint64_t diff = i ^ idx;
        uint64_t mask = ((diff | (~diff + 1)) >> 63) - 1;
        result |= table[i] & mask;
    }
    return result;
}

/* ============================================================
 * Implementation
 * ============================================================ */

uint64_t approx_exp(uint64_t x_q60) {
    /* Step 1: Transform to positive domain. u = N*ln2 - x */
    uint64_t u_q60 = N_LN2_Q60 - x_q60;

    /* Step 2: nh = floor(u / ln2) via multiply-by-reciprocal */
    uint64_t nh = mulh64(u_q60, RECIP_LN2_Q63) >> 59;

    /* Remainder r = u - nh*ln2, with overshoot correction */
    uint64_t r_q60 = u_q60 - nh * LN2_Q60;
    uint64_t overshoot = r_q60 >> 63;
    nh -= overshoot;
    r_q60 += overshoot * LN2_Q60;

    /* Standard adjustment: if r >= ln2, increment nh */
    uint64_t adj = ct_ge(r_q60, LN2_Q60);
    nh += adj;
    r_q60 -= adj * LN2_Q60;

    /* Step 3: Sub-interval decomposition.
     * rh = floor(2*r/ln2), rl = r - rh*(ln2/2) */
    uint64_t rh = mulh64(r_q60, RECIP_LN2_2_Q62) >> 58;
    uint64_t rl_q60 = r_q60 - rh * LN2_2_Q60;
    uint64_t adj2 = ct_ge(rl_q60, LN2_2_Q60);
    rh += adj2;
    rl_q60 -= adj2 * LN2_2_Q60;

    /* Step 4: Convert rl to Q63 and evaluate degree-9 Horner polynomial */
    uint64_t rl_q63 = rl_q60 << 3;
    uint64_t X = rl_q63 << 1;  /* Q64 */

    uint64_t acc = V2_C9;
    acc = mulh64(acc, X) + V2_C8;
    acc = mulh64(acc, X) + V2_C7;
    acc = mulh64(acc, X) + V2_C6;
    acc = mulh64(acc, X) + V2_C5;
    acc = mulh64(acc, X) + V2_C4;
    acc = mulh64(acc, X) + V2_C3;
    acc = mulh64(acc, X) + V2_C2;
    acc = mulh64(acc, X) + V2_C1;
    acc = mulh64(acc, X) + (UINT64_C(1) << 63);  /* C0 = 1.0 in Q63 */

    /* Step 5: Table lookup 2^{rh/2} (constant-time) */
    uint64_t table_val = ct_lookup2(TABLE_2, rh);

    /* Step 6: Combine exp(rl) * 2^{rh/2}.
     * mulh64(Q63, Q63) -> Q62. Saturate on overflow. */
    uint64_t val_q62 = mulh64(acc, table_val);
    uint64_t ovf_mask = -(val_q62 >> 63);
    val_q62 = (val_q62 & ~ovf_mask) | (UINT64_C(0x7FFFFFFFFFFFFFFF) & ovf_mask);

    /* Step 7: Apply 2^{nh - N_SHIFT} via right-shift.
     * result_q63 = (val_q62 << 1) >> (N_SHIFT - nh) */
    uint64_t shifted = val_q62 << 1;
    uint64_t shift = N_SHIFT - nh;
    uint64_t shift_mask = ((uint64_t)(63 - shift) >> 63) - 1;
    return (shifted >> (shift & 63)) & shift_mask;
}
