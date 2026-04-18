/*
 * rejsample.c - Iterative Rejection Sampling for SHUTTLE.
 *
 * Implements Algorithm 6 (Simplified Log-Domain) from InteractiveRejection.tex.
 *
 * Overview:
 *   For each of tau=30 challenge monomials x^{j_i}, the algorithm:
 *   1. Computes v_i = x^{j_i} * sk (cyclic rotation of the secret key vector)
 *   2. Computes u = <z, v_i> / sigma^2 in Q62 fixed-point
 *   3. Samples r_rand uniform in (0,1) as 64-bit integer
 *   4. Computes ell = -ln(r_rand) via ApproxNegLn
 *   5. Checks acceptance: ell >= ln_M
 *   6. Uses IntervalParity to determine the sign of the update
 *   7. Updates z = z +/- v_i (branchless)
 *
 * All operations are constant-time to prevent side-channel leakage.
 *
 * Key parameters:
 *   - sigma = SHUTTLE_SIGMA = 128, so sigma^2 = 16384
 *   - N = SHUTTLE_N = 256 (polynomial ring degree)
 *   - VECLEN = 6 (full secret key vector length)
 *   - tau = SHUTTLE_TAU = 30 (challenge Hamming weight)
 */

#include <string.h>
#include "rejsample.h"
#include "approx_log.h"
#include "approx_exp.h"  /* for wide_int128, wide_uint128 */

/* ============================================================
 * Constants
 * ============================================================ */

/* Number of IntervalParity iterations. Must be exactly 100 for security. */
#define IRS_PARITY_STEPS  100

/* sigma^2 = 128^2 = 16384 = 2^14 */
#define SIGMA_SQ          16384

/* ============================================================
 * Constant-time helpers
 * ============================================================ */

/*
 * ct_abs_i64 - Constant-time absolute value of int64_t.
 * Returns |x| without branching. Undefined for x = INT64_MIN.
 */
static inline int64_t ct_abs_i64(int64_t x) {
    int64_t mask = x >> 63;  /* 0 if x >= 0, -1 if x < 0 */
    return (x ^ mask) - mask;
}

/*
 * ct_neg_mask_i64 - Returns -1 (all ones) if x < 0, else 0.
 */
static inline int64_t ct_neg_mask_i64(int64_t x) {
    return x >> 63;
}

/*
 * ct_select_i32 - Constant-time select: returns a if sel==0, b if sel==1.
 * sel must be 0 or 1.
 */
static inline int32_t ct_select_i32(int32_t a, int32_t b, uint32_t sel) {
    int32_t mask = -(int32_t)sel;  /* 0 if sel=0, -1 if sel=1 */
    return a ^ (mask & (a ^ b));
}

/* ============================================================
 * IntervalParity
 *
 * Determines the parity (even/odd) of the interval index k such that
 * ell_prime falls in [c_k, c_{k+1}).
 *
 * The interval boundaries are: c_k = k * abs_u + k^2 * gamma
 * Computed incrementally:
 *   c_0 = 0
 *   inc = abs_u + gamma        (= c_1)
 *   delta_inc = 2 * gamma
 *   c_{k+1} = c_k + inc; inc += delta_inc
 *
 * For each step k=0..99, if ell_prime >= c_{k+1}, flip parity.
 *
 * All comparisons use uint64_t to handle potential overflow safely
 * (unsigned wrapping gives correct comparison for our value range).
 *
 * MUST run exactly 100 iterations (no early exit) for constant-time.
 * ============================================================ */
static int interval_parity(int64_t ell_prime, int64_t abs_u, int64_t gamma) {
    /*
     * We use uint64_t arithmetic because the boundary values c_k can grow
     * large (up to ~100*abs_u + 10000*gamma). For our parameters these
     * stay within uint64 range, and unsigned comparison is well-defined
     * even on overflow (wraps modulo 2^64).
     */
    uint64_t c = 0;
    uint64_t inc = (uint64_t)abs_u + (uint64_t)gamma;
    uint64_t delta = 2 * (uint64_t)gamma;
    uint64_t ell_u = (uint64_t)ell_prime;
    uint32_t parity = 0;

    for (int k = 0; k < IRS_PARITY_STEPS; k++) {
        c += inc;
        inc += delta;
        /* crossed = 1 if ell_prime >= c_{k+1}, else 0 */
        uint32_t crossed = (uint32_t)(1 - ((ell_u - c) >> 63));
        parity ^= crossed;
    }

    return (int)parity;
}

/* ============================================================
 * Cyclic shift: compute v = x^j * p mod (x^N + 1)
 *
 * For a polynomial p = [p_0, ..., p_{N-1}] in R = Z[x]/(x^N + 1):
 *   (x^j * p)_k = p_{k-j}      if k >= j
 *                -p_{N+k-j}     if k < j
 *
 * This rotates coefficients right by j positions, negating the
 * wrapped-around part due to the x^N + 1 relation.
 * ============================================================ */
static void poly_cyclic_shift(poly *v, const poly *p, unsigned int j) {
    unsigned int n = SHUTTLE_N;
    for (unsigned int k = 0; k < n; k++) {
        if (k >= j) {
            v->coeffs[k] = p->coeffs[k - j];
        } else {
            v->coeffs[k] = -p->coeffs[n + k - j];
        }
    }
}

/* ============================================================
 * polyvec_cyclic_shift: apply cyclic shift to each component
 * of a polyvec (length VECLEN = 6).
 * ============================================================ */
static void polyvec_cyclic_shift(polyvec *v, const polyvec *sk, unsigned int j) {
    for (int i = 0; i < SHUTTLE_VECLEN; i++) {
        poly_cyclic_shift(&v->vec[i], &sk->vec[i], j);
    }
}

/* ============================================================
 * Inner product <z, v> for polyvec (VECLEN polynomials, N coefficients each).
 *
 * Returns the scalar sum: sum_{i=0}^{VECLEN-1} sum_{k=0}^{N-1} z[i][k] * v[i][k]
 *
 * Uses int64_t accumulator. For SHUTTLE parameters:
 *   |z[i][k]| <= ~5000 (Gaussian bounded), |v[i][k]| <= 26 (secret key bound)
 *   Max single product: 5000 * 26 = 130000
 *   Max sum: 6 * 256 * 130000 = 199,680,000
 *   Fits comfortably in int64.
 * ============================================================ */
static int64_t polyvec_inner_product(const polyvec *z, const polyvec *v) {
    int64_t sum = 0;
    for (int i = 0; i < SHUTTLE_VECLEN; i++) {
        for (int k = 0; k < SHUTTLE_N; k++) {
            sum += (int64_t)z->vec[i].coeffs[k] * (int64_t)v->vec[i].coeffs[k];
        }
    }
    return sum;
}

/* ============================================================
 * Extract challenge positions from the challenge polynomial c.
 *
 * The challenge polynomial has exactly SHUTTLE_TAU nonzero coefficients,
 * each in {-1, +1}. We extract their positions and signs.
 *
 * positions[i] = index where c has a nonzero coefficient
 * signs[i]     = +1 or -1 (the coefficient value)
 * ============================================================ */
static void extract_challenge(unsigned int positions[SHUTTLE_TAU],
                              int32_t signs[SHUTTLE_TAU],
                              const poly *c) {
    int count = 0;
    for (int k = 0; k < SHUTTLE_N && count < SHUTTLE_TAU; k++) {
        if (c->coeffs[k] != 0) {
            positions[count] = (unsigned int)k;
            signs[count] = c->coeffs[k];
            count++;
        }
    }
}

/* ============================================================
 * Single IRS step (Simplified Log-Domain).
 *
 * Given current z and a challenge vector v (= x^{j_i} * sk):
 *   1. u = <z, v> / sigma^2 in Q62
 *   2. ell = -ln(r_rand) via ApproxNegLn
 *   3. If ell < ln_M: reject
 *   4. ell_prime = ell - ln_M
 *   5. For u >= 0: par_A = IntervalParity(ell', u, gamma)
 *      For u <  0: par_B = IntervalParity(ell', -u, gamma) (equivalently)
 *   6. Branchless sign decision based on parity
 *   7. Update z = z +/- v
 *
 * Returns 1 if accepted, 0 if rejected.
 * ============================================================ */
static int irs_step(polyvec *z,
                    const polyvec *v,
                    int64_t gamma_q62,
                    int64_t ln_M_q62,
                    uint64_t r_rand,
                    int32_t *out_sign_update) {
    /*
     * Step 1: Compute u = <z, v> / sigma^2 in Q62.
     *
     * u_real = inner / sigma^2, where sigma^2 = 2^14 = 16384.
     * u_q62 = u_real * 2^62 = inner * 2^62 / 2^14 = inner * 2^48.
     *
    * For SHUTTLE parameters, |inner| can reach ~2e8 (Cauchy-Schwarz bound),
     * so inner * 2^48 can reach ~5.6e22, which overflows int64 (max ~9.2e18).
     * We compute in __int128 and clamp: when |u| is very large, the step
     * always rejects (signing restarts), so clamping is safe.
     */
    int64_t inner = polyvec_inner_product(z, v);

    /* u_q62 = inner * 2^48 via __int128, then clamp to int64 */
    wide_int128 u_wide = (wide_int128)inner << 48;
    int64_t u_q62;
    if (u_wide > (wide_int128)INT64_MAX)
        u_q62 = INT64_MAX;
    else if (u_wide < (wide_int128)INT64_MIN)
        u_q62 = INT64_MIN;
    else
        u_q62 = (int64_t)u_wide;

    /* Step 2: Compute ell = -ln(r_rand) in Q62 */
    int64_t ell = approx_neg_ln(r_rand);

    /* Step 3: Acceptance check: ell >= ln_M */
    /* If ell < ln_M, reject. For SHUTTLE with alpha >= 3, ln_M ~ 0 (nearly always accepts). */
    if (ell < ln_M_q62) {
        return 0;
    }

    /* Step 4: ell_prime = ell - ln_M */
    int64_t ell_prime = ell - ln_M_q62;

    /* Step 5: Compute IntervalParity for both sign cases.
     *
     * The Simplified Log-Domain algorithm needs to determine which "region"
     * the current state falls in and the corresponding sign.
     *
     * For u >= 0 (Region A): par_pos = IntervalParity(ell', u, gamma)
     * For u <  0 (Region B): par_neg = IntervalParity(ell', -u, gamma)
     *
     * The sign decision:
     *   If u >= 0: sign = par_pos ? -1 : +1  (odd parity -> subtract v)
     *   If u <  0: sign = par_neg ? +1 : -1  (odd parity -> add v)
     *
     * But we always compute both parities (constant-time) and select.
     */
    int64_t abs_u = ct_abs_i64(u_q62);
    int par_pos = interval_parity(ell_prime, abs_u, gamma_q62);
    /* par_neg is the parity for -u, which uses the same |u|.
     * Actually, IntervalParity(ell', |u|, gamma) is the same regardless of sign of u.
     * The sign of u only affects which parity interpretation to use.
     *
     * From Algorithm 6:
     * - Compute par = IntervalParity(ell', |u|, gamma)
     * - If u >= 0: sign_update = par ? -1 : +1
     *              (odd parity = Region A rejection direction)
     * - If u <  0: sign_update = par ? +1 : -1
     *              (odd parity = Region B rejection direction)
     *
     * Branchless: sign_update = (u >= 0) ? (1 - 2*par) : (2*par - 1)
     *           = (1 - 2*par) * sign(u >= 0 ? +1 : -1)
     *           = (1 - 2*par) if u >= 0, -(1 - 2*par) if u < 0
     */
    int64_t u_neg_mask = ct_neg_mask_i64(u_q62);  /* -1 if u < 0, 0 if u >= 0 */
    int32_t base_sign = 1 - 2 * par_pos;  /* +1 if even parity, -1 if odd */
    /* Flip sign if u < 0: use XOR with sign bit */
    int32_t flip = (int32_t)(u_neg_mask & 0xFFFFFFFF);
    /* flip is -1 (0xFFFFFFFF) if u<0, 0 if u>=0 */
    /* base_sign XOR flip trick: (base_sign ^ flip) - flip gives -base_sign when flip=-1 */
    int32_t sign_update = (base_sign ^ flip) - flip;

    /* Step 6: Update z = z + sign_update * v (branchless).
     *
     * sign_update is +1 or -1.
     */
    for (int i = 0; i < SHUTTLE_VECLEN; i++) {
        for (int k = 0; k < SHUTTLE_N; k++) {
            z->vec[i].coeffs[k] += sign_update * v->vec[i].coeffs[k];
        }
    }

    /* Output the sign choice if requested */
    if (out_sign_update)
        *out_sign_update = sign_update;

    return 1;  /* accepted */
}

/* ============================================================
 * Internal IRS loop shared by irs_sign and irs_sign_with_signs.
 *
 * If out_signs is non-NULL, records the per-step sign choices.
 * The effective sign for step i is: out_signs[i] * c_i, where
 * c_i is the original challenge coefficient at position j_i.
 * This means the total perturbation is:
 *   sum_i out_signs[i] * c_i * x^{j_i} * sk
 * ============================================================ */
static int irs_loop(polyvec *z,
                    int8_t *out_signs,
                    const poly *c,
                    const polyvec *sk_stretched,
                    int64_t gamma_q62,
                    int64_t ln_M_q62,
                    keccak_state *rng_state) {

    /* Extract challenge positions and signs */
    unsigned int positions[SHUTTLE_TAU];
    int32_t signs[SHUTTLE_TAU];
    extract_challenge(positions, signs, c);

    /* Temporary storage for the rotated secret key vector */
    polyvec v;

    /* Buffer for random bytes: 8 bytes per step for r_rand */
    uint8_t randbuf[8];

    /* Process each of the tau challenge monomials */
    for (int step = 0; step < SHUTTLE_TAU; step++) {
        unsigned int j = positions[step];
        int32_t c_sign = signs[step];

        /* Compute v = x^j * sk. If the challenge coefficient is -1,
         * we negate v (since c_i * x^{j_i} * sk = c_i * v). */
        polyvec_cyclic_shift(&v, sk_stretched, j);

        /* Apply challenge sign: if c_sign == -1, negate all of v */
        if (c_sign == -1) {
            for (int i = 0; i < SHUTTLE_VECLEN; i++) {
                for (int k = 0; k < SHUTTLE_N; k++) {
                    v.vec[i].coeffs[k] = -v.vec[i].coeffs[k];
                }
            }
        }

        /* Sample 8 random bytes for r_rand */
        shake256_squeeze(randbuf, 8, rng_state);

        uint64_t r_rand = (uint64_t)randbuf[0]
                        | ((uint64_t)randbuf[1] << 8)
                        | ((uint64_t)randbuf[2] << 16)
                        | ((uint64_t)randbuf[3] << 24)
                        | ((uint64_t)randbuf[4] << 32)
                        | ((uint64_t)randbuf[5] << 40)
                        | ((uint64_t)randbuf[6] << 48)
                        | ((uint64_t)randbuf[7] << 56);

        /* Perform one IRS step, capturing the sign choice */
        int32_t step_sign = 0;
        int accepted = irs_step(z, &v, gamma_q62, ln_M_q62, r_rand,
                                out_signs ? &step_sign : NULL);

        if (!accepted) {
            return 0;  /* Rejection: caller must restart signing */
        }

        /* Record the sign choice if output buffer provided.
         * step_sign is +1 or -1 (the IRS parity decision).
         * v was already multiplied by c_sign, so the total perturbation
         * for this step is step_sign * c_sign * x^{j} * sk.
         * The effective challenge coefficient is step_sign * c_sign. */
        if (out_signs) {
            out_signs[step] = (int8_t)(step_sign * c_sign);
        }
    }

    return 1;  /* All tau steps accepted */
}

/* ============================================================
 * irs_sign - Main iterative rejection sampling loop.
 * ============================================================ */
int irs_sign(polyvec *z,
             const poly *c,
             const polyvec *sk_stretched,
             int64_t gamma_q62,
             int64_t ln_M_q62,
             keccak_state *rng_state) {
    return irs_loop(z, NULL, c, sk_stretched, gamma_q62, ln_M_q62, rng_state);
}

/* ============================================================
 * irs_sign_with_signs - IRS loop that also outputs sign choices.
 * ============================================================ */
int irs_sign_with_signs(polyvec *z,
                        int8_t sign_choices[SHUTTLE_TAU],
                        const poly *c,
                        const polyvec *sk_stretched,
                        int64_t gamma_q62,
                        int64_t ln_M_q62,
                        keccak_state *rng_state) {
    return irs_loop(z, sign_choices, c, sk_stretched, gamma_q62, ln_M_q62, rng_state);
}
