/*
 * irs.c -- SHUTTLE RejectSample / R reference (scalar, integer-only, CT).
 *
 * The rejection-free inner masking transition (Algorithms alg:RejectSample
 * and alg:Ryv).  See irs.h for the full contract: the fresh 0x09||seed_y
 * context (K2), the ascending-j traversal (K5), the base-2->ln multiply by
 * 2 r^2 ln2 (K12 / MS-C2), the N=29 / 15-boundary-pair interval test, and
 * the isochrony / secret-vs-public leakage argument (the P13 anchor).
 *
 * Constant-time discipline (irs.h "ISOCHRONY / LEAKAGE"):
 *   - the ONLY branch on the IRS path is the ascending-j  if (c[j])  gate,
 *     whitelisted because c is PUBLIC (recomputed by the verifier);
 *   - the sign-normalize, the 15-pair interval test, and the z += flag*v
 *     update are branchless (two's-complement masks);
 *   - no division/modulo; the u-vs-boundary comparison is exact in
 * __int128.
 *
 * The inner products run on the RAW signed scheme-domain coefficients (NOT
 * mod-q-centered): sk_tilde and y are small signed integers, so we do NOT
 * use poly_sqnorm (which centers mod q and would corrupt them) -- we
 * accumulate the products directly in int64.
 */
#include "irs.h"

#include <stdint.h>
#include <string.h>

#include "params.h"
#include "poly.h"
#include "sampler_u.h"
#include "symmetric.h"

/* The u-vs-boundary comparison and the 2 r^2 ln2 multiply use GNU __int128
 * (a GCC/Clang extension ISO C does not define; -Wpedantic flags it).
 * Localize the suppression to this TU's __int128 use, exactly as P06 did.
 */
#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wpedantic"
#endif

/* ===================================================================== *
 *  IRS R-transition fixed-point constant (R3 / MS-C2)                    *
 *  Owned by tools/gen_irs_consts.py; re-derived by `make check-consts`.  *
 * ===================================================================== */
/* @@AUTOGEN:irs_consts@@ BEGIN */
/* IRS R-transition fixed-point constant (R3 / MS-C2; gen_irs_consts.py).
 *
 *   2 r^2 ln 2 = 2*825^2*ln2 = 943546.5995372255524442...
 *   R2LN2_QSHIFT = 44   (the Q-scale F; u is carried at Q44)
 *   R2LN2_QF     = round(2 r^2 ln2 * 2^44) = 16599047320634951608
 *                = 0xE65BA997B45887B8  (fits uint64_t:
 * 16599047320634951608 < 2^64)
 *
 * relative error of the rounded constant = 2^-65.1577; amplified additive
 * natural-log error |ln U|*relerr <= 50.53*2^-65.1577 ~ 2^-59.4986
 * (binding |ln U|); total delta_log ~ 2^-58.6572, accumulated delta_tau ~
 * 2^-46.8043
 * (< 2^-45 budget).  See gen_irs_consts.py +
 * log/irs_consts_derivation.txt.
 */
#define R2LN2_QSHIFT 44
#define R2LN2_QF UINT64_C(16599047320634951608)
/* @@AUTOGEN:irs_consts@@ END */

_Static_assert(TWO_RSQ == 1361250L, "2 r^2 = 2*825^2 = 1361250");

/* ===================================================================== *
 *  Exact integer inner products (raw signed coeffs, int64 accum)        *
 * ===================================================================== */

/*
 * V = <sk_tilde, sk_tilde> = sum over KVEC polys, n coeffs each.  sk_tilde
 * = StretchS(sk) is norm-bounded (||sk_tilde|| <= B_k ~ 296, so V ~ 88000
 * < 2^17); the int64 accumulation cannot overflow.  Computed ONCE per
 * RejectSample and reused for every shift (isometry,
 * K5/Description.tex:1461).
 */
static int64_t sk_tilde_norm2(const poly sk_tilde[KVEC])
{
    int64_t acc = 0;
    unsigned i, k;
    for (i = 0; i < KVEC; ++i)
        for (k = 0; k < N; ++k) {
            int64_t c = (int64_t)sk_tilde[i].coeffs[k];
            acc += c * c;
        }
    return acc;
}

/*
 * t = <y, v> = sum over KVEC polys, n coeffs each.  |y| <= ~5000 (wide
 * Gaussian ~11 sigma of r=825), |v| <= ||sk_tilde||_inf <= ~26 per coeff,
 * so |t| <= KVEC*N*5000*26 ~ 2^28; int64 is safe (the V-tail term -m^2 V
 * with m<=28 dominates the boundary magnitude, still well within int64).
 */
static int64_t inner_y_v(const poly y[KVEC], const poly v[KVEC])
{
    int64_t acc = 0;
    unsigned i, k;
    for (i = 0; i < KVEC; ++i)
        for (k = 0; k < N; ++k)
            acc += (int64_t)y[i].coeffs[k] * (int64_t)v[i].coeffs[k];
    return acc;
}

/*
 * v <- sk_tilde . X^j  in R = Z[X]/(X^n+1): a negacyclic right-rotation by
 * j, negating the part that wraps past degree n.  Per-component over all
 * KVEC polys: v[i].coeffs[k] = (k >= j) ?  sk_tilde[i].coeffs[k-j] :
 * -sk_tilde[i].coeffs[n + k - j]. (j is PUBLIC -- it is a challenge index
 * -- so the index arithmetic is fine; the from-scratch rotation is the KAT
 * reference, S10.)
 */
static void poly_shift_negacyclic(poly v[KVEC], const poly src[KVEC],
                                  unsigned int j)
{
    unsigned i, k;
    for (i = 0; i < KVEC; ++i) {
        for (k = 0; k < j; ++k)
            v[i].coeffs[k] = -src[i].coeffs[N + k - j];
        for (k = j; k < N; ++k)
            v[i].coeffs[k] = src[i].coeffs[k - j];
    }
}

/* z[i] += flag * v[i]  for every coeff, flag in {-1,+1}, branchless. */
static void poly_axpy_flag(poly z[KVEC], const poly v[KVEC], int64_t flag)
{
    unsigned i, k;
    int32_t f = (int32_t)flag; /* -1 or +1 */
    for (i = 0; i < KVEC; ++i)
        for (k = 0; k < N; ++k)
            z[i].coeffs[k] += f * v[i].coeffs[k];
}

/* ===================================================================== *
 *  The u fixed-point form (K12 / MS-C2)                                  *
 * ===================================================================== */

/*
 * u_q44 = (2 r^2 ln2) * log2(U) at the pinned scale Q44, in __int128.
 * log2(U) = frac_q62/2^62 - a (the unfolded SamplerU pair).  Computed as:
 *   u_frac = round(R2LN2_QF * frac_q62 / 2^62) = (R2LN2_QF*frac + 2^61) >>
 * 62 u_a    = a * R2LN2_QF u      = u_frac - u_a frac_q62 is in [0, 2^62)
 * (nonnegative log2(b)), so R2LN2_QF*frac fits an unsigned __int128 (<
 * 2^126).  The result is a signed Q44 value.
 */
static __int128 sampler_u_to_u_q44(sampler_u_res ell)
{
    unsigned __int128 prod = (unsigned __int128)R2LN2_QF *
                             (unsigned __int128)(uint64_t)ell.frac_q62;
    /* round-to-nearest >> 62 (unbiased), then to signed Q44. */
    unsigned __int128 u_frac = (prod + ((unsigned __int128)1 << 61)) >> 62;
    __int128 u_a = (__int128)ell.a * (__int128)R2LN2_QF;
    return (__int128)u_frac - u_a;
}

/* ===================================================================== *
 *  One R transition (Algorithm alg:Ryv)                                 *
 * ===================================================================== */

/*
 * R_transition: draw ell from SamplerU(ctx) (already drawn and passed in
 * by the caller so the x2 batch can be wired -- see reject_sample), form
 * u, sign-normalize v so t = <y,v> > 0, run the 15 boundary-pair interval
 * tests, and apply z <- z + flag*v.  `v` is the shift sk_tilde.X^j; it is
 * mutated (sign-normalized) in place.  z is mutated.
 *
 * Branchless throughout: the t<=0 flip and the flag selection use
 * two's-complement masks; the loop runs all 15 pairs with no early exit.
 */
static void R_transition(poly z[KVEC], poly v[KVEC], int64_t V,
                         sampler_u_res ell)
{
    __int128 u = sampler_u_to_u_q44(ell);
    int64_t t = inner_y_v(z, v);
    int64_t flagmask; /* sign-normalize: flip iff t <= 0 (spec uses <=) */
    int64_t flag64;
    unsigned i;

    /*
     * Sign-normalize so t > 0.  The spec condition is t <= 0 (inclusive at
     * 0).  flipmask = ((t - 1) >> 63) is all-ones iff (t-1) < 0 iff t <= 0
     * (for in-range t); it MUST flip at t == 0 to match the spec's <=.
     * Negate per-coeff via two's complement (x ^ mask) - mask, and t
     * likewise.
     */
    flagmask = (t - 1) >> 63; /* -1 iff t <= 0, else 0 */
    {
        int32_t m32 = (int32_t)flagmask; /* 0 or -1 */
        unsigned k, ii;
        for (ii = 0; ii < KVEC; ++ii)
            for (k = 0; k < N; ++k) {
                int32_t c = v[ii].coeffs[k];
                v[ii].coeffs[k] = (c ^ m32) - m32; /* c, or -c iff mask */
            }
    }
    t = (t ^ flagmask) - flagmask; /* |t| (now t > 0) */

    /*
     * 15 boundary pairs cover the N=29 truncated terms m=0..28.  Pair i
     * tests lo = -2(2i+1)t - (2i+1)^2 V   (the m=2i+1 side) hi = -4 i t -
     * 4 i^2 V        (the m=2i   side) and sets flag=1 iff (lo<<F) < u <=
     * (hi<<F).  At most one i matches, but we run ALL 15 and OR the
     * condition (CT).  flag starts at -1.
     */
    flag64 = -1;
    for (i = 0; i < IRS_BDRY; ++i) {
        int64_t two_i1 = (int64_t)(2u * i + 1u); /* 2i+1 */
        int64_t four_i = (int64_t)(4u * i);      /* 4i   */
        int64_t lo = -2 * two_i1 * t - two_i1 * two_i1 * V;
        int64_t hi = -four_i * t - (int64_t)(4u * i * i) * V;
        __int128 lo_s = (__int128)lo << R2LN2_QSHIFT;
        __int128 hi_s = (__int128)hi << R2LN2_QSHIFT;
        /* cond = (lo_s < u) && (u <= hi_s), as a 0/-1 mask. */
        int64_t c_lo = (int64_t)((lo_s < u) ? 1 : 0);
        int64_t c_hi = (int64_t)((u <= hi_s) ? 1 : 0);
        int64_t cond =
            -(c_lo & c_hi); /* -1 iff in the half-open interval */
        /* flag = cond ? 1 : flag  (branchless select). */
        flag64 = (flag64 & ~cond) | ((int64_t)1 & cond);
    }

    poly_axpy_flag(z, v, flag64);
}

/* ===================================================================== *
 *  RejectSample (Algorithm alg:RejectSample)                            *
 * ===================================================================== */

void reject_sample(xof_ctx *ctx, poly z[KVEC], const poly y[KVEC],
                   const poly *c, const poly sk_tilde[KVEC])
{
    poly v[KVEC];
    int64_t V;
    unsigned j;

    /* z <- y  (copy; y may alias z safely after this). */
    memcpy(z, y, KVEC * sizeof(poly));

    /* V = <sk_tilde, sk_tilde>, computed ONCE (isometry, K5). */
    V = sk_tilde_norm2(sk_tilde);

    /*
     * Strict ascending traversal j = 0..n-1, transition iff c[j] == 1
     * (K5). c is PUBLIC, so this branch is a whitelisted public-data
     * branch (irs.h leakage note).  Each matching j draws ONE SamplerU
     * value off ctx, in ascending-j order, so the squeeze schedule is
     * pinned by ascending j.
     *
     * The SamplerU draw of transition k is independent of the SECOND
     * draw's value, so consecutive matching j's could be x2-batched
     * (sampler_u_x2, S9); the from-scratch scalar path here is the KAT
     * reference -- the x2 pairing is wired and proven bit-identical in
     * test_irs.  We keep the scalar per-j path here for the reference
     * oracle.
     */
    for (j = 0; j < N; ++j) {
        if (c->coeffs[j] == 1) {
            sampler_u_res ell = sampler_u(ctx);
            poly_shift_negacyclic(v, sk_tilde, j); /* v = sk_tilde . X^j */
            R_transition(z, v, V, ell);
        }
    }
}

#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC diagnostic pop
#endif
