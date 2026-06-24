/*
 * avx512/poly_ntt.c -- AVX-512 fork of the NTT shim (P03), behind
 * USE_AVX512_NTT.  Opt-in (built only by the AVX512 targets).
 *
 * The AVX-512 family is UNIFORM across all three configs:
 *   ntt/invntt: (poly, qdata, ztab, scale);  pointwise: (c,a,b,qdata).
 * For q59393 (n=1024) this is the 2x512-coeff SUPERBLOCK kernel; the shim
 * does not see that -- the entry signature is identical.
 *
 * SIGNED mod-2q canonicalization (q15361 only): the AVX-512 invntt has NO
 * reduce_avx export, so the centered->[0,q) step is open-coded here.  The
 * signed AVX-512 invntt finishes with red16+montmul per ZMM, so its output
 * is a centered montmul result (|x| < q); a single conditional +q for
 * negative lanes maps it to [0,q).  This MUST be bit-identical to the
 * scalar smod and the AVX2 reduce_avx+cond-add-q path (K10), enforced by
 * test_freeze_avx.c / t_canon2q.c.
 */
#include "poly_ntt.h"

#if SHUTTLE_NTT_SIGNED
/* Centered signed int16 (|x| < q) -> [0,q): add q to negative lanes
 * (constant-time sign-mask).  Same final step as the AVX2 fork. */
static void canon_signed_to_unsigned(poly16 *a)
{
    int16_t *p = (int16_t *)a->coeffs;
    int i;
    for (i = 0; i < N; i++) {
        int16_t x = p[i];
        x = (int16_t)(x + ((x >> 15) & (int16_t)Q)); /* x<0 ? x+q : x */
        p[i] = x;
    }
}
#endif

void poly_ntt(poly16 *a)
{
#if SHUTTLE_MODE == 128
    s256_ntt_avx512((int16_t *)a->coeffs, s256_ntt512_qdata,
                    s256_ntt512_zetas_fwd, s256_ntt512_scale);
#elif SHUTTLE_MODE == 256
    s512_ntt_avx512((int16_t *)a->coeffs, s512_ntt512_qdata,
                    s512_ntt512_zetas_fwd, s512_ntt512_scale);
#elif SHUTTLE_MODE == 512
    s1024_ntt_avx512((int16_t *)a->coeffs, s1024_ntt512_qdata,
                     s1024_ntt512_zetas_fwd, s1024_ntt512_scale);
#endif
}

void poly_invntt_tomont(poly16 *a)
{
#if SHUTTLE_MODE == 128
    s256_invntt_tomont_avx512((int16_t *)a->coeffs, s256_ntt512_qdata,
                              s256_ntt512_zetas_inv, s256_ntt512_scale);
    canon_signed_to_unsigned(a); /* centered (|x|<q) -> [0,q) */
#elif SHUTTLE_MODE == 256
    s512_invntt_tomont_avx512((int16_t *)a->coeffs, s512_ntt512_qdata,
                              s512_ntt512_zetas_inv, s512_ntt512_scale);
#elif SHUTTLE_MODE == 512
    s1024_invntt_tomont_avx512((int16_t *)a->coeffs, s1024_ntt512_qdata,
                               s1024_ntt512_zetas_inv, s1024_ntt512_scale);
#endif
}

void poly_pointwise_montgomery(poly16 *c, const poly16 *a, const poly16 *b)
{
#if SHUTTLE_MODE == 128
    s256_pointwise_avx512((int16_t *)c->coeffs, (const int16_t *)a->coeffs,
                          (const int16_t *)b->coeffs, s256_ntt512_qdata);
#elif SHUTTLE_MODE == 256
    s512_pointwise_avx512((int16_t *)c->coeffs, (const int16_t *)a->coeffs,
                          (const int16_t *)b->coeffs, s512_ntt512_qdata);
#elif SHUTTLE_MODE == 512
    s1024_pointwise_avx512((int16_t *)c->coeffs,
                           (const int16_t *)a->coeffs,
                           (const int16_t *)b->coeffs, s1024_ntt512_qdata);
#endif
}

void poly_ntt_canonical(poly16 *a)
{
    /* Canonical (ref bit-reversed) wire order via the scalar oracle (K1).
     */
    static int inited = 0;
    if (!inited) {
        ntt_ref_init();
        inited = 1;
    }
    ntt_ref(a->coeffs);
}

void poly_ntt_import(poly16 *a)
{
    /* Canonical order -> AVX-512 backend-native slot order via nttunpack
     * (replays the in-superblock shuffle ladder for q59393). */
    int32_t src[N];
    int i;
    for (i = 0; i < N; i++)
        src[i] = (int32_t)a->coeffs[i];
#if SHUTTLE_MODE == 128
    s256_nttunpack_avx512((int16_t *)a->coeffs, src);
#elif SHUTTLE_MODE == 256
    s512_nttunpack_avx512((int16_t *)a->coeffs, src);
#elif SHUTTLE_MODE == 512
    s1024_nttunpack_avx512((int16_t *)a->coeffs, src);
#endif
}
