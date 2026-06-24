/*
 * avx2/poly_ntt.c -- AVX2 fork of the NTT shim (P03), behind USE_AVX2_NTT.
 *
 * Hard-codes the TWO AVX2 signature families so the upward poly_ntt
 * contract is uniform:
 *   - SIGNED q15361 (SHUTTLE_MODE==128): ntt/invntt take (poly, qdata,
 * ztab, z0/z0inv, scale); pointwise takes (c,a,b,qdata); plus
 * s256_reduce_avx.
 *   - UNSIGNED q61441/q59393: ntt/invntt take (poly, qdata, ztab, cross,
 * ninv); pointwise takes (c,a,b,qdata). AVX2 NTT output is in the
 * backend-native shuffle-network permutation (NOT the scalar bit-reversed
 * order).
 *
 * poly16 (uint16_t[N]) is reinterpret-cast to int16_t* at each backend
 * call: the asm's signed-vs-unsigned view is only its interpretation of
 * the same 16 bits; the cast is sound because sizeof(poly16)==2*N
 * (asserted in poly.h).
 *
 * SIGNED mod-2q canonicalization (q15361 only; feeds P09 LiftToModTwoQ,
 * K13): s256_invntt_tomont_avx leaves a SIGNED value up to ~2q in
 * magnitude.  Before the output is [0,q) (the input convention
 * LiftToModTwoQ requires), we run s256_reduce_avx (one red16 pass ->
 * centered |x| <~ 0.5002q) THEN a final conditional +q for negative lanes.
 * red16 alone is NOT enough -- the +q for negative lanes is mandatory.
 * This canonicalization MUST be bit-identical scalar==avx2==avx512 (K10),
 * enforced by t_canon2q.c / test_freeze_avx.c. The unsigned configs
 * already yield [0,q) directly, no extra step.
 */
#include "poly_ntt.h"

#if SHUTTLE_NTT_SIGNED
/* Map back the centered red16 output of s256_reduce_avx into [0,q): add q
 * to every negative int16 lane (constant-time sign-mask).  Operates in
 * place on the poly16 viewed as signed int16. */
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
    s256_ntt_avx((int16_t *)a->coeffs, s256_ntt_qdata, s256_ntt_zetas_fwd,
                 s256_ntt_z0, s256_ntt_scale);
#elif SHUTTLE_MODE == 256
    s512_ntt_avx((int16_t *)a->coeffs, s512_ntt_qdata, s512_ntt_zetas_fwd,
                 s512_ntt_cross_fwd, s512_ntt_ninv);
#elif SHUTTLE_MODE == 512
    s1024_ntt_avx((int16_t *)a->coeffs, s1024_ntt_qdata,
                  s1024_ntt_zetas_fwd, s1024_ntt_cross_fwd,
                  s1024_ntt_ninv);
#endif
}

void poly_invntt_tomont(poly16 *a)
{
#if SHUTTLE_MODE == 128
    s256_invntt_tomont_avx((int16_t *)a->coeffs, s256_ntt_qdata,
                           s256_ntt_zetas_inv, s256_ntt_z0inv,
                           s256_ntt_scale);
    s256_reduce_avx(
        (int16_t *)a->coeffs);   /* lazy ~2q -> centered |x|<~q/2 */
    canon_signed_to_unsigned(a); /* centered -> [0,q) */
#elif SHUTTLE_MODE == 256
    s512_invntt_tomont_avx((int16_t *)a->coeffs, s512_ntt_qdata,
                           s512_ntt_zetas_inv, s512_ntt_cross_inv,
                           s512_ntt_ninv); /* already [0,q) */
#elif SHUTTLE_MODE == 512
    s1024_invntt_tomont_avx((int16_t *)a->coeffs, s1024_ntt_qdata,
                            s1024_ntt_zetas_inv, s1024_ntt_cross_inv,
                            s1024_ntt_ninv); /* already [0,q) */
#endif
}

void poly_pointwise_montgomery(poly16 *c, const poly16 *a, const poly16 *b)
{
#if SHUTTLE_MODE == 128
    s256_pointwise_avx((int16_t *)c->coeffs, (const int16_t *)a->coeffs,
                       (const int16_t *)b->coeffs, s256_ntt_qdata);
#elif SHUTTLE_MODE == 256
    s512_pointwise_avx((int16_t *)c->coeffs, (const int16_t *)a->coeffs,
                       (const int16_t *)b->coeffs, s512_ntt_qdata);
#elif SHUTTLE_MODE == 512
    s1024_pointwise_avx((int16_t *)c->coeffs, (const int16_t *)a->coeffs,
                        (const int16_t *)b->coeffs, s1024_ntt_qdata);
#endif
}

void poly_ntt_canonical(poly16 *a)
{
    /* Wire-byte (canonical, ref bit-reversed) order.  The scalar ntt_ref
     * order is the canonical reference; the AVX2 backend's own permutation
     * is NOT it. So canonical NTT for an AVX2 build uses the scalar oracle
     * to keep wire bytes byte-exact with the reference KAT (K1). */
    static int inited = 0;
    if (!inited) {
        ntt_ref_init();
        inited = 1;
    }
    ntt_ref(a->coeffs);
}

void poly_ntt_import(poly16 *a)
{
    /* Canonical (ref bit-reversed) order -> AVX2 backend-native slot order
     * via nttunpack.  nttunpack consumes int32 standard-order [0,q)
     * samples; here the input is already a canonical-order [0,q) poly16,
     * so widen each lane to int32 and replay the forward shuffle ladder
     * (no butterflies). */
    int32_t src[N];
    int i;
    for (i = 0; i < N; i++)
        src[i] = (int32_t)a->coeffs[i];
#if SHUTTLE_MODE == 128
    s256_nttunpack_avx((int16_t *)a->coeffs, src);
#elif SHUTTLE_MODE == 256
    s512_nttunpack_avx((int16_t *)a->coeffs, src);
#elif SHUTTLE_MODE == 512
    s1024_nttunpack_avx((int16_t *)a->coeffs, src);
#endif
}
