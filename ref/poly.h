#ifndef SHUTTLE_POLY_H
#define SHUTTLE_POLY_H

#include <stdint.h>
#include "params.h"

typedef struct {
  int32_t coeffs[SHUTTLE_N];
} poly;

#define poly_reduce SHUTTLE_NAMESPACE(poly_reduce)
void poly_reduce(poly *a);
#define poly_caddq SHUTTLE_NAMESPACE(poly_caddq)
void poly_caddq(poly *a);

#define poly_add SHUTTLE_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub SHUTTLE_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

#define poly_ntt SHUTTLE_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont SHUTTLE_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);
#define poly_pointwise_montgomery SHUTTLE_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_chknorm SHUTTLE_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);

#define poly_sq_norm SHUTTLE_NAMESPACE(poly_sq_norm)
int64_t poly_sq_norm(const poly *a);

#define poly_uniform SHUTTLE_NAMESPACE(poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SHUTTLE_SEEDBYTES],
                  uint16_t nonce);

#define poly_uniform_eta SHUTTLE_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[SHUTTLE_CRHBYTES],
                      uint16_t nonce);

#define poly_challenge SHUTTLE_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[SHUTTLE_CTILDEBYTES]);

#define poly_decompose SHUTTLE_NAMESPACE(poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);
#define poly_make_hint SHUTTLE_NAMESPACE(poly_make_hint)
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
#define poly_use_hint SHUTTLE_NAMESPACE(poly_use_hint)
void poly_use_hint(poly *b, const poly *a, const poly *h);

#define polyeta_pack SHUTTLE_NAMESPACE(polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_unpack SHUTTLE_NAMESPACE(polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);

#define polypk_pack SHUTTLE_NAMESPACE(polypk_pack)
void polypk_pack(uint8_t *r, const poly *a);
#define polypk_unpack SHUTTLE_NAMESPACE(polypk_unpack)
void polypk_unpack(poly *r, const uint8_t *a);

#define polyz_pack SHUTTLE_NAMESPACE(polyz_pack)
void polyz_pack(uint8_t *r, const poly *a);
#define polyz_unpack SHUTTLE_NAMESPACE(polyz_unpack)
void polyz_unpack(poly *r, const uint8_t *a);

#define polyz0_pack SHUTTLE_NAMESPACE(polyz0_pack)
void polyz0_pack(uint8_t *r, const poly *a);
#define polyz0_unpack SHUTTLE_NAMESPACE(polyz0_unpack)
void polyz0_unpack(poly *r, const uint8_t *a);

#define polyw1_pack SHUTTLE_NAMESPACE(polyw1_pack)
void polyw1_pack(uint8_t *r, const poly *a);

/* ============================================================
 * z-hi split / combine helpers (rANS encoding path).
 *
 * The signature packer splits each z polynomial into HighBits + LowBits.
 * HighBits feed the unified z-hi rANS context; LowBits are bit-packed
 * uniformly. See agent/rANS/SHUTTLE_rANS.tex §2.5 (z^(0) special case)
 * and §3.4 (z^(1..lenS) main case).
 *
 *   polyz1_split: step alpha_r = 64 (cross-mode constant). Used for
 *                 z^(1..lenS). HighBits land at scale r/alpha_r,
 *                 LowBits in [-alpha_r/2, alpha_r/2).
 *   polyz0_split: step alpha_0' = alpha_r/alpha_1. Used for z^(0).
 *                 The choice alpha_0' = alpha_r/alpha_1 makes
 *                 HighBits_{alpha_0'}(z^(0)) share the same effective
 *                 scale r/alpha_r as HighBits_{alpha_r}(z^(i)), so all
 *                 hi coefficients can share one rANS frequency table.
 *
 * Round-half-up convention: ties at -step/2 keep lo = -step/2 (hi = 0).
 * Matches rounding.c::highbits_mod_2q for consistency.
 *
 * Combine routines are exact inverses (step is a power of two, so the
 * split is round-trip perfect).
 *
 * No vocabulary check is needed here: with the 11*sigma SampleY
 * truncation, HighBits stays inside the rANS vocabulary tight bound
 * derived in SHUTTLE_rANS.tex §3.
 * ============================================================ */
#define polyz1_split SHUTTLE_NAMESPACE(polyz1_split)
void polyz1_split(int32_t *hi, int32_t *lo, const poly *a);

#define polyz1_combine SHUTTLE_NAMESPACE(polyz1_combine)
void polyz1_combine(poly *a, const int32_t *hi, const int32_t *lo);

#define polyz0_split SHUTTLE_NAMESPACE(polyz0_split)
void polyz0_split(int32_t *hi, int32_t *lo, const poly *a);

#define polyz0_combine SHUTTLE_NAMESPACE(polyz0_combine)
void polyz0_combine(poly *a, const int32_t *hi, const int32_t *lo);

/* Bit-pack / unpack the LowBits parts.
 *
 *   polyz1_lo_pack / unpack: ALPHA_R_BITS = 6 bits/coef.
 *   polyz0_lo_pack / unpack: ALPHA_0P_BITS = 3 (mode-128) or 2 bits/coef.
 *
 * The bijection [-step/2, step/2) -> [0, step-1] is the low-bits mask:
 *   packed = lo & (step - 1)
 * with the inverse
 *   lo = (packed >= step/2) ? packed - step : packed.
 *
 * Buffer sizes: SHUTTLE_POLYZ1_LO_PACKEDBYTES / SHUTTLE_POLYZ0_LO_PACKEDBYTES
 * are defined in params.h so that SHUTTLE_BYTES can reference them.
 */
#define polyz1_lo_pack SHUTTLE_NAMESPACE(polyz1_lo_pack)
void polyz1_lo_pack(uint8_t *r, const int32_t *lo);

#define polyz1_lo_unpack SHUTTLE_NAMESPACE(polyz1_lo_unpack)
void polyz1_lo_unpack(int32_t *lo, const uint8_t *r);

#define polyz0_lo_pack SHUTTLE_NAMESPACE(polyz0_lo_pack)
void polyz0_lo_pack(uint8_t *r, const int32_t *lo);

#define polyz0_lo_unpack SHUTTLE_NAMESPACE(polyz0_lo_unpack)
void polyz0_lo_unpack(int32_t *lo, const uint8_t *r);

/* ============================================================
 * mod 2q helpers (Phase 6b-1)
 * ============================================================ */

/* Per-coefficient HighBits for polynomial with coeffs in [0, 2q). */
#define poly_highbits_mod_2q SHUTTLE_NAMESPACE(poly_highbits_mod_2q)
void poly_highbits_mod_2q(poly *out, const poly *w);

/* Per-coefficient LSB extraction into a bitmap (N/8 bytes, bit-packed). */
#define poly_lsb_extract SHUTTLE_NAMESPACE(poly_lsb_extract)
void poly_lsb_extract(uint8_t *bitmap, const poly *w);

/* Per-coefficient lift: out[j] = (2 * u[j] + q * (Y0[j] & 1)) mod 2q. */
#define poly_lift_to_2q SHUTTLE_NAMESPACE(poly_lift_to_2q)
void poly_lift_to_2q(poly *out, const poly *u_mod_q, const poly *Y0);

/* In-place w[j] <- (w[j] - 2 * z2[j]) mod 2q. */
#define poly_sub_2z2_mod2q SHUTTLE_NAMESPACE(poly_sub_2z2_mod2q)
void poly_sub_2z2_mod2q(poly *w, const poly *z2);

/* CompressY on the first slot of y (NGCC-Signature Alg 8):
 *   out[j] = round_half_up(in[j] / alpha_1)
 *
 * Since alpha_1 is a power of 2, implemented as arithmetic right shift
 * with a +alpha_1/2 bias. Produces integer coefficients (the "compressed"
 * Y_0). The remaining slots of y are left untouched by CompressY and do
 * not need this function. */
#define poly_compress_y_slot0 SHUTTLE_NAMESPACE(poly_compress_y_slot0)
void poly_compress_y_slot0(poly *out, const poly *in);

/* ============================================================
 * Per-poly MakeHint / UseHint (Phase 6b-3).
 * ============================================================ */

/* poly_make_hint_mod2q: per-coefficient MakeHint. See rounding.h. */
#define poly_make_hint_mod2q SHUTTLE_NAMESPACE(poly_make_hint_mod2q)
void poly_make_hint_mod2q(poly *h, const poly *w, const poly *z2);

/* poly_use_hint_wh_mod2q: per-coefficient recovery of signer's w_h. */
#define poly_use_hint_wh_mod2q SHUTTLE_NAMESPACE(poly_use_hint_wh_mod2q)
void poly_use_hint_wh_mod2q(poly *w_h, const poly *tilde_w, const poly *h);

/* poly_recover_z2_mod2q: per-coefficient recovery of z_2 from (w_h, w_0, tilde_w). */
#define poly_recover_z2_mod2q SHUTTLE_NAMESPACE(poly_recover_z2_mod2q)
void poly_recover_z2_mod2q(poly *z2_out,
                           const poly *w_h,
                           const poly *w_0,
                           const poly *tilde_w);

#endif
