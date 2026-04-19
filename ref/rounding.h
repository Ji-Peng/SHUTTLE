#ifndef SHUTTLE_ROUNDING_H
#define SHUTTLE_ROUNDING_H

#include <stdint.h>
#include "params.h"

#define decompose SHUTTLE_NAMESPACE(decompose)
int32_t decompose(int32_t *a0, int32_t a);

#define make_hint SHUTTLE_NAMESPACE(make_hint)
unsigned int make_hint(int32_t z2, int32_t comY);

#define use_hint SHUTTLE_NAMESPACE(use_hint)
int32_t use_hint(int32_t a, unsigned int hint);

/* ============================================================
 * mod 2q helpers (Phase 6b). See docs/NGCC_Sign/SHUTTLE_draft.md
 * sections 6, 13, 15 for background.
 * ============================================================ */

/* highbits_mod_2q: compute HighBits(w) with round-half-up toward +inf,
 * wrapped to [0, HINT_MAX) per spec Alg 9 / HAETAE's decompose_hint.
 *
 * Input: w in [0, 2q).
 * Output: integer bucket index in [0, HINT_MAX). Constant-time. */
#define highbits_mod_2q SHUTTLE_NAMESPACE(highbits_mod_2q)
int32_t highbits_mod_2q(int32_t w);

/* lift_to_2q: compute (2 * u + q * parity) mod 2q, where u is the
 * canonical mod q residue ([0, q)) and parity is 0 or 1 (lift bit).
 *
 * Used to turn the mod q output of NTT into a mod 2q value matching the
 * derived public matrix's factor-2 + q*j structure. Constant-time. */
#define lift_to_2q SHUTTLE_NAMESPACE(lift_to_2q)
int32_t lift_to_2q(int32_t u_mod_q, int32_t parity);

/* ============================================================
 * MakeHint / UseHint (Phase 6b-3). See SHUTTLE_draft.md section 6 for
 * the algebraic identity, and section 15 for the HighBits convention.
 * ============================================================ */

/* make_hint_mod2q: single-coefficient MakeHint (spec Alg 9).
 *   w_h         = round(w / alpha_h)                     [in 0, HINT_MAX)]
 *   tilde_w_h   = round((w - 2*z2) mod 2q) / alpha_h      [in 0, HINT_MAX)]
 *   returns      w_h - tilde_w_h   (signed, in (-HINT_MAX, HINT_MAX)).
 *
 * Input: w in [0, 2q), z2 any int32_t (bounded small).
 * Output: small signed integer (typically within +/- 20 per coefficient). */
#define make_hint_mod2q SHUTTLE_NAMESPACE(make_hint_mod2q)
int32_t make_hint_mod2q(int32_t w, int32_t z2_coef);

/* use_hint_wh_mod2q: recover signer's w_h given verifier's tilde_w and hint h.
 *   tilde_w_h = round(tilde_w / alpha_h)
 *   w_h       = (tilde_w_h + h) mod HINT_MAX
 *
 * Returns w_h in [0, HINT_MAX). Constant-time. */
#define use_hint_wh_mod2q SHUTTLE_NAMESPACE(use_hint_wh_mod2q)
int32_t use_hint_wh_mod2q(int32_t tilde_w, int32_t h);

/* recover_z2_coef_mod2q: recover z_2 coefficient given signer's (w_h, w_0)
 * and verifier's tilde_w.
 *   w_app = alpha_h * w_h + w_0             (approximate comY)
 *   diff  = (w_app - tilde_w)               (in (-2q, 2q))
 *   centered = diff canonicalized to (-q, q]
 *   z_2_recov = centered / 2                 (signed integer div)
 *
 * The recovered z_2 may differ from the true z_2 by up to +/- alpha_h/2
 * per coefficient due to lost mid-bits; this error is absorbed by the
 * verifier's B_v bound. */
#define recover_z2_coef_mod2q SHUTTLE_NAMESPACE(recover_z2_coef_mod2q)
int32_t recover_z2_coef_mod2q(int32_t w_h, int32_t w_0, int32_t tilde_w);

#endif
