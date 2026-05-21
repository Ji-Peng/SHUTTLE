/*
 * packing.h - Serialization of keys and signatures for SHUTTLE
 * (NGCC-Signature Alg 2 compressed form, two-stream rANS entropy coding).
 *
 * Two rANS streams, one shared frequency table for z-hi + (when scales
 * coincide) for hint as well. See SHUTTLE_rANS.tex §3.3 for the design.
 *   z-hi : HighBits_{alpha_0'}(z^(0)) ⨁ HighBits_{alpha_r}(z^(1..lenS))
 *          via shuttle_rans_encode_zhi
 *   hint : MakeHint output, via shuttle_rans_encode_hint
 *
 * LowBits are bit-packed uniformly:
 *   - z^(0)        : ALPHA_0P_BITS/coef via polyz0_lo_pack
 *   - z^(1..lenS)  : ALPHA_R_BITS = 6/coef via polyz1_lo_pack
 *
 * Byte layout (see packing.c::OFF_* for the offsets):
 *   seedC || irs_signs
 *   || uint16 zhi_rans_len  || rANS(z-hi)      + pad to ZHI_RESERVED
 *   || polyz0_lo_pack(lo(z^(0)))
 *   || L * polyz1_lo_pack(lo(z^(1..lenS)))
 *   || uint16 hint_rans_len || rANS(hint)      + pad to HINT_RESERVED
 *
 * SHUTTLE_BYTES = CTILDEBYTES + IRS_SIGNBYTES + ZHI_BLOCK_BYTES
 *               + POLYZ0_LO_PACKEDBYTES + L * POLYZ1_LO_PACKEDBYTES
 *               + HINT_BLOCK_BYTES.
 *
 * Signer-side reject conditions (caller restarts IRS):
 *   - Any rANS-encoded stream exceeds its reservation (rc = -2).
 *   (OOV cannot happen: the rANS vocabulary covers the 11*sigma tight
 *    bound, so every legal signature coefficient is in-vocabulary.)
 *
 * Verifier-side rejects:
 *   - Length field exceeds reservation.
 *   - rANS decode underflows or fails the final-state check (x != L).
 */

#ifndef SHUTTLE_PACKING_H
#define SHUTTLE_PACKING_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "polyvec.h"

#define pack_pk SHUTTLE_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[SHUTTLE_PUBLICKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const polyveck *b);

#define unpack_pk SHUTTLE_NAMESPACE(unpack_pk)
void unpack_pk(uint8_t rho[SHUTTLE_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[SHUTTLE_PUBLICKEYBYTES]);

#define pack_sk SHUTTLE_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[SHUTTLE_SECRETKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const uint8_t tr[SHUTTLE_TRBYTES],
             const uint8_t key[SHUTTLE_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e);

#define unpack_sk SHUTTLE_NAMESPACE(unpack_sk)
void unpack_sk(uint8_t rho[SHUTTLE_SEEDBYTES],
               uint8_t tr[SHUTTLE_TRBYTES],
               uint8_t key[SHUTTLE_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[SHUTTLE_SECRETKEYBYTES]);

/* pack_sig
 *
 * Inputs:
 *   sig        : output buffer, SHUTTLE_BYTES bytes.
 *   c_tilde    : CTILDEBYTES challenge seed.
 *   irs_signs  : TAU-element +/-1 sign vector (IRS output).
 *   z_1        : (1 + L) polys. z_1[0] is the CompressY'd Z_0, z_1[1..L] are
 *                the full-range z[1..L]. Passed as a pointer to the first
 *                polynomial of the array (caller supplies an array of length
 *                1 + L; we access z_1[0] ... z_1[L]).
 *   h          : polyveck of M hint polynomials (MakeHint output).
 *
 * Returns 0 on success; a negative code on rANS overflow / OOV. */
#define pack_sig SHUTTLE_NAMESPACE(pack_sig)
int pack_sig(uint8_t sig[SHUTTLE_BYTES],
             const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
             const int8_t irs_signs[SHUTTLE_TAU],
             const poly *z_1,               /* array of 1 + L polys */
             const polyveck *h);

/* unpack_sig
 *
 * Outputs:
 *   c_tilde    : seedC.
 *   irs_signs  : TAU-element +/-1.
 *   z_1        : array of 1 + L polys.
 *   h          : recovered hint.
 * Inputs:
 *   sig        : byte stream of length SHUTTLE_BYTES.
 *
 * Returns 0 on success; nonzero on malformed stream. */
#define unpack_sig SHUTTLE_NAMESPACE(unpack_sig)
int unpack_sig(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
               int8_t irs_signs[SHUTTLE_TAU],
               poly *z_1,                   /* array of 1 + L polys */
               polyveck *h,
               const uint8_t sig[SHUTTLE_BYTES]);

#endif /* SHUTTLE_PACKING_H */
