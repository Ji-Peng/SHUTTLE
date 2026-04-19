/*
 * packing_v3.h - NGCC-Signature Alg 2 compressed signature format, v3.
 *
 * v3 entropy-codes THREE distributions with rANS:
 *   1. Z_0 = CompressY(z^(0)) as a single Gaussian-like distribution.
 *      Replaces the Phase 6c-era 11-bit fixed polyz0_pack; see
 *      `shuttle_rans_encode_z0`.
 *   2. HighBits(z^{(1..L)}) via shuttle_rans_encode_z1; the matching
 *      LowBits part is bit-packed uniformly at ALPHA_H_BITS bits/coef
 *      (polyz1_lo_pack).
 *   3. Hint h via shuttle_rans_encode (hint table).
 *
 * Byte layout (mode-128, SHUTTLE_BYTES_V3 = 1450 B):
 *   [  0 ..  32)          seedC (c_tilde)
 *   [ 32 ..  36)          irs_signs bitmap (ceil(TAU/8))
 *   [ 36 ..  38)          uint16 LE: Z_0 rANS length
 *   [ 38 .. 278)          Z_0 rANS stream + zero pad          (240 B reserved)
 *   [278 .. 950)          L * polyz1_lo_pack(lo[i])           (3 * 224 = 672 B)
 *   [950 .. 952)          uint16 LE: z1 rANS length
 *   [952 ..1192)          z1 rANS high stream + zero pad      (240 B reserved)
 *   [1192..1194)          uint16 LE: hint rANS length
 *   [1194..1450)          hint rANS + zero pad                (256 B reserved)
 *
 * SHUTTLE_BYTES_V3 = CTILDEBYTES + IRS_SIGNBYTES
 *                  + Z0_RANS_BLOCK_BYTES
 *                  + L * POLYZ1_LO_PACKEDBYTES
 *                  + Z1_RANS_BLOCK_BYTES
 *                  + HINT_BLOCK_BYTES.
 *
 * Signer-side reject conditions (caller must restart IRS):
 *   - Any Z_0 coefficient falls outside the z0 rANS vocabulary.
 *   - Any HighBits value of z^{(1..L)} falls outside the z1 vocabulary.
 *   - Any rANS-encoded stream exceeds its reservation.
 *
 * Verifier-side rejects:
 *   - Length field exceeds reservation.
 *   - rANS decode underflows.
 */

#ifndef SHUTTLE_PACKING_V3_H
#define SHUTTLE_PACKING_V3_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "polyvec.h"

/* Total v3 signature size. */
#define SHUTTLE_BYTES_V3  ( SHUTTLE_CTILDEBYTES \
                          + SHUTTLE_IRS_SIGNBYTES \
                          + SHUTTLE_Z0_RANS_BLOCK_BYTES \
                          + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES \
                          + SHUTTLE_Z1_RANS_BLOCK_BYTES \
                          + SHUTTLE_HINT_BLOCK_BYTES )

/* pack_sig_v3
 *
 * Inputs mirror pack_sig_v2. Returns 0 on success; negative on
 * z1/hint rANS OOV or overflow (caller rejects the signing round). */
#define pack_sig_v3 SHUTTLE_NAMESPACE(pack_sig_v3)
int pack_sig_v3(uint8_t sig[SHUTTLE_BYTES_V3],
                const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                const int8_t irs_signs[SHUTTLE_TAU],
                const poly *z_1,               /* array of 1 + L polys */
                const polyveck *h);

/* unpack_sig_v3
 *
 * Outputs mirror unpack_sig_v2. Returns 0 on success; nonzero on malformed
 * stream. */
#define unpack_sig_v3 SHUTTLE_NAMESPACE(unpack_sig_v3)
int unpack_sig_v3(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                  int8_t irs_signs[SHUTTLE_TAU],
                  poly *z_1,                   /* array of 1 + L polys */
                  polyveck *h,
                  const uint8_t sig[SHUTTLE_BYTES_V3]);

#endif /* SHUTTLE_PACKING_V3_H */
