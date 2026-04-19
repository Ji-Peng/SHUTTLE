/*
 * packing_v3.h - NGCC-Signature Alg 2 compressed signature format, v3.
 *
 * v3 extends v2 by additionally entropy-coding the HighBits of z[1..L]
 * via the Phase 6c z1 rANS table. For each z_1[i] (i in 1..L) we:
 *   1. Split z_1[i] into (hi[i], lo[i]) using polyz1_split (round-half-up).
 *   2. Bit-pack lo[i] at ALPHA_H_BITS bits/coef (fixed width, uniform).
 *   3. Feed all L*N hi values through shuttle_rans_encode_z1 into a single
 *      reserved block, prefixed by a 2-byte little-endian length.
 *
 * Byte layout (mode-128 example, target SHUTTLE_BYTES_V3 = 1560 B):
 *   [  0 ..  32)   seedC (c_tilde)
 *   [ 32 ..  36)   irs_signs bitmap (ceil(TAU/8))
 *   [ 36 .. 388)   polyz0_pack(Z_0)                  (352 B)
 *   [388 ..1060)   L * polyz1_lo_pack(lo[i])         (3 * 224 = 672 B)
 *   [1060..1062)   uint16 LE: actual z1 rANS length
 *   [1062..1302)   z1 rANS high stream + zero pad    (240 B reserved)
 *   [1302..1304)   uint16 LE: actual hint rANS length
 *   [1304..1560)   hint rANS + zero pad              (256 B reserved)
 *
 * SHUTTLE_BYTES_V3 = CTILDEBYTES + IRS_SIGNBYTES + POLYZ0_PACKEDBYTES
 *                  + L * POLYZ1_LO_PACKEDBYTES
 *                  + Z1_RANS_BLOCK_BYTES + HINT_BLOCK_BYTES.
 *
 * Signer-side reject conditions (caller must restart IRS):
 *   - Any HighBits value of z_1[1..L] falls outside the z1 rANS vocabulary.
 *   - The encoded z1 rANS or hint rANS stream exceeds its reservation.
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
                          + SHUTTLE_POLYZ0_PACKEDBYTES \
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
