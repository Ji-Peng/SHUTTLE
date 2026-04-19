/*
 * packing_v2.h - NGCC-Signature Alg 2 compressed signature format.
 *
 * Layout (mode-128 example, 1990 B total):
 *   [  0 ..  32)   seedC (c_tilde)
 *   [ 32 ..  36)   irs_signs bitmap (TAU bits, packed to ceil(TAU/8) bytes)
 *   [ 36 .. 388)   polyz0_pack(Z_0)          (352 B @ 11 bits/coef, N coef)
 *   [388 .. 1732)  polyz_pack(z[1..L])       (L=3 polys * 448 B)
 *   [1732..1734)   uint16 LE: actual rANS hint length
 *   [1734..1990)   rANS-encoded hint + zero padding (reserved 256 B)
 *
 * Total SHUTTLE_BYTES_V2 = CTILDEBYTES + IRS_SIGNBYTES + POLYZ0_PACKEDBYTES
 *                        + L * POLYZ_PACKEDBYTES + 2 + HINT_RESERVED_BYTES.
 *
 * pack_sig_v2 may fail if:
 *   - A hint coefficient is outside the rANS vocabulary (OOV).
 *   - The rANS-encoded byte stream exceeds the reserved space.
 * In both cases the signer must reject this round and restart the IRS
 * loop with fresh randomness.
 *
 * unpack_sig_v2 may fail if:
 *   - rANS length field exceeds the reserved space.
 *   - rANS decode underflows.
 */

#ifndef SHUTTLE_PACKING_V2_H
#define SHUTTLE_PACKING_V2_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "polyvec.h"

/* Mode-dependent hint reservation lives in params.h (shared with v3). */

/* Total v2 signature size. */
#define SHUTTLE_BYTES_V2  ( SHUTTLE_CTILDEBYTES \
                          + SHUTTLE_IRS_SIGNBYTES \
                          + SHUTTLE_POLYZ0_PACKEDBYTES \
                          + SHUTTLE_L * SHUTTLE_POLYZ_PACKEDBYTES \
                          + SHUTTLE_HINT_BLOCK_BYTES )

/* ============================================================
 * Pack / unpack
 * ============================================================ */

/* pack_sig_v2
 *
 * Inputs:
 *   sig        : output buffer, SHUTTLE_BYTES_V2 bytes.
 *   c_tilde    : CTILDEBYTES challenge seed.
 *   irs_signs  : TAU-element +/-1 sign vector (IRS output).
 *   z_1        : (1 + L) polys. z_1[0] is the CompressY'd Z_0, z_1[1..L] are
 *                the full-range z[1..L]. Passed as a pointer to the first
 *                polynomial of the array (caller supplies an array of length
 *                1 + L; we access z_1[0] ... z_1[L]).
 *   h          : polyveck of M hint polynomials (MakeHint output).
 *
 * Returns 0 on success; a negative code on rANS overflow / OOV. */
#define pack_sig_v2 SHUTTLE_NAMESPACE(pack_sig_v2)
int pack_sig_v2(uint8_t sig[SHUTTLE_BYTES_V2],
                const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                const int8_t irs_signs[SHUTTLE_TAU],
                const poly *z_1,               /* array of 1 + L polys */
                const polyveck *h);

/* unpack_sig_v2
 *
 * Outputs:
 *   c_tilde    : seedC.
 *   irs_signs  : TAU-element +/-1.
 *   z_1        : array of 1 + L polys.
 *   h          : recovered hint.
 * Inputs:
 *   sig        : byte stream of length SHUTTLE_BYTES_V2.
 *
 * Returns 0 on success; nonzero on malformed stream. */
#define unpack_sig_v2 SHUTTLE_NAMESPACE(unpack_sig_v2)
int unpack_sig_v2(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
                  int8_t irs_signs[SHUTTLE_TAU],
                  poly *z_1,                   /* array of 1 + L polys */
                  polyveck *h,
                  const uint8_t sig[SHUTTLE_BYTES_V2]);

#endif /* SHUTTLE_PACKING_V2_H */
