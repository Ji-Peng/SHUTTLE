/*
 * packing.h - Serialization of keys and signatures for SHUTTLE (NGCC-Signature
 * Alg 2 compressed form with full rANS entropy coding).
 *
 * Three distributions are entropy-coded with independent rANS tables
 * (see rans_tables.h):
 *   1. Z_0 = CompressY(z^(0))                    via shuttle_rans_encode_z0
 *   2. HighBits(z^{(1..L)})                      via shuttle_rans_encode_z1
 *      (LowBits bit-packed uniformly at ALPHA_H_BITS/coef by polyz1_lo_pack)
 *   3. Hint h                                     via shuttle_rans_encode
 *
 * Byte layout (mode-128, SHUTTLE_BYTES = 1450 B):
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
 * SHUTTLE_BYTES = CTILDEBYTES + IRS_SIGNBYTES + Z0_RANS_BLOCK_BYTES
 *               + L * POLYZ1_LO_PACKEDBYTES + Z1_RANS_BLOCK_BYTES
 *               + HINT_BLOCK_BYTES.
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
