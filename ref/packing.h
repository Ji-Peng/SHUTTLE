/*
 * packing.h -- byte primitives + scheme (de)serialization for SHUTTLE
 * (P04).
 *
 * This is the wire-format contract every other component serializes
 * through.  Three layers:
 *
 *   1. Byte primitives (LITTLE-ENDIAN; data-independent schedule):
 *        integer_to_bytes / bytes_to_integer  -- LE multi-byte integers.
 *        poly_to_bytes / bytes_to_poly        -- LSB-first d-bit coeff
 * pack. These are K10-load-bearing: their byte schedule depends ONLY on
 * the loop index, never on the coefficient value, so the pk/sk/com bytes
 *      are KAT-stable and constant-time across ref/avx2/avx512.
 *
 *   2. Scheme serializers:
 *        pack_pk   / unpack_pk    (pkEncode  / pkDecode)
 *        pack_sk   / unpack_sk    (skEncode  / skDecode)
 *        pack_com  / unpack_com   (EncodeCom / decode-with-range-check)
 *
 *   3. ct_range_reject -- the shared branchless range guard used by the
 *      non-power-of-2 EncodeCom / hint decode (K14) and the secret-range
 *      checks in unpack_sk.  NEVER reduce a non-canonical value mod H_h /
 *      the bound -- range-CHECK then REJECT (injectivity / SUF-CMA).
 *
 * BIT/BYTE ORDER (binding, Description.tex 116-126, 691-774):
 *   - integer_to_bytes / bytes_to_integer: byte i = floor(x / 256^i) mod
 *     256 (little-endian).
 *   - poly_to_bytes / bytes_to_poly: each coeff's d bits are packed
 *     LSB-FIRST, concatenated across coeffs, the final byte zero-padded.
 *   (BytesToBits -- the MSB-first-within-a-byte convention used ONLY by
 *    SamplerU in P08 -- is the OPPOSITE order and is explicitly OUT OF
 *    SCOPE here; do not conflate it with this LSB-first packing.)
 */
#ifndef SHUTTLE_PACKING_H
#define SHUTTLE_PACKING_H

#include <stdint.h>

#include "params.h"
#include "poly.h"

/* ---------------------------------------------------------------------- *
 *  Byte primitives (little-endian; data-independent schedule)            *
 * ----------------------------------------------------------------------
 */

/* integer_to_bytes: write `len` little-endian bytes of x.  Returns 0 on
 * success, -1 if x does not fit in `len` bytes (bottom check; only the
 * small fixed lengths 2/4 are ever requested by SHUTTLE callers).  The
 * output schedule (out[i] = (x >> 8*i) & 0xFF) is data-independent. */
int integer_to_bytes(uint8_t *out, uint64_t x, unsigned len);

/* bytes_to_integer: sum_i in[i] * 256^i (little-endian). */
uint64_t bytes_to_integer(const uint8_t *in, unsigned len);

/* poly_to_bytes: pack N coeffs, each in [0, 2^d), d bits LSB-first, final
 * byte zero-padded.  Writes ceil(N*d/8) bytes.  Precondition: every coeff
 * 0 <= w[i] < 2^d (caller guarantees; the d-bit mask is applied anyway so
 * any high bits are silently dropped, keeping the schedule fixed). */
void poly_to_bytes(uint8_t *out, const poly *w, unsigned d);

/* bytes_to_poly: inverse of poly_to_bytes; regroups d-bit LSB-first chunks
 * into N coeffs in [0, 2^d).  No range check (per-field callers add it).
 */
void bytes_to_poly(poly *w, const uint8_t *in, unsigned d);

/* ---------------------------------------------------------------------- *
 *  ct_range_reject -- shared branchless range guard (K14)                *
 *                                                                        *
 *  Returns 0 if lo <= v <= hi, else a non-zero "fail" accumulator.  The *
 *  caller OR-accumulates the result across a fixed-length loop and *
 *  branches only ONCE at the end (no mid-loop branch on the value), so *
 *  the reject decision does not leak per-coefficient timing.  Used by *
 *  unpack_com (comY_h in [0,H_h), the non-power-of-2 gate), unpack_sk *
 *  (secret ranges), and shared with the rANS hint decode in P10. *
 *                                                                        *
 *  This is a CHECK-then-REJECT, NOT a mod reduction: reducing an *
 *  out-of-range value mod H_h would make two distinct byte strings *
 *  decode to the same value, breaking injectivity -> breaking SUF-CMA. *
 * ----------------------------------------------------------------------
 */
uint32_t ct_range_reject(int32_t v, int32_t lo, int32_t hi);

/* ---------------------------------------------------------------------- *
 *  Public key (pkEncode / pkDecode)                                      *
 * ----------------------------------------------------------------------
 */

/* pk = seedA || PolyToBytes(b[0]/alpha_b, d_b) || ... (EM polys).  `b` is
 * coefficient-domain in [0,q), each coeff an exact multiple of alpha_b
 * (guaranteed by RoundB in P09).  The /alpha_b is a power-of-two right
 * shift (public data). */
void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES],
             const uint8_t seedA[SEEDBYTES], const poly b[EM]);

/* pkDecode: read seedA, then EM blocks; RANGE-CHECK each b1 field <
 * ceil(q/alpha_b) (2^d_b can exceed ceil(q/alpha_b), so a malformed pk can
 * carry an over-range b1) and reject; else rescale *alpha_b into [0,q).
 * Returns 0 on success, -1 on an out-of-range field. */
int unpack_pk(uint8_t seedA[SEEDBYTES], poly b[EM],
              const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);

/* ---------------------------------------------------------------------- *
 *  Secret key (skEncode / skDecode)                                      *
 * ----------------------------------------------------------------------
 */

/* sk = seedA || PolyToBytes(b/alpha_b, d_b)[EM]
 *         || masterSeed (K) || tr
 *         || PolyToBytes(s + BS_ENC, d_s)[ELL]
 *         || PolyToBytes(e' + BE_ENC, d_e)[EM].
 * The two secret segments are UNEQUAL length when ELL != EM (256/512). */
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const uint8_t seedA[SEEDBYTES], const poly b[EM],
             const uint8_t masterSeed[CHALLENGESEEDBYTES],
             const uint8_t tr[CHALLENGESEEDBYTES], const poly s[ELL],
             const poly ep[EM]);

/* skDecode: parse fixed-length segments L->R; rescale b *alpha_b; apply
 * the -BS_ENC / -BE_ENC unshifts; RANGE-CHECK every decoded secret coeff
 * into
 * [-BS_ENC, BS_ENC] / [-BE_ENC, BE_ENC] and reject if any is out of range.
 * Returns 0 on success, -1 on a malformed sk. */
int unpack_sk(uint8_t seedA[SEEDBYTES], poly b[EM],
              uint8_t masterSeed[CHALLENGESEEDBYTES],
              uint8_t tr[CHALLENGESEEDBYTES], poly s[ELL], poly ep[EM],
              const uint8_t sk[CRYPTO_SECRETKEYBYTES]);

/* ---------------------------------------------------------------------- *
 *  Commitment (EncodeCom / decode) -- used by HashCh in Sign/Verify      *
 * ----------------------------------------------------------------------
 */

/* EncodeCom = PolyToBytes(comY_h, d_h) || PolyToBytes(comY_0, 1).
 * Precondition (caller-guaranteed, NO encode-side check): comY_h[k] in
 * [0,H_h), comY_0[k] in {0,1}. */
#define ENCODECOM_BYTES (POLYWH_PACKEDBYTES + POLYW0_PACKEDBYTES)
void pack_com(uint8_t out[ENCODECOM_BYTES], const poly *comY_h,
              const poly *comY_0);

/* RANGE-CHECK decode (K14): rejects (returns -1) if ANY comY_h coeff >=
 * H_h (the gap [H_h, 2^d_h) is non-empty since H_h=30/120/58 is not a
 * power of two and d_h=5/7/6 over-covers it).  comY_0 is a 1-bit field so
 * every decoded value is trivially in {0,1}.  NEVER reduces mod H_h.
 * Returns 0 on success.  This is the canonical-encoding / injectivity
 * gate; it is primarily a negative-test / robustness surface (the live
 * Verify flow re-encodes its own reconstructed commitment rather than
 * round-tripping through unpack_com -- see 04-Poly-Packing.md). */
int unpack_com(poly *comY_h, poly *comY_0,
               const uint8_t in[ENCODECOM_BYTES]);

#endif /* SHUTTLE_PACKING_H */
