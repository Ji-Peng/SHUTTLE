#include <stdint.h>
#include <string.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "rounding.h"
#include "fips202.h"

/*************************************************
* Name:        poly_reduce
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in approximately [-Q/2, Q/2].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    a->coeffs[i] = reduce32(a->coeffs[i]);
}

/*************************************************
* Name:        poly_caddq
*
* Description: For all coefficients of in/out polynomial add Q if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_caddq(poly *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    a->coeffs[i] = caddq(a->coeffs[i]);
}

/*************************************************
* Name:        poly_add
*
* Description: Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
**************************************************/
void poly_add(poly *c, const poly *a, const poly *b) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    c->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtracted from first input polynomial
**************************************************/
void poly_sub(poly *c, const poly *a, const poly *b) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    c->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

/*************************************************
* Name:        poly_ntt
*
* Description: Inplace forward NTT. Coefficients can grow by
*              up to 8*Q in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt(poly *a) {
  ntt(a->coeffs);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Inplace inverse NTT and multiplication by 2^{32}.
*              Input coefficients need to be less than Q in absolute
*              value and output coefficients are again bounded by Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_tomont(poly *a) {
  invntt_tomont(a->coeffs);
}

/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    c->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);
}

/*************************************************
* Name:        poly_chknorm
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients were reduced by reduce32().
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B and 1 otherwise.
**************************************************/
int poly_chknorm(const poly *a, int32_t B) {
  unsigned int i;
  int32_t t;

  if(B > (SHUTTLE_Q - 1) / 8)
    return 1;

  for(i = 0; i < SHUTTLE_N; ++i) {
    /* Absolute value */
    t = a->coeffs[i] >> 31;
    t = a->coeffs[i] - (t & 2 * a->coeffs[i]);

    if(t >= B)
      return 1;
  }

  return 0;
}

/*************************************************
* Name:        poly_sq_norm
*
* Description: Compute the squared L2 norm of the polynomial,
*              i.e., sum of coeffs[i]^2.
*
* Arguments:   - const poly *a: pointer to polynomial
*
* Returns sum of squares as int64_t.
**************************************************/
int64_t poly_sq_norm(const poly *a) {
  unsigned int i;
  int64_t norm = 0;

  for(i = 0; i < SHUTTLE_N; ++i)
    norm += (int64_t)a->coeffs[i] * a->coeffs[i];

  return norm;
}

/**************************************************************/
/*********** Sampling routines ********************************/
/**************************************************************/

/*************************************************
* Name:        rej_uniform
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*              Uses 2-byte sampling with 14-bit mask for q=15361.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients.
**************************************************/
static unsigned int rej_uniform(int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t;

  ctr = pos = 0;
  while(ctr < len && pos + 2 <= buflen) {
    t  = buf[pos++];
    t |= (uint32_t)buf[pos++] << 8;
    t &= 0x3FFF; /* mask to 14 bits */

    if(t < SHUTTLE_Q)
      a[ctr++] = t;
  }

  return ctr;
}

/*************************************************
* Name:        poly_uniform
*
* Description: Sample polynomial with uniformly random coefficients
*              in [0,Q-1] by performing rejection sampling on the
*              output stream of SHAKE128(seed|nonce).
*              Acceptance rate: 15361/16384 ~ 93.8%.
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
/* Need at least 256 * 2 / 0.938 ~ 546 bytes; use 4 SHAKE128 blocks = 672 */
#define POLY_UNIFORM_NBLOCKS ((546 + SHAKE128_RATE - 1) / SHAKE128_RATE)
void poly_uniform(poly *a,
                  const uint8_t seed[SHUTTLE_SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS * SHAKE128_RATE;
  uint8_t buf[POLY_UNIFORM_NBLOCKS * SHAKE128_RATE + 2];
  keccak_state state;

  /* Absorb seed || nonce (LE 2 bytes) */
  uint8_t inbuf[SHUTTLE_SEEDBYTES + 2];
  memcpy(inbuf, seed, SHUTTLE_SEEDBYTES);
  inbuf[SHUTTLE_SEEDBYTES + 0] = (uint8_t)(nonce);
  inbuf[SHUTTLE_SEEDBYTES + 1] = (uint8_t)(nonce >> 8);

  shake128_absorb_once(&state, inbuf, SHUTTLE_SEEDBYTES + 2);
  shake128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  ctr = rej_uniform(a->coeffs, SHUTTLE_N, buf, buflen);

  while(ctr < SHUTTLE_N) {
    off = buflen % 2;
    for(i = 0; i < off; ++i)
      buf[i] = buf[buflen - off + i];

    shake128_squeezeblocks(buf + off, 1, &state);
    buflen = SHAKE128_RATE + off;
    ctr += rej_uniform(a->coeffs + ctr, SHUTTLE_N - ctr, buf, buflen);
  }
}

/*************************************************
* Name:        poly_uniform_eta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-ETA,ETA] using CBD(eta=1).
*              Each coefficient needs 2 bits: a = bit0, b = bit1, coeff = a - b.
*              Total: 256 * 2 / 8 = 64 bytes needed.
*              Uses SHAKE256(seed || nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
void poly_uniform_eta(poly *a,
                      const uint8_t seed[SHUTTLE_CRHBYTES],
                      uint16_t nonce)
{
  uint8_t buf[SHAKE256_RATE]; /* 136 bytes, >= 64 needed */
  keccak_state state;
  unsigned int i;

  /* Absorb seed || nonce (LE 2 bytes) */
  uint8_t inbuf[SHUTTLE_CRHBYTES + 2];
  memcpy(inbuf, seed, SHUTTLE_CRHBYTES);
  inbuf[SHUTTLE_CRHBYTES + 0] = (uint8_t)(nonce);
  inbuf[SHUTTLE_CRHBYTES + 1] = (uint8_t)(nonce >> 8);

  shake256_absorb_once(&state, inbuf, SHUTTLE_CRHBYTES + 2);
  shake256_squeezeblocks(buf, 1, &state);

  /* CBD(eta=1): 2 bits per coefficient, 4 coefficients per byte */
  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    uint8_t byte = buf[i];
    int32_t a0, b0, a1, b1, a2, b2, a3, b3;

    a0 = (byte >> 0) & 1;
    b0 = (byte >> 1) & 1;
    a->coeffs[4 * i + 0] = a0 - b0;

    a1 = (byte >> 2) & 1;
    b1 = (byte >> 3) & 1;
    a->coeffs[4 * i + 1] = a1 - b1;

    a2 = (byte >> 4) & 1;
    b2 = (byte >> 5) & 1;
    a->coeffs[4 * i + 2] = a2 - b2;

    a3 = (byte >> 6) & 1;
    b3 = (byte >> 7) & 1;
    a->coeffs[4 * i + 3] = a3 - b3;
  }
}

/*************************************************
* Name:        poly_challenge
*
* Description: Implementation of SampleC (NGCC-Signature Alg, sign line 181).
*              Samples polynomial with exactly TAU coefficients set to +1
*              and the remaining (N - TAU) coefficients set to 0, using the
*              SHAKE-256 output stream of seed as pseudorandom choice for
*              the Fisher-Yates shuffle.
*
*              Unlike Dilithium's challenge which places tau nonzeros drawn
*              from {-1, +1}, SHUTTLE restricts the challenge to {0, +1}.
*              Sign randomization is handled separately by the IRS routine
*              through per-monomial sign bits (irs_signs), not by the
*              challenge polynomial itself.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const uint8_t seed[]: byte array containing seed of length
*                                      CTILDEBYTES
**************************************************/
void poly_challenge(poly *c, const uint8_t seed[SHUTTLE_CTILDEBYTES]) {
  unsigned int i, b, pos;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, seed, SHUTTLE_CTILDEBYTES);
  shake256_finalize(&state);
  shake256_squeezeblocks(buf, 1, &state);

  pos = 0;

  for(i = 0; i < SHUTTLE_N; ++i)
    c->coeffs[i] = 0;

  /* Fisher-Yates: place TAU ones at fresh random positions in [0, N). */
  for(i = SHUTTLE_N - SHUTTLE_TAU; i < SHUTTLE_N; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }

      b = buf[pos++];
    } while(b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1;
  }
}

/**************************************************************/
/*********** Decomposition wrappers ***************************/
/**************************************************************/

/*************************************************
* Name:        poly_decompose
*
* Description: For all coefficients c of the input polynomial,
*              compute high and low bits c1, c0 such that
*              c mod Q = c1 * 2*ALPHA_H + c0.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_decompose(poly *a1, poly *a0, const poly *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    a1->coeffs[i] = decompose(&a0->coeffs[i], a->coeffs[i]);
}

/*************************************************
* Name:        poly_make_hint
*
* Description: Compute hint polynomial.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
**************************************************/
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1) {
  unsigned int i, s = 0;

  for(i = 0; i < SHUTTLE_N; ++i) {
    h->coeffs[i] = make_hint(a0->coeffs[i], a1->coeffs[i]);
    s += h->coeffs[i];
  }

  return s;
}

/*************************************************
* Name:        poly_use_hint
*
* Description: Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/
void poly_use_hint(poly *b, const poly *a, const poly *h) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N; ++i)
    b->coeffs[i] = use_hint(a->coeffs[i], h->coeffs[i]);
}

/**************************************************************/
/*********** Packing routines *********************************/
/**************************************************************/

/*************************************************
* Name:        polyeta_pack
*
* Description: Bit-pack polynomial with coefficients in [-ETA,ETA].
*              For eta=1: coefficients in {-1,0,1} mapped to {2,1,0}
*              via t = ETA - coeff, stored as 2 bits. 4 coefficients per byte.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            SHUTTLE_POLYETA_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyeta_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint8_t t[4];

  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    t[0] = (uint8_t)(SHUTTLE_ETA - a->coeffs[4 * i + 0]);
    t[1] = (uint8_t)(SHUTTLE_ETA - a->coeffs[4 * i + 1]);
    t[2] = (uint8_t)(SHUTTLE_ETA - a->coeffs[4 * i + 2]);
    t[3] = (uint8_t)(SHUTTLE_ETA - a->coeffs[4 * i + 3]);

    r[i] = (t[0] & 0x03)
         | ((t[1] & 0x03) << 2)
         | ((t[2] & 0x03) << 4)
         | ((t[3] & 0x03) << 6);
  }
}

/*************************************************
* Name:        polyeta_unpack
*
* Description: Unpack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyeta_unpack(poly *r, const uint8_t *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    r->coeffs[4 * i + 0] = SHUTTLE_ETA - ((a[i] >> 0) & 0x03);
    r->coeffs[4 * i + 1] = SHUTTLE_ETA - ((a[i] >> 2) & 0x03);
    r->coeffs[4 * i + 2] = SHUTTLE_ETA - ((a[i] >> 4) & 0x03);
    r->coeffs[4 * i + 3] = SHUTTLE_ETA - ((a[i] >> 6) & 0x03);
  }
}

/*************************************************
* Name:        polypk_pack
*
* Description: Bit-pack polynomial with 14-bit unsigned coefficients [0, Q-1].
*              7 bytes per 4 coefficients (14 * 4 = 56 bits = 7 bytes).
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            SHUTTLE_POLYPK_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polypk_pack(uint8_t *r, const poly *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    uint32_t c0 = (uint32_t)a->coeffs[4 * i + 0];
    uint32_t c1 = (uint32_t)a->coeffs[4 * i + 1];
    uint32_t c2 = (uint32_t)a->coeffs[4 * i + 2];
    uint32_t c3 = (uint32_t)a->coeffs[4 * i + 3];

    r[7 * i + 0] = (uint8_t)(c0);
    r[7 * i + 1] = (uint8_t)(c0 >> 8) | (uint8_t)(c1 << 6);
    r[7 * i + 2] = (uint8_t)(c1 >> 2);
    r[7 * i + 3] = (uint8_t)(c1 >> 10) | (uint8_t)(c2 << 4);
    r[7 * i + 4] = (uint8_t)(c2 >> 4);
    r[7 * i + 5] = (uint8_t)(c2 >> 12) | (uint8_t)(c3 << 2);
    r[7 * i + 6] = (uint8_t)(c3 >> 6);
  }
}

/*************************************************
* Name:        polypk_unpack
*
* Description: Unpack polynomial with 14-bit unsigned coefficients.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polypk_unpack(poly *r, const uint8_t *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    r->coeffs[4 * i + 0] = ((uint32_t)a[7 * i + 0]
                           | ((uint32_t)a[7 * i + 1] << 8)) & 0x3FFF;
    r->coeffs[4 * i + 1] = (((uint32_t)a[7 * i + 1] >> 6)
                           | ((uint32_t)a[7 * i + 2] << 2)
                           | ((uint32_t)a[7 * i + 3] << 10)) & 0x3FFF;
    r->coeffs[4 * i + 2] = (((uint32_t)a[7 * i + 3] >> 4)
                           | ((uint32_t)a[7 * i + 4] << 4)
                           | ((uint32_t)a[7 * i + 5] << 12)) & 0x3FFF;
    r->coeffs[4 * i + 3] = (((uint32_t)a[7 * i + 5] >> 2)
                           | ((uint32_t)a[7 * i + 6] << 6)) & 0x3FFF;
  }
}

/*************************************************
* Name:        polyz_pack
*
* Description: Bit-pack polynomial with signed coefficients.
*              After the alpha_1 stretch, z-slot coefficients live in
*              [-Z_BOUND, Z_BOUND] where Z_BOUND = 11*sigma + alpha_1*tau.
*              Both modes satisfy 2*Z_BOUND < 2^14 so the 14-bit offset-
*              binary encoding (t = Z_BOUND - coeff) still fits.
*              7 bytes per 4 coefficients.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            SHUTTLE_POLYZ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyz_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[4];

  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    /* Map from [-Z_BOUND, Z_BOUND] to [0, 2*Z_BOUND] */
    t[0] = (uint32_t)(SHUTTLE_Z_BOUND - a->coeffs[4 * i + 0]);
    t[1] = (uint32_t)(SHUTTLE_Z_BOUND - a->coeffs[4 * i + 1]);
    t[2] = (uint32_t)(SHUTTLE_Z_BOUND - a->coeffs[4 * i + 2]);
    t[3] = (uint32_t)(SHUTTLE_Z_BOUND - a->coeffs[4 * i + 3]);

    r[7 * i + 0] = (uint8_t)(t[0]);
    r[7 * i + 1] = (uint8_t)(t[0] >> 8) | (uint8_t)(t[1] << 6);
    r[7 * i + 2] = (uint8_t)(t[1] >> 2);
    r[7 * i + 3] = (uint8_t)(t[1] >> 10) | (uint8_t)(t[2] << 4);
    r[7 * i + 4] = (uint8_t)(t[2] >> 4);
    r[7 * i + 5] = (uint8_t)(t[2] >> 12) | (uint8_t)(t[3] << 2);
    r[7 * i + 6] = (uint8_t)(t[3] >> 6);
  }
}

/*************************************************
* Name:        polyz_unpack
*
* Description: Unpack polynomial with signed coefficients from packed
*              14-bit representation.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyz_unpack(poly *r, const uint8_t *a) {
  unsigned int i;

  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    uint32_t t0, t1, t2, t3;

    t0 = ((uint32_t)a[7 * i + 0]
        | ((uint32_t)a[7 * i + 1] << 8)) & 0x3FFF;
    t1 = (((uint32_t)a[7 * i + 1] >> 6)
        | ((uint32_t)a[7 * i + 2] << 2)
        | ((uint32_t)a[7 * i + 3] << 10)) & 0x3FFF;
    t2 = (((uint32_t)a[7 * i + 3] >> 4)
        | ((uint32_t)a[7 * i + 4] << 4)
        | ((uint32_t)a[7 * i + 5] << 12)) & 0x3FFF;
    t3 = (((uint32_t)a[7 * i + 5] >> 2)
        | ((uint32_t)a[7 * i + 6] << 6)) & 0x3FFF;

    /* Map back from [0, 2*Z_BOUND] to [-Z_BOUND, Z_BOUND] */
    r->coeffs[4 * i + 0] = SHUTTLE_Z_BOUND - (int32_t)t0;
    r->coeffs[4 * i + 1] = SHUTTLE_Z_BOUND - (int32_t)t1;
    r->coeffs[4 * i + 2] = SHUTTLE_Z_BOUND - (int32_t)t2;
    r->coeffs[4 * i + 3] = SHUTTLE_Z_BOUND - (int32_t)t3;
  }
}

/*************************************************
* Name:        polyw1_pack
*
* Description: Bit-pack polynomial w1 with W1_BITS-bit coefficients.
*              Values are in [0, W1_MAX]; W1_BITS = 6 for SHUTTLE-128
*              (values fit in [0,52]) and W1_BITS = 5 for SHUTTLE-256
*              (values fit in [0,26]).
*              6-bit encoding: 3 bytes per 4 coefficients.
*              5-bit encoding: 5 bytes per 8 coefficients.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            SHUTTLE_POLYW1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyw1_pack(uint8_t *r, const poly *a) {
  unsigned int i;

#if SHUTTLE_W1_BITS == 6
  for(i = 0; i < SHUTTLE_N / 4; ++i) {
    r[3 * i + 0] = (uint8_t)(a->coeffs[4 * i + 0])
                 | (uint8_t)(a->coeffs[4 * i + 1] << 6);
    r[3 * i + 1] = (uint8_t)(a->coeffs[4 * i + 1] >> 2)
                 | (uint8_t)(a->coeffs[4 * i + 2] << 4);
    r[3 * i + 2] = (uint8_t)(a->coeffs[4 * i + 2] >> 4)
                 | (uint8_t)(a->coeffs[4 * i + 3] << 2);
  }
#elif SHUTTLE_W1_BITS == 5
  for(i = 0; i < SHUTTLE_N / 8; ++i) {
    const int32_t *c = &a->coeffs[8 * i];
    r[5 * i + 0] = (uint8_t)(c[0])
                 | (uint8_t)(c[1] << 5);
    r[5 * i + 1] = (uint8_t)(c[1] >> 3)
                 | (uint8_t)(c[2] << 2)
                 | (uint8_t)(c[3] << 7);
    r[5 * i + 2] = (uint8_t)(c[3] >> 1)
                 | (uint8_t)(c[4] << 4);
    r[5 * i + 3] = (uint8_t)(c[4] >> 4)
                 | (uint8_t)(c[5] << 1)
                 | (uint8_t)(c[6] << 6);
    r[5 * i + 4] = (uint8_t)(c[6] >> 2)
                 | (uint8_t)(c[7] << 3);
  }
#else
#  error "Unsupported SHUTTLE_W1_BITS (expected 5 or 6)"
#endif
}

/*************************************************
* Name:        polyz0_pack
*
* Description: Bit-pack the first response slot z[0] (compressed by alpha_1).
*              Coefficients live in a signed range whose absolute value fits
*              in SHUTTLE_Z0_BITS-1 bits; we encode them with SHUTTLE_Z0_BITS
*              (= 11 here) signed bits via offset binary (t = 2^{Z0_BITS-1} - c).
*              11 bits/coeff packs as 11 bytes per 8 coefficients.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            SHUTTLE_POLYZ0_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyz0_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[8];
  const int32_t center = 1 << (SHUTTLE_Z0_BITS - 1);

  for(i = 0; i < SHUTTLE_N / 8; ++i) {
    unsigned int j;
    for(j = 0; j < 8; ++j)
      t[j] = (uint32_t)(center - a->coeffs[8 * i + j]);

    /* 11 bits/coeff, 88 bits / 8 bytes per 8 coefficients -> 11 bytes */
    r[11 * i +  0] = (uint8_t)(t[0]);
    r[11 * i +  1] = (uint8_t)(t[0] >> 8) | (uint8_t)(t[1] << 3);
    r[11 * i +  2] = (uint8_t)(t[1] >> 5) | (uint8_t)(t[2] << 6);
    r[11 * i +  3] = (uint8_t)(t[2] >> 2);
    r[11 * i +  4] = (uint8_t)(t[2] >> 10) | (uint8_t)(t[3] << 1);
    r[11 * i +  5] = (uint8_t)(t[3] >> 7) | (uint8_t)(t[4] << 4);
    r[11 * i +  6] = (uint8_t)(t[4] >> 4) | (uint8_t)(t[5] << 7);
    r[11 * i +  7] = (uint8_t)(t[5] >> 1);
    r[11 * i +  8] = (uint8_t)(t[5] >> 9) | (uint8_t)(t[6] << 2);
    r[11 * i +  9] = (uint8_t)(t[6] >> 6) | (uint8_t)(t[7] << 5);
    r[11 * i + 10] = (uint8_t)(t[7] >> 3);
  }
}

/*************************************************
* Name:        polyz0_unpack
*
* Description: Inverse of polyz0_pack.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyz0_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  const int32_t center = 1 << (SHUTTLE_Z0_BITS - 1);
  const uint32_t mask = ((uint32_t)1 << SHUTTLE_Z0_BITS) - 1;

  for(i = 0; i < SHUTTLE_N / 8; ++i) {
    uint32_t t[8];

    t[0] = ((uint32_t)a[11 * i + 0]
          | ((uint32_t)a[11 * i + 1] << 8)) & mask;
    t[1] = (((uint32_t)a[11 * i + 1] >> 3)
          | ((uint32_t)a[11 * i + 2] << 5)) & mask;
    t[2] = (((uint32_t)a[11 * i + 2] >> 6)
          | ((uint32_t)a[11 * i + 3] << 2)
          | ((uint32_t)a[11 * i + 4] << 10)) & mask;
    t[3] = (((uint32_t)a[11 * i + 4] >> 1)
          | ((uint32_t)a[11 * i + 5] << 7)) & mask;
    t[4] = (((uint32_t)a[11 * i + 5] >> 4)
          | ((uint32_t)a[11 * i + 6] << 4)) & mask;
    t[5] = (((uint32_t)a[11 * i + 6] >> 7)
          | ((uint32_t)a[11 * i + 7] << 1)
          | ((uint32_t)a[11 * i + 8] << 9)) & mask;
    t[6] = (((uint32_t)a[11 * i + 8] >> 2)
          | ((uint32_t)a[11 * i + 9] << 6)) & mask;
    t[7] = (((uint32_t)a[11 * i + 9] >> 5)
          | ((uint32_t)a[11 * i + 10] << 3)) & mask;

    for(unsigned int j = 0; j < 8; ++j)
      r->coeffs[8 * i + j] = center - (int32_t)t[j];
  }
}

/**************************************************************/
/*********** z[1..L] split/combine (Phase 6c) *****************/
/**************************************************************/

/*************************************************
* Name:        polyz1_split
*
* Description: Per-coefficient round-half-up split of a z[1..L] polynomial
*              into HighBits and LowBits parts. The HighBits land in the
*              mode's z1 rANS vocabulary; the LowBits fit in ALPHA_H_BITS
*              bits using the (-alpha_h/2, alpha_h/2] -> [0, alpha_h-1]
*              bijection implemented by polyz1_lo_pack.
*
*              Uses arithmetic right shift (sign-preserving) on signed int32_t
*              which all SHUTTLE reference targets satisfy (GCC, Clang).
*              Matches the semantics used by rounding.c::highbits_mod_2q.
*
* Arguments:   - int32_t *hi: output HighBits, length N
*              - int32_t *lo: output LowBits (in (-alpha_h/2, alpha_h/2]), length N
*              - const poly *a: input polynomial, coeffs in [-Z_BOUND, Z_BOUND]
**************************************************/
void polyz1_split(int32_t *hi, int32_t *lo, const poly *a) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i) {
    int32_t z = a->coeffs[i];
    int32_t h = (z + SHUTTLE_HALF_ALPHA_H) >> SHUTTLE_ALPHA_H_BITS;
    hi[i] = h;
    lo[i] = z - (h << SHUTTLE_ALPHA_H_BITS);
  }
}

/*************************************************
* Name:        polyz1_combine
*
* Description: Inverse of polyz1_split: a[i] = hi[i] * alpha_h + lo[i].
*              Because alpha_h is a power of two, equivalent to a left shift.
**************************************************/
void polyz1_combine(poly *a, const int32_t *hi, const int32_t *lo) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    a->coeffs[i] = (hi[i] << SHUTTLE_ALPHA_H_BITS) + lo[i];
}

/*************************************************
* Name:        polyz1_lo_pack / polyz1_lo_unpack
*
* Description: Bit-pack N low-coefficients at ALPHA_H_BITS bits each.
*
*              For mode-128 (ALPHA_H_BITS=7): 8 coef -> 7 bytes.
*              For mode-256 (ALPHA_H_BITS=8): trivial byte copy.
*
*              Packed representation is `lo & (alpha_h - 1)` so the wrap
*              around zero is lossless and the unpacker recovers lo by
*              treating the high half of [0, alpha_h-1] as negative.
**************************************************/
#if SHUTTLE_ALPHA_H_BITS == 7
void polyz1_lo_pack(uint8_t *r, const int32_t *lo) {
  unsigned int i;
  uint32_t c[8];
  for(i = 0; i < SHUTTLE_N / 8; ++i) {
    unsigned int j;
    for(j = 0; j < 8; ++j)
      c[j] = (uint32_t)(lo[8 * i + j] & (SHUTTLE_ALPHA_H - 1));

    r[7 * i + 0] = (uint8_t)(c[0]      | (c[1] << 7));
    r[7 * i + 1] = (uint8_t)((c[1] >> 1) | (c[2] << 6));
    r[7 * i + 2] = (uint8_t)((c[2] >> 2) | (c[3] << 5));
    r[7 * i + 3] = (uint8_t)((c[3] >> 3) | (c[4] << 4));
    r[7 * i + 4] = (uint8_t)((c[4] >> 4) | (c[5] << 3));
    r[7 * i + 5] = (uint8_t)((c[5] >> 5) | (c[6] << 2));
    r[7 * i + 6] = (uint8_t)((c[6] >> 6) | (c[7] << 1));
  }
}

void polyz1_lo_unpack(int32_t *lo, const uint8_t *r) {
  unsigned int i;
  uint32_t c[8];
  const int32_t half = SHUTTLE_HALF_ALPHA_H;
  const int32_t alpha = SHUTTLE_ALPHA_H;
  for(i = 0; i < SHUTTLE_N / 8; ++i) {
    unsigned int j;
    c[0] = ((uint32_t)r[7 * i + 0])                                 & 0x7F;
    c[1] = (((uint32_t)r[7 * i + 0] >> 7) | ((uint32_t)r[7 * i + 1] << 1)) & 0x7F;
    c[2] = (((uint32_t)r[7 * i + 1] >> 6) | ((uint32_t)r[7 * i + 2] << 2)) & 0x7F;
    c[3] = (((uint32_t)r[7 * i + 2] >> 5) | ((uint32_t)r[7 * i + 3] << 3)) & 0x7F;
    c[4] = (((uint32_t)r[7 * i + 3] >> 4) | ((uint32_t)r[7 * i + 4] << 4)) & 0x7F;
    c[5] = (((uint32_t)r[7 * i + 4] >> 3) | ((uint32_t)r[7 * i + 5] << 5)) & 0x7F;
    c[6] = (((uint32_t)r[7 * i + 5] >> 2) | ((uint32_t)r[7 * i + 6] << 6)) & 0x7F;
    c[7] = ( (uint32_t)r[7 * i + 6] >> 1)                                  & 0x7F;
    for(j = 0; j < 8; ++j) {
      int32_t v = (int32_t)c[j];
      lo[8 * i + j] = (v >= half) ? v - alpha : v;
    }
  }
}
#elif SHUTTLE_ALPHA_H_BITS == 8
void polyz1_lo_pack(uint8_t *r, const int32_t *lo) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    r[i] = (uint8_t)(lo[i] & 0xFF);
}

void polyz1_lo_unpack(int32_t *lo, const uint8_t *r) {
  unsigned int i;
  const int32_t half = SHUTTLE_HALF_ALPHA_H;
  const int32_t alpha = SHUTTLE_ALPHA_H;
  for(i = 0; i < SHUTTLE_N; ++i) {
    int32_t v = (int32_t)r[i];
    lo[i] = (v >= half) ? v - alpha : v;
  }
}
#else
#  error "polyz1_lo_pack/unpack unsupported ALPHA_H_BITS"
#endif

/**************************************************************/
/*********** mod 2q helpers (Phase 6b-1) **********************/
/**************************************************************/

/*************************************************
* Name:        poly_highbits_mod_2q
*
* Description: Apply highbits_mod_2q coefficient-wise. Input coefficients
*              must be in [0, 2q); outputs are bucket indices in [0, HINT_MAX).
**************************************************/
void poly_highbits_mod_2q(poly *out, const poly *w) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    out->coeffs[i] = highbits_mod_2q(w->coeffs[i]);
}

/*************************************************
* Name:        poly_lsb_extract
*
* Description: Extract the least significant bit of each coefficient of w
*              (expected in [0, 2q)) into a bitmap of N/8 bytes. Bit i is
*              the LSB of coefficient i. Bits are little-endian within each
*              byte (bit 0 of byte 0 = coefficient 0, bit 7 of byte 0 =
*              coefficient 7, etc.).
**************************************************/
void poly_lsb_extract(uint8_t *bitmap, const poly *w) {
  unsigned int i;
  memset(bitmap, 0, SHUTTLE_N / 8);
  for(i = 0; i < SHUTTLE_N; ++i) {
    uint8_t bit = (uint8_t)(w->coeffs[i] & 1);
    bitmap[i >> 3] |= (uint8_t)(bit << (i & 7));
  }
}

/*************************************************
* Name:        poly_lift_to_2q
*
* Description: Lift a mod q polynomial U (coefficients in [0, q)) to mod 2q
*              using the per-coefficient parity of Y0: out[j] = 2*U[j] +
*              q * (Y0[j] & 1), reduced to [0, 2q). Only the LSB of Y0[j]
*              is used; higher bits are ignored.
**************************************************/
void poly_lift_to_2q(poly *out, const poly *u_mod_q, const poly *Y0) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    out->coeffs[i] = lift_to_2q(u_mod_q->coeffs[i], Y0->coeffs[i]);
}

/*************************************************
* Name:        poly_sub_2z2_mod2q
*
* Description: In-place w <- (w - 2*z2) mod 2q, coefficient-wise. Input w
*              must already be in [0, 2q); z2 may have arbitrary int32_t
*              coefficients. Result is canonicalized to [0, 2q).
**************************************************/
void poly_sub_2z2_mod2q(poly *w, const poly *z2) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    w->coeffs[i] = reduce_mod_2q(w->coeffs[i] - 2 * z2->coeffs[i]);
}

/*************************************************
* Name:        poly_compress_y_slot0
*
* Description: CompressY on the first slot. Implements
*              out[j] = round_half_up(in[j] / alpha_1)
*              using arithmetic right shift with +alpha_1/2 bias.
*
*              Equivalent to HAETAE's decompose_z1 highbits for our
*              alpha_1 factor. Works on full int32_t range (positive and
*              negative) because arithmetic right shift on two's complement
*              performs floor division toward -infinity, which combined
*              with the +alpha_1/2 bias gives the desired round-half-up.
*
* Arguments:   - poly *out: output (compressed) polynomial
*              - const poly *in: input polynomial (unbounded int32_t)
**************************************************/
void poly_compress_y_slot0(poly *out, const poly *in) {
  unsigned int i;
  const int32_t bias = SHUTTLE_ALPHA_1 >> 1;
  for(i = 0; i < SHUTTLE_N; ++i)
    out->coeffs[i] = (in->coeffs[i] + bias) >> SHUTTLE_ALPHA_1_BITS;
}

/*************************************************
* Name:        poly_make_hint_mod2q
**************************************************/
void poly_make_hint_mod2q(poly *h, const poly *w, const poly *z2) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    h->coeffs[i] = make_hint_mod2q(w->coeffs[i], z2->coeffs[i]);
}

/*************************************************
* Name:        poly_use_hint_wh_mod2q
**************************************************/
void poly_use_hint_wh_mod2q(poly *w_h, const poly *tilde_w, const poly *h) {
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    w_h->coeffs[i] = use_hint_wh_mod2q(tilde_w->coeffs[i], h->coeffs[i]);
}

/*************************************************
* Name:        poly_recover_z2_mod2q
**************************************************/
void poly_recover_z2_mod2q(poly *z2_out,
                           const poly *w_h,
                           const poly *w_0,
                           const poly *tilde_w)
{
  unsigned int i;
  for(i = 0; i < SHUTTLE_N; ++i)
    z2_out->coeffs[i] = recover_z2_coef_mod2q(w_h->coeffs[i],
                                               w_0->coeffs[i],
                                               tilde_w->coeffs[i]);
}
