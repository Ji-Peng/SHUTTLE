/*
 * sign.c - Key generation, signing, and verification for SHUTTLE.
 *
 * SHUTTLE lattice-based signature scheme overview:
 *
 * Key structure:
 *   Secret key: sk = [1, s, e] (VECLEN=6 polynomials)
 *   Public key: pk = (rho, b) where b = B * sk (mod q)
 *   B = [a_gen | A_gen | I_m], an m x VECLEN public matrix
 *
 * Signing:
 *   1. Sample Gaussian y (VECLEN=6 polynomials, sigma=128)
 *   2. Compute w = B*y, decompose to (w1, w0)
 *   3. seed_c = H(w1_packed || mu), challenge c = SampleC(seed_c)
 *   4. IRS: z = y + c_eff * sk (per-monomial signs chosen by IRS)
 *   5. Norm check: ||z1||^2 < BS_SQ
 *   6. Output (seed_c, irs_signs, z)
 *
 * Verification:
 *   1. Build c_eff from (seed_c, irs_signs)
 *   2. Compute w = B*z - c_eff * b
 *   3. Check H(HighBits(w) || mu) == seed_c and ||z1||^2 <= BV_SQ
 */

#include <stdint.h>
#include <string.h>
#include "sign.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "rounding.h"
#include "fips202.h"
#include "symmetric.h"
#include "sampler.h"
#include "rejsample.h"
#include "randombytes.h"
#include "reduce.h"

/* ============================================================
 * Helper: Build the full secret key vector [alpha_1, s, e].
 *
 * The constant slot holds alpha_1 (the "StretchS" operation from
 * NGCC-Signature). The scaling is absorbed symmetrically by signer and
 * verifier: b = B * [alpha_1, s, e] at keygen, and the IRS update
 * z[0] = y[0] + c_eff * alpha_1 matches the B*z - c_eff*b algebra used
 * at verify time. polyz_pack's bound was widened to SHUTTLE_Z_BOUND to
 * accommodate the expanded |z[0]| range.
 *
 * Lossy CompressY / decomposed polyz0 packing + MakeHint/UseHint come
 * in a follow-up phase (hint + rANS) to actually shrink the signature.
 * ============================================================ */
static void build_sk_full(polyvec *sk_full,
                          const polyvecl *s,
                          const polyveck *e)
{
  unsigned int i, j;

  for(j = 0; j < SHUTTLE_N; ++j)
    sk_full->vec[0].coeffs[j] = 0;
  sk_full->vec[0].coeffs[0] = SHUTTLE_ALPHA_1;

  for(i = 0; i < SHUTTLE_L; ++i)
    sk_full->vec[1 + i] = s->vec[i];

  for(i = 0; i < SHUTTLE_M; ++i)
    sk_full->vec[1 + SHUTTLE_L + i] = e->vec[i];
}

/* ============================================================
 * Helper: Compute commitment w = INTT(B_hat * NTT(v)).
 * ============================================================ */
static void compute_commitment(polyveck *w,
                               const polyveck *a_gen_hat,
                               const polyveck A_gen_hat[SHUTTLE_L],
                               const polyvec *v)
{
  polyvec v_hat;

  v_hat = *v;
  polyvec_ntt(&v_hat);
  polyvec_matrix_pointwise_montgomery(w, a_gen_hat, A_gen_hat, &v_hat);
  polyveck_invntt_tomont(w);
  polyveck_reduce(w);
  polyveck_caddq(w);
}

/* ============================================================
 * Helper: Compute mu = SHAKE-256(tr || m, 64 bytes).
 * ============================================================ */
static void compute_mu(uint8_t mu[SHUTTLE_CRHBYTES],
                       const uint8_t tr[SHUTTLE_TRBYTES],
                       const uint8_t *m, size_t mlen)
{
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, tr, SHUTTLE_TRBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, SHUTTLE_CRHBYTES, &state);
}

/* ============================================================
 * Helper: Compute seed_c = SHAKE-256(w1_packed || mu, CTILDEBYTES).
 * ============================================================ */
static void compute_seed_c(uint8_t seed_c[SHUTTLE_CTILDEBYTES],
                           const polyveck *w1,
                           const uint8_t mu[SHUTTLE_CRHBYTES])
{
  keccak_state state;
  uint8_t w1_packed[SHUTTLE_M * SHUTTLE_POLYW1_PACKEDBYTES];

  polyveck_pack_w1(w1_packed, w1);

  shake256_init(&state);
  shake256_absorb(&state, w1_packed,
                  SHUTTLE_M * SHUTTLE_POLYW1_PACKEDBYTES);
  shake256_absorb(&state, mu, SHUTTLE_CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(seed_c, SHUTTLE_CTILDEBYTES, &state);
}

/* ============================================================
 * Helper: Build the effective challenge polynomial c_eff.
 *
 * Given the challenge c (from seed_c) with TAU nonzero terms and
 * the IRS sign choices, construct c_eff where each monomial's
 * coefficient is multiplied by the corresponding IRS sign.
 *
 * c has coefficients c_i in {-1, 0, +1}.
 * irs_signs[j] is the effective sign for the j-th nonzero term.
 * c_eff = sum_j irs_signs[j] * x^{pos_j} where pos_j is the
 * position of the j-th nonzero coefficient in c.
 * ============================================================ */
static void build_c_eff(poly *c_eff,
                        const uint8_t seed_c[SHUTTLE_CTILDEBYTES],
                        const int8_t irs_signs[SHUTTLE_TAU])
{
  poly c;
  unsigned int i, j;

  /* Reconstruct the challenge c from seed_c */
  poly_challenge(&c, seed_c);

  /* Zero out c_eff */
  for(i = 0; i < SHUTTLE_N; ++i)
    c_eff->coeffs[i] = 0;

  /* Walk through c, replacing each nonzero with the IRS sign */
  j = 0;
  for(i = 0; i < SHUTTLE_N && j < SHUTTLE_TAU; ++i) {
    if(c.coeffs[i] != 0) {
      c_eff->coeffs[i] = (int32_t)irs_signs[j];
      j++;
    }
  }
}

/*************************************************
* Name:        crypto_sign_keypair
**************************************************/
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
  /* ExpandSeeds per NGCC-Signature Alg KeyGen, lines 115-116:
   *   input  = xi || I2B(lenE, 1)         // seedBytes + 1 bytes
   *   output = 4 * seedBytes total, split as:
   *            seedA      : seedBytes
   *            seedsk     : 2 * seedBytes
   *            masterSeed : seedBytes
   */
  uint8_t inbuf[SHUTTLE_SEEDBYTES + 1];
  uint8_t outbuf[4 * SHUTTLE_SEEDBYTES];
  uint8_t tr[SHUTTLE_TRBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyveck a_gen, b;
  polyveck A_gen[SHUTTLE_L];
  polyvecl s;
  polyveck e;
  int64_t norm_sq;

  for(;;) {
    randombytes(inbuf, SHUTTLE_SEEDBYTES);
    inbuf[SHUTTLE_SEEDBYTES] = (uint8_t)SHUTTLE_M;      /* I2B(lenE, 1) */

    shake256(outbuf, 4 * SHUTTLE_SEEDBYTES,
             inbuf, SHUTTLE_SEEDBYTES + 1);

    rho      = outbuf;                                  /* seedA      : seedBytes */
    rhoprime = rho + SHUTTLE_SEEDBYTES;                 /* seedsk     : 2*seedBytes */
    key      = rhoprime + SHUTTLE_SKSEEDBYTES;          /* masterSeed : seedBytes */

    polyvec_matrix_expand(&a_gen, A_gen, rho);
    polyvecl_uniform_eta(&s, rhoprime, 0);
    polyveck_uniform_eta(&e, rhoprime, SHUTTLE_L);

    norm_sq = 1;
    norm_sq += polyvecl_sq_norm(&s);
    norm_sq += polyveck_sq_norm(&e);

    if(norm_sq < (int64_t)SHUTTLE_BK * SHUTTLE_BK
       && norm_sq > (int64_t)SHUTTLE_BK_LOW * SHUTTLE_BK_LOW)
      break;
  }

  {
    polyveck a_gen_hat;
    polyveck A_gen_hat[SHUTTLE_L];
    polyvec sk_full;
    unsigned int j;

    a_gen_hat = a_gen;
    polyveck_ntt(&a_gen_hat);
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }

    build_sk_full(&sk_full, &s, &e);
    compute_commitment(&b, &a_gen_hat, A_gen_hat, &sk_full);
  }

  pack_pk(pk, rho, &b);
  shake256(tr, SHUTTLE_TRBYTES, pk, SHUTTLE_PUBLICKEYBYTES);
  pack_sk(sk, rho, tr, key, &s, &e);

  return 0;
}

/*************************************************
* Name:        crypto_sign_signature
**************************************************/
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk)
{
  uint8_t rho[SHUTTLE_SEEDBYTES];
  uint8_t tr[SHUTTLE_TRBYTES];
  uint8_t key[SHUTTLE_SEEDBYTES];
  uint8_t mu[SHUTTLE_CRHBYTES];
  uint8_t rhoprime[SHUTTLE_CRHBYTES];
  uint8_t seed_c[SHUTTLE_CTILDEBYTES];
  uint8_t rnd[SHUTTLE_RNDBYTES];
  polyvecl s;
  polyveck e;
  polyveck a_gen_hat_store;
  polyveck A_gen_hat[SHUTTLE_L];
  polyvec sk_full, y, z;
  polyveck w, w1, w0;
  poly c_poly;
  int8_t irs_signs[SHUTTLE_TAU];
  int64_t gamma_q62, norm_sq;
  unsigned int i, j, n;
  unsigned int nonce = 0;
  int irs_ok;
  int16_t y_buf0[SHUTTLE_N], y_buf1[SHUTTLE_N];
  int16_t y_buf2[SHUTTLE_N], y_buf3[SHUTTLE_N];

  unpack_sk(rho, tr, key, &s, &e, sk);
  build_sk_full(&sk_full, &s, &e);

  norm_sq = polyvec_sq_norm(&sk_full);
  gamma_q62 = norm_sq * ((int64_t)1 << 47);
  int64_t ln_M_q62 = 0;

  compute_mu(mu, tr, m, mlen);

  {
    keccak_state state;

    randombytes(rnd, SHUTTLE_RNDBYTES);
    shake256_init(&state);
    shake256_absorb(&state, key, SHUTTLE_SEEDBYTES);
    shake256_absorb(&state, rnd, SHUTTLE_RNDBYTES);
    shake256_absorb(&state, mu, SHUTTLE_CRHBYTES);
    shake256_finalize(&state);
    shake256_squeeze(rhoprime, SHUTTLE_CRHBYTES, &state);
  }

  {
    polyveck a_gen;
    polyveck A_gen[SHUTTLE_L];

    polyvec_matrix_expand(&a_gen, A_gen, rho);
    a_gen_hat_store = a_gen;
    polyveck_ntt(&a_gen_hat_store);
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
  }

  for(;;) {
    /* Sample y */
    sample_gauss_N_4x(y_buf0, y_buf1, y_buf2, y_buf3,
                       rhoprime, nonce, nonce + 1, nonce + 2, nonce + 3,
                       SHUTTLE_N, SHUTTLE_N, SHUTTLE_N, SHUTTLE_N);
    nonce += 4;

    for(n = 0; n < SHUTTLE_N; ++n) {
      y.vec[0].coeffs[n] = (int32_t)y_buf0[n];
      y.vec[1].coeffs[n] = (int32_t)y_buf1[n];
      y.vec[2].coeffs[n] = (int32_t)y_buf2[n];
      y.vec[3].coeffs[n] = (int32_t)y_buf3[n];
    }

    sample_gauss_N_4x(y_buf0, y_buf1, y_buf2, y_buf3,
                       rhoprime, nonce, nonce + 1, 0, 0,
                       SHUTTLE_N, SHUTTLE_N, 0, 0);
    nonce += 2;

    for(n = 0; n < SHUTTLE_N; ++n) {
      y.vec[4].coeffs[n] = (int32_t)y_buf0[n];
      y.vec[5].coeffs[n] = (int32_t)y_buf1[n];
    }

    /* Commitment */
    compute_commitment(&w, &a_gen_hat_store, A_gen_hat, &y);
    polyveck_decompose(&w1, &w0, &w);
    compute_seed_c(seed_c, &w1, mu);
    poly_challenge(&c_poly, seed_c);

    /* IRS with sign output */
    z = y;
    {
      keccak_state rng_state;
      uint8_t irs_seed[SHUTTLE_SEEDBYTES + SHUTTLE_CRHBYTES +
                        SHUTTLE_CTILDEBYTES];

      memcpy(irs_seed, key, SHUTTLE_SEEDBYTES);
      memcpy(irs_seed + SHUTTLE_SEEDBYTES, mu, SHUTTLE_CRHBYTES);
      memcpy(irs_seed + SHUTTLE_SEEDBYTES + SHUTTLE_CRHBYTES,
             seed_c, SHUTTLE_CTILDEBYTES);
      shake256_init(&rng_state);
      shake256_absorb(&rng_state, irs_seed, sizeof(irs_seed));
      shake256_finalize(&rng_state);

      irs_ok = irs_sign_with_signs(&z, irs_signs, &c_poly, &sk_full,
                                    gamma_q62, ln_M_q62, &rng_state);
    }

    if(!irs_ok)
      continue;

    /* Norm check: full ||z||_2 < B_s per NGCC-Signature Alg Sign, line 187.
     *
     * z is split as (z_1, z_2) per Sign line 191, where z_1 is the first
     * (L+1) polynomials and z_2 is the subsequent M polynomials. Both
     * halves participate in the rejection bound. */
    norm_sq = polyvec_sq_norm(&z);
    if(norm_sq >= (int64_t)SHUTTLE_BS_SQ)
      continue;

    /* Pack signature */
    pack_sig(sig, seed_c, irs_signs, &z);
    *siglen = SHUTTLE_BYTES;

    return 0;
  }
}

/*************************************************
* Name:        crypto_sign_verify
**************************************************/
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk)
{
  uint8_t rho[SHUTTLE_SEEDBYTES];
  uint8_t tr[SHUTTLE_TRBYTES];
  uint8_t mu[SHUTTLE_CRHBYTES];
  uint8_t seed_c[SHUTTLE_CTILDEBYTES];
  uint8_t seed_c_prime[SHUTTLE_CTILDEBYTES];
  int8_t irs_signs[SHUTTLE_TAU];
  polyveck b;
  polyveck a_gen_hat;
  polyveck A_gen_hat[SHUTTLE_L];
  polyvec z;
  polyveck w, w1, w0;
  poly c_eff, c_eff_hat;
  unsigned int i, j;

  if(siglen != SHUTTLE_BYTES)
    return -1;

  unpack_pk(rho, &b, pk);

  if(unpack_sig(seed_c, irs_signs, &z, sig))
    return -1;

  /* Verify-side norm check: full ||(z_1, z_2)||_2 <= B_v per NGCC-Signature
   * Alg Verify, line 246. With the current (uncompressed) signature format
   * the entire z vector is published so polyvec_sq_norm gives ||(z_1, z_2)||^2
   * directly. */
  if(polyvec_sq_norm(&z) > (int64_t)SHUTTLE_BV_SQ)
    return -2;

  /* mu */
  shake256(tr, SHUTTLE_TRBYTES, pk, SHUTTLE_PUBLICKEYBYTES);
  compute_mu(mu, tr, m, mlen);

  /* Expand matrix */
  {
    polyveck a_gen;
    polyveck A_gen[SHUTTLE_L];

    polyvec_matrix_expand(&a_gen, A_gen, rho);
    a_gen_hat = a_gen;
    polyveck_ntt(&a_gen_hat);
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
  }

  /* Compute w = B*z */
  compute_commitment(&w, &a_gen_hat, A_gen_hat, &z);

  /* Build c_eff from seed_c and irs_signs, then subtract c_eff * b */
  build_c_eff(&c_eff, seed_c, irs_signs);

  {
    polyveck b_hat, ct;

    c_eff_hat = c_eff;
    poly_ntt(&c_eff_hat);

    b_hat = b;
    polyveck_ntt(&b_hat);

    polyveck_pointwise_poly_montgomery(&ct, &c_eff_hat, &b_hat);
    polyveck_invntt_tomont(&ct);

    polyveck_sub(&w, &w, &ct);
  }

  polyveck_reduce(&w);
  polyveck_caddq(&w);

  /* Decompose and check seed */
  polyveck_decompose(&w1, &w0, &w);
  compute_seed_c(seed_c_prime, &w1, mu);

  for(i = 0; i < SHUTTLE_CTILDEBYTES; ++i) {
    if(seed_c[i] != seed_c_prime[i])
      return -3;
  }

  return 0;
}

/*************************************************
* Name:        crypto_sign
**************************************************/
int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk)
{
  size_t siglen;
  int ret;

  memmove(sm + SHUTTLE_BYTES, m, mlen);
  ret = crypto_sign_signature(sm, &siglen, sm + SHUTTLE_BYTES, mlen, sk);
  *smlen = siglen + mlen;

  return ret;
}

/*************************************************
* Name:        crypto_sign_open
**************************************************/
int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk)
{
  int ret;

  if(smlen < SHUTTLE_BYTES)
    return -1;

  *mlen = smlen - SHUTTLE_BYTES;

  ret = crypto_sign_verify(sm, SHUTTLE_BYTES,
                           sm + SHUTTLE_BYTES, *mlen,
                           pk);

  if(ret) {
    *mlen = 0;
    return -1;
  }

  memmove(m, sm + SHUTTLE_BYTES, *mlen);
  return 0;
}
