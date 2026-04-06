/*
 * sign.c - Key generation, signing, and verification for NGCC_SIGN.
 *
 * NGCC_SIGN lattice-based signature scheme overview:
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
 * Helper: Build the full secret key vector [1, s, e].
 * ============================================================ */
static void build_sk_full(polyvec *sk_full,
                          const polyvecl *s,
                          const polyveck *e)
{
  unsigned int i, j;

  for(j = 0; j < NGCC_SIGN_N; ++j)
    sk_full->vec[0].coeffs[j] = 0;
  sk_full->vec[0].coeffs[0] = 1;

  for(i = 0; i < NGCC_SIGN_L; ++i)
    sk_full->vec[1 + i] = s->vec[i];

  for(i = 0; i < NGCC_SIGN_M; ++i)
    sk_full->vec[1 + NGCC_SIGN_L + i] = e->vec[i];
}

/* ============================================================
 * Helper: Compute commitment w = INTT(B_hat * NTT(v)).
 * ============================================================ */
static void compute_commitment(polyveck *w,
                               const polyveck *a_gen_hat,
                               const polyveck A_gen_hat[NGCC_SIGN_L],
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
static void compute_mu(uint8_t mu[NGCC_SIGN_CRHBYTES],
                       const uint8_t tr[NGCC_SIGN_TRBYTES],
                       const uint8_t *m, size_t mlen)
{
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, tr, NGCC_SIGN_TRBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, NGCC_SIGN_CRHBYTES, &state);
}

/* ============================================================
 * Helper: Compute seed_c = SHAKE-256(w1_packed || mu, CTILDEBYTES).
 * ============================================================ */
static void compute_seed_c(uint8_t seed_c[NGCC_SIGN_CTILDEBYTES],
                           const polyveck *w1,
                           const uint8_t mu[NGCC_SIGN_CRHBYTES])
{
  keccak_state state;
  uint8_t w1_packed[NGCC_SIGN_M * NGCC_SIGN_POLYW1_PACKEDBYTES];

  polyveck_pack_w1(w1_packed, w1);

  shake256_init(&state);
  shake256_absorb(&state, w1_packed,
                  NGCC_SIGN_M * NGCC_SIGN_POLYW1_PACKEDBYTES);
  shake256_absorb(&state, mu, NGCC_SIGN_CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(seed_c, NGCC_SIGN_CTILDEBYTES, &state);
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
                        const uint8_t seed_c[NGCC_SIGN_CTILDEBYTES],
                        const int8_t irs_signs[NGCC_SIGN_TAU])
{
  poly c;
  unsigned int i, j;

  /* Reconstruct the challenge c from seed_c */
  poly_challenge(&c, seed_c);

  /* Zero out c_eff */
  for(i = 0; i < NGCC_SIGN_N; ++i)
    c_eff->coeffs[i] = 0;

  /* Walk through c, replacing each nonzero with the IRS sign */
  j = 0;
  for(i = 0; i < NGCC_SIGN_N && j < NGCC_SIGN_TAU; ++i) {
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
  uint8_t seedbuf[2 * NGCC_SIGN_SEEDBYTES + NGCC_SIGN_CRHBYTES];
  uint8_t tr[NGCC_SIGN_TRBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyveck a_gen, b;
  polyveck A_gen[NGCC_SIGN_L];
  polyvecl s;
  polyveck e;
  int64_t norm_sq;

  for(;;) {
    randombytes(seedbuf, NGCC_SIGN_SEEDBYTES);
    seedbuf[NGCC_SIGN_SEEDBYTES] = NGCC_SIGN_L;
    seedbuf[NGCC_SIGN_SEEDBYTES + 1] = NGCC_SIGN_M;
    shake256(seedbuf,
             2 * NGCC_SIGN_SEEDBYTES + NGCC_SIGN_CRHBYTES,
             seedbuf,
             NGCC_SIGN_SEEDBYTES + 2);
    rho      = seedbuf;
    rhoprime = rho + NGCC_SIGN_SEEDBYTES;
    key      = rhoprime + NGCC_SIGN_CRHBYTES;

    polyvec_matrix_expand(&a_gen, A_gen, rho);
    polyvecl_uniform_eta(&s, rhoprime, 0);
    polyveck_uniform_eta(&e, rhoprime, NGCC_SIGN_L);

    norm_sq = 1;
    norm_sq += polyvecl_sq_norm(&s);
    norm_sq += polyveck_sq_norm(&e);

    if(norm_sq < (int64_t)NGCC_SIGN_BK * NGCC_SIGN_BK
       && norm_sq > (int64_t)NGCC_SIGN_BK_LOW * NGCC_SIGN_BK_LOW)
      break;
  }

  {
    polyveck a_gen_hat;
    polyveck A_gen_hat[NGCC_SIGN_L];
    polyvec sk_full;
    unsigned int j;

    a_gen_hat = a_gen;
    polyveck_ntt(&a_gen_hat);
    for(j = 0; j < NGCC_SIGN_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }

    build_sk_full(&sk_full, &s, &e);
    compute_commitment(&b, &a_gen_hat, A_gen_hat, &sk_full);
  }

  pack_pk(pk, rho, &b);
  shake256(tr, NGCC_SIGN_TRBYTES, pk, NGCC_SIGN_PUBLICKEYBYTES);
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
  uint8_t rho[NGCC_SIGN_SEEDBYTES];
  uint8_t tr[NGCC_SIGN_TRBYTES];
  uint8_t key[NGCC_SIGN_SEEDBYTES];
  uint8_t mu[NGCC_SIGN_CRHBYTES];
  uint8_t rhoprime[NGCC_SIGN_CRHBYTES];
  uint8_t seed_c[NGCC_SIGN_CTILDEBYTES];
  uint8_t rnd[NGCC_SIGN_RNDBYTES];
  polyvecl s;
  polyveck e;
  polyveck a_gen_hat_store;
  polyveck A_gen_hat[NGCC_SIGN_L];
  polyvec sk_full, y, z;
  polyveck w, w1, w0;
  poly c_poly;
  int8_t irs_signs[NGCC_SIGN_TAU];
  int64_t gamma_q62, norm_sq;
  unsigned int i, j, n;
  unsigned int nonce = 0;
  int irs_ok;
  int16_t y_buf0[NGCC_SIGN_N], y_buf1[NGCC_SIGN_N];
  int16_t y_buf2[NGCC_SIGN_N], y_buf3[NGCC_SIGN_N];

  unpack_sk(rho, tr, key, &s, &e, sk);
  build_sk_full(&sk_full, &s, &e);

  norm_sq = polyvec_sq_norm(&sk_full);
  gamma_q62 = norm_sq * ((int64_t)1 << 47);
  int64_t ln_M_q62 = 0;

  compute_mu(mu, tr, m, mlen);

  {
    keccak_state state;

    randombytes(rnd, NGCC_SIGN_RNDBYTES);
    shake256_init(&state);
    shake256_absorb(&state, key, NGCC_SIGN_SEEDBYTES);
    shake256_absorb(&state, rnd, NGCC_SIGN_RNDBYTES);
    shake256_absorb(&state, mu, NGCC_SIGN_CRHBYTES);
    shake256_finalize(&state);
    shake256_squeeze(rhoprime, NGCC_SIGN_CRHBYTES, &state);
  }

  {
    polyveck a_gen;
    polyveck A_gen[NGCC_SIGN_L];

    polyvec_matrix_expand(&a_gen, A_gen, rho);
    a_gen_hat_store = a_gen;
    polyveck_ntt(&a_gen_hat_store);
    for(j = 0; j < NGCC_SIGN_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
  }

  for(;;) {
    /* Sample y */
    sample_gauss_N_4x(y_buf0, y_buf1, y_buf2, y_buf3,
                       rhoprime, nonce, nonce + 1, nonce + 2, nonce + 3,
                       NGCC_SIGN_N, NGCC_SIGN_N, NGCC_SIGN_N, NGCC_SIGN_N);
    nonce += 4;

    for(n = 0; n < NGCC_SIGN_N; ++n) {
      y.vec[0].coeffs[n] = (int32_t)y_buf0[n];
      y.vec[1].coeffs[n] = (int32_t)y_buf1[n];
      y.vec[2].coeffs[n] = (int32_t)y_buf2[n];
      y.vec[3].coeffs[n] = (int32_t)y_buf3[n];
    }

    sample_gauss_N_4x(y_buf0, y_buf1, y_buf2, y_buf3,
                       rhoprime, nonce, nonce + 1, 0, 0,
                       NGCC_SIGN_N, NGCC_SIGN_N, 0, 0);
    nonce += 2;

    for(n = 0; n < NGCC_SIGN_N; ++n) {
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
      uint8_t irs_seed[NGCC_SIGN_SEEDBYTES + NGCC_SIGN_CRHBYTES +
                        NGCC_SIGN_CTILDEBYTES];

      memcpy(irs_seed, key, NGCC_SIGN_SEEDBYTES);
      memcpy(irs_seed + NGCC_SIGN_SEEDBYTES, mu, NGCC_SIGN_CRHBYTES);
      memcpy(irs_seed + NGCC_SIGN_SEEDBYTES + NGCC_SIGN_CRHBYTES,
             seed_c, NGCC_SIGN_CTILDEBYTES);
      shake256_init(&rng_state);
      shake256_absorb(&rng_state, irs_seed, sizeof(irs_seed));
      shake256_finalize(&rng_state);

      irs_ok = irs_sign_with_signs(&z, irs_signs, &c_poly, &sk_full,
                                    gamma_q62, ln_M_q62, &rng_state);
    }

    if(!irs_ok)
      continue;

    /* Norm check on z1 */
    {
      polyvecl z1;
      for(i = 0; i < NGCC_SIGN_L; ++i)
        z1.vec[i] = z.vec[1 + i];
      norm_sq = polyvecl_sq_norm(&z1);
      if(norm_sq >= (int64_t)NGCC_SIGN_BS_SQ)
        continue;
    }

    /* Pack signature */
    pack_sig(sig, seed_c, irs_signs, &z);
    *siglen = NGCC_SIGN_BYTES;

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
  uint8_t rho[NGCC_SIGN_SEEDBYTES];
  uint8_t tr[NGCC_SIGN_TRBYTES];
  uint8_t mu[NGCC_SIGN_CRHBYTES];
  uint8_t seed_c[NGCC_SIGN_CTILDEBYTES];
  uint8_t seed_c_prime[NGCC_SIGN_CTILDEBYTES];
  int8_t irs_signs[NGCC_SIGN_TAU];
  polyveck b;
  polyveck a_gen_hat;
  polyveck A_gen_hat[NGCC_SIGN_L];
  polyvec z;
  polyveck w, w1, w0;
  poly c_eff, c_eff_hat;
  unsigned int i, j;

  if(siglen != NGCC_SIGN_BYTES)
    return -1;

  unpack_pk(rho, &b, pk);

  if(unpack_sig(seed_c, irs_signs, &z, sig))
    return -1;

  /* Norm check on z1 */
  {
    polyvecl z1;
    for(i = 0; i < NGCC_SIGN_L; ++i)
      z1.vec[i] = z.vec[1 + i];
    if(polyvecl_sq_norm(&z1) > (int64_t)NGCC_SIGN_BV_SQ)
      return -1;
  }

  /* mu */
  shake256(tr, NGCC_SIGN_TRBYTES, pk, NGCC_SIGN_PUBLICKEYBYTES);
  compute_mu(mu, tr, m, mlen);

  /* Expand matrix */
  {
    polyveck a_gen;
    polyveck A_gen[NGCC_SIGN_L];

    polyvec_matrix_expand(&a_gen, A_gen, rho);
    a_gen_hat = a_gen;
    polyveck_ntt(&a_gen_hat);
    for(j = 0; j < NGCC_SIGN_L; ++j) {
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

  for(i = 0; i < NGCC_SIGN_CTILDEBYTES; ++i) {
    if(seed_c[i] != seed_c_prime[i])
      return -1;
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

  memmove(sm + NGCC_SIGN_BYTES, m, mlen);
  ret = crypto_sign_signature(sm, &siglen, sm + NGCC_SIGN_BYTES, mlen, sk);
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

  if(smlen < NGCC_SIGN_BYTES)
    return -1;

  *mlen = smlen - NGCC_SIGN_BYTES;

  ret = crypto_sign_verify(sm, NGCC_SIGN_BYTES,
                           sm + NGCC_SIGN_BYTES, *mlen,
                           pk);

  if(ret) {
    *mlen = 0;
    return -1;
  }

  memmove(m, sm + NGCC_SIGN_BYTES, *mlen);
  return 0;
}
