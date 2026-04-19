/*
 * sign.c - NGCC-Signature Alg 2 sign / verify.
 *
 * Mod 2q commitment + CompressY on slot 0 (alpha_1 stretch) + MakeHint /
 * UseHint for the z_2 portion. The on-wire signature is produced by
 * pack_sig (packing.c), which entropy-codes Z_0, HighBits(z[1..L]) and
 * hint h with three independent rANS tables.
 */

#include <stdint.h>
#include <stddef.h>
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
 * Helpers
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

static void pack_w_h_for_hash(uint8_t out[SHUTTLE_M * SHUTTLE_N],
                              const polyveck *w_h)
{
  unsigned int i, j;
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      out[i * SHUTTLE_N + j] = (uint8_t)w_h->vec[i].coeffs[j];
}

static void compute_seed_c(uint8_t seed_c[SHUTTLE_CTILDEBYTES],
                              const polyveck *w_h,
                              const uint8_t w_0_bitmap[SHUTTLE_M * SHUTTLE_N / 8],
                              const uint8_t mu[SHUTTLE_CRHBYTES],
                              const uint8_t rho[SHUTTLE_SEEDBYTES])
{
  keccak_state state;
  uint8_t w_h_packed[SHUTTLE_M * SHUTTLE_N];
  pack_w_h_for_hash(w_h_packed, w_h);

  shake256_init(&state);
  shake256_absorb(&state, w_h_packed, sizeof w_h_packed);
  shake256_absorb(&state, w_0_bitmap, SHUTTLE_M * SHUTTLE_N / 8);
  shake256_absorb(&state, mu, SHUTTLE_CRHBYTES);
  shake256_absorb(&state, rho, SHUTTLE_SEEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(seed_c, SHUTTLE_CTILDEBYTES, &state);
}

static void build_c_eff(poly *c_eff,
                           const uint8_t seed_c[SHUTTLE_CTILDEBYTES],
                           const int8_t irs_signs[SHUTTLE_TAU])
{
  poly c;
  unsigned int i, j;
  poly_challenge(&c, seed_c);
  for(i = 0; i < SHUTTLE_N; ++i) c_eff->coeffs[i] = 0;
  j = 0;
  for(i = 0; i < SHUTTLE_N && j < SHUTTLE_TAU; ++i) {
    if(c.coeffs[i] != 0) {
      c_eff->coeffs[i] = (int32_t)irs_signs[j];
      j++;
    }
  }
}

/* ============================================================
 * Keypair.
 * ============================================================ */
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
  uint8_t inbuf[SHUTTLE_SEEDBYTES + 1];
  uint8_t outbuf[4 * SHUTTLE_SEEDBYTES];
  uint8_t tr[SHUTTLE_TRBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyveck a_gen, b, e;
  polyveck A_gen[SHUTTLE_L];
  polyvecl s;
  int64_t norm_sq;

  for(;;) {
    randombytes(inbuf, SHUTTLE_SEEDBYTES);
    inbuf[SHUTTLE_SEEDBYTES] = (uint8_t)SHUTTLE_M;

    shake256(outbuf, 4 * SHUTTLE_SEEDBYTES,
             inbuf, SHUTTLE_SEEDBYTES + 1);

    rho      = outbuf;
    rhoprime = rho + SHUTTLE_SEEDBYTES;
    key      = rhoprime + SHUTTLE_SKSEEDBYTES;

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
    polyveck A_gen_hat[SHUTTLE_L];
    unsigned int j;
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
    compute_b(&b, &a_gen, A_gen_hat, &s, &e);
  }

  pack_pk(pk, rho, &b);
  shake256(tr, SHUTTLE_TRBYTES, pk, SHUTTLE_PUBLICKEYBYTES);
  pack_sk(sk, rho, tr, key, &s, &e);

  return 0;
}

/* ============================================================
 * Sign: Alg 2 flow + rANS-compressed packing.
 * ============================================================ */
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                             const uint8_t *m, size_t mlen,
                             const uint8_t *sk)
{
  uint8_t rho[SHUTTLE_SEEDBYTES];
  uint8_t tr[SHUTTLE_TRBYTES];
  uint8_t key[SHUTTLE_SEEDBYTES];
  uint8_t mu[SHUTTLE_CRHBYTES];
  uint8_t rhoprime[SHUTTLE_CRHBYTES];
  uint8_t rnd[SHUTTLE_RNDBYTES];
  uint8_t seed_c[SHUTTLE_CTILDEBYTES];
  polyvecl s;
  polyveck e;
  polyveck a_gen_hat, b_hat, b;
  polyveck A_gen_hat[SHUTTLE_L];
  polyvec sk_full, y, z;
  polyveck comY, w_h, z_2;
  poly c_poly;
  int8_t irs_signs[SHUTTLE_TAU];
  int64_t norm_sq;
  unsigned int i, j, n;
  unsigned int nonce = 0;
  int irs_ok;
  int16_t y_buf0[SHUTTLE_N], y_buf1[SHUTTLE_N];
  int16_t y_buf2[SHUTTLE_N], y_buf3[SHUTTLE_N];
  uint8_t w_0_bitmap[SHUTTLE_M * SHUTTLE_N / 8];

  unpack_sk(rho, tr, key, &s, &e, sk);
  build_sk_full(&sk_full, &s, &e);

  int64_t gamma_q62 = polyvec_sq_norm(&sk_full) * ((int64_t)1 << 47);
  int64_t ln_M_q62  = 0;

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
    polyveck a_gen, A_gen_raw[SHUTTLE_L];
    polyvec_matrix_expand(&a_gen, A_gen_raw, rho);
    for(i = 0; i < SHUTTLE_M; ++i) {
      a_gen_hat.vec[i] = a_gen.vec[i];
      poly_ntt(&a_gen_hat.vec[i]);
    }
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen_raw[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
    compute_b(&b, &a_gen, A_gen_hat, &s, &e);
    b_hat = b;
    polyveck_ntt(&b_hat);
  }

  for(;;) {
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

    polyvec Y;
    poly_compress_y_slot0(&Y.vec[0], &y.vec[0]);
    for(i = 1; i < SHUTTLE_VECLEN; ++i) Y.vec[i] = y.vec[i];

    compute_commitment_mod2q(&comY, &a_gen_hat, A_gen_hat, &b_hat, &Y);

    polyveck_highbits_mod_2q(&w_h, &comY);
    polyveck_lsb_extract(w_0_bitmap, &comY);

    compute_seed_c(seed_c, &w_h, w_0_bitmap, mu, rho);
    poly_challenge(&c_poly, seed_c);

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
    if(!irs_ok) continue;

    poly z0_comp;
    poly_compress_y_slot0(&z0_comp, &z.vec[0]);
    z.vec[0] = z0_comp;

    norm_sq = polyvec_sq_norm(&z);
    if(norm_sq >= (int64_t)SHUTTLE_BS_SQ)
      continue;

    poly z_1[1 + SHUTTLE_L];
    for(i = 0; i < 1 + SHUTTLE_L; ++i)
      z_1[i] = z.vec[i];
    for(i = 0; i < SHUTTLE_M; ++i)
      z_2.vec[i] = z.vec[1 + SHUTTLE_L + i];

    polyveck h;
    polyveck_make_hint_mod2q(&h, &comY, &z_2);

    /* Reject on rANS OOV / overflow (signer restarts IRS with fresh y). */
    int rc = pack_sig(sig, seed_c, irs_signs, z_1, &h);
    if(rc != 0) continue;

    *siglen = SHUTTLE_BYTES;
    return 0;
  }
}

/* ============================================================
 * Verify: Alg 3 flow.
 * ============================================================ */
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *pk)
{
  uint8_t rho[SHUTTLE_SEEDBYTES];
  uint8_t tr[SHUTTLE_TRBYTES];
  uint8_t mu[SHUTTLE_CRHBYTES];
  uint8_t seed_c[SHUTTLE_CTILDEBYTES];
  uint8_t seed_c_prime[SHUTTLE_CTILDEBYTES];
  uint8_t w_0_prime_bitmap[SHUTTLE_M * SHUTTLE_N / 8];
  int8_t irs_signs[SHUTTLE_TAU];
  polyveck b, a_gen_hat, b_hat;
  polyveck A_gen_hat[SHUTTLE_L];
  poly z_1[1 + SHUTTLE_L];
  polyveck h;
  polyveck tilde_w, w_h;
  polyveck w_0_prime_poly;
  poly c_eff;
  unsigned int i, j;

  if(siglen != SHUTTLE_BYTES)
    return -1;

  unpack_pk(rho, &b, pk);
  if(unpack_sig(seed_c, irs_signs, z_1, &h, sig))
    return -1;

  shake256(tr, SHUTTLE_TRBYTES, pk, SHUTTLE_PUBLICKEYBYTES);
  compute_mu(mu, tr, m, mlen);

  {
    polyveck a_gen, A_gen_raw[SHUTTLE_L];
    polyvec_matrix_expand(&a_gen, A_gen_raw, rho);
    for(i = 0; i < SHUTTLE_M; ++i) {
      a_gen_hat.vec[i] = a_gen.vec[i];
      poly_ntt(&a_gen_hat.vec[i]);
    }
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen_raw[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
    b_hat = b;
    polyveck_ntt(&b_hat);
  }

  build_c_eff(&c_eff, seed_c, irs_signs);

  compute_tilde_w_mod2q(&tilde_w, &a_gen_hat, A_gen_hat, &b_hat, z_1, &c_eff);

  polyveck_lsb_extract(w_0_prime_bitmap, &tilde_w);
  polyveck_lsb_from_bitmap(&w_0_prime_poly, w_0_prime_bitmap);

  polyveck_use_hint_wh_mod2q(&w_h, &tilde_w, &h);

  compute_seed_c(seed_c_prime, &w_h, w_0_prime_bitmap, mu, rho);
  for(i = 0; i < SHUTTLE_CTILDEBYTES; ++i)
    if(seed_c[i] != seed_c_prime[i])
      return -2;

  polyveck z_2_recov;
  polyveck_recover_z2_mod2q(&z_2_recov, &w_h, &w_0_prime_poly, &tilde_w);

  int64_t norm_sq = 0;
  for(i = 0; i < 1 + SHUTTLE_L; ++i)
    norm_sq += poly_sq_norm(&z_1[i]);
  for(i = 0; i < SHUTTLE_M; ++i)
    norm_sq += poly_sq_norm(&z_2_recov.vec[i]);

  if(norm_sq > (int64_t)SHUTTLE_BV_SQ)
    return -3;

  return 0;
}

/* ============================================================
 * Combined API.
 * ============================================================ */
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

int crypto_sign_open(uint8_t *m, size_t *mlen,
                        const uint8_t *sm, size_t smlen,
                        const uint8_t *pk)
{
  if(smlen < SHUTTLE_BYTES) return -1;
  *mlen = smlen - SHUTTLE_BYTES;
  int ret = crypto_sign_verify(sm, SHUTTLE_BYTES,
                                   sm + SHUTTLE_BYTES, *mlen, pk);
  if(ret) { *mlen = 0; return -1; }
  memmove(m, sm + SHUTTLE_BYTES, *mlen);
  return 0;
}
