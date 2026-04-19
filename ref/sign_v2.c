/*
 * sign_v2.c - NGCC-Signature Alg 2 sign/verify (Phase 6b-5c).
 *
 * Parallel to sign.c (Phase 6a). Uses mod 2q commitment (via
 * compute_commitment_mod2q), CompressY on the first slot of y and z,
 * MakeHint/UseHint for the z_2 portion, and pack_sig_v2 for the signature
 * (rANS-compressed hint).
 *
 * NOTE ON KEYS
 * ============
 * v2 keys have a DIFFERENT b value than Phase 6a: v2 stores
 *   b = a_gen + A_gen*s + e (mod q)             [spec Alg KeyGen line 124]
 * whereas Phase 6a stored
 *   b = alpha_1*a_gen + A_gen*s + e + (MONT-artifact)
 *
 * Concretely: keys generated via crypto_sign_keypair_v2 are VERIFY-compatible
 * only with crypto_sign_signature_v2 / crypto_sign_verify_v2. Phase 6a keys
 * are verify-compatible only with Phase 6a sign/verify.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#include "sign_v2.h"
#include "packing.h"
#include "packing_v2.h"
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

/* sk_full_v2 = [alpha_1, s, e] (stretched; only used inside IRS). */
static void build_sk_full_v2(polyvec *sk_full,
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

/* mu = SHAKE256(tr || m), CRHBYTES output. */
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

/* Pack HighBits(comY) for seedC hashing. 1 byte per coefficient (HINT_MAX < 256). */
static void pack_w_h_for_hash(uint8_t out[SHUTTLE_M * SHUTTLE_N],
                              const polyveck *w_h)
{
  unsigned int i, j;
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      out[i * SHUTTLE_N + j] = (uint8_t)w_h->vec[i].coeffs[j];
}

/* seedC = SHAKE256(HighBits || LSB || mu || rho). Alg 2 L178. */
static void compute_seed_c_v2(uint8_t seed_c[SHUTTLE_CTILDEBYTES],
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

/* Build c_eff from c_poly (the +1-only challenge) and per-monomial IRS signs. */
static void build_c_eff_v2(poly *c_eff,
                           const uint8_t seed_c[SHUTTLE_CTILDEBYTES],
                           const int8_t irs_signs[SHUTTLE_TAU])
{
  poly c;
  unsigned int i, j;
  poly_challenge(&c, seed_c);           /* c has +1 at TAU positions */
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
 * Keypair (v2): b = a_gen + A_gen*s + e (mod q), unstretched.
 * Otherwise matches Phase 6a structure.
 * ============================================================ */
int crypto_sign_keypair_v2(uint8_t *pk, uint8_t *sk)
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

  /* b = a_gen + A_gen*s + e (mod q), cleanly computed. */
  {
    polyveck A_gen_hat[SHUTTLE_L];
    unsigned int j;
    for(j = 0; j < SHUTTLE_L; ++j) {
      A_gen_hat[j] = A_gen[j];
      polyveck_ntt(&A_gen_hat[j]);
    }
    compute_b_v2(&b, &a_gen, A_gen_hat, &s, &e);
  }

  pack_pk(pk, rho, &b);
  shake256(tr, SHUTTLE_TRBYTES, pk, SHUTTLE_PUBLICKEYBYTES);
  pack_sk(sk, rho, tr, key, &s, &e);

  return 0;
}

/* ============================================================
 * Sign (v2): Alg 2 flow.
 * ============================================================ */
int crypto_sign_signature_v2(uint8_t *sig, size_t *siglen,
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
  build_sk_full_v2(&sk_full, &s, &e);

  /* Precompute sk_full norm (for IRS internals). */
  {
    int64_t ns = polyvec_sq_norm(&sk_full);
    (void)ns;     /* gamma_q62 computed below */
  }
  int64_t gamma_q62 = polyvec_sq_norm(&sk_full) * ((int64_t)1 << 47);
  int64_t ln_M_q62  = 0;

  compute_mu(mu, tr, m, mlen);

  /* rhoprime derivation (same scheme as Phase 6a). */
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

  /* Expand public matrix, compute b = a + A*s + e (needed as b_hat by the
   * commitment machinery), then NTT-transform. */
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
    /* b = a + A*s + e (clean, via compute_b_v2). */
    compute_b_v2(&b, &a_gen, A_gen_hat, &s, &e);
    b_hat = b;
    polyveck_ntt(&b_hat);
  }

  /* Outer rejection loop. */
  for(;;) {
    /* Sample y (VECLEN = 1 + L + M polynomials from Gaussian). */
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

    /* CompressY on slot 0: build the vector Y = (Y_0, y[1..L+M]) for commit. */
    polyvec Y;
    poly_compress_y_slot0(&Y.vec[0], &y.vec[0]);
    for(i = 1; i < SHUTTLE_VECLEN; ++i) Y.vec[i] = y.vec[i];

    /* Compute comY = hat_A * Y (mod 2q). */
    compute_commitment_mod2q(&comY, &a_gen_hat, A_gen_hat, &b_hat, &Y);

    /* HighBits and LSB. */
    polyveck_highbits_mod_2q(&w_h, &comY);
    polyveck_lsb_extract(w_0_bitmap, &comY);

    /* seedC = H(w_h || w_0 || mu || rho). */
    compute_seed_c_v2(seed_c, &w_h, w_0_bitmap, mu, rho);
    poly_challenge(&c_poly, seed_c);

    /* IRS on y with stretched sk (alpha_1 in slot 0). */
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

    /* CompressY on IRS output's slot 0:  z[0] = round(y'[0]/alpha_1). */
    poly z0_comp;
    poly_compress_y_slot0(&z0_comp, &z.vec[0]);
    z.vec[0] = z0_comp;

    /* Norm check: full ||z||^2 < B_s^2. */
    norm_sq = polyvec_sq_norm(&z);
    if(norm_sq >= (int64_t)SHUTTLE_BS_SQ)
      continue;

    /* Split (z_1, z_2): z_1 = first 1+L polys, z_2 = last M polys. */
    poly z_1[1 + SHUTTLE_L];
    for(i = 0; i < 1 + SHUTTLE_L; ++i)
      z_1[i] = z.vec[i];
    for(i = 0; i < SHUTTLE_M; ++i)
      z_2.vec[i] = z.vec[1 + SHUTTLE_L + i];

    /* Build c_eff for MakeHint: MakeHint uses tilde_w = comY - 2*z_2, and
     * z_2 already carries the c_eff-scaled e contribution (via IRS). So we
     * just call polyveck_make_hint_mod2q directly with comY and z_2. */
    polyveck h;
    polyveck_make_hint_mod2q(&h, &comY, &z_2);

    /* Pack signature; reject if rANS overflows or OOV. */
    int rc = pack_sig_v2(sig, seed_c, irs_signs, z_1, &h);
    if(rc != 0) continue;

    *siglen = SHUTTLE_BYTES_V2;
    return 0;
  }
}

/* ============================================================
 * Verify (v2): Alg 3 flow.
 * ============================================================ */
int crypto_sign_verify_v2(const uint8_t *sig, size_t siglen,
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

  if(siglen != SHUTTLE_BYTES_V2)
    return -1;

  unpack_pk(rho, &b, pk);
  if(unpack_sig_v2(seed_c, irs_signs, z_1, &h, sig))
    return -1;

  /* mu = H(tr || m), tr = H(pk). */
  shake256(tr, SHUTTLE_TRBYTES, pk, SHUTTLE_PUBLICKEYBYTES);
  compute_mu(mu, tr, m, mlen);

  /* Expand public matrix and NTT-transform a_gen, A_gen, b. */
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

  /* Build c_eff = SampleC(seed_c) with IRS signs applied. */
  build_c_eff_v2(&c_eff, seed_c, irs_signs);

  /* tilde_w = hat_A_1 * z_1 - q * c_eff * j (mod 2q). */
  compute_tilde_w_mod2q(&tilde_w, &a_gen_hat, A_gen_hat, &b_hat, z_1, &c_eff);

  /* Extract LSB bitmap and expanded LSB polyveck (for z_2 recovery). */
  polyveck_lsb_extract(w_0_prime_bitmap, &tilde_w);
  polyveck_lsb_from_bitmap(&w_0_prime_poly, w_0_prime_bitmap);

  /* Recover w_h = UseHint(tilde_w, h). */
  polyveck_use_hint_wh_mod2q(&w_h, &tilde_w, &h);

  /* seedC' = H(w_h || w_0' || mu || rho); compare with seed_c. */
  compute_seed_c_v2(seed_c_prime, &w_h, w_0_prime_bitmap, mu, rho);
  for(i = 0; i < SHUTTLE_CTILDEBYTES; ++i)
    if(seed_c[i] != seed_c_prime[i])
      return -2;

  /* Recover z_2 and do verification norm check ||(z_1, z_2)||^2 <= B_v^2. */
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
int crypto_sign_v2(uint8_t *sm, size_t *smlen,
                   const uint8_t *m, size_t mlen,
                   const uint8_t *sk)
{
  size_t siglen;
  int ret;
  memmove(sm + SHUTTLE_BYTES_V2, m, mlen);
  ret = crypto_sign_signature_v2(sm, &siglen, sm + SHUTTLE_BYTES_V2, mlen, sk);
  *smlen = siglen + mlen;
  return ret;
}

int crypto_sign_open_v2(uint8_t *m, size_t *mlen,
                        const uint8_t *sm, size_t smlen,
                        const uint8_t *pk)
{
  if(smlen < SHUTTLE_BYTES_V2) return -1;
  *mlen = smlen - SHUTTLE_BYTES_V2;
  int ret = crypto_sign_verify_v2(sm, SHUTTLE_BYTES_V2,
                                   sm + SHUTTLE_BYTES_V2, *mlen, pk);
  if(ret) { *mlen = 0; return -1; }
  memmove(m, sm + SHUTTLE_BYTES_V2, *mlen);
  return 0;
}
