/*
 * profile_sign.c - Per-component cycle breakdown of the SHUTTLE signing flow.
 *
 * Mirrors crypto_sign_signature() from sign.c and pack_sig() from packing.c
 * inline, with cpucycles() probes around every component. Implementation-
 * agnostic: prints SHUTTLE_NAMESPACETOP so both ref/ and avx2/ builds
 * emit a self-identifying header.
 *
 * Only cycles from the FIRST successful iteration are recorded per run;
 * rejected iterations (IRS reject, norm reject, rANS OOV/overflow) loop
 * back and do not contribute to stage totals.
 *
 * The "Pack" step is split into three rows so the rANS encoding cost can
 * be read directly:
 *   6. Finalize       (post-IRS, pre-pack arithmetic)
 *   7. Pack plumbing  (header bytes, z[1..L] lo split+pack, length/pad)
 *   8. rANS encode    (Z_0 + z1 hi + hint, with per-stream sub-rows)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../params.h"
#include "../sign.h"
#include "../packing.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../rounding.h"
#include "../fips202.h"
#include "../symmetric.h"
#include "../sampler.h"
#include "../rejsample.h"
#include "../randombytes.h"
#include "../reduce.h"
#include "../shuttle_rans.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NRUNS 1000

#define SHUTTLE_STR_(x) #x
#define SHUTTLE_STR(x)  SHUTTLE_STR_(x)

/* Signature-layout offsets, identical to pack_sig's file-local defines. */
#define OFF_C_TILDE     0
#define OFF_IRS_SIGNS   (OFF_C_TILDE + SHUTTLE_CTILDEBYTES)
#define OFF_Z0_LEN      (OFF_IRS_SIGNS + SHUTTLE_IRS_SIGNBYTES)
#define OFF_Z0_DATA     (OFF_Z0_LEN + 2)
#define OFF_Z1_LO       (OFF_Z0_DATA + SHUTTLE_Z0_RANS_RESERVED_BYTES)
#define OFF_Z1_HI_LEN   (OFF_Z1_LO + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES)
#define OFF_Z1_HI_DATA  (OFF_Z1_HI_LEN + 2)
#define OFF_HINT_LEN    (OFF_Z1_HI_DATA + SHUTTLE_Z1_RANS_RESERVED_BYTES)
#define OFF_HINT_DATA   (OFF_HINT_LEN + 2)

/* ============================================================
 * Local helpers duplicated from sign.c (verbatim).
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

static void compute_mu_local(uint8_t mu[SHUTTLE_CRHBYTES],
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

static void compute_seed_c_local(uint8_t seed_c[SHUTTLE_CTILDEBYTES],
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

/* ============================================================
 * Stats helpers.
 * ============================================================ */
static int cmp_u64(const void *a, const void *b) {
  uint64_t va = *(const uint64_t *)a;
  uint64_t vb = *(const uint64_t *)b;
  if(va < vb) return -1;
  if(va > vb) return 1;
  return 0;
}
static uint64_t median_of(uint64_t *t, size_t n) {
  qsort(t, n, sizeof(uint64_t), cmp_u64);
  return t[n / 2];
}
static uint64_t average_of(const uint64_t *t, size_t n) {
  uint64_t s = 0;
  for(size_t i = 0; i < n; ++i) s += t[i];
  return s / n;
}

/* Per-component cycle buffers. */
static uint64_t t_precomp[NRUNS];
static uint64_t t_gauss[NRUNS];
static uint64_t t_commit[NRUNS];
static uint64_t t_challenge[NRUNS];
static uint64_t t_irs[NRUNS];
static uint64_t t_finalize[NRUNS];       /* CompressY z[0], norm, MakeHint */
static uint64_t t_pack_plumb[NRUNS];     /* header, lo split+pack, prefix, pad */
static uint64_t t_rans_z0[NRUNS];
static uint64_t t_rans_z1[NRUNS];
static uint64_t t_rans_hint[NRUNS];
static uint64_t t_rans_total[NRUNS];
static uint64_t t_total[NRUNS];

/* ============================================================
 * Legend printed once at startup so the cycle table is self-documenting.
 * ============================================================ */
static void print_legend(void)
{
  printf("Stage legend:\n");
  printf("  1. Precomputation\n"
         "       unpack_sk, build_sk_full (slot 0 = ALPHA_1)\n"
         "       polyvec_sq_norm -> gamma_q62, ln_M_q62 = 0\n"
         "       compute_mu (SHAKE-256 of tr || msg)\n"
         "       rhoprime = SHAKE-256(key || rnd || mu)\n"
         "       polyvec_matrix_expand + per-row poly_ntt / polyveck_ntt\n"
         "       compute_b, polyveck_ntt -> b_hat\n");
  printf("  2. Gaussian sampling (y)\n"
         "       sample_gauss_N_4x (lanes 0..3, SHUTTLE_N each)\n"
         "       sample_gauss_N_4x (lanes 4..5 + 2 idle, SHUTTLE_N each)\n"
         "       y-buf int16 -> int32 widening copy (6 polys x SHUTTLE_N)\n");
  printf("  3. Commitment (mod-2q)\n"
         "       poly_compress_y_slot0 (alpha_1 stretch on Y.vec[0])\n"
         "       compute_commitment_mod2q (NTT matvec over modulus 2q)\n");
  printf("  4. Challenge derivation\n"
         "       polyveck_highbits_mod_2q -> w_h\n"
         "       polyveck_lsb_extract     -> w_0 bitmap\n"
         "       compute_seed_c = SHAKE-256(w_h || w_0_bitmap || mu || rho)\n"
         "       poly_challenge (SampleC from seed_c -> tau monomials)\n");
  printf("  5. IRS (iterative rejection sampling)\n"
         "       irs_sign_with_signs: tau inner steps, each does\n"
         "         - approx_neg_ln over the running log-domain score\n"
         "         - IntervalParity constant-time check\n"
         "         - inner-product update on sk_stretched\n"
         "         - branchless sign decision + vector update\n");
  printf("  6. Finalize (post-IRS, pre-pack)\n"
         "       poly_compress_y_slot0 on z.vec[0]\n"
         "       polyvec_sq_norm + bound check (SHUTTLE_BS_SQ)\n"
         "       polyveck_make_hint_mod2q -> h\n");
  printf("  7. Pack plumbing (non-rANS parts of pack_sig)\n"
         "       seed_c memcpy into sig[OFF_C_TILDE]\n"
         "       irs_signs -> bit-packed bitmap\n"
         "       z_1[0] -> int32 flat buffer\n"
         "       z_1[1..L] polyz1_split + polyz1_lo_pack (L polys)\n"
         "       h -> int32 flat buffer\n"
         "       3x length prefix + zero padding to reserved budget\n");
  printf("  8. rANS encode (isolated)\n"
         "       shuttle_rans_encode_z0 (Z_0 stream, 1 poly)\n"
         "       shuttle_rans_encode_z1 (HighBits z[1..L], L*N syms)\n"
         "       shuttle_rans_encode    (hint h,         M*N syms)\n\n");
}

/* ============================================================
 * Main: replicate crypto_sign_signature + pack_sig with per-stage timers.
 * ============================================================ */
int main(void)
{
  uint8_t pk[SHUTTLE_PUBLICKEYBYTES];
  uint8_t sk[SHUTTLE_SECRETKEYBYTES];
  uint8_t sig[SHUTTLE_BYTES];
  uint8_t msg[32];
  unsigned int run;

  printf("=== %s Signing Profiler (%s) ===\n",
         CRYPTO_ALGNAME, SHUTTLE_STR(SHUTTLE_NAMESPACETOP));
  printf("Parameters: MODE=%d, N=%d, Q=%d, L=%d, M=%d, SIGMA=%d, TAU=%d\n",
         SHUTTLE_MODE, SHUTTLE_N, SHUTTLE_Q, SHUTTLE_L, SHUTTLE_M,
         SHUTTLE_SIGMA, SHUTTLE_TAU);
  printf("PK=%d bytes, SK=%d bytes, SIG=%d bytes\n",
         SHUTTLE_PUBLICKEYBYTES, SHUTTLE_SECRETKEYBYTES, SHUTTLE_BYTES);
  printf("rANS reserved budgets: Z_0=%d, z1_hi=%d, hint=%d\n",
         SHUTTLE_Z0_RANS_RESERVED_BYTES,
         SHUTTLE_Z1_RANS_RESERVED_BYTES,
         SHUTTLE_HINT_RESERVED_BYTES);
  printf("Runs: %d\n\n", NRUNS);

  print_legend();

  if(crypto_sign_keypair(pk, sk) != 0) {
    printf("keypair failed\n");
    return 1;
  }

  for(run = 0; run < NRUNS; run++) {
    randombytes(msg, 32);

    uint64_t c_start, c_pre, c_gauss, c_commit, c_challenge, c_irs;
    uint64_t c_fin, c_plumb_hdr, c_rans_z0_end, c_plumb_lo;
    uint64_t c_rans_z1_end, c_plumb_hflat, c_rans_hint_end;

    uint64_t run_gauss = 0, run_commit = 0, run_challenge = 0;
    uint64_t run_irs = 0, run_fin = 0;
    uint64_t run_plumb = 0;
    uint64_t run_rans_z0 = 0, run_rans_z1 = 0, run_rans_hint = 0;

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

    /* ---- Stage 1: Precomputation ---- */
    c_start = cpucycles();
    unpack_sk(rho, tr, key, &s, &e, sk);
    build_sk_full(&sk_full, &s, &e);
    int64_t gamma_q62 = polyvec_sq_norm(&sk_full) * ((int64_t)1 << 47);
    int64_t ln_M_q62  = 0;

    compute_mu_local(mu, tr, msg, sizeof msg);

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
    c_pre = cpucycles();

    /* ---- Outer loop: keep timing from the final successful iteration. ---- */
    for(;;) {
      /* --- Stage 2: Gaussian sampling --- */
      uint64_t s2_start = cpucycles();
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
      c_gauss = cpucycles();
      run_gauss = c_gauss - s2_start;

      /* --- Stage 3: Commitment (CompressY slot0 + mod-2q commit) --- */
      polyvec Y;
      poly_compress_y_slot0(&Y.vec[0], &y.vec[0]);
      for(i = 1; i < SHUTTLE_VECLEN; ++i) Y.vec[i] = y.vec[i];
      compute_commitment_mod2q(&comY, &a_gen_hat, A_gen_hat, &b_hat, &Y);
      c_commit = cpucycles();
      run_commit = c_commit - c_gauss;

      /* --- Stage 4: Challenge derivation --- */
      polyveck_highbits_mod_2q(&w_h, &comY);
      polyveck_lsb_extract(w_0_bitmap, &comY);
      compute_seed_c_local(seed_c, &w_h, w_0_bitmap, mu, rho);
      poly_challenge(&c_poly, seed_c);
      c_challenge = cpucycles();
      run_challenge = c_challenge - c_commit;

      /* --- Stage 5: IRS --- */
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
      c_irs = cpucycles();
      run_irs = c_irs - c_challenge;
      if(!irs_ok) continue;

      /* --- Stage 6: Finalize (CompressY z[0] + norm + MakeHint) --- */
      {
        poly z0_comp;
        poly_compress_y_slot0(&z0_comp, &z.vec[0]);
        z.vec[0] = z0_comp;
      }
      norm_sq = polyvec_sq_norm(&z);
      if(norm_sq >= (int64_t)SHUTTLE_BS_SQ) continue;

      poly z_1[1 + SHUTTLE_L];
      for(i = 0; i < 1 + SHUTTLE_L; ++i) z_1[i] = z.vec[i];
      for(i = 0; i < SHUTTLE_M; ++i) z_2.vec[i] = z.vec[1 + SHUTTLE_L + i];

      polyveck h;
      polyveck_make_hint_mod2q(&h, &comY, &z_2);
      c_fin = cpucycles();
      run_fin = c_fin - c_irs;

      /* --- Stages 7 & 8: Pack sig, inlined from packing.c ---
       *
       * Plumbing time accumulates into run_plumb; each rANS encode is
       * measured separately so the isolated rANS cost is visible.
       * If any rANS call returns OOV (rc != 0) we fall back to continue
       * just like pack_sig in sign.c.
       */

      /* 7a. header: c_tilde memcpy + irs_signs bitmap */
      uint64_t plumb0 = cpucycles();
      memcpy(&sig[OFF_C_TILDE], seed_c, SHUTTLE_CTILDEBYTES);
      memset(&sig[OFF_IRS_SIGNS], 0, SHUTTLE_IRS_SIGNBYTES);
      for(i = 0; i < SHUTTLE_TAU; ++i) {
        if(irs_signs[i] > 0)
          sig[OFF_IRS_SIGNS + (i >> 3)] |= (uint8_t)(1u << (i & 7));
      }
      /* 7b. Z_0 flat buffer */
      int32_t z0_flat[SHUTTLE_N];
      for(i = 0; i < SHUTTLE_N; ++i) z0_flat[i] = z_1[0].coeffs[i];
      c_plumb_hdr = cpucycles();
      run_plumb += c_plumb_hdr - plumb0;

      /* 8a. Z_0 rANS */
      size_t z0_rans_len = 0;
      int rc = shuttle_rans_encode_z0(&sig[OFF_Z0_DATA], &z0_rans_len,
                                       SHUTTLE_Z0_RANS_RESERVED_BYTES,
                                       z0_flat, SHUTTLE_N);
      c_rans_z0_end = cpucycles();
      run_rans_z0 = c_rans_z0_end - c_plumb_hdr;
      if(rc != 0) continue;

      /* 7c. Z_0 length prefix + pad, then z_1[1..L] split + lo pack */
      sig[OFF_Z0_LEN + 0] = (uint8_t)(z0_rans_len & 0xFF);
      sig[OFF_Z0_LEN + 1] = (uint8_t)((z0_rans_len >> 8) & 0xFF);
      if(z0_rans_len < SHUTTLE_Z0_RANS_RESERVED_BYTES)
        memset(&sig[OFF_Z0_DATA + z0_rans_len], 0,
               SHUTTLE_Z0_RANS_RESERVED_BYTES - z0_rans_len);

      int32_t z1_hi_flat[SHUTTLE_L * SHUTTLE_N];
      {
        int32_t lo_scratch[SHUTTLE_N];
        for(i = 0; i < SHUTTLE_L; ++i) {
          polyz1_split(&z1_hi_flat[i * SHUTTLE_N], lo_scratch, &z_1[1 + i]);
          polyz1_lo_pack(&sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES],
                          lo_scratch);
        }
      }
      c_plumb_lo = cpucycles();
      run_plumb += c_plumb_lo - c_rans_z0_end;

      /* 8b. z1 hi rANS */
      size_t z1_rans_len = 0;
      rc = shuttle_rans_encode_z1(&sig[OFF_Z1_HI_DATA], &z1_rans_len,
                                   SHUTTLE_Z1_RANS_RESERVED_BYTES,
                                   z1_hi_flat, SHUTTLE_L * SHUTTLE_N);
      c_rans_z1_end = cpucycles();
      run_rans_z1 = c_rans_z1_end - c_plumb_lo;
      if(rc != 0) continue;

      /* 7d. z1 hi length prefix + pad, then h flat buffer */
      sig[OFF_Z1_HI_LEN + 0] = (uint8_t)(z1_rans_len & 0xFF);
      sig[OFF_Z1_HI_LEN + 1] = (uint8_t)((z1_rans_len >> 8) & 0xFF);
      if(z1_rans_len < SHUTTLE_Z1_RANS_RESERVED_BYTES)
        memset(&sig[OFF_Z1_HI_DATA + z1_rans_len], 0,
               SHUTTLE_Z1_RANS_RESERVED_BYTES - z1_rans_len);

      int32_t h_flat[SHUTTLE_M * SHUTTLE_N];
      for(i = 0; i < SHUTTLE_M; ++i)
        for(j = 0; j < SHUTTLE_N; ++j)
          h_flat[i * SHUTTLE_N + j] = h.vec[i].coeffs[j];
      c_plumb_hflat = cpucycles();
      run_plumb += c_plumb_hflat - c_rans_z1_end;

      /* 8c. hint rANS */
      size_t hint_rans_len = 0;
      rc = shuttle_rans_encode(&sig[OFF_HINT_DATA], &hint_rans_len,
                                SHUTTLE_HINT_RESERVED_BYTES,
                                h_flat, SHUTTLE_M * SHUTTLE_N);
      c_rans_hint_end = cpucycles();
      run_rans_hint = c_rans_hint_end - c_plumb_hflat;
      if(rc != 0) continue;

      /* 7e. hint length prefix + pad */
      uint64_t plumb_end_start = cpucycles();
      sig[OFF_HINT_LEN + 0] = (uint8_t)(hint_rans_len & 0xFF);
      sig[OFF_HINT_LEN + 1] = (uint8_t)((hint_rans_len >> 8) & 0xFF);
      if(hint_rans_len < SHUTTLE_HINT_RESERVED_BYTES)
        memset(&sig[OFF_HINT_DATA + hint_rans_len], 0,
               SHUTTLE_HINT_RESERVED_BYTES - hint_rans_len);
      run_plumb += cpucycles() - plumb_end_start;

      break;
    }

    t_precomp[run]    = c_pre - c_start;
    t_gauss[run]      = run_gauss;
    t_commit[run]     = run_commit;
    t_challenge[run]  = run_challenge;
    t_irs[run]        = run_irs;
    t_finalize[run]   = run_fin;
    t_pack_plumb[run] = run_plumb;
    t_rans_z0[run]    = run_rans_z0;
    t_rans_z1[run]    = run_rans_z1;
    t_rans_hint[run]  = run_rans_hint;
    t_rans_total[run] = run_rans_z0 + run_rans_z1 + run_rans_hint;
    t_total[run]      = t_precomp[run] + t_gauss[run] + t_commit[run]
                        + t_challenge[run] + t_irs[run] + t_finalize[run]
                        + t_pack_plumb[run] + t_rans_total[run];
  }

  /* ---- Report ---- */
  uint64_t m_pre = median_of(t_precomp, NRUNS);
  uint64_t m_gau = median_of(t_gauss, NRUNS);
  uint64_t m_com = median_of(t_commit, NRUNS);
  uint64_t m_cha = median_of(t_challenge, NRUNS);
  uint64_t m_irs = median_of(t_irs, NRUNS);
  uint64_t m_fin = median_of(t_finalize, NRUNS);
  uint64_t m_pak = median_of(t_pack_plumb, NRUNS);
  uint64_t m_rz0 = median_of(t_rans_z0, NRUNS);
  uint64_t m_rz1 = median_of(t_rans_z1, NRUNS);
  uint64_t m_rhn = median_of(t_rans_hint, NRUNS);
  uint64_t m_rto = median_of(t_rans_total, NRUNS);
  uint64_t m_tot = median_of(t_total, NRUNS);

  uint64_t a_pre = average_of(t_precomp, NRUNS);
  uint64_t a_gau = average_of(t_gauss, NRUNS);
  uint64_t a_com = average_of(t_commit, NRUNS);
  uint64_t a_cha = average_of(t_challenge, NRUNS);
  uint64_t a_irs = average_of(t_irs, NRUNS);
  uint64_t a_fin = average_of(t_finalize, NRUNS);
  uint64_t a_pak = average_of(t_pack_plumb, NRUNS);
  uint64_t a_rz0 = average_of(t_rans_z0, NRUNS);
  uint64_t a_rz1 = average_of(t_rans_z1, NRUNS);
  uint64_t a_rhn = average_of(t_rans_hint, NRUNS);
  uint64_t a_rto = average_of(t_rans_total, NRUNS);
  uint64_t a_tot = average_of(t_total, NRUNS);

  uint64_t sum_med = m_pre + m_gau + m_com + m_cha + m_irs + m_fin + m_pak + m_rto;
  double denom = (double)sum_med;
  double rans_denom = (double)m_rto;

  printf("Component                             median    average pct(med)\n");
  printf("-------------------------------------------------------------------\n");
  printf("1. Precomputation                  %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_pre, (unsigned long)a_pre,
         100.0 * (double)m_pre / denom);
  printf("2. Gaussian sampling (y)           %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_gau, (unsigned long)a_gau,
         100.0 * (double)m_gau / denom);
  printf("3. Commitment (mod-2q)             %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_com, (unsigned long)a_com,
         100.0 * (double)m_com / denom);
  printf("4. Challenge derivation            %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_cha, (unsigned long)a_cha,
         100.0 * (double)m_cha / denom);
  printf("5. IRS (rejection sampling)        %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_irs, (unsigned long)a_irs,
         100.0 * (double)m_irs / denom);
  printf("6. Finalize (norm + MakeHint)      %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_fin, (unsigned long)a_fin,
         100.0 * (double)m_fin / denom);
  printf("7. Pack plumbing                   %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_pak, (unsigned long)a_pak,
         100.0 * (double)m_pak / denom);
  printf("8. rANS encode (total)             %9lu  %9lu   %5.1f%%\n",
         (unsigned long)m_rto, (unsigned long)a_rto,
         100.0 * (double)m_rto / denom);
  printf("     Z_0  (shuttle_rans_encode_z0) %9lu  %9lu   (%4.1f%% of rANS)\n",
         (unsigned long)m_rz0, (unsigned long)a_rz0,
         100.0 * (double)m_rz0 / rans_denom);
  printf("     z1_hi(shuttle_rans_encode_z1) %9lu  %9lu   (%4.1f%% of rANS)\n",
         (unsigned long)m_rz1, (unsigned long)a_rz1,
         100.0 * (double)m_rz1 / rans_denom);
  printf("     hint (shuttle_rans_encode)    %9lu  %9lu   (%4.1f%% of rANS)\n",
         (unsigned long)m_rhn, (unsigned long)a_rhn,
         100.0 * (double)m_rhn / rans_denom);
  printf("-------------------------------------------------------------------\n");
  printf("Sum of components                  %9lu  %9lu\n",
         (unsigned long)sum_med,
         (unsigned long)(a_pre + a_gau + a_com + a_cha + a_irs + a_fin
                         + a_pak + a_rto));
  printf("Total (end-to-end)                 %9lu  %9lu\n",
         (unsigned long)m_tot, (unsigned long)a_tot);

  return 0;
}
