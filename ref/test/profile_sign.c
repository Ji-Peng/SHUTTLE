/*
 * profile_sign.c - Profiling experiment for NGCC_SIGN Sign().
 *
 * Measures cycle counts of each component inside the signing flow:
 *   1. SK decode + precomputation (matrix expand, NTT, mu, rhoprime)
 *   2. Gaussian sampling of y (6 polynomials)
 *   3. Commitment computation (NTT, matrix-vector multiply, INTT)
 *   4. Challenge derivation (decompose, hash, SampleC)
 *   5. IRS (iterative rejection sampling with log-domain)
 *   6. Norm check + pack signature
 *
 * Runs NRUNS full signing iterations and reports median/average for
 * each component, plus the breakdown as percentage of total.
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
#include "cpucycles.h"
#include "speed_print.h"

#define NRUNS 1000

/* Sort uint64_t array for median computation */
static int cmp_uint64(const void *a, const void *b) {
  uint64_t va = *(const uint64_t *)a;
  uint64_t vb = *(const uint64_t *)b;
  if(va < vb) return -1;
  if(va > vb) return 1;
  return 0;
}

static uint64_t median(uint64_t *t, size_t n) {
  qsort(t, n, sizeof(uint64_t), cmp_uint64);
  return t[n / 2];
}

static uint64_t average(uint64_t *t, size_t n) {
  uint64_t s = 0;
  for(size_t i = 0; i < n; i++) s += t[i];
  return s / n;
}

/* ============================================================
 * Replicated signing logic with per-component timing.
 * This mirrors crypto_sign_signature() from sign.c but adds
 * cycle counting at each stage.
 * ============================================================ */

/* Helper: same as in sign.c */
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

static void compute_seed_c(uint8_t seed_c[NGCC_SIGN_CTILDEBYTES],
                           const polyveck *w1,
                           const uint8_t mu[NGCC_SIGN_CRHBYTES])
{
  keccak_state state;
  uint8_t w1_packed[NGCC_SIGN_M * NGCC_SIGN_POLYW1_PACKEDBYTES];
  polyveck_pack_w1(w1_packed, w1);
  shake256_init(&state);
  shake256_absorb(&state, w1_packed, NGCC_SIGN_M * NGCC_SIGN_POLYW1_PACKEDBYTES);
  shake256_absorb(&state, mu, NGCC_SIGN_CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(seed_c, NGCC_SIGN_CTILDEBYTES, &state);
}

/* Component cycle arrays */
static uint64_t t_precomp[NRUNS];
static uint64_t t_gauss[NRUNS];
static uint64_t t_commit[NRUNS];
static uint64_t t_challenge[NRUNS];
static uint64_t t_irs[NRUNS];
static uint64_t t_norm_pack[NRUNS];
static uint64_t t_total[NRUNS];

int main(void)
{
  uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES];
  uint8_t sk[NGCC_SIGN_SECRETKEYBYTES];
  uint8_t sig[NGCC_SIGN_BYTES];
  size_t siglen;
  uint8_t msg[32];
  unsigned int run;
  uint64_t c0, c1, c2, c3, c4, c5, c6;

  printf("=== NGCC_SIGN Signing Profiler ===\n");
  printf("Parameters: N=%d, Q=%d, L=%d, M=%d, SIGMA=%d, TAU=%d\n",
         NGCC_SIGN_N, NGCC_SIGN_Q, NGCC_SIGN_L, NGCC_SIGN_M,
         NGCC_SIGN_SIGMA, NGCC_SIGN_TAU);
  printf("Runs: %d\n\n", NRUNS);

  /* Generate a key pair to use */
  crypto_sign_keypair(pk, sk);

  for(run = 0; run < NRUNS; run++) {
    /* Change message each run to get variety */
    randombytes(msg, 32);

    /* ---- Inline sign with timing ---- */
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
    unsigned int j, n;
    unsigned int nonce = 0;
    int irs_ok;
    int16_t y_buf0[NGCC_SIGN_N], y_buf1[NGCC_SIGN_N];
    int16_t y_buf2[NGCC_SIGN_N], y_buf3[NGCC_SIGN_N];

    /* === PHASE 1: Precomputation === */
    c0 = cpucycles();

    unpack_sk(rho, tr, key, &s, &e, sk);
    build_sk_full(&sk_full, &s, &e);
    norm_sq = polyvec_sq_norm(&sk_full);
    gamma_q62 = norm_sq * ((int64_t)1 << 47);
    int64_t ln_M_q62 = 0;

    compute_mu(mu, tr, msg, 32);

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
      polyveck A_gen_local[NGCC_SIGN_L];
      polyvec_matrix_expand(&a_gen, A_gen_local, rho);
      a_gen_hat_store = a_gen;
      polyveck_ntt(&a_gen_hat_store);
      for(j = 0; j < NGCC_SIGN_L; ++j) {
        A_gen_hat[j] = A_gen_local[j];
        polyveck_ntt(&A_gen_hat[j]);
      }
    }

    c1 = cpucycles();
    t_precomp[run] = c1 - c0;

    /* Signing loop (we only time one successful attempt) */
    for(;;) {
      /* === PHASE 2: Gaussian sampling === */
      c1 = cpucycles();

      sample_gauss_N_4x(y_buf0, y_buf1, y_buf2, y_buf3,
                         rhoprime, nonce, nonce+1, nonce+2, nonce+3,
                         NGCC_SIGN_N, NGCC_SIGN_N, NGCC_SIGN_N, NGCC_SIGN_N);
      nonce += 4;
      for(n = 0; n < NGCC_SIGN_N; ++n) {
        y.vec[0].coeffs[n] = (int32_t)y_buf0[n];
        y.vec[1].coeffs[n] = (int32_t)y_buf1[n];
        y.vec[2].coeffs[n] = (int32_t)y_buf2[n];
        y.vec[3].coeffs[n] = (int32_t)y_buf3[n];
      }
      sample_gauss_N_4x(y_buf0, y_buf1, y_buf2, y_buf3,
                         rhoprime, nonce, nonce+1, 0, 0,
                         NGCC_SIGN_N, NGCC_SIGN_N, 0, 0);
      nonce += 2;
      for(n = 0; n < NGCC_SIGN_N; ++n) {
        y.vec[4].coeffs[n] = (int32_t)y_buf0[n];
        y.vec[5].coeffs[n] = (int32_t)y_buf1[n];
      }

      c2 = cpucycles();
      t_gauss[run] = c2 - c1;

      /* === PHASE 3: Commitment === */
      compute_commitment(&w, &a_gen_hat_store, A_gen_hat, &y);

      c3 = cpucycles();
      t_commit[run] = c3 - c2;

      /* === PHASE 4: Challenge derivation === */
      polyveck_decompose(&w1, &w0, &w);
      compute_seed_c(seed_c, &w1, mu);
      poly_challenge(&c_poly, seed_c);

      c4 = cpucycles();
      t_challenge[run] = c4 - c3;

      /* === PHASE 5: IRS === */
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

      c5 = cpucycles();
      t_irs[run] = c5 - c4;

      if(!irs_ok) {
        /* On rejection, restart from Gaussian; re-measure only the successful attempt */
        continue;
      }

      /* === PHASE 6: Norm check + pack === */
      {
        polyvecl z1;
        unsigned int i;
        for(i = 0; i < NGCC_SIGN_L; ++i)
          z1.vec[i] = z.vec[1 + i];
        norm_sq = polyvecl_sq_norm(&z1);
        if(norm_sq >= (int64_t)NGCC_SIGN_BS_SQ)
          continue;
      }

      pack_sig(sig, seed_c, irs_signs, &z);
      siglen = NGCC_SIGN_BYTES;

      c6 = cpucycles();
      t_norm_pack[run] = c6 - c5;
      t_total[run] = c6 - c0;

      break;  /* Success — exit the signing loop */
    }
  }

  /* ---- Print results ---- */
  printf("Component breakdown (cycles, median of %d runs):\n\n", NRUNS);

  uint64_t med_precomp   = median(t_precomp, NRUNS);
  uint64_t med_gauss     = median(t_gauss, NRUNS);
  uint64_t med_commit    = median(t_commit, NRUNS);
  uint64_t med_challenge = median(t_challenge, NRUNS);
  uint64_t med_irs       = median(t_irs, NRUNS);
  uint64_t med_norm_pack = median(t_norm_pack, NRUNS);
  uint64_t med_total     = median(t_total, NRUNS);

  uint64_t avg_precomp   = average(t_precomp, NRUNS);
  uint64_t avg_gauss     = average(t_gauss, NRUNS);
  uint64_t avg_commit    = average(t_commit, NRUNS);
  uint64_t avg_challenge = average(t_challenge, NRUNS);
  uint64_t avg_irs       = average(t_irs, NRUNS);
  uint64_t avg_norm_pack = average(t_norm_pack, NRUNS);
  uint64_t avg_total     = average(t_total, NRUNS);

  /* Sum of components (median) */
  uint64_t sum_med = med_precomp + med_gauss + med_commit
                   + med_challenge + med_irs + med_norm_pack;

  printf("%-28s %10s %10s %8s\n",
         "Component", "median", "average", "pct(med)");
  printf("--------------------------------------------------------------\n");
  printf("%-28s %10lu %10lu %7.1f%%\n",
         "1. Precomputation",
         (unsigned long)med_precomp, (unsigned long)avg_precomp,
         100.0 * med_precomp / sum_med);
  printf("   (sk decode, ExpandA,\n");
  printf("    NTT, mu, rhoprime)\n");
  printf("%-28s %10lu %10lu %7.1f%%\n",
         "2. Gaussian sampling (y)",
         (unsigned long)med_gauss, (unsigned long)avg_gauss,
         100.0 * med_gauss / sum_med);
  printf("   (sample_gauss_N_4x x2)\n");
  printf("%-28s %10lu %10lu %7.1f%%\n",
         "3. Commitment (B*y)",
         (unsigned long)med_commit, (unsigned long)avg_commit,
         100.0 * med_commit / sum_med);
  printf("   (NTT, mat-vec mul, INTT)\n");
  printf("%-28s %10lu %10lu %7.1f%%\n",
         "4. Challenge derivation",
         (unsigned long)med_challenge, (unsigned long)avg_challenge,
         100.0 * med_challenge / sum_med);
  printf("   (decompose, hash, SampleC)\n");
  printf("%-28s %10lu %10lu %7.1f%%\n",
         "5. IRS (rejection sampling)",
         (unsigned long)med_irs, (unsigned long)avg_irs,
         100.0 * med_irs / sum_med);
  printf("   (30 steps: ApproxLn,\n");
  printf("    IntervalParity, inner prod)\n");
  printf("%-28s %10lu %10lu %7.1f%%\n",
         "6. Norm check + pack sig",
         (unsigned long)med_norm_pack, (unsigned long)avg_norm_pack,
         100.0 * med_norm_pack / sum_med);
  printf("--------------------------------------------------------------\n");
  printf("%-28s %10lu %10lu\n",
         "Sum of components",
         (unsigned long)sum_med, (unsigned long)(avg_precomp + avg_gauss +
         avg_commit + avg_challenge + avg_irs + avg_norm_pack));
  printf("%-28s %10lu %10lu\n",
         "Total (end-to-end)",
         (unsigned long)med_total, (unsigned long)avg_total);

  printf("\n");

  /* Sub-component micro-benchmarks.
   * For very fast operations, batch BATCH_INNER iterations per timing sample
   * to amortize rdtsc/lfence overhead. */
  printf("=== Individual operation benchmarks ===\n\n");

  #define OP_NRUNS 10000
  #define BATCH_INNER 100
  uint64_t t_op[OP_NRUNS];

  #define CYC cpucycles

  /* NTT */
  {
    poly a;
    for(unsigned int i = 0; i < NGCC_SIGN_N; i++)
      a.coeffs[i] = i % NGCC_SIGN_Q;
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      poly tmp = a;
      uint64_t s = CYC();
      for(int b = 0; b < BATCH_INNER; b++) { poly_ntt(&tmp); poly_invntt_tomont(&tmp); }
      t_op[r] = (CYC() - s) / BATCH_INNER / 2;
    }
    print_results("poly_ntt (approx):", t_op, OP_NRUNS);
  }

  /* Pointwise mul */
  {
    poly a, b, c;
    for(unsigned int i = 0; i < NGCC_SIGN_N; i++) {
      a.coeffs[i] = i % NGCC_SIGN_Q;
      b.coeffs[i] = (i * 7) % NGCC_SIGN_Q;
    }
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      uint64_t s = CYC();
      for(int bt = 0; bt < BATCH_INNER; bt++) poly_pointwise_montgomery(&c, &a, &b);
      t_op[r] = (CYC() - s) / BATCH_INNER;
    }
    print_results("poly_pointwise_mont:", t_op, OP_NRUNS);
  }

  /* sample_gauss_N (single, N=256) */
  {
    int16_t buf[NGCC_SIGN_N];
    uint8_t seed[NGCC_SIGN_SEEDBYTES];
    randombytes(seed, NGCC_SIGN_SEEDBYTES);
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      uint64_t s = CYC();
      sample_gauss_N(buf, seed, r, NGCC_SIGN_N);
      t_op[r] = CYC() - s;
    }
    print_results("sample_gauss_N (N=256):", t_op, OP_NRUNS);
  }

  /* sample_gauss_N_4x (4 x N=256) */
  {
    int16_t b0[NGCC_SIGN_N], b1[NGCC_SIGN_N];
    int16_t b2[NGCC_SIGN_N], b3[NGCC_SIGN_N];
    uint8_t seed[NGCC_SIGN_SEEDBYTES];
    randombytes(seed, NGCC_SIGN_SEEDBYTES);
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      uint64_t s = CYC();
      sample_gauss_N_4x(b0, b1, b2, b3, seed,
                         4*r, 4*r+1, 4*r+2, 4*r+3,
                         NGCC_SIGN_N, NGCC_SIGN_N,
                         NGCC_SIGN_N, NGCC_SIGN_N);
      t_op[r] = CYC() - s;
    }
    print_results("sample_gauss_N_4x (4xN):", t_op, OP_NRUNS);
  }

  /* ExpandA (matrix expansion) */
  {
    polyveck a_gen;
    polyveck A_gen_local[NGCC_SIGN_L];
    uint8_t rho_tmp[NGCC_SIGN_SEEDBYTES];
    randombytes(rho_tmp, NGCC_SIGN_SEEDBYTES);
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      uint64_t s = CYC();
      polyvec_matrix_expand(&a_gen, A_gen_local, rho_tmp);
      t_op[r] = CYC() - s;
    }
    print_results("ExpandA (matrix expand):", t_op, OP_NRUNS);
  }

  /* ApproxNegLn */
  {
    extern int64_t approx_neg_ln(uint64_t);
    volatile int64_t sink;
    uint64_t rv = 0x123456789ABCDEF0ULL;
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      uint64_t s = CYC();
      for(int b = 0; b < BATCH_INNER; b++) sink = approx_neg_ln(rv ^ (uint64_t)(r*BATCH_INNER+b));
      t_op[r] = (CYC() - s) / BATCH_INNER;
    }
    (void)sink;
    print_results("approx_neg_ln (single):", t_op, OP_NRUNS);
  }

  /* poly_uniform (single poly from SHAKE128) */
  {
    poly a;
    uint8_t seed[NGCC_SIGN_SEEDBYTES];
    randombytes(seed, NGCC_SIGN_SEEDBYTES);
    for(unsigned int r = 0; r < OP_NRUNS; r++) {
      uint64_t s = CYC();
      poly_uniform(&a, seed, (uint16_t)r);
      t_op[r] = CYC() - s;
    }
    print_results("poly_uniform (SHAKE128):", t_op, OP_NRUNS);
  }

  return 0;
}
