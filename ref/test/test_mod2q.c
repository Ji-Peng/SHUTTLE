/*
 * test_mod2q.c - unit tests for Phase 6b-1 mod 2q infrastructure.
 *
 * Tests:
 *   1. Scalar lift_to_2q correctness (exhaustive over (u, parity)).
 *   2. Scalar highbits_mod_2q boundary behavior (round + wrap).
 *   3. Per-coefficient poly_lift_to_2q / poly_sub_2z2_mod2q round-trip.
 *   4. Signer / verifier algebraic identity:
 *      hat_A_1 . z_1 - q*c*j == comY - 2*z_2   (mod 2q)
 *      computed using schoolbook polynomial multiplication as ground truth.
 *   5. LSB invariants:
 *      LSB(comY) == Y_0 mod 2    (per coefficient)
 *      LSB(tilde_w) == LSB(comY) (since 2*z_2 is even)
 *
 * These tests DO NOT exercise NTT-based fast multiplication; they use an
 * O(N^2) schoolbook multiply so that the mod 2q math is validated in
 * isolation from any NTT bugs. Phase 6b-2 will re-verify with the NTT path.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../params.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../ntt.h"
#include "../rounding.h"
#include "../reduce.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

/* ------------------------------------------------------------
 * Helpers: direct (schoolbook) polynomial arithmetic in Z_q[x]/(x^N + 1).
 * O(N^2); used as ground truth, NOT performance-critical.
 * ------------------------------------------------------------ */

/* Reduce arbitrary int32_t to canonical [0, q). */
static int32_t cmod_q(int64_t a) {
  int64_t r = a % SHUTTLE_Q;
  if (r < 0) r += SHUTTLE_Q;
  return (int32_t)r;
}

/* c <- a * b in R_q (negacyclic convolution). */
static void poly_mul_school_modq(poly *c, const poly *a, const poly *b) {
  int64_t tmp[SHUTTLE_N];
  unsigned int i, j;

  for (i = 0; i < SHUTTLE_N; ++i) tmp[i] = 0;

  for (i = 0; i < SHUTTLE_N; ++i) {
    for (j = 0; j < SHUTTLE_N; ++j) {
      int64_t prod = (int64_t)a->coeffs[i] * b->coeffs[j];
      unsigned int k = i + j;
      if (k < SHUTTLE_N) {
        tmp[k] += prod;
      } else {
        tmp[k - SHUTTLE_N] -= prod;        /* x^N = -1 */
      }
    }
    /* Periodically reduce to avoid int64 overflow on bigger N. */
    if ((i & 0x1F) == 0x1F) {
      for (unsigned int k = 0; k < SHUTTLE_N; ++k) tmp[k] = tmp[k] % SHUTTLE_Q;
    }
  }
  for (i = 0; i < SHUTTLE_N; ++i) c->coeffs[i] = cmod_q(tmp[i]);
}

/* c <- a + b in R_q. */
static void poly_add_modq(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  for (i = 0; i < SHUTTLE_N; ++i)
    c->coeffs[i] = cmod_q((int64_t)a->coeffs[i] + b->coeffs[i]);
}

/* c <- a - b in R_q. */
static void poly_sub_modq(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  for (i = 0; i < SHUTTLE_N; ++i)
    c->coeffs[i] = cmod_q((int64_t)a->coeffs[i] - b->coeffs[i]);
}

/* Generate a random poly with coefficients in [0, q). */
static void poly_random_modq(poly *a) {
  unsigned int i;
  for (i = 0; i < SHUTTLE_N; ++i)
    a->coeffs[i] = rand() % SHUTTLE_Q;
}

/* Generate a small poly with coefficients in [-bound, bound]. */
static void poly_random_small(poly *a, int32_t bound) {
  unsigned int i;
  for (i = 0; i < SHUTTLE_N; ++i) {
    int32_t x = (rand() % (2 * bound + 1)) - bound;
    a->coeffs[i] = x;
  }
}

/* Make a sparse challenge poly: tau positions = +1, rest = 0. */
static void poly_random_challenge(poly *c) {
  unsigned int i, idx;
  memset(c->coeffs, 0, sizeof c->coeffs);
  for (i = 0; i < SHUTTLE_TAU; ++i) {
    do { idx = rand() % SHUTTLE_N; } while (c->coeffs[idx] != 0);
    c->coeffs[idx] = 1;
  }
}

/* Compare two polyvecks mod 2q. */
static int polyveck_eq_mod_2q(const polyveck *a, const polyveck *b) {
  unsigned int i, j;
  for (i = 0; i < SHUTTLE_M; ++i) {
    for (j = 0; j < SHUTTLE_N; ++j) {
      int32_t ai = reduce_mod_2q(a->vec[i].coeffs[j]);
      int32_t bi = reduce_mod_2q(b->vec[i].coeffs[j]);
      if (ai != bi) {
        printf("    mismatch at [%u][%u]: %d vs %d\n", i, j, ai, bi);
        return 0;
      }
    }
  }
  return 1;
}

/* ------------------------------------------------------------
 * Test 1: scalar lift_to_2q identity.
 * ------------------------------------------------------------ */
static int test_scalar_lift(void) {
  printf("Test 1: lift_to_2q identity for u in [0,q), parity in {0,1}\n");
  /* Sample across the range. */
  int32_t test_us[] = {0, 1, 17, SHUTTLE_Q/2, SHUTTLE_Q-1};
  for (unsigned k = 0; k < sizeof test_us / sizeof test_us[0]; ++k) {
    int32_t u = test_us[k];
    for (int p = 0; p <= 1; ++p) {
      int32_t got = lift_to_2q(u, p);
      int32_t want = (2 * u + SHUTTLE_Q * p) % SHUTTLE_DQ;
      if (got != want)
        FAIL("lift_to_2q(%d,%d) = %d, want %d", u, p, got, want);
    }
  }
  /* Add 200 random cases. */
  for (int k = 0; k < 200; ++k) {
    int32_t u = rand() % SHUTTLE_Q;
    int p = rand() & 1;
    int32_t got = lift_to_2q(u, p);
    int32_t want = (2 * u + SHUTTLE_Q * p) % SHUTTLE_DQ;
    if (got != want) FAIL("random lift_to_2q(%d,%d) = %d, want %d", u, p, got, want);
  }
  PASS("scalar lift_to_2q matches (2*u + q*p) mod 2q");
  return 0;
}

/* ------------------------------------------------------------
 * Test 2: scalar highbits_mod_2q boundary behavior.
 * ------------------------------------------------------------ */
static int test_scalar_highbits(void) {
  printf("Test 2: highbits_mod_2q round + wrap\n");

  /* w = 0 -> hb = 0 */
  if (highbits_mod_2q(0) != 0) FAIL("highbits(0) != 0, got %d", highbits_mod_2q(0));

  /* w = alpha_h/2 - 1 -> hb = 0 (round down) */
  if (highbits_mod_2q(SHUTTLE_ALPHA_H/2 - 1) != 0)
    FAIL("highbits(alpha_h/2 - 1) != 0, got %d", highbits_mod_2q(SHUTTLE_ALPHA_H/2 - 1));

  /* w = alpha_h/2 -> hb = 1 (round up) */
  if (highbits_mod_2q(SHUTTLE_ALPHA_H/2) != 1)
    FAIL("highbits(alpha_h/2) != 1, got %d", highbits_mod_2q(SHUTTLE_ALPHA_H/2));

  /* w = alpha_h -> hb = 1 */
  if (highbits_mod_2q(SHUTTLE_ALPHA_H) != 1)
    FAIL("highbits(alpha_h) != 1, got %d", highbits_mod_2q(SHUTTLE_ALPHA_H));

  /* w = HINT_MOD - 1 = 2q - 3 -> hb might wrap */
  {
    int32_t w = SHUTTLE_HINT_MOD - 1;
    int32_t hb = highbits_mod_2q(w);
    int32_t expected = ((w + SHUTTLE_HALF_ALPHA_H) / SHUTTLE_ALPHA_H) % SHUTTLE_HINT_MAX;
    if (hb != expected)
      FAIL("highbits(2q-3) = %d, expected %d", hb, expected);
  }

  /* Sweep a few random w in [0, 2q). */
  for (int k = 0; k < 500; ++k) {
    int32_t w = rand() % SHUTTLE_DQ;
    int32_t hb = highbits_mod_2q(w);
    int32_t expected = (w + SHUTTLE_HALF_ALPHA_H) / SHUTTLE_ALPHA_H;
    if (expected >= SHUTTLE_HINT_MAX) expected -= SHUTTLE_HINT_MAX;
    if (hb != expected)
      FAIL("highbits(%d) = %d, expected %d", w, hb, expected);
    if (hb < 0 || hb >= SHUTTLE_HINT_MAX)
      FAIL("highbits(%d) = %d out of [0, HINT_MAX=%d)", w, hb, SHUTTLE_HINT_MAX);
  }
  PASS("highbits_mod_2q matches round(w / alpha_h) wrapped to [0, HINT_MAX)");
  return 0;
}

/* ------------------------------------------------------------
 * Test 3: per-coefficient poly helpers.
 * ------------------------------------------------------------ */
static int test_poly_helpers(void) {
  printf("Test 3: poly_lift / poly_sub_2z2 round-trip\n");

  poly u_mod_q, Y0, lifted, z2, w_sub;
  poly_random_modq(&u_mod_q);
  poly_random_modq(&Y0);

  poly_lift_to_2q(&lifted, &u_mod_q, &Y0);
  /* Validate each coeff. */
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    int32_t want = (2 * u_mod_q.coeffs[i] + (Y0.coeffs[i] & 1) * SHUTTLE_Q) % SHUTTLE_DQ;
    if (lifted.coeffs[i] != want)
      FAIL("lift coeff %u: got %d, want %d", i, lifted.coeffs[i], want);
    if (lifted.coeffs[i] < 0 || lifted.coeffs[i] >= SHUTTLE_DQ)
      FAIL("lift coeff %u out of [0,2q): %d", i, lifted.coeffs[i]);
  }

  /* sub_2z2: start from lifted, subtract 2*z2, check per-coeff. */
  poly_random_small(&z2, SHUTTLE_BOUND);   /* z2 small */
  w_sub = lifted;
  poly_sub_2z2_mod2q(&w_sub, &z2);
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    int64_t want = (int64_t)lifted.coeffs[i] - 2 * (int64_t)z2.coeffs[i];
    int32_t want_m = (int32_t)(want % SHUTTLE_DQ);
    if (want_m < 0) want_m += SHUTTLE_DQ;
    if (w_sub.coeffs[i] != want_m)
      FAIL("sub_2z2 coeff %u: got %d, want %d", i, w_sub.coeffs[i], want_m);
  }
  PASS("poly_lift_to_2q and poly_sub_2z2_mod2q match element-wise formulas");
  return 0;
}

/* ------------------------------------------------------------
 * Test 4: signer / verifier algebraic identity.
 * ------------------------------------------------------------ */
static int test_signer_verifier_identity(void) {
  printf("Test 4: hat_A_1 . z_1 - q*c*j == comY - 2*z_2 (mod 2q)\n");

  /* Generate key material. */
  polyveck a_gen, b, e;
  polyveck A_gen[SHUTTLE_L];
  polyvecl s;
  polyvec Y;                     /* (Y_0, y[1..L+M]) — the input to hat_A */
  poly c;                        /* challenge */
  unsigned i, j;

  /* Random a_gen. */
  for (i = 0; i < SHUTTLE_M; ++i) poly_random_modq(&a_gen.vec[i]);

  /* Random A_gen (L columns, each a polyveck). */
  for (j = 0; j < SHUTTLE_L; ++j)
    for (i = 0; i < SHUTTLE_M; ++i)
      poly_random_modq(&A_gen[j].vec[i]);

  /* Small s and e. */
  for (j = 0; j < SHUTTLE_L; ++j) poly_random_small(&s.vec[j], 1);
  for (i = 0; i < SHUTTLE_M; ++i) poly_random_small(&e.vec[i], 1);

  /* b = a + A*s + e (mod q). */
  for (i = 0; i < SHUTTLE_M; ++i) {
    poly tmp, acc;
    memset(acc.coeffs, 0, sizeof acc.coeffs);
    poly_mul_school_modq(&tmp, &A_gen[0].vec[i], &s.vec[0]);
    acc = tmp;
    for (j = 1; j < SHUTTLE_L; ++j) {
      poly_mul_school_modq(&tmp, &A_gen[j].vec[i], &s.vec[j]);
      poly_add_modq(&acc, &acc, &tmp);
    }
    poly_add_modq(&b.vec[i], &a_gen.vec[i], &acc);
    poly_add_modq(&b.vec[i], &b.vec[i], &e.vec[i]);
  }

  /* Random Y (the input to hat_A). */
  for (i = 0; i < SHUTTLE_VECLEN; ++i) poly_random_modq(&Y.vec[i]);

  /* Random challenge c. */
  poly_random_challenge(&c);

  /* --- Signer path: comY = hat_A . Y (mod 2q). --- */
  /* U_sign[i] = (a[i] - b[i]) * Y_0 + sum_j A_gen[j][i] * Y[1+j] + Y[1+L+i]  (mod q) */
  polyveck U_sign, comY;
  for (i = 0; i < SHUTTLE_M; ++i) {
    poly tmp, acc, ab;
    poly_sub_modq(&ab, &a_gen.vec[i], &b.vec[i]);
    poly_mul_school_modq(&acc, &ab, &Y.vec[0]);
    for (j = 0; j < SHUTTLE_L; ++j) {
      poly_mul_school_modq(&tmp, &A_gen[j].vec[i], &Y.vec[1 + j]);
      poly_add_modq(&acc, &acc, &tmp);
    }
    poly_add_modq(&U_sign.vec[i], &acc, &Y.vec[1 + SHUTTLE_L + i]);
  }
  /* comY = lift(U_sign, Y_0). */
  polyveck_lift_to_2q(&comY, &U_sign, &Y.vec[0]);

  /* --- Verifier: build z_1, z_2 from Y and c using idealized IRS
     (z = Y + c*sk_stretched), then compute tilde_w and compare. --- */
  /* z_1 = (Z_0, y[1..L] + c*s)   where Z_0 = Y_0 + c.
   * z_2 = (y[L+1..L+M] + c*e). */
  poly Z_0;
  poly_add_modq(&Z_0, &Y.vec[0], &c);     /* Z_0 = Y_0 + c (mod q) */

  polyvecl z1_mid;                         /* z[1..L] */
  for (j = 0; j < SHUTTLE_L; ++j) {
    poly tmp;
    poly_mul_school_modq(&tmp, &c, &s.vec[j]);
    poly_add_modq(&z1_mid.vec[j], &Y.vec[1 + j], &tmp);
  }

  polyveck z2;
  for (i = 0; i < SHUTTLE_M; ++i) {
    poly tmp;
    poly_mul_school_modq(&tmp, &c, &e.vec[i]);
    poly_add_modq(&z2.vec[i], &Y.vec[1 + SHUTTLE_L + i], &tmp);
  }

  /* --- Verifier path: tilde_w = hat_A_1 . z_1 - q*c*j (mod 2q). --- */
  /* U_ver[i] = (a[i] - b[i]) * Z_0 + sum_j A_gen[j][i] * z1_mid[j]  (mod q) */
  polyveck U_ver;
  for (i = 0; i < SHUTTLE_M; ++i) {
    poly tmp, acc, ab;
    poly_sub_modq(&ab, &a_gen.vec[i], &b.vec[i]);
    poly_mul_school_modq(&acc, &ab, &Z_0);
    for (j = 0; j < SHUTTLE_L; ++j) {
      poly_mul_school_modq(&tmp, &A_gen[j].vec[i], &z1_mid.vec[j]);
      poly_add_modq(&acc, &acc, &tmp);
    }
    U_ver.vec[i] = acc;
  }
  /* Lift parity is (Z_0 - c) mod 2 = Y_0 mod 2, same as signer. */
  polyveck tilde_w;
  polyveck_lift_to_2q(&tilde_w, &U_ver, &Y.vec[0]);

  /* Check: tilde_w == comY - 2*z_2 (mod 2q). */
  polyveck expected;
  expected = comY;
  polyveck_sub_2z2_mod2q(&expected, &z2);
  if (!polyveck_eq_mod_2q(&tilde_w, &expected))
    FAIL("signer/verifier tilde_w mismatch (see first bad coeff above)%s", "");

  PASS("tilde_w identity holds mod 2q");
  return 0;
}

/* ------------------------------------------------------------
 * Test 5: LSB invariants.
 * ------------------------------------------------------------ */
static int test_lsb_invariants(void) {
  printf("Test 5: LSB(comY) and LSB(tilde_w) invariants\n");

  /* Quick check: after lift_to_2q with arbitrary u and Y0,
     LSB of the lifted value equals LSB of Y0 (since q is odd and
     2*u is even, so r = 2*u + q*(Y0&1) has LSB = q*(Y0&1) mod 2 = Y0&1). */
  poly u, Y0, lifted;
  uint8_t bitmap[SHUTTLE_N / 8];
  poly_random_modq(&u);
  poly_random_modq(&Y0);
  poly_lift_to_2q(&lifted, &u, &Y0);
  poly_lsb_extract(bitmap, &lifted);
  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    uint8_t got = (bitmap[i >> 3] >> (i & 7)) & 1;
    uint8_t want = (uint8_t)(Y0.coeffs[i] & 1);
    if (got != want)
      FAIL("LSB mismatch at %u: got %u, want %u (Y0=%d)",
           i, got, want, Y0.coeffs[i]);
  }
  PASS("LSB(lift(u, Y0)) == Y0 mod 2 per coefficient");
  return 0;
}

/* ------------------------------------------------------------
 * Test 6: poly_compress_y_slot0.
 * ------------------------------------------------------------ */
static int test_compress_y_slot0(void) {
  printf("Test 6: poly_compress_y_slot0 (round_half_up(in / alpha_1))\n");

  poly in, out;
  /* Random small-to-medium input to exercise both signs. */
  for (unsigned i = 0; i < SHUTTLE_N; ++i)
    in.coeffs[i] = (rand() % (2 * 1500)) - 1500;

  poly_compress_y_slot0(&out, &in);

  for (unsigned i = 0; i < SHUTTLE_N; ++i) {
    /* Reference: floor((x + alpha_1/2) / alpha_1) */
    int32_t bias = SHUTTLE_ALPHA_1 / 2;
    int32_t x = in.coeffs[i] + bias;
    /* Floor division toward -infinity. */
    int32_t want;
    if (x >= 0)
      want = x / SHUTTLE_ALPHA_1;
    else
      want = -(((-x) + SHUTTLE_ALPHA_1 - 1) / SHUTTLE_ALPHA_1);

    if (out.coeffs[i] != want)
      FAIL("compress coeff %u: in=%d got=%d want=%d",
           i, in.coeffs[i], out.coeffs[i], want);
  }

  /* Spot-check known values. */
  struct { int32_t in; int32_t want; } cases[] = {
    { 0, 0 },
    { SHUTTLE_ALPHA_1 / 2 - 1, 0 },         /* below half */
    { SHUTTLE_ALPHA_1 / 2,     1 },         /* at half, round up */
    { SHUTTLE_ALPHA_1 - 1,     1 },
    { SHUTTLE_ALPHA_1,         1 },
    { -SHUTTLE_ALPHA_1 / 2,    0 },         /* -half -> round to 0 */
    { -SHUTTLE_ALPHA_1 / 2 - 1, -1 },       /* just below -half */
    { -SHUTTLE_ALPHA_1,        -1 },
    { 3 * SHUTTLE_ALPHA_1 / 2, 2 },
    { -3 * SHUTTLE_ALPHA_1 / 2, -1 },       /* half at -3/2 rounds to -1 (half-up toward +inf) */
  };
  for (unsigned k = 0; k < sizeof cases / sizeof cases[0]; ++k) {
    poly a, b;
    memset(a.coeffs, 0, sizeof a.coeffs);
    a.coeffs[0] = cases[k].in;
    poly_compress_y_slot0(&b, &a);
    if (b.coeffs[0] != cases[k].want)
      FAIL("compress spot-check in=%d got=%d want=%d",
           cases[k].in, b.coeffs[0], cases[k].want);
  }
  PASS("poly_compress_y_slot0 matches round_half_up(in / alpha_1)");
  return 0;
}

/* ------------------------------------------------------------
 * Test 7: compute_commitment_mod2q (NTT path) matches schoolbook.
 *
 * Mirrors test 4's setup but also calls compute_commitment_mod2q and
 * verifies the NTT-based result equals the schoolbook-computed comY.
 * ------------------------------------------------------------ */
static int test_commitment_ntt_matches_schoolbook(void) {
  printf("Test 7: compute_commitment_mod2q (NTT) == schoolbook comY\n");

  polyveck a_gen, b, e;
  polyveck A_gen[SHUTTLE_L];
  polyvecl s;
  polyvec v;
  unsigned i, j;

  /* Same generation path as test 4. */
  for (i = 0; i < SHUTTLE_M; ++i) poly_random_modq(&a_gen.vec[i]);
  for (j = 0; j < SHUTTLE_L; ++j)
    for (i = 0; i < SHUTTLE_M; ++i)
      poly_random_modq(&A_gen[j].vec[i]);
  for (j = 0; j < SHUTTLE_L; ++j) poly_random_small(&s.vec[j], 1);
  for (i = 0; i < SHUTTLE_M; ++i) poly_random_small(&e.vec[i], 1);

  /* b = a + A*s + e (mod q). */
  for (i = 0; i < SHUTTLE_M; ++i) {
    poly tmp, acc;
    memset(acc.coeffs, 0, sizeof acc.coeffs);
    poly_mul_school_modq(&tmp, &A_gen[0].vec[i], &s.vec[0]);
    acc = tmp;
    for (j = 1; j < SHUTTLE_L; ++j) {
      poly_mul_school_modq(&tmp, &A_gen[j].vec[i], &s.vec[j]);
      poly_add_modq(&acc, &acc, &tmp);
    }
    poly_add_modq(&b.vec[i], &a_gen.vec[i], &acc);
    poly_add_modq(&b.vec[i], &b.vec[i], &e.vec[i]);
  }

  /* Random v (VECLEN). */
  for (i = 0; i < SHUTTLE_VECLEN; ++i) poly_random_modq(&v.vec[i]);

  /* --- Schoolbook comY. --- */
  polyveck U_sign, comY_ref;
  for (i = 0; i < SHUTTLE_M; ++i) {
    poly tmp, acc, ab;
    poly_sub_modq(&ab, &a_gen.vec[i], &b.vec[i]);
    poly_mul_school_modq(&acc, &ab, &v.vec[0]);
    for (j = 0; j < SHUTTLE_L; ++j) {
      poly_mul_school_modq(&tmp, &A_gen[j].vec[i], &v.vec[1 + j]);
      poly_add_modq(&acc, &acc, &tmp);
    }
    poly_add_modq(&U_sign.vec[i], &acc, &v.vec[1 + SHUTTLE_L + i]);
  }
  polyveck_lift_to_2q(&comY_ref, &U_sign, &v.vec[0]);

  /* --- NTT-based comY via compute_commitment_mod2q. --- */
  polyveck a_gen_hat, b_hat, comY_ntt;
  polyveck A_gen_hat[SHUTTLE_L];

  /* NTT of a_gen, b, A_gen (inputs in mod q [0, q) since schoolbook uses cmod_q). */
  for (i = 0; i < SHUTTLE_M; ++i) {
    a_gen_hat.vec[i] = a_gen.vec[i];
    b_hat.vec[i]     = b.vec[i];
    poly_ntt(&a_gen_hat.vec[i]);
    poly_ntt(&b_hat.vec[i]);
  }
  for (j = 0; j < SHUTTLE_L; ++j)
    for (i = 0; i < SHUTTLE_M; ++i) {
      A_gen_hat[j].vec[i] = A_gen[j].vec[i];
      poly_ntt(&A_gen_hat[j].vec[i]);
    }

  compute_commitment_mod2q(&comY_ntt, &a_gen_hat, A_gen_hat, &b_hat, &v);

  /* Compare. */
  if (!polyveck_eq_mod_2q(&comY_ntt, &comY_ref))
    FAIL("NTT comY != schoolbook comY (see mismatch above)%s", "");
  PASS("compute_commitment_mod2q (NTT) matches schoolbook reference");
  return 0;
}

/* ------------------------------------------------------------
 * Test 8: MakeHint / UseHint w_h round-trip.
 *
 * Pick random w in [0, 2q) and small z_2. Compute h = MakeHint(z_2, w),
 * tilde_w = w - 2*z_2 (mod 2q), and verify UseHint(tilde_w, h) recovers
 * signer's HighBits(w) per coefficient.
 * ------------------------------------------------------------ */
static int test_hint_wh_roundtrip(void) {
  printf("Test 8: MakeHint / UseHint w_h round-trip\n");

  polyveck w, z2, h, tilde_w, w_h_ref, w_h_recov;
  unsigned int i, j;

  /* Random w in [0, 2q). */
  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j)
      w.vec[i].coeffs[j] = rand() % SHUTTLE_DQ;

  /* Random z_2 with small coefficients (|z_2| <= 11*sigma + tau*eta). */
  for (i = 0; i < SHUTTLE_M; ++i)
    poly_random_small(&z2.vec[i], SHUTTLE_BOUND + SHUTTLE_TAU);

  /* Signer: h = MakeHint(z_2, w). */
  polyveck_make_hint_mod2q(&h, &w, &z2);

  /* Signer: tilde_w = w - 2*z_2 (mod 2q). */
  tilde_w = w;
  polyveck_sub_2z2_mod2q(&tilde_w, &z2);

  /* Signer's reference w_h. */
  polyveck_highbits_mod_2q(&w_h_ref, &w);

  /* Verifier: recover w_h via UseHint(tilde_w, h). */
  polyveck_use_hint_wh_mod2q(&w_h_recov, &tilde_w, &h);

  /* Compare. */
  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j)
      if (w_h_recov.vec[i].coeffs[j] != w_h_ref.vec[i].coeffs[j])
        FAIL("w_h mismatch at [%u][%u]: recov=%d ref=%d (h=%d tilde_w=%d w=%d)",
             i, j, w_h_recov.vec[i].coeffs[j], w_h_ref.vec[i].coeffs[j],
             h.vec[i].coeffs[j], tilde_w.vec[i].coeffs[j], w.vec[i].coeffs[j]);

  PASS("UseHint(tilde_w, h) recovers HighBits(w) exactly");
  return 0;
}

/* ------------------------------------------------------------
 * Test 9: z_2 recovery approximate identity.
 *
 * Same setup as test 8. Additionally verify recover_z2 gives a result
 * within +/- alpha_h/2 of the true z_2 per coefficient.
 * ------------------------------------------------------------ */
static int test_z2_recovery(void) {
  printf("Test 9: z_2 recovery within +/- alpha_h/2\n");

  polyveck w, z2, h, tilde_w, w_h_recov, w_0, z2_recov;
  unsigned int i, j;

  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j)
      w.vec[i].coeffs[j] = rand() % SHUTTLE_DQ;

  for (i = 0; i < SHUTTLE_M; ++i)
    poly_random_small(&z2.vec[i], SHUTTLE_BOUND + SHUTTLE_TAU);

  polyveck_make_hint_mod2q(&h, &w, &z2);
  tilde_w = w;
  polyveck_sub_2z2_mod2q(&tilde_w, &z2);
  polyveck_use_hint_wh_mod2q(&w_h_recov, &tilde_w, &h);

  /* w_0 = LSB(w) expanded as polyveck. */
  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j)
      w_0.vec[i].coeffs[j] = w.vec[i].coeffs[j] & 1;

  /* Verifier reconstructs z_2 approximately. */
  polyveck_recover_z2_mod2q(&z2_recov, &w_h_recov, &w_0, &tilde_w);

  /* Check approximate equality. */
  int32_t max_err = 0;
  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j) {
      int32_t err = z2_recov.vec[i].coeffs[j] - z2.vec[i].coeffs[j];
      if (err < 0) err = -err;
      if (err > max_err) max_err = err;
      if (err > SHUTTLE_ALPHA_H / 2)
        FAIL("z_2 recovery error too big at [%u][%u]: recov=%d true=%d err=%d (bound alpha_h/2=%d)",
             i, j, z2_recov.vec[i].coeffs[j], z2.vec[i].coeffs[j], err, SHUTTLE_ALPHA_H / 2);
    }
  printf("  max |z2_recov - z2|_inf = %d (bound alpha_h/2 = %d)\n",
         max_err, SHUTTLE_ALPHA_H / 2);
  PASS("z_2 recovery within alpha_h/2 per coefficient");
  return 0;
}

/* ------------------------------------------------------------
 * Test 10: hint pack/unpack round-trip.
 * ------------------------------------------------------------ */
static int test_hint_pack_basic(void) {
  printf("Test 10: basic hint pack/unpack round-trip\n");

  polyveck h, h2;
  uint8_t buf[SHUTTLE_HINT_PACKEDBYTES_BASIC];
  unsigned int i, j;

  /* Hint values typically in +/- 20 range. Generate some. */
  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j) {
      int32_t v = (rand() % 41) - 20;   /* [-20, 20] */
      h.vec[i].coeffs[j] = v;
    }

  polyveck_hint_pack_basic(buf, &h);
  polyveck_hint_unpack_basic(&h2, buf);

  for (i = 0; i < SHUTTLE_M; ++i)
    for (j = 0; j < SHUTTLE_N; ++j)
      if (h.vec[i].coeffs[j] != h2.vec[i].coeffs[j])
        FAIL("hint pack round-trip mismatch at [%u][%u]: %d vs %d",
             i, j, h.vec[i].coeffs[j], h2.vec[i].coeffs[j]);

  printf("  sig_bytes (basic hint only) = %d B\n", SHUTTLE_HINT_PACKEDBYTES_BASIC);
  PASS("polyveck_hint_pack_basic round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(42);

  printf("=== test_mod2q (MODE=%d, N=%d, Q=%d, L=%d, M=%d, ALPHA_H=%d, TAU=%d) ===\n",
         SHUTTLE_MODE, SHUTTLE_N, SHUTTLE_Q, SHUTTLE_L, SHUTTLE_M,
         SHUTTLE_ALPHA_H, SHUTTLE_TAU);
  printf("  DQ=%d, HALF_ALPHA_H=%d, HINT_MOD=%d, HINT_MAX=%d, ALPHA_1=%d\n",
         SHUTTLE_DQ, SHUTTLE_HALF_ALPHA_H, SHUTTLE_HINT_MOD, SHUTTLE_HINT_MAX,
         SHUTTLE_ALPHA_1);

  ret |= test_scalar_lift();
  ret |= test_scalar_highbits();
  ret |= test_poly_helpers();
  ret |= test_signer_verifier_identity();
  ret |= test_lsb_invariants();
  ret |= test_compress_y_slot0();
  ret |= test_commitment_ntt_matches_schoolbook();
  ret |= test_hint_wh_roundtrip();
  ret |= test_z2_recovery();
  ret |= test_hint_pack_basic();

  if (ret == 0) {
    printf("\n=== All mod 2q tests PASSED ===\n");
  } else {
    printf("\n=== Some mod 2q tests FAILED ===\n");
  }
  return ret;
}
