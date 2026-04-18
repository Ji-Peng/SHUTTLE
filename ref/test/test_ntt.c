/*
 * test_ntt.c - NTT correctness tests for SHUTTLE.
 *
 * Tests:
 *  1) NTT -> INTT roundtrip
 *  2) NTT convolution matches schoolbook multiplication mod x^n+1
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../params.h"
#include "../ntt.h"
#include "../reduce.h"
#include "../poly.h"
#include "../randombytes.h"

#define N SHUTTLE_N
#define Q SHUTTLE_Q

/* Schoolbook polynomial multiplication mod x^n+1, mod q */
static void poly_schoolbook(int32_t c[N], const int32_t a[N], const int32_t b[N]) {
  int64_t t;
  for (int i = 0; i < N; i++) {
    t = 0;
    for (int j = 0; j <= i; j++)
      t += (int64_t)a[j] * b[i - j];
    for (int j = i + 1; j < N; j++)
      t -= (int64_t)a[j] * b[N + i - j];  /* x^n = -1 */
    /* Reduce mod q */
    t %= Q;
    if (t < 0) t += Q;
    c[i] = (int32_t)t;
  }
}

/* Test 1: NTT -> INTT roundtrip */
static int test_ntt_roundtrip(void) {
  int32_t a[N], a_orig[N];

  printf("Test NTT roundtrip... ");

  /* Initialize with small random values */
  for (int i = 0; i < N; i++) {
    a[i] = (rand() % (2 * Q)) - Q;  /* values in [-Q, Q] */
    a_orig[i] = a[i];
  }

  /* Forward NTT */
  ntt(a);

  /* Inverse NTT (returns in Montgomery form: a * 2^32 mod q) */
  invntt_tomont(a);

  /* The result should be a_orig * MONT mod q */
  for (int i = 0; i < N; i++) {
    /* Reduce a[i] to [0, q) */
    int32_t got = freeze(a[i]);
    /* Expected: a_orig[i] * MONT mod q */
    int32_t expect = (int32_t)(((int64_t)a_orig[i] * SHUTTLE_MONT % Q + Q) % Q);
    if (got != expect) {
      printf("FAIL at index %d: got %d, expected %d (orig=%d)\n",
             i, got, expect, a_orig[i]);
      return 1;
    }
  }
  printf("PASS\n");
  return 0;
}

/* Test 2: NTT convolution vs schoolbook */
static int test_ntt_convolution(void) {
  int32_t a[N], b[N], c_ntt[N], c_school[N];

  printf("Test NTT convolution... ");

  /* Small random polynomials in [0, q) */
  for (int i = 0; i < N; i++) {
    a[i] = rand() % Q;
    b[i] = rand() % Q;
  }

  /* Schoolbook multiplication */
  poly_schoolbook(c_school, a, b);

  /* NTT multiplication */
  int32_t a_ntt[N], b_ntt[N];
  memcpy(a_ntt, a, sizeof(a));
  memcpy(b_ntt, b, sizeof(b));
  ntt(a_ntt);
  ntt(b_ntt);

  /* Pointwise multiply (in Montgomery domain) */
  for (int i = 0; i < N; i++)
    c_ntt[i] = montgomery_reduce((int64_t)a_ntt[i] * b_ntt[i]);

  /* Inverse NTT */
  invntt_tomont(c_ntt);

  /* invntt_tomont: result = a*b mod q (the Montgomery factors cancel) */
  for (int i = 0; i < N; i++) {
    c_ntt[i] = freeze(c_ntt[i]);
  }

  /* Compare: c_ntt[i] should be c_school[i] directly */
  for (int i = 0; i < N; i++) {
    if (c_ntt[i] != c_school[i]) {
      printf("FAIL at index %d: got %d, expected %d\n", i, c_ntt[i], c_school[i]);
      return 1;
    }
  }
  printf("PASS\n");
  return 0;
}

/* Test 3: Using poly abstraction for multiply */
static int test_poly_multiply(void) {
  poly a, b, c;

  printf("Test poly NTT multiply... ");

  /* Small random polynomials */
  for (int i = 0; i < N; i++) {
    a.coeffs[i] = rand() % Q;
    b.coeffs[i] = rand() % Q;
  }

  /* Schoolbook for reference */
  int32_t c_school[N];
  poly_schoolbook(c_school, a.coeffs, b.coeffs);

  /* NTT-based multiply */
  poly a_hat, b_hat;
  memcpy(&a_hat, &a, sizeof(poly));
  memcpy(&b_hat, &b, sizeof(poly));
  poly_ntt(&a_hat);
  poly_ntt(&b_hat);
  poly_pointwise_montgomery(&c, &a_hat, &b_hat);
  poly_invntt_tomont(&c);
  poly_caddq(&c);

  /* Compare mod q: result should be a*b directly */
  for (int i = 0; i < N; i++) {
    int32_t got = freeze(c.coeffs[i]);
    if (got != c_school[i]) {
      printf("FAIL at index %d: got %d, expected %d\n", i, got, c_school[i]);
      return 1;
    }
  }
  printf("PASS\n");
  return 0;
}

int main(void) {
  int ret = 0;

  srand(42);  /* Deterministic for reproducibility */

  ret |= test_ntt_roundtrip();
  ret |= test_ntt_convolution();
  ret |= test_poly_multiply();

  if (ret == 0)
    printf("\nAll NTT tests passed!\n");
  else
    printf("\nSome NTT tests FAILED!\n");

  return ret;
}
