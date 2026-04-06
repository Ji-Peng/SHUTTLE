/*
 * test_sign.c - Basic functional tests for the NGCC_SIGN signature scheme.
 *
 * Tests:
 *   1. Key generation: generates a keypair, prints sizes
 *   2. Sign + verify: signs a test message and verifies it
 *   3. Forgery detection: flips a bit in the signature and checks rejection
 *   4. Stress test: 100 rounds of sign-verify with random messages
 *   5. crypto_sign / crypto_sign_open (combined API)
 *
 * Compile:
 *   gcc -O3 -std=c99 -o test_sign test/test_sign.c sign.c packing.c \
 *       rejsample.c approx_log.c poly.c polyvec.c ntt.c reduce.c rounding.c \
 *       sampler.c sampler_4x.c approx_exp.c fips202.c symmetric-shake.c \
 *       randombytes.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../sign.h"
#include "../params.h"
#include "../randombytes.h"

#define MLEN_MAX 256
#define NROUNDS  100

static void print_hex(const char *label, const uint8_t *data, size_t len)
{
  printf("%s (%zu bytes): ", label, len);
  for(size_t i = 0; i < len && i < 32; ++i)
    printf("%02x", data[i]);
  if(len > 32)
    printf("...");
  printf("\n");
}

int main(void)
{
  uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES];
  uint8_t sk[NGCC_SIGN_SECRETKEYBYTES];
  uint8_t sig[NGCC_SIGN_BYTES];
  uint8_t msg[MLEN_MAX];
  size_t siglen, mlen;
  int ret;
  unsigned int i;

  printf("=== NGCC_SIGN Functional Test ===\n\n");
  printf("Parameters:\n");
  printf("  N = %d, Q = %d, L = %d, M = %d\n",
         NGCC_SIGN_N, NGCC_SIGN_Q, NGCC_SIGN_L, NGCC_SIGN_M);
  printf("  ETA = %d, TAU = %d, SIGMA = %d\n",
         NGCC_SIGN_ETA, NGCC_SIGN_TAU, NGCC_SIGN_SIGMA);
  printf("  Public key size:  %d bytes\n", NGCC_SIGN_PUBLICKEYBYTES);
  printf("  Secret key size:  %d bytes\n", NGCC_SIGN_SECRETKEYBYTES);
  printf("  Signature size:   %d bytes\n", NGCC_SIGN_BYTES);
  printf("  BS_SQ = %d, BV_SQ = %d\n\n",
         NGCC_SIGN_BS_SQ, NGCC_SIGN_BV_SQ);

  /* ---- Test 1: Key Generation ---- */
  printf("[Test 1] Key generation...\n");
  ret = crypto_sign_keypair(pk, sk);
  if(ret != 0) {
    printf("  FAIL: keypair returned %d\n", ret);
    return 1;
  }
  print_hex("  pk", pk, NGCC_SIGN_PUBLICKEYBYTES);
  print_hex("  sk", sk, NGCC_SIGN_SECRETKEYBYTES);
  printf("  PASS\n\n");

  /* ---- Test 2: Sign + Verify ---- */
  printf("[Test 2] Sign and verify a test message...\n");
  mlen = 33;
  memcpy(msg, "NGCC_SIGN test message, round 0!", mlen);
  msg[mlen - 1] = '\0';

  ret = crypto_sign_signature(sig, &siglen, msg, mlen, sk);
  if(ret != 0) {
    printf("  FAIL: signature returned %d\n", ret);
    return 1;
  }
  printf("  Signature length: %zu bytes\n", siglen);
  print_hex("  sig", sig, siglen);

  ret = crypto_sign_verify(sig, siglen, msg, mlen, pk);
  if(ret != 0) {
    printf("  FAIL: valid signature rejected (ret=%d)\n", ret);
    return 1;
  }
  printf("  Verification: PASS\n\n");

  /* ---- Test 3: Forgery Detection ---- */
  printf("[Test 3] Forgery detection (bit flip in signature)...\n");
  {
    uint8_t sig_bad[NGCC_SIGN_BYTES];
    memcpy(sig_bad, sig, NGCC_SIGN_BYTES);

    /* Flip a bit in the z1 portion of the signature */
    sig_bad[NGCC_SIGN_CTILDEBYTES + 10] ^= 0x01;

    ret = crypto_sign_verify(sig_bad, NGCC_SIGN_BYTES, msg, mlen, pk);
    if(ret == 0) {
      printf("  FAIL: forged signature accepted!\n");
      return 1;
    }
    printf("  Forged sig correctly rejected (ret=%d)\n", ret);
  }

  /* Flip a bit in the message */
  {
    uint8_t msg_bad[MLEN_MAX];
    memcpy(msg_bad, msg, mlen);
    msg_bad[0] ^= 0x80;

    ret = crypto_sign_verify(sig, siglen, msg_bad, mlen, pk);
    if(ret == 0) {
      printf("  FAIL: wrong-message signature accepted!\n");
      return 1;
    }
    printf("  Wrong-message sig correctly rejected (ret=%d)\n", ret);
  }

  /* Flip a bit near the end of the z portion */
  {
    uint8_t sig_bad[NGCC_SIGN_BYTES];
    memcpy(sig_bad, sig, NGCC_SIGN_BYTES);

    /* Flip a bit in the last polynomial of z (near end of signature) */
    sig_bad[NGCC_SIGN_BYTES - 10] ^= 0x04;

    ret = crypto_sign_verify(sig_bad, NGCC_SIGN_BYTES, msg, mlen, pk);
    if(ret == 0) {
      printf("  FAIL: tail-flipped signature accepted!\n");
      return 1;
    }
    printf("  Tail-flipped sig correctly rejected (ret=%d)\n", ret);
  }
  printf("  PASS\n\n");

  /* ---- Test 4: Stress Test ---- */
  printf("[Test 4] Stress test: %d rounds of sign+verify...\n", NROUNDS);
  for(i = 0; i < NROUNDS; ++i) {
    /* Random message of random length [1, MLEN_MAX] */
    uint8_t rmsg[MLEN_MAX];
    uint8_t rlen_buf;
    size_t rmlen;

    randombytes(&rlen_buf, 1);
    rmlen = (rlen_buf % MLEN_MAX) + 1;
    randombytes(rmsg, rmlen);

    ret = crypto_sign_signature(sig, &siglen, rmsg, rmlen, sk);
    if(ret != 0) {
      printf("  FAIL at round %u: sign returned %d\n", i, ret);
      return 1;
    }

    ret = crypto_sign_verify(sig, siglen, rmsg, rmlen, pk);
    if(ret != 0) {
      printf("  FAIL at round %u: verify returned %d\n", i, ret);
      return 1;
    }

    if((i + 1) % 10 == 0)
      printf("  Completed %u/%d rounds\n", i + 1, NROUNDS);
  }
  printf("  All %d rounds PASS\n\n", NROUNDS);

  /* ---- Test 5: Combined sign/open API ---- */
  printf("[Test 5] Combined crypto_sign / crypto_sign_open...\n");
  {
    uint8_t sm[NGCC_SIGN_BYTES + MLEN_MAX];
    uint8_t m2[MLEN_MAX];
    size_t smlen, m2len;

    mlen = 16;
    memcpy(msg, "Combined API tst", mlen);

    ret = crypto_sign(sm, &smlen, msg, mlen, sk);
    if(ret != 0) {
      printf("  FAIL: crypto_sign returned %d\n", ret);
      return 1;
    }
    printf("  Signed message length: %zu bytes\n", smlen);

    ret = crypto_sign_open(m2, &m2len, sm, smlen, pk);
    if(ret != 0) {
      printf("  FAIL: crypto_sign_open returned %d\n", ret);
      return 1;
    }

    if(m2len != mlen || memcmp(m2, msg, mlen) != 0) {
      printf("  FAIL: recovered message differs\n");
      return 1;
    }
    printf("  Recovered message matches original\n");
    printf("  PASS\n\n");
  }

  printf("=== All tests PASSED ===\n");
  return 0;
}
