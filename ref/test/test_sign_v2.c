/*
 * test_sign_v2.c - end-to-end tests for Phase 6b-5c (Alg 2 compressed sig).
 *
 * Tests:
 *   1. Key generation.
 *   2. Sign + verify roundtrip on a test message.
 *   3. Forgery detection (bit-flip in signature + message tampering).
 *   4. 100-round stress test with random messages.
 *   5. Combined crypto_sign_v2 / crypto_sign_open_v2.
 *
 * Key difference from test_sign.c:
 *   - Uses _v2 public entry points (Alg 2 flow with mod 2q, CompressY,
 *     MakeHint/UseHint, rANS hint packing).
 *   - Signature size is SHUTTLE_BYTES_V2 (smaller than Phase 6a).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../sign_v2.h"
#include "../packing_v2.h"
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
  uint8_t pk[SHUTTLE_PUBLICKEYBYTES];
  uint8_t sk[SHUTTLE_SECRETKEYBYTES];
  uint8_t sig[SHUTTLE_BYTES_V2];
  uint8_t msg[MLEN_MAX];
  size_t siglen, mlen;
  int ret;
  unsigned int i;

  printf("=== SHUTTLE v2 Functional Test ===\n\n");
  printf("Parameters:\n");
  printf("  MODE = %d, N = %d, Q = %d, L = %d, M = %d\n",
         SHUTTLE_MODE, SHUTTLE_N, SHUTTLE_Q, SHUTTLE_L, SHUTTLE_M);
  printf("  ETA = %d, TAU = %d, SIGMA = %d, ALPHA_1 = %d, ALPHA_H = %d\n",
         SHUTTLE_ETA, SHUTTLE_TAU, SHUTTLE_SIGMA, SHUTTLE_ALPHA_1, SHUTTLE_ALPHA_H);
  printf("  Public key:  %d bytes\n", SHUTTLE_PUBLICKEYBYTES);
  printf("  Secret key:  %d bytes\n", SHUTTLE_SECRETKEYBYTES);
  printf("  Signature:   %d bytes (v2, incl. rANS-compressed hint)\n", SHUTTLE_BYTES_V2);
  printf("  BS_SQ = %lu, BV_SQ = %lu\n\n",
         (unsigned long)SHUTTLE_BS_SQ, (unsigned long)SHUTTLE_BV_SQ);

  /* ---- Test 1: Key Generation ---- */
  printf("[Test 1] Key generation...\n");
  ret = crypto_sign_keypair_v2(pk, sk);
  if(ret != 0) { printf("  FAIL: keypair returned %d\n", ret); return 1; }
  print_hex("  pk", pk, SHUTTLE_PUBLICKEYBYTES);
  print_hex("  sk", sk, SHUTTLE_SECRETKEYBYTES);
  printf("  PASS\n\n");

  /* ---- Test 2: Sign + Verify ---- */
  printf("[Test 2] Sign and verify...\n");
  mlen = 33;
  memcpy(msg, "SHUTTLE v2 test message, round 0", mlen);

  ret = crypto_sign_signature_v2(sig, &siglen, msg, mlen, sk);
  if(ret != 0) { printf("  FAIL: sign returned %d\n", ret); return 1; }
  printf("  Signature length: %zu bytes\n", siglen);
  print_hex("  sig", sig, siglen);

  ret = crypto_sign_verify_v2(sig, siglen, msg, mlen, pk);
  if(ret != 0) { printf("  FAIL: valid signature rejected (ret=%d)\n", ret); return 1; }
  printf("  Verification: PASS\n\n");

  /* ---- Test 3: Forgery Detection ---- */
  printf("[Test 3] Forgery detection...\n");
  {
    uint8_t sig_bad[SHUTTLE_BYTES_V2];
    memcpy(sig_bad, sig, SHUTTLE_BYTES_V2);
    sig_bad[SHUTTLE_CTILDEBYTES + 10] ^= 0x01;
    ret = crypto_sign_verify_v2(sig_bad, SHUTTLE_BYTES_V2, msg, mlen, pk);
    if(ret == 0) { printf("  FAIL: forged signature accepted!\n"); return 1; }
    printf("  Forged sig rejected (ret=%d)\n", ret);
  }
  {
    uint8_t msg_bad[MLEN_MAX];
    memcpy(msg_bad, msg, mlen);
    msg_bad[0] ^= 0x80;
    ret = crypto_sign_verify_v2(sig, siglen, msg_bad, mlen, pk);
    if(ret == 0) { printf("  FAIL: tampered message accepted!\n"); return 1; }
    printf("  Tampered message rejected (ret=%d)\n", ret);
  }
  printf("  PASS\n\n");

  /* ---- Test 4: Stress ---- */
  printf("[Test 4] Stress test (%d rounds)...\n", NROUNDS);
  for(i = 0; i < NROUNDS; ++i) {
    uint8_t rmsg[MLEN_MAX];
    uint8_t rlen_buf;
    size_t rmlen;
    randombytes(&rlen_buf, 1);
    rmlen = (rlen_buf % MLEN_MAX) + 1;
    randombytes(rmsg, rmlen);

    ret = crypto_sign_signature_v2(sig, &siglen, rmsg, rmlen, sk);
    if(ret != 0) { printf("  FAIL at round %u: sign returned %d\n", i, ret); return 1; }
    ret = crypto_sign_verify_v2(sig, siglen, rmsg, rmlen, pk);
    if(ret != 0) { printf("  FAIL at round %u: verify returned %d\n", i, ret); return 1; }
    if((i + 1) % 10 == 0)
      printf("  Completed %u/%d rounds\n", i + 1, NROUNDS);
  }
  printf("  All %d rounds PASS\n\n", NROUNDS);

  /* ---- Test 5: Combined API ---- */
  printf("[Test 5] Combined crypto_sign_v2 / crypto_sign_open_v2...\n");
  {
    uint8_t sm[SHUTTLE_BYTES_V2 + MLEN_MAX];
    uint8_t m2[MLEN_MAX];
    size_t smlen, m2len;

    mlen = 16;
    memcpy(msg, "Combined API tst", mlen);

    ret = crypto_sign_v2(sm, &smlen, msg, mlen, sk);
    if(ret != 0) { printf("  FAIL: crypto_sign_v2 returned %d\n", ret); return 1; }
    printf("  Signed message length: %zu bytes\n", smlen);

    ret = crypto_sign_open_v2(m2, &m2len, sm, smlen, pk);
    if(ret != 0) { printf("  FAIL: crypto_sign_open_v2 returned %d\n", ret); return 1; }
    if(m2len != mlen || memcmp(m2, msg, mlen) != 0) {
      printf("  FAIL: recovered message differs\n"); return 1;
    }
    printf("  Recovered message matches\n  PASS\n\n");
  }

  printf("=== All v2 tests PASSED ===\n");
  return 0;
}
