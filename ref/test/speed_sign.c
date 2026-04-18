#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "../params.h"
#include "../sign.h"
#include "../randombytes.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 1000

static uint64_t t[NTESTS];

int main(void) {
  uint8_t pk[SHUTTLE_PUBLICKEYBYTES];
  uint8_t sk[SHUTTLE_SECRETKEYBYTES];
  uint8_t sig[SHUTTLE_BYTES];
  size_t siglen;
  uint8_t msg[32] = {0};

  printf("=== SHUTTLE Benchmarks - Ref ===\n");
  printf("Parameters: N=%d, Q=%d, L=%d, M=%d, SIGMA=%d, TAU=%d\n",
         SHUTTLE_N, SHUTTLE_Q, SHUTTLE_L, SHUTTLE_M,
         SHUTTLE_SIGMA, SHUTTLE_TAU);
  printf("PK=%d bytes, SK=%d bytes, SIG=%d bytes\n\n",
         SHUTTLE_PUBLICKEYBYTES, SHUTTLE_SECRETKEYBYTES, SHUTTLE_BYTES);

  /* KeyGen */
  for (int i = 0; i < NTESTS; i++) {
    t[i] = cpucycles();
    crypto_sign_keypair(pk, sk);
  }
  print_results("crypto_sign_keypair:", t, NTESTS);

  /* Sign */
  crypto_sign_keypair(pk, sk);
  for (int i = 0; i < NTESTS; i++) {
    msg[0] = i & 0xFF;
    t[i] = cpucycles();
    crypto_sign_signature(sig, &siglen, msg, 32, sk);
  }
  print_results("crypto_sign_signature:", t, NTESTS);

  /* Verify */
  crypto_sign_signature(sig, &siglen, msg, 32, sk);
  for (int i = 0; i < NTESTS; i++) {
    t[i] = cpucycles();
    crypto_sign_verify(sig, siglen, msg, 32, pk);
  }
  print_results("crypto_sign_verify:", t, NTESTS);

  return 0;
}
