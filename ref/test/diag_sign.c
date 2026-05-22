#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "../sign.h"
#include "../packing.h"
#include "../shuttle_rans.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../params.h"
#include "../randombytes.h"
#include "../rejsample.h"
#include "../poly.h"

extern int pack_sig(uint8_t *, const uint8_t *, const int8_t *,
                    const poly *, const polyveck *);

int main(void) {
  uint8_t pk[SHUTTLE_PUBLICKEYBYTES];
  uint8_t sk[SHUTTLE_SECRETKEYBYTES];
  crypto_sign_keypair(pk, sk);
  printf("keypair OK\n");

  uint8_t sig[SHUTTLE_BYTES];
  size_t siglen;
  uint8_t msg[16] = "test";
  int ret = crypto_sign_signature(sig, &siglen, msg, 16, sk);
  printf("sign ret=%d siglen=%zu\n", ret, siglen);
  return 0;
}
