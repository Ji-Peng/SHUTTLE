/*
 * speed_sign.c - Performance benchmarks for SHUTTLE top-level API.
 *
 * Implementation-agnostic: prints the CRYPTO_ALGNAME (including the
 * _ref / _avx2 suffix via SHUTTLE_NAMESPACETOP) so the same file can be
 * linked from both ref/ and avx2/ and still emit a self-identifying
 * header.
 *
 * Benchmarks:
 *   1) crypto_sign_keypair
 *   2) crypto_sign_signature
 *   3) crypto_sign_verify
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "../sign.h"
#include "../params.h"
#include "../randombytes.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 1000
#define MLEN   59

#define SHUTTLE_STR_(x) #x
#define SHUTTLE_STR(x)  SHUTTLE_STR_(x)

int main(void) {
    uint8_t pk[SHUTTLE_PUBLICKEYBYTES];
    uint8_t sk[SHUTTLE_SECRETKEYBYTES];
    uint8_t sig[SHUTTLE_BYTES];
    uint8_t msg[MLEN];
    size_t siglen;
    uint64_t t[NTESTS];
    unsigned int i;

    randombytes(msg, MLEN);

    printf("=== %s Benchmarks (%s) ===\n",
           CRYPTO_ALGNAME, SHUTTLE_STR(SHUTTLE_NAMESPACETOP));
    printf("Parameters: N=%d, Q=%d, L=%d, M=%d, SIGMA=%d, TAU=%d\n",
           SHUTTLE_N, SHUTTLE_Q, SHUTTLE_L, SHUTTLE_M,
           SHUTTLE_SIGMA, SHUTTLE_TAU);
    printf("PK=%d bytes, SK=%d bytes, SIG=%d bytes\n\n",
           SHUTTLE_PUBLICKEYBYTES, SHUTTLE_SECRETKEYBYTES,
           SHUTTLE_BYTES);

    /* ---- 1. KeyGen ---- */
    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_sign_keypair(pk, sk);
    }
    print_results("crypto_sign_keypair:", t, NTESTS);

    /* ---- 2. Sign ---- */
    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_sign_signature(sig, &siglen, msg, MLEN, sk);
    }
    print_results("crypto_sign_signature:", t, NTESTS);

    /* ---- 3. Verify ---- */
    crypto_sign_signature(sig, &siglen, msg, MLEN, sk);
    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_sign_verify(sig, siglen, msg, MLEN, pk);
    }
    print_results("crypto_sign_verify:", t, NTESTS);

    return 0;
}
