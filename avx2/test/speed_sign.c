/*
 * speed_sign.c - Performance benchmarks for NGCC_SIGN (AVX2).
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

int main(void) {
    uint8_t pk[NGCC_SIGN_PUBLICKEYBYTES];
    uint8_t sk[NGCC_SIGN_SECRETKEYBYTES];
    uint8_t sig[NGCC_SIGN_BYTES];
    uint8_t msg[MLEN];
    size_t siglen;
    uint64_t t[NTESTS];
    unsigned int i;

    randombytes(msg, MLEN);

    printf("=== NGCC_SIGN Benchmarks - AVX2 ===\n");
    printf("Parameters: N=%d, Q=%d, L=%d, M=%d, SIGMA=%d, TAU=%d\n",
           NGCC_SIGN_N, NGCC_SIGN_Q, NGCC_SIGN_L, NGCC_SIGN_M,
           NGCC_SIGN_SIGMA, NGCC_SIGN_TAU);
    printf("PK=%d bytes, SK=%d bytes, SIG=%d bytes\n\n",
           NGCC_SIGN_PUBLICKEYBYTES, NGCC_SIGN_SECRETKEYBYTES,
           NGCC_SIGN_BYTES);

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
