/*
 * kat_dump.c - Emit a deterministic signature dump for byte-level
 * comparison between SHUTTLE ref/ and avx2/ implementations.
 *
 * Protocol (one record per iteration):
 *   count=<i> seed=<hex48> mlen=<mlen> msg=<hex> pk=<hex> sk=<hex> sm=<hex>
 *
 * Each iteration:
 *   1. Derive a per-iteration 48-byte seed from a top-level 48-byte seed
 *      via SHAKE-256 (so ref and avx2 agree on the seed schedule without
 *      having to consume randombytes() in a specific pre-sign pattern).
 *   2. Re-seed the DRBG from that per-iteration seed.
 *   3. Sample a random message of length mlen bytes (mlen cycles through
 *      a small set).
 *   4. Run keypair -> sign_combined (`crypto_sign`) -> open. Abort on
 *      open failure (implementations that disagree on correctness
 *      would also fail the make check tests upstream).
 *   5. Print the record.
 *
 * Invoked as:
 *   kat_dump <nrounds>
 *
 * nrounds defaults to 20.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "rng.h"
#include "fips202.h"
#include "sign.h"
#include "params.h"

#define MLEN_MAX 128

static void hex_print(const uint8_t *buf, size_t len)
{
    for (size_t i = 0; i < len; ++i) printf("%02x", buf[i]);
}

/* Derive per-iteration seed i from master seed via independent SHAKE-256
 * stream. Using a separate SHAKE instance leaves the main DRBG stream
 * untouched. */
static void derive_iter_seed(uint8_t out[48],
                             const uint8_t master[48],
                             uint32_t i)
{
    keccak_state st;
    uint8_t nonce[4];
    nonce[0] = (uint8_t)(i >>  0);
    nonce[1] = (uint8_t)(i >>  8);
    nonce[2] = (uint8_t)(i >> 16);
    nonce[3] = (uint8_t)(i >> 24);
    shake256_init(&st);
    shake256_absorb(&st, master, 48);
    shake256_absorb(&st, nonce, 4);
    shake256_finalize(&st);
    shake256_squeeze(out, 48, &st);
}

int main(int argc, char **argv)
{
    unsigned nrounds = 20;
    if (argc >= 2) nrounds = (unsigned)atoi(argv[1]);

    /* Fixed master seed so runs are reproducible. */
    uint8_t master[48];
    for (unsigned i = 0; i < 48; ++i) master[i] = (uint8_t)(0xA0 + i);

    printf("algorithm=%s mode=%d n=%d l=%d m=%d pk=%d sk=%d sig=%d\n",
           CRYPTO_ALGNAME, SHUTTLE_MODE, SHUTTLE_N, SHUTTLE_L, SHUTTLE_M,
           SHUTTLE_PUBLICKEYBYTES, SHUTTLE_SECRETKEYBYTES, SHUTTLE_BYTES);
    printf("master=");
    hex_print(master, 48);
    printf("\n");

    uint8_t pk[SHUTTLE_PUBLICKEYBYTES];
    uint8_t sk[SHUTTLE_SECRETKEYBYTES];
    uint8_t sm[SHUTTLE_BYTES + MLEN_MAX];
    uint8_t m2[MLEN_MAX];
    uint8_t msg[MLEN_MAX];

    for (unsigned i = 0; i < nrounds; ++i) {
        uint8_t iter_seed[48];
        derive_iter_seed(iter_seed, master, i);
        randombytes_seed(iter_seed);

        /* Cycle mlen through a few values so we exercise different msg
         * lengths without spending an extra randombytes call on length
         * selection (which would shift the randomness schedule). */
        size_t mlen = ((size_t)i * 7 + 3) % MLEN_MAX + 1;
        randombytes(msg, mlen);

        if (crypto_sign_keypair(pk, sk) != 0) {
            fprintf(stderr, "keypair failed at round %u\n", i);
            return 1;
        }

        size_t smlen;
        if (crypto_sign(sm, &smlen, msg, mlen, sk) != 0) {
            fprintf(stderr, "sign failed at round %u\n", i);
            return 1;
        }

        size_t m2len;
        if (crypto_sign_open(m2, &m2len, sm, smlen, pk) != 0) {
            fprintf(stderr, "open failed at round %u\n", i);
            return 1;
        }
        if (m2len != mlen || memcmp(m2, msg, mlen) != 0) {
            fprintf(stderr, "recovered message mismatch at round %u\n", i);
            return 1;
        }

        printf("count=%u seed=", i);
        hex_print(iter_seed, 48);
        printf(" mlen=%zu msg=", mlen);
        hex_print(msg, mlen);
        printf(" pk=");
        hex_print(pk, SHUTTLE_PUBLICKEYBYTES);
        printf(" sk=");
        hex_print(sk, SHUTTLE_SECRETKEYBYTES);
        printf(" smlen=%zu sm=", smlen);
        hex_print(sm, smlen);
        printf("\n");
    }

    return 0;
}
