/*
 * gen_vectors.c - Emit deterministic test vectors for the SHUTTLE
 * signature scheme, split into three phase files:
 *
 *   keygen_mode<MODE>.kat   : (seed) -> (pk, sk)
 *   sign_mode<MODE>.kat     : (sk, msg, sign_seed) -> (sig)
 *   verify_mode<MODE>.kat   : (pk, msg, sig, expected) records, including
 *                             positive and negative cases (sig/msg tampered).
 *
 * The DRBG is the SHAKE-256-based randombytes() in rng.c, so every record
 * is reproducible from the seeds shown in the file.
 *
 * Master seed is fixed (0xA0..0xCF). Per-iteration seeds are derived from
 * (master || domain || i_LE32) via an independent SHAKE-256 stream so the
 * keygen / sign randomness schedules are independent.
 *
 * Usage:
 *   gen_vectors <out_dir> [<nrounds>]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "rng.h"
#include "fips202.h"
#include "sign.h"
#include "params.h"

#define MAX_NROUNDS  64
#define MAX_MLEN     256

static const size_t MLEN_CYCLE[] = { 1, 16, 33, 64, 100, 128, 200, 256 };
#define MLEN_CYCLE_LEN (sizeof(MLEN_CYCLE) / sizeof(MLEN_CYCLE[0]))

/* Fixed master seed so vectors are reproducible. */
static const uint8_t MASTER_SEED[48] = {
    0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7,
    0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF,
    0xB0, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7,
    0xB8, 0xB9, 0xBA, 0xBB, 0xBC, 0xBD, 0xBE, 0xBF,
    0xC0, 0xC1, 0xC2, 0xC3, 0xC4, 0xC5, 0xC6, 0xC7,
    0xC8, 0xC9, 0xCA, 0xCB, 0xCC, 0xCD, 0xCE, 0xCF
};

static void hex_fprint(FILE *f, const uint8_t *buf, size_t len)
{
    for (size_t i = 0; i < len; ++i) fprintf(f, "%02x", buf[i]);
}

/* Derive a 48-byte seed via independent SHAKE-256 stream:
 *   out = SHAKE256(master || domain || I2B(i, 4)). */
static void derive_seed(uint8_t out[48],
                        const uint8_t master[48],
                        const char *domain,
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
    shake256_absorb(&st, (const uint8_t *)domain, strlen(domain));
    shake256_absorb(&st, nonce, 4);
    shake256_finalize(&st);
    shake256_squeeze(out, 48, &st);
}

/* Derive a deterministic message of mlen bytes for round i.
 * Independent of the keygen/sign DRBG so the message is fixed regardless
 * of how many DRBG bytes those phases consume. */
static void derive_msg(uint8_t *msg, size_t mlen,
                       const uint8_t master[48], uint32_t i)
{
    keccak_state st;
    uint8_t nonce[4];
    nonce[0] = (uint8_t)(i >>  0);
    nonce[1] = (uint8_t)(i >>  8);
    nonce[2] = (uint8_t)(i >> 16);
    nonce[3] = (uint8_t)(i >> 24);
    shake256_init(&st);
    shake256_absorb(&st, master, 48);
    shake256_absorb(&st, (const uint8_t *)"msg", 3);
    shake256_absorb(&st, nonce, 4);
    shake256_finalize(&st);
    shake256_squeeze(msg, mlen, &st);
}

static FILE *open_kat(const char *out_dir, const char *phase)
{
    char path[512];
    int n = snprintf(path, sizeof path,
                     "%s/%s_mode%d.kat", out_dir, phase, SHUTTLE_MODE);
    if (n <= 0 || (size_t)n >= sizeof path) {
        fprintf(stderr, "output path too long\n");
        exit(1);
    }
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    fprintf(stderr, "  writing %s\n", path);
    return f;
}

static void write_keygen_header(FILE *f)
{
    fprintf(f, "# %s KeyGen test vectors\n", CRYPTO_ALGNAME);
    fprintf(f, "# algorithm = %s\n", CRYPTO_ALGNAME);
    fprintf(f, "# mode = %d\n", SHUTTLE_MODE);
    fprintf(f, "# n = %d, l = %d, m = %d\n",
            SHUTTLE_N, SHUTTLE_L, SHUTTLE_M);
    fprintf(f, "# pk_bytes = %d\n", SHUTTLE_PUBLICKEYBYTES);
    fprintf(f, "# sk_bytes = %d\n", SHUTTLE_SECRETKEYBYTES);
    fprintf(f, "# Each record: feed `seed` (48 B) into the SHAKE-256 DRBG\n");
    fprintf(f, "# (rng.c) and call crypto_sign_keypair to reproduce pk,sk.\n\n");
}

static void write_sign_header(FILE *f)
{
    fprintf(f, "# %s Sign test vectors\n", CRYPTO_ALGNAME);
    fprintf(f, "# algorithm = %s\n", CRYPTO_ALGNAME);
    fprintf(f, "# mode = %d\n", SHUTTLE_MODE);
    fprintf(f, "# sig_bytes = %d (fixed; rANS streams are zero-padded\n"
               "#              to reserved slots, so siglen == sig_bytes)\n",
            SHUTTLE_BYTES);
    fprintf(f, "# Each record: re-seed the DRBG with `sign_seed`, then call\n");
    fprintf(f, "# crypto_sign_signature(sig,&siglen,msg,mlen,sk) with the sk\n");
    fprintf(f, "# from the matching keygen[count] record. The signer uses one\n");
    fprintf(f, "# 32-byte randombytes() draw for `rnd`; the rest of the\n");
    fprintf(f, "# randomness is derived from sk and msg.\n\n");
}

static void write_verify_header(FILE *f)
{
    fprintf(f, "# %s Verify test vectors\n", CRYPTO_ALGNAME);
    fprintf(f, "# algorithm = %s\n", CRYPTO_ALGNAME);
    fprintf(f, "# mode = %d\n", SHUTTLE_MODE);
    fprintf(f, "# Each record: call crypto_sign_verify(sig,siglen,msg,mlen,pk).\n");
    fprintf(f, "# `expected` is `accept` (rc == 0) or `reject` (rc != 0).\n");
    fprintf(f, "# Subcases per round:\n");
    fprintf(f, "#   original     : untampered (sig,msg) -> accept\n");
    fprintf(f, "#   sig_tampered : one bit flipped in the rANS Z_0 region -> reject\n");
    fprintf(f, "#   msg_tampered : MSB flipped on msg[0]                 -> reject\n\n");
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr,
                "usage: %s <out_dir> [<nrounds>]\n", argv[0]);
        return 1;
    }
    const char *out_dir = argv[1];
    unsigned nrounds = 10;
    if (argc >= 3) {
        long v = strtol(argv[2], NULL, 10);
        if (v <= 0 || v > MAX_NROUNDS) {
            fprintf(stderr, "nrounds out of range (1..%d)\n", MAX_NROUNDS);
            return 1;
        }
        nrounds = (unsigned)v;
    }

    fprintf(stderr, "[gen_vectors] %s, mode=%d, nrounds=%u\n",
            CRYPTO_ALGNAME, SHUTTLE_MODE, nrounds);
    fprintf(stderr, "[gen_vectors] master seed = ");
    for (unsigned i = 0; i < 48; ++i) fprintf(stderr, "%02x", MASTER_SEED[i]);
    fprintf(stderr, "\n");

    FILE *f_kg  = open_kat(out_dir, "keygen");
    FILE *f_sg  = open_kat(out_dir, "sign");
    FILE *f_vf  = open_kat(out_dir, "verify");

    write_keygen_header(f_kg);
    write_sign_header(f_sg);
    write_verify_header(f_vf);

    /* Heap-allocate to avoid blowing up the stack at MODE=256. */
    uint8_t *pk  = malloc(SHUTTLE_PUBLICKEYBYTES);
    uint8_t *sk  = malloc(SHUTTLE_SECRETKEYBYTES);
    uint8_t *sig = malloc(SHUTTLE_BYTES);
    uint8_t *sig_bad = malloc(SHUTTLE_BYTES);
    uint8_t *msg = malloc(MAX_MLEN);
    uint8_t *msg_bad = malloc(MAX_MLEN);
    if (!pk || !sk || !sig || !sig_bad || !msg || !msg_bad) {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    for (unsigned i = 0; i < nrounds; ++i) {
        uint8_t kg_seed[48], sg_seed[48];
        derive_seed(kg_seed, MASTER_SEED, "keygen", i);
        derive_seed(sg_seed, MASTER_SEED, "sign",   i);

        size_t mlen = MLEN_CYCLE[i % MLEN_CYCLE_LEN];
        derive_msg(msg, mlen, MASTER_SEED, i);

        /* ---- KeyGen ---- */
        randombytes_seed(kg_seed);
        if (crypto_sign_keypair(pk, sk) != 0) {
            fprintf(stderr, "keypair failed at round %u\n", i);
            return 1;
        }

        fprintf(f_kg, "count = %u\n", i);
        fprintf(f_kg, "seed = ");  hex_fprint(f_kg, kg_seed, 48);  fprintf(f_kg, "\n");
        fprintf(f_kg, "pk = ");    hex_fprint(f_kg, pk, SHUTTLE_PUBLICKEYBYTES); fprintf(f_kg, "\n");
        fprintf(f_kg, "sk = ");    hex_fprint(f_kg, sk, SHUTTLE_SECRETKEYBYTES); fprintf(f_kg, "\n\n");

        /* ---- Sign ---- */
        size_t siglen = 0;
        randombytes_seed(sg_seed);
        if (crypto_sign_signature(sig, &siglen, msg, mlen, sk) != 0) {
            fprintf(stderr, "sign failed at round %u\n", i);
            return 1;
        }

        fprintf(f_sg, "count = %u\n", i);
        fprintf(f_sg, "sk = ");        hex_fprint(f_sg, sk, SHUTTLE_SECRETKEYBYTES); fprintf(f_sg, "\n");
        fprintf(f_sg, "mlen = %zu\n",  mlen);
        fprintf(f_sg, "msg = ");       hex_fprint(f_sg, msg, mlen); fprintf(f_sg, "\n");
        fprintf(f_sg, "sign_seed = "); hex_fprint(f_sg, sg_seed, 48); fprintf(f_sg, "\n");
        fprintf(f_sg, "siglen = %zu\n", siglen);
        fprintf(f_sg, "sig = ");       hex_fprint(f_sg, sig, siglen); fprintf(f_sg, "\n\n");

        /* ---- Verify (3 subcases) ---- */
        /* a) original -> accept */
        if (crypto_sign_verify(sig, siglen, msg, mlen, pk) != 0) {
            fprintf(stderr, "verify(original) rejected at round %u\n", i);
            return 1;
        }
        fprintf(f_vf, "count = %u\n",  i);
        fprintf(f_vf, "subcase = original\n");
        fprintf(f_vf, "pk = ");     hex_fprint(f_vf, pk, SHUTTLE_PUBLICKEYBYTES); fprintf(f_vf, "\n");
        fprintf(f_vf, "mlen = %zu\n", mlen);
        fprintf(f_vf, "msg = ");    hex_fprint(f_vf, msg, mlen); fprintf(f_vf, "\n");
        fprintf(f_vf, "siglen = %zu\n", siglen);
        fprintf(f_vf, "sig = ");    hex_fprint(f_vf, sig, siglen); fprintf(f_vf, "\n");
        fprintf(f_vf, "expected = accept\n\n");

        /* b) bit-flip inside the rANS Z_0 area (offset chosen per
         *    test_sign.c Test 3: c_tilde + irs_signs + 10) -> reject */
        memcpy(sig_bad, sig, siglen);
        size_t flip_off = (size_t)SHUTTLE_CTILDEBYTES + SHUTTLE_IRS_SIGNBYTES + 10;
        if (flip_off >= siglen) flip_off = siglen / 2;
        sig_bad[flip_off] ^= 0x01;
        if (crypto_sign_verify(sig_bad, siglen, msg, mlen, pk) == 0) {
            fprintf(stderr, "verify(sig_tampered) accepted at round %u\n", i);
            return 1;
        }
        fprintf(f_vf, "count = %u\n",  i);
        fprintf(f_vf, "subcase = sig_tampered\n");
        fprintf(f_vf, "pk = ");     hex_fprint(f_vf, pk, SHUTTLE_PUBLICKEYBYTES); fprintf(f_vf, "\n");
        fprintf(f_vf, "mlen = %zu\n", mlen);
        fprintf(f_vf, "msg = ");    hex_fprint(f_vf, msg, mlen); fprintf(f_vf, "\n");
        fprintf(f_vf, "siglen = %zu\n", siglen);
        fprintf(f_vf, "sig = ");    hex_fprint(f_vf, sig_bad, siglen); fprintf(f_vf, "\n");
        fprintf(f_vf, "expected = reject\n\n");

        /* c) MSB-flip on msg[0] -> reject */
        memcpy(msg_bad, msg, mlen);
        msg_bad[0] ^= 0x80;
        if (crypto_sign_verify(sig, siglen, msg_bad, mlen, pk) == 0) {
            fprintf(stderr, "verify(msg_tampered) accepted at round %u\n", i);
            return 1;
        }
        fprintf(f_vf, "count = %u\n",  i);
        fprintf(f_vf, "subcase = msg_tampered\n");
        fprintf(f_vf, "pk = ");     hex_fprint(f_vf, pk, SHUTTLE_PUBLICKEYBYTES); fprintf(f_vf, "\n");
        fprintf(f_vf, "mlen = %zu\n", mlen);
        fprintf(f_vf, "msg = ");    hex_fprint(f_vf, msg_bad, mlen); fprintf(f_vf, "\n");
        fprintf(f_vf, "siglen = %zu\n", siglen);
        fprintf(f_vf, "sig = ");    hex_fprint(f_vf, sig, siglen); fprintf(f_vf, "\n");
        fprintf(f_vf, "expected = reject\n\n");
    }

    fclose(f_kg);
    fclose(f_sg);
    fclose(f_vf);

    free(pk); free(sk); free(sig); free(sig_bad); free(msg); free(msg_bad);

    fprintf(stderr, "[gen_vectors] OK (%u rounds)\n", nrounds);
    return 0;
}
