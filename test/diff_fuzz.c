/*
 * diff_fuzz.c -- cross-backend differential fuzz driver for SHUTTLE (P13 FUZZ-1).
 *
 * For a deterministic schedule of DIFF_FUZZ_N (seed, msg) cases per mode, run
 *   keygen(xi) -> sign(sk, msg, rnd) -> verify(pk, sig, msg)
 * and emit a NORMALIZED, backend-independent transcript: one `CASE` line per
 * case carrying FNV-1a hashes of pk / sk / sig plus siglen and the verify code.
 * `integration_audit_matrix.sh` runs this on ref/avx2/avx512, normalizes the
 * `CASE` / `PASS diff_fuzz` lines, and `cmp`s every backend against the first:
 * any mismatch is the ML-DSA-bug class (functional-clean but backend-divergent
 * on some input, the AABBCC escape). Catches both one-backend-only and
 * all-backend-shared divergences when paired with the Python reference KAT.
 *
 * It ALSO cross-verifies: sign on THIS backend, verify on THIS backend, AND a
 * tampered-sig must reject -- so the transcript also encodes the reject path.
 *
 * Build (per mode m), NGCC default, RAW or rANS:
 *   gcc -O2 -std=c99 -I. -Itest -I../tools -Intt/<qset> -DDISABLE_NAMESPACE=1 \
 *       -DSHUTTLE_MODE=<m> -DDIFF_FUZZ_N=32 \
 *       test/diff_fuzz.c <full sign source list> -o out/diff_fuzz_<m>
 *   (the avx2/avx512 backends compile the SAME driver via -I../ref/test.)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "api.h"
#include "params.h"
#include "rans.h"

int crypto_sign_keypair_xi(uint8_t *pk, uint8_t *sk, const uint8_t xi[SEEDBYTES]);
int crypto_sign_signature_rnd(uint8_t *sig, size_t *siglen, const uint8_t *m,
                              size_t mlen, const uint8_t *sk,
                              const uint8_t rnd[SEEDBYTES]);

#ifndef DIFF_FUZZ_N
#    define DIFF_FUZZ_N 32
#endif

#if defined(SIG_RAW)
#    define SIG_BYTES SIG_RAW_PACKED_BYTES
#else
#    define SIG_BYTES CRYPTO_BYTES
#endif

/* FNV-1a 64-bit over a byte buffer (backend-independent). */
static uint64_t fnv1a(const uint8_t *p, size_t n)
{
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

/* Deterministic per-case seed schedule: SplitMix64 from (mode, case). The seed
 * is the SAME on every backend, so a divergence is a real backend bug. */
static uint64_t splitmix64(uint64_t *s)
{
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void fill_det(uint64_t *s, uint8_t *p, size_t n)
{
    for (size_t i = 0; i < n; i++)
        p[i] = (uint8_t)(splitmix64(s) & 0xFF);
}

int main(void)
{
    uint8_t *pk = malloc(CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk = malloc(CRYPTO_SECRETKEYBYTES);
    uint8_t *sig = malloc(SIG_BYTES);
    uint8_t xi[SEEDBYTES], rnd[SEEDBYTES];
    uint8_t msg[96];
    if (!pk || !sk || !sig) {
        puts("diff_fuzz: OOM");
        return 2;
    }

    printf("== SHUTTLE-%d diff_fuzz (DIFF_FUZZ_N=%d) ==\n", (int)LAMBDA,
           (int)DIFF_FUZZ_N);

    int fails = 0;
    for (int caseno = 0; caseno < DIFF_FUZZ_N; caseno++) {
        /* deterministic seed = (mode<<16) ^ caseno -- identical per backend. */
        uint64_t s = ((uint64_t)LAMBDA << 32) ^ (uint64_t)(0xABCDu) ^
                     (uint64_t)caseno;
        fill_det(&s, xi, sizeof xi);
        fill_det(&s, rnd, sizeof rnd);
        /* variable message length in [1, sizeof msg], deterministic. */
        size_t mlen = 1 + (size_t)(splitmix64(&s) % sizeof msg);
        fill_det(&s, msg, mlen);

        if (crypto_sign_keypair_xi(pk, sk, xi) != 0) {
            /* keygen norm-window can reject this xi; emit a deterministic
             * KEYGEN_REJECT line (identical across backends) and continue. */
            printf("CASE %d KEYGEN_REJECT\n", caseno);
            continue;
        }
        size_t sl = 0;
        int rc = crypto_sign_signature_rnd(sig, &sl, msg, mlen, sk, rnd);
        if (rc != 0) {
            printf("CASE %d SIGN_FAIL rc=%d\n", caseno, rc);
            fails++;
            continue;
        }
        int v = crypto_sign_verify(sig, sl, msg, mlen, pk);
        /* tamper: flip a byte in the z1 region; must reject. */
        uint8_t saved = sig[sl / 2];
        sig[sl / 2] ^= 0x40;
        int vt = crypto_sign_verify(sig, sl, msg, mlen, pk);
        sig[sl / 2] = saved;

        uint64_t hpk = fnv1a(pk, CRYPTO_PUBLICKEYBYTES);
        uint64_t hsk = fnv1a(sk, CRYPTO_SECRETKEYBYTES);
        uint64_t hsig = fnv1a(sig, sl);
        printf("CASE %d mlen=%zu siglen=%zu pk=%016llx sk=%016llx sig=%016llx "
               "verify=%d tamper=%d\n",
               caseno, mlen, sl, (unsigned long long)hpk,
               (unsigned long long)hsk, (unsigned long long)hsig, v, vt);
        if (v != 0) {
            printf("  FAIL: genuine sig rejected (v=%d)\n", v);
            fails++;
        }
        if (vt == 0) {
            printf("  FAIL: tampered sig ACCEPTED (vt=%d)\n", vt);
            fails++;
        }
    }

    free(pk);
    free(sk);
    free(sig);
    if (fails == 0)
        printf("PASS diff_fuzz mode=%d cases=%d\n", (int)LAMBDA,
               (int)DIFF_FUZZ_N);
    else
        printf("FAIL diff_fuzz mode=%d fails=%d\n", (int)LAMBDA, fails);
    return fails == 0 ? 0 : 1;
}
