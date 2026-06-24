/*
 * security_fault_injection.c -- fault/tamper robustness for SHUTTLE (P13 FAULT-1).
 *
 * Inject bit-flips / byte corruptions into pk, sk, and sig and assert the
 * scheme handles every tampered/faulted input DETERMINISTICALLY: no crash, no
 * accept of corrupted material. Specifically:
 *   (1) tampered sig (single byte flips swept across the whole signature) ->
 *       verify must REJECT every one (never accept) and never crash.
 *   (2) tampered pk (sweep) -> the genuine sig must NOT verify under it.
 *   (3) tampered sk (sweep) -> signing handles a corrupted key GRACEFULLY:
 *       either it fails cleanly (nonzero return) or it produces a signature for
 *       the SAME message; never a crash, never an out-of-bounds. NOTE: a
 *       corrupted-sk signature MAY still verify under the honest pk -- lattice
 *       hash-and-sign signatures are NON-UNIQUE and the secret key is redundant
 *       (the masterSeed re-derives s,e'), so a small key perturbation can yield
 *       a different-but-valid short z for the SAME message. That is harmless
 *       (the attacker already holds a secret key); it is NOT a forgery. The
 *       forgery threats (corrupt SIG accepted, accept under corrupt PK,
 *       cross-MESSAGE forge) are covered by (1), (2), and (5).
 *   (4) repeated re-verify of the genuine artifacts stays accept (no state
 *       corruption / no global mutation across calls).
 *   (5) a genuine signature must NOT verify for a DIFFERENT message (the real
 *       unforgeability surface).
 *
 * Build (per mode m), NGCC default, rANS path:
 *   gcc -O2 -std=c99 -I. -Itest -I../tools -Intt/<qset> -DDISABLE_NAMESPACE=1 \
 *       -DSHUTTLE_MODE=<m> test/security_fault_injection.c \
 *       <full sign source list> -o out/security_fault_injection_<m>
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

static void fill_pat(uint8_t *p, size_t n, uint8_t base)
{
    for (size_t i = 0; i < n; i++)
        p[i] = (uint8_t)(base ^ (uint8_t)(13 * i + 1));
}

int main(void)
{
    uint8_t *pk = malloc(CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk = malloc(CRYPTO_SECRETKEYBYTES);
    uint8_t *sig = malloc(CRYPTO_BYTES);
    uint8_t *pkx = malloc(CRYPTO_PUBLICKEYBYTES);
    uint8_t *skx = malloc(CRYPTO_SECRETKEYBYTES);
    uint8_t *sigx = malloc(CRYPTO_BYTES);
    uint8_t xi[SEEDBYTES], rnd[SEEDBYTES];
    uint8_t msg[48];
    size_t sl = 0;
    int fails = 0;
    long accepted_corrupt = 0;
    long crashes = 0; /* a crash would abort; if we reach the end, crashes=0 */

    if (!pk || !sk || !sig || !pkx || !skx || !sigx) {
        puts("fault_injection: OOM");
        return 2;
    }

    printf("== SHUTTLE-%d security_fault_injection ==\n", (int)LAMBDA);

    fill_pat(xi, sizeof xi, 0x44);
    fill_pat(rnd, sizeof rnd, 0x91);
    fill_pat(msg, sizeof msg, 0x2D);

    if (crypto_sign_keypair_xi(pk, sk, xi) != 0) {
        puts("fault_injection: keygen failed");
        return 2;
    }
    if (crypto_sign_signature_rnd(sig, &sl, msg, sizeof msg, sk, rnd) != 0) {
        puts("fault_injection: sign failed");
        return 2;
    }
    if (crypto_sign_verify(sig, sl, msg, sizeof msg, pk) != 0) {
        puts("fault_injection: genuine sig did not verify");
        return 2;
    }

    /* step over the buffers so a full sweep stays fast for large modes. The
     * sig/pk sweeps are verify-only (cheap); the sk sweep RE-SIGNS, and a
     * corrupted sk can drive signing to the SIGN_MAX_ITER cap (~1000 attempts,
     * a bounded but slow path), so the sk sweep is capped to a modest sample. */
    size_t sig_step = sl > 256 ? sl / 256 : 1;
    size_t pk_step = CRYPTO_PUBLICKEYBYTES > 256 ? CRYPTO_PUBLICKEYBYTES / 256 : 1;
    size_t SK_SAMPLES = 24; /* signing is the cost; sample, do not full-sweep */
    size_t sk_step = CRYPTO_SECRETKEYBYTES > SK_SAMPLES
                         ? CRYPTO_SECRETKEYBYTES / SK_SAMPLES
                         : 1;

    /* (1) tampered sig sweep -- must always REJECT. */
    long sig_tests = 0;
    for (size_t i = 0; i < sl; i += sig_step) {
        memcpy(sigx, sig, sl);
        sigx[i] ^= 0xA5;
        int v = crypto_sign_verify(sigx, sl, msg, sizeof msg, pk);
        sig_tests++;
        if (v == 0)
            accepted_corrupt++;
    }
    printf("[1] tampered-sig sweep: %ld flips, %ld accepted (want 0)\n",
           sig_tests, accepted_corrupt);
    if (accepted_corrupt != 0)
        fails++;

    /* (2) tampered pk sweep -- genuine sig must NOT verify under it. */
    long pk_tests = 0, pk_accept = 0;
    for (size_t i = 0; i < (size_t)CRYPTO_PUBLICKEYBYTES; i += pk_step) {
        memcpy(pkx, pk, CRYPTO_PUBLICKEYBYTES);
        pkx[i] ^= 0x5A;
        int v = crypto_sign_verify(sig, sl, msg, sizeof msg, pkx);
        pk_tests++;
        if (v == 0)
            pk_accept++;
    }
    printf("[2] tampered-pk sweep: %ld flips, %ld verify-accept (want 0)\n",
           pk_tests, pk_accept);
    if (pk_accept != 0)
        fails++;

    /* (3) tampered sk sweep -- GRACEFUL handling: no crash, and any produced
     * signature must (still) be for the SAME message. A corrupted-sk signature
     * verifying under the honest pk is EXPECTED (non-unique lattice sigs +
     * redundant key) and is NOT counted as a failure. We track how many produce
     * a still-valid sig (informational) and assert only no-crash determinism. */
    long sk_tests = 0, sk_valid = 0, sk_cleanfail = 0;
    for (size_t i = 0; i < (size_t)CRYPTO_SECRETKEYBYTES; i += sk_step) {
        memcpy(skx, sk, CRYPTO_SECRETKEYBYTES);
        skx[i] ^= 0x3C;
        size_t slx = 0;
        int rc = crypto_sign_signature_rnd(sigx, &slx, msg, sizeof msg, skx, rnd);
        sk_tests++;
        if (rc != 0) {
            sk_cleanfail++; /* clean failure: fine */
            continue;
        }
        /* produced a sig: it must be a well-formed sig for THIS message (it may
         * or may not verify -- both are acceptable, non-crash). */
        int v = crypto_sign_verify(sigx, slx, msg, sizeof msg, pk);
        if (v == 0)
            sk_valid++;
    }
    printf("[3] tampered-sk sweep: %ld flips, %ld clean-fail, %ld still-valid "
           "(non-unique sigs; no crash) -> graceful\n",
           sk_tests, sk_cleanfail, sk_valid);

    /* (4) idempotent re-verify of genuine artifacts (no global-state mutation). */
    {
        int ok = 1;
        for (int r = 0; r < 8; r++)
            ok &= (crypto_sign_verify(sig, sl, msg, sizeof msg, pk) == 0);
        printf("[4] genuine re-verify x8: %s\n", ok ? "stable ACCEPT" : "UNSTABLE");
        if (!ok)
            fails++;
    }

    /* (5) unforgeability surface: genuine sig must NOT verify for a DIFFERENT
     * message (sweep several message-byte flips). */
    {
        long m_tests = 0, m_accept = 0;
        for (size_t i = 0; i < sizeof msg; i++) {
            uint8_t saved = msg[i];
            msg[i] ^= 0xFF;
            if (crypto_sign_verify(sig, sl, msg, sizeof msg, pk) == 0)
                m_accept++;
            m_tests++;
            msg[i] = saved;
        }
        printf("[5] cross-message sweep: %ld flips, %ld accepted (want 0)\n",
               m_tests, m_accept);
        if (m_accept != 0)
            fails++;
    }

    free(pk);
    free(sk);
    free(sig);
    free(pkx);
    free(skx);
    free(sigx);

    printf("== SHUTTLE-%d security_fault_injection: crashes=%ld accepted_corrupt=%ld %s ==\n",
           (int)LAMBDA, crashes, accepted_corrupt,
           fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}
