/*
 * parser_negative.c -- negative parser corpus for SHUTTLE verify (P13 PARSE-1).
 *
 * The verifier on a malformed signature MUST be a DETERMINISTIC PUBLIC REJECT
 * (crypto_sign_verify != 0), never a crash and never an accept of corrupted
 * material. This cleanly separates verify-side parser branches from secret
 * branches (SECRET_PUBLIC_AUDIT Verify/rANS/Packing). The corpus covers:
 *   - wrong length (too short / too long / off-by-one) -> usage error (-2)
 *   - flipped seedC bytes -> challenge mismatch -> reject (-1)
 *   - zeroed / all-0xFF signature body -> reject
 *   - rANS-region mutations (K15: container length byte, reserve/padding, the
 *     interior state bytes -- at least 6 distinct mutations) -> reject
 *   - hint-region mutations targeting the non-power-of-2 hint range (K14: a
 *     decoded hint >= H_h must REJECT, NEVER be reduced mod H_h)
 * Every case asserts a nonzero return and (under no sanitizer crash) survival.
 *
 * Build (per mode m), default rANS path (NOT -DSIG_RAW: the rANS K15 cases
 * need the real wire format):
 *   gcc -O2 -std=c99 -I. -Itest -I../tools -Intt/<qset> -DDISABLE_NAMESPACE=1 \
 *       -DSHUTTLE_MODE=<m> test/parser_negative.c <full sign source list> \
 *       -o out/parser_negative_<m>
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

static int g_fails;
static int g_cases;

static void fill_pat(uint8_t *p, size_t n, uint8_t base)
{
    for (size_t i = 0; i < n; i++)
        p[i] = (uint8_t)(base ^ (uint8_t)(7 * i));
}

/* Assert: verifying (sig[0..len)) against (pk, msg) REJECTS (nonzero). */
static void expect_reject(const char *name, const uint8_t *sig, size_t len,
                          const uint8_t *msg, size_t mlen, const uint8_t *pk)
{
    g_cases++;
    int v = crypto_sign_verify(sig, len, msg, mlen, pk);
    if (v == 0) {
        printf("  FAIL[%s]: malformed sig ACCEPTED (v=0)\n", name);
        g_fails++;
    } else {
        printf("  ok[%s]: reject v=%d\n", name, v);
    }
}

int main(void)
{
    uint8_t *pk = malloc(CRYPTO_PUBLICKEYBYTES);
    uint8_t *sk = malloc(CRYPTO_SECRETKEYBYTES);
    uint8_t *sig = malloc(CRYPTO_BYTES);
    uint8_t *mut = malloc(CRYPTO_BYTES + 16);
    uint8_t xi[SEEDBYTES], rnd[SEEDBYTES];
    uint8_t msg[40];
    size_t sl = 0;
    if (!pk || !sk || !sig || !mut) {
        puts("parser_negative: OOM");
        return 2;
    }

    printf("== SHUTTLE-%d parser_negative ==\n", (int)LAMBDA);

    fill_pat(xi, sizeof xi, 0x31);
    fill_pat(rnd, sizeof rnd, 0x57);
    fill_pat(msg, sizeof msg, 0xC3);

    /* Get one genuine (pk, sig, msg) over the rANS wire format. */
    if (crypto_sign_keypair_xi(pk, sk, xi) != 0) {
        puts("parser_negative: keygen failed");
        return 2;
    }
    if (crypto_sign_signature_rnd(sig, &sl, msg, sizeof msg, sk, rnd) != 0) {
        puts("parser_negative: sign failed");
        return 2;
    }
    /* sanity: the genuine sig verifies (positive control). */
    if (crypto_sign_verify(sig, sl, msg, sizeof msg, pk) != 0) {
        puts("parser_negative: genuine sig did NOT verify (cannot test)");
        return 2;
    }
    printf("positive control: genuine sig (len=%zu) verifies\n", sl);

    /* ---- length mutations ---- */
    expect_reject("len-1", sig, sl - 1, msg, sizeof msg, pk);
    expect_reject("len+1", sig, sl + 1, msg, sizeof msg, pk);
    expect_reject("len/2", sig, sl / 2, msg, sizeof msg, pk);
    expect_reject("len0", sig, 0, msg, sizeof msg, pk);
    expect_reject("len=CRYPTO_BYTES", sig, CRYPTO_BYTES, msg, sizeof msg, pk);

    /* ---- seedC (challenge) mutations: byte flips in [0, CHALLENGESEEDBYTES) ---- */
    memcpy(mut, sig, sl);
    mut[0] ^= 0x01;
    expect_reject("seedC[0]^1", mut, sl, msg, sizeof msg, pk);
    memcpy(mut, sig, sl);
    mut[CHALLENGESEEDBYTES - 1] ^= 0x80;
    expect_reject("seedC[last]^0x80", mut, sl, msg, sizeof msg, pk);

    /* ---- whole-body zero / all-ones ---- */
    memset(mut, 0x00, sl);
    expect_reject("all-zero", mut, sl, msg, sizeof msg, pk);
    memset(mut, 0xFF, sl);
    expect_reject("all-ones", mut, sl, msg, sizeof msg, pk);

    /* ---- rANS-region mutations (K15): the bytes AFTER seedC are the rANS
     * container (rlen + state + reserve/padding). At least 6 distinct
     * single/multi-byte mutations into the body -> decode-fail / non-canonical
     * state / support / re-encode mismatch -> REJECT. ---- */
    {
        size_t body_off = CHALLENGESEEDBYTES; /* first rANS container byte */
        size_t span = sl - body_off;
        size_t pts[6];
        pts[0] = body_off;                 /* container length/header */
        pts[1] = body_off + 1;
        pts[2] = body_off + span / 4;
        pts[3] = body_off + span / 2;      /* interior state */
        pts[4] = body_off + (3 * span) / 4;
        pts[5] = sl - 1;                   /* trailing reserve/padding */
        for (int i = 0; i < 6; i++) {
            char nm[32];
            memcpy(mut, sig, sl);
            mut[pts[i]] ^= 0x55; /* flip several bits in one body byte */
            snprintf(nm, sizeof nm, "rANS-body@%zu", pts[i]);
            expect_reject(nm, mut, sl, msg, sizeof msg, pk);
        }
        /* zero the entire rANS container (decode must fail, never wrap). */
        memcpy(mut, sig, sl);
        memset(mut + body_off, 0x00, span);
        expect_reject("rANS-body-zeroed", mut, sl, msg, sizeof msg, pk);
        /* set the trailing reserve/padding nonzero (non-canonical). */
        memcpy(mut, sig, sl);
        mut[sl - 1] |= 0x7F;
        mut[sl - 2] |= 0x7F;
        expect_reject("rANS-padding-nonzero", mut, sl, msg, sizeof msg, pk);
    }

    /* ---- hint-range mutation (K14): drive the decoded hint out of [0,H_h).
     * We flip high bits across the body region that holds the hint encoding so
     * a decoded value lands in the non-empty gap [H_h, 2^d_h); the parser must
     * RANGE-CHECK and REJECT, never reduce mod H_h. (The exact byte is format-
     * dependent; we sweep a window and assert at least the verify rejects.) ---- */
    {
        size_t body_off = CHALLENGESEEDBYTES;
        for (size_t off = body_off; off < sl; off += (sl - body_off) / 4 + 1) {
            char nm[40];
            memcpy(mut, sig, sl);
            mut[off] = 0xFE; /* push the field high (toward >= H_h) */
            snprintf(nm, sizeof nm, "hint-range@%zu(K14)", off);
            expect_reject(nm, mut, sl, msg, sizeof msg, pk);
        }
    }

    /* ---- wrong message (sig genuine but bound to a different msg) ---- */
    {
        uint8_t msg2[40];
        memcpy(msg2, msg, sizeof msg2);
        msg2[0] ^= 0xFF;
        expect_reject("wrong-msg", sig, sl, msg2, sizeof msg2, pk);
    }

    free(pk);
    free(sk);
    free(sig);
    free(mut);

    printf("== SHUTTLE-%d parser_negative: %d cases, %s ==\n", (int)LAMBDA,
           g_cases, g_fails == 0 ? "PASS" : "FAIL");
    return g_fails == 0 ? 0 : 1;
}
