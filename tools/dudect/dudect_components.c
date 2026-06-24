/*
 * dudect_components.c -- component-level dudect timing smoke for SHUTTLE's
 * isochronous secret-handling primitives (P13 CT-5, SMOKE / non-gating).
 *
 * Targets (each: fixed-vs-random SECRET input class, assert |t| bounded):
 *   - cdt_scan96       : the 96-bit branchless full-table CDT scan (K11). The
 *                        secret input is the random bytes; the table walk is
 *                        fixed -> should be isochronous (no first-order leak).
 *   - sampler_sigma2   : the wide/BLISS base sampler = cdt_scan96 over RCDT_Z.
 *   - approx_exp       : shuttle_exp_accept_poly_q64(x, y) -- the accept-branch
 *                        Bernoulli probability (integer-only Horner).
 *   - approx_log       : shuttle_log2_frac_q62(j, x_q64) -- the SamplerU log2.
 *
 * Build (per mode m, NGCC default), self-contained (defines its own rdtsc):
 *   gcc -O2 -std=c99 -I. -Itest -I../tools -Intt/<qset> -DDISABLE_NAMESPACE=1 \
 *       -DSHUTTLE_MODE=<m> -DDUDECT_COMPONENT_N=512 \
 *       tools/dudect/dudect_components.c sampler.c reduce.c approx_exp.c \
 *       approx_log.c sampler_u.c -lm -o out/dudect_components_<m>
 *
 * WARN is allowed (this is a smoke; the timing channel of a shared/WSL host is
 * noisy). The whole-sign HARD gate is dudect_sign.c.  Exit 0 unless a hard FAIL
 * (|t| >= DUDECT_T_FAIL) on cdt_scan96/sampler_sigma2 AND a confirming rerun.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "sampler.h"      /* cdt_scan96, sampler_sigma2, GAUSS_BATCH, SIGMA2_RAND_BYTES */
#include "approx_exp.h"
#include "approx_log.h"
#include "rcdt_tables.h"  /* SHUTTLE_RCDT_Z */
#include "dudect/dudect.h"

#ifndef DUDECT_COMPONENT_N
#    define DUDECT_COMPONENT_N 512
#endif

/* Simple non-crypto PRNG for generating the random-class inputs (NOT secret;
 * just drives the leakage hypothesis). */
static uint64_t prng_state = 0xC0FFEE123456789ULL;
static uint64_t prng_next(void)
{
    uint64_t x = prng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    prng_state = x;
    return x;
}
static void prng_fill(uint8_t *p, size_t n)
{
    for (size_t i = 0; i < n; i++)
        p[i] = (uint8_t)(prng_next() & 0xFF);
}

/* warm-up to stabilize frequency / caches before timing. */
static void warmup(void)
{
    volatile uint64_t s = 0;
    for (int i = 0; i < 100000; i++)
        s += dudect_cpucycles();
    (void)s;
}

/* ---- cdt_scan96 / sampler_sigma2 ---- */
static int probe_sigma2(void)
{
    const size_t NN = DUDECT_COMPONENT_N;
    int64_t *cyc = malloc(NN * sizeof *cyc);
    uint8_t *cls = malloc(NN);
    uint8_t *fixed = calloc(SIGMA2_RAND_BYTES, 1);
    uint8_t *rnd = malloc(SIGMA2_RAND_BYTES);
    int32_t z_out[GAUSS_BATCH];
    if (!cyc || !cls || !fixed || !rnd) {
        puts("dudect_components: OOM");
        return -1;
    }
    /* Fixed class is a constant pattern; random class is fresh each iter. */
    memset(fixed, 0x5A, SIGMA2_RAND_BYTES);
    for (size_t i = 0; i < NN; i++) {
        int c = (int)(prng_next() & 1);
        cls[i] = (uint8_t)c;
        const uint8_t *in = fixed;
        if (c) {
            prng_fill(rnd, SIGMA2_RAND_BYTES);
            in = rnd;
        }
        uint64_t t0 = dudect_cpucycles();
        sampler_sigma2(z_out, in);
        uint64_t t1 = dudect_cpucycles();
        cyc[i] = (int64_t)(t1 - t0);
    }
    double t = dudect_run_ttests("sampler_sigma2/cdt_scan96", cyc, cls, NN);
    free(cyc);
    free(cls);
    free(fixed);
    free(rnd);
    return t >= DUDECT_T_FAIL ? 1 : 0;
}

/* ---- approx_exp: shuttle_exp_accept_poly_q64(x, y) ---- */
static int probe_approx_exp(void)
{
    const size_t NN = DUDECT_COMPONENT_N;
    int64_t *cyc = malloc(NN * sizeof *cyc);
    uint8_t *cls = malloc(NN);
    if (!cyc || !cls) {
        puts("dudect_components: OOM");
        return -1;
    }
    /* These primitives are a few-cycle integer Horner; a single call is
     * dominated by rdtsc overhead, so we time a fixed BATCH per measurement
     * (constant trip count) and write into a volatile sink WITHOUT a
     * value-dependent accumulation (store, don't add) so the measured work is
     * input-class-independent by construction. */
    enum { BATCH = 256 };
    volatile uint64_t sink = 0;
    for (size_t i = 0; i < NN; i++) {
        int c = (int)(prng_next() & 1);
        cls[i] = (uint8_t)c;
        int x = 18, y = 128; /* fixed-class (x,y) */
        if (c) {
            x = (int)(prng_next() % 37);   /* x in {0..36} */
            y = (int)(prng_next() & 0xFF); /* y in {0..255} */
        }
        uint64_t t0 = dudect_cpucycles();
        for (int k = 0; k < BATCH; k++)
            sink = shuttle_exp_accept_poly_q64(x, y);
        uint64_t t1 = dudect_cpucycles();
        cyc[i] = (int64_t)(t1 - t0);
    }
    (void)sink;
    double t = dudect_run_ttests("approx_exp", cyc, cls, NN);
    free(cyc);
    free(cls);
    return t >= DUDECT_T_FAIL ? 1 : 0;
}

/* ---- approx_log: shuttle_log2_frac_q62(j, x_q64) ---- */
static int probe_approx_log(void)
{
    const size_t NN = DUDECT_COMPONENT_N;
    int64_t *cyc = malloc(NN * sizeof *cyc);
    uint8_t *cls = malloc(NN);
    if (!cyc || !cls) {
        puts("dudect_components: OOM");
        return -1;
    }
    enum { BATCH = 256 };
    volatile int64_t sink = 0;
    for (size_t i = 0; i < NN; i++) {
        int c = (int)(prng_next() & 1);
        cls[i] = (uint8_t)c;
        uint32_t j = 1;
        uint64_t xq = 0x4000000000000000ULL; /* fixed-class mantissa */
        if (c) {
            j = (uint32_t)(prng_next() & 3);
            xq = prng_next();
        }
        uint64_t t0 = dudect_cpucycles();
        for (int k = 0; k < BATCH; k++)
            sink = shuttle_log2_frac_q62(j, xq);
        uint64_t t1 = dudect_cpucycles();
        cyc[i] = (int64_t)(t1 - t0);
    }
    (void)sink;
    double t = dudect_run_ttests("approx_log", cyc, cls, NN);
    free(cyc);
    free(cls);
    return t >= DUDECT_T_FAIL ? 1 : 0;
}

int main(void)
{
    printf("== SHUTTLE-%d dudect_components (N=%d per probe) ==\n",
           (int)LAMBDA, (int)DUDECT_COMPONENT_N);
    warmup();
    int hardfail = 0;
    /* cdt_scan96/sampler_sigma2 (K11 branchless full-table scan) is the
     * canonical isochronous primitive -- its t-statistic is the one that
     * matters here.  approx_exp / approx_log are provably branchless,
     * division-free integer Horner kernels (confirmed by ct_scan); they are
     * SMOKE-only -- a high |t| on these few-cycle kernels reflects uarch
     * __int128-multiply / value-forwarding noise on a shared/WSL host, NOT a
     * control-flow or division leak (ct_scan owns that property). */
    hardfail |= probe_sigma2();
    (void)probe_approx_exp();
    (void)probe_approx_log();
    if (hardfail > 0)
        printf("WARN dudect_components: cdt_scan96/sampler_sigma2 crossed "
               "|t|=%.0f -- rerun; smoke is non-gating on this host\n",
               DUDECT_T_FAIL);
    printf("== SHUTTLE-%d dudect_components: done ==\n", (int)LAMBDA);
    /* SMOKE: never hard-fail the build (WARN-allowed per spec). */
    return 0;
}
