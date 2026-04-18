/*
 * sampler_4x.c - 4-way parallel Gaussian sampling using shake256x4.
 *
 * Strategy: large initial 4-way squeeze (enough for ~all 256 samples),
 * then sequential mini-batch processing per lane. Rare continuation
 * squeezes if rejection causes a lane to need more bytes.
 */

#include "sampler.h"
#include "approx_exp.h"
#include <string.h>

/*
 * Buffer per lane: enough for ~308 attempts (256 accepted at ~83%).
 * Per-mini-batch bytes (16 attempts): SIGMA2_RAND_BYTES + Y_RAND_BYTES +
 *   16 * GAUSS_RAND_BYTES = 192 + 12 + 128 = 332  (~20.75 bytes / attempt)
 * sign_bytes(32) + ceil(308/16) * 332 = 32 + 20 * 332 = 6672 bytes.
 * ceil(6672 / 136) = 50 blocks. Use 56 blocks = 7616 bytes for margin.
 */
#define INIT_4X_NBLOCKS 56
#define BUF_4X_SIZE     ((INIT_4X_NBLOCKS + SQUEEZE_NBLOCKS) * STREAM256_BLOCKBYTES + 32)

/* ============================================================
 * Batched uniform y sampler (same as sampler.c; duplicated static to
 * match existing per-file pattern and keep linkage minimal).
 * See sampler.c for the AVX2 analysis comment on why this is scalar.
 * ============================================================ */
static void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand) {
    const uint32_t y_mask = (1U << SHUTTLE_Y_BITS) - 1U;
    uint64_t lo = load_le64(rand);
    uint32_t hi = load_le32(rand + 8);

    for (int i = 0; i < 10; i++)
        y_out[i] = (uint8_t)((lo >> (SHUTTLE_Y_BITS * i)) & y_mask);

    y_out[10] = (uint8_t)(((lo >> 60) | ((uint64_t)hi << 4)) & y_mask);

    for (int i = 11; i < 16; i++) {
        unsigned shift = (unsigned)(SHUTTLE_Y_BITS * i - 64);
        y_out[i] = (uint8_t)((hi >> shift) & y_mask);
    }
}

/* ============================================================
 * Single Gaussian attempt (same math as sampler.c).
 * See sampler.c (ref or avx2) for the full error analysis of the
 * non-shift rejection path.
 * ============================================================ */
static int sample_gauss_attempt(int16_t *r,
                                uint32_t x,
                                uint32_t y,
                                const uint8_t *rej_rand,
                                const uint8_t *signs,
                                size_t idx) {
    int32_t candidate = ((int32_t)x << SHUTTLE_K_BITS) | (int32_t)y;

    uint64_t rand_tail    = load_le64(rej_rand);
    uint8_t  sign_r0      = (uint8_t)(rand_tail & 1U);
    uint64_t rand_rej_63  = rand_tail >> 1;

    uint32_t t   = (uint32_t)y + ((uint32_t)x << SHUTTLE_TWO_K_BITS);
    uint64_t num = (uint64_t)y * (uint64_t)t;

#if SHUTTLE_SIGMA == 128
    uint64_t exp_q60 = num << 45;
#else
    uint64_t exp_q60 = mulh64((uint64_t)num << 44, SHUTTLE_INV_2SIGMA2_Q64) << 16;
#endif

    uint64_t exp_val_q63 = approx_exp(exp_q60);
    int accepted = (rand_rej_63 < exp_val_q63) ? 1 : 0;

    uint64_t r_is_zero = (candidate == 0) ? 1ULL : 0ULL;
    accepted &= (int)(1 - (r_is_zero & (1 - (uint64_t)sign_r0)));

    uint8_t sign = (signs[idx / 8] >> (idx % 8)) & 1;
    int16_t cand16 = (int16_t)candidate;
    *r = sign ? (int16_t)(-cand16) : cand16;

    return accepted;
}

void sample_gauss_N_4x(int16_t *r0, int16_t *r1,
                        int16_t *r2, int16_t *r3,
                        const uint8_t seed[SHUTTLE_SEEDBYTES],
                        uint64_t nonce0, uint64_t nonce1,
                        uint64_t nonce2, uint64_t nonce3,
                        size_t len0, size_t len1,
                        size_t len2, size_t len3) {
    /* Prepare 4 input buffers: seed || nonce (LE 8 bytes) */
    uint8_t inbuf[4][SHUTTLE_SEEDBYTES + 8];
    uint64_t nonces[4] = {nonce0, nonce1, nonce2, nonce3};
    size_t len[4] = {len0, len1, len2, len3};

    for (int j = 0; j < 4; j++) {
        memcpy(inbuf[j], seed, SHUTTLE_SEEDBYTES);
        for (int i = 0; i < 8; i++)
            inbuf[j][SHUTTLE_SEEDBYTES + i] = (uint8_t)(nonces[j] >> (8 * i));
    }

    /* 4-way SHAKE-256 absorb */
    keccakx4_state state4x;
    shake256x4_absorb_once(&state4x,
                           inbuf[0], inbuf[1], inbuf[2], inbuf[3],
                           SHUTTLE_SEEDBYTES + 8);

    /* Per-lane buffers: large initial squeeze */
    static uint8_t buf[4][BUF_4X_SIZE];
    shake256x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3],
                             INIT_4X_NBLOCKS, &state4x);

    int16_t *r[4]        = {r0, r1, r2, r3};
    size_t pos[4], avail[4], coefcnt[4];

    /* Extract signs and initialize per-lane state */
    uint8_t signs[4][SHUTTLE_N / 8];
    for (int j = 0; j < 4; j++) {
        size_t sb = len[j] / 8;
        if (len[j] > 0)
            memcpy(signs[j], buf[j], sb);
        pos[j] = sb;
        avail[j] = INIT_4X_NBLOCKS * STREAM256_BLOCKBYTES - sb;
        coefcnt[j] = 0;
    }

    /* Mini-batch processing */
    while (coefcnt[0] < len0 || coefcnt[1] < len1 ||
           coefcnt[2] < len2 || coefcnt[3] < len3) {

        /* Continuation squeeze if any active lane runs low (rare) */
        int need_refill = 0;
        for (int j = 0; j < 4; j++) {
            if (coefcnt[j] < len[j] &&
                avail[j] < (size_t)MINIBATCH_RAND_BYTES) {
                need_refill = 1;
                break;
            }
        }
        if (need_refill) {
            for (int j = 0; j < 4; j++) {
                if (avail[j] > 0 && pos[j] > 0)
                    memmove(buf[j], buf[j] + pos[j], avail[j]);
                pos[j] = 0;
            }
            shake256x4_squeezeblocks(
                buf[0] + avail[0], buf[1] + avail[1],
                buf[2] + avail[2], buf[3] + avail[3],
                SQUEEZE_NBLOCKS, &state4x);
            for (int j = 0; j < 4; j++)
                avail[j] += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }

        /* Process one mini-batch per active lane */
        for (int j = 0; j < 4; j++) {
            if (coefcnt[j] >= len[j]) continue;
            if (avail[j] < (size_t)MINIBATCH_RAND_BYTES) continue;

            int16_t z[GAUSS_BATCH];
            uint8_t y[GAUSS_BATCH];

            sampler_sigma2(z, buf[j] + pos[j]);
            pos[j]   += SIGMA2_RAND_BYTES;
            avail[j] -= SIGMA2_RAND_BYTES;

            sampler_y(y, buf[j] + pos[j]);
            pos[j]   += Y_RAND_BYTES;
            avail[j] -= Y_RAND_BYTES;

            for (int k = 0; k < GAUSS_BATCH && coefcnt[j] < len[j]; k++) {
                int accepted = sample_gauss_attempt(
                    &r[j][coefcnt[j]],
                    (uint32_t)z[k],
                    (uint32_t)y[k],
                    buf[j] + pos[j],
                    signs[j], coefcnt[j]);
                pos[j]   += GAUSS_RAND_BYTES;
                avail[j] -= GAUSS_RAND_BYTES;

                if (accepted) {
                    coefcnt[j]++;
                }
            }
        }
    }
}
