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
 * sign_bytes(32) + 308 * (192/16 + 9) = 32 + 308*21 = 6500 bytes
 * ceil(6500/136) = 48 blocks. Use 56 blocks = 7616 bytes for margin.
 *
 * Note: NGCC has smaller per-attempt bytes (9 vs 20) than SM_Sign,
 * so we need fewer blocks.
 */
#define INIT_4X_NBLOCKS 56
#define BUF_4X_SIZE     ((INIT_4X_NBLOCKS + SQUEEZE_NBLOCKS) * STREAM256_BLOCKBYTES + 32)

/* ============================================================
 * Single Gaussian attempt (same as sampler.c)
 * ============================================================ */
static int sample_gauss_sigma128(int16_t *r, const uint8_t *rand,
                                 int16_t x, const uint8_t *signs,
                                 size_t idx) {
    uint32_t y = rand[0] & 0x3FU;
    int32_t candidate = (int32_t)x * NGCC_SIGN_K + (int32_t)y;

    uint64_t rand_tail = load_le64(rand + 1);
    uint8_t sign_r0 = rand_tail & 1;
    uint64_t rand_rej_63 = rand_tail >> 1;

    uint32_t t = y + (uint32_t)(2 * NGCC_SIGN_K) * (uint32_t)x;
    uint64_t num = (uint64_t)y * (uint64_t)t;
    uint64_t exp_q60 = num << 45;

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
                        const uint8_t seed[NGCC_SIGN_SEEDBYTES],
                        uint64_t nonce0, uint64_t nonce1,
                        uint64_t nonce2, uint64_t nonce3,
                        size_t len0, size_t len1,
                        size_t len2, size_t len3) {
    /* Prepare 4 input buffers: seed || nonce (LE 8 bytes) */
    uint8_t inbuf[4][NGCC_SIGN_SEEDBYTES + 8];
    uint64_t nonces[4] = {nonce0, nonce1, nonce2, nonce3};
    size_t len[4] = {len0, len1, len2, len3};

    for (int j = 0; j < 4; j++) {
        memcpy(inbuf[j], seed, NGCC_SIGN_SEEDBYTES);
        for (int i = 0; i < 8; i++)
            inbuf[j][NGCC_SIGN_SEEDBYTES + i] = (uint8_t)(nonces[j] >> (8 * i));
    }

    /* 4-way SHAKE-256 absorb */
    keccakx4_state state4x;
    shake256x4_absorb_once(&state4x,
                           inbuf[0], inbuf[1], inbuf[2], inbuf[3],
                           NGCC_SIGN_SEEDBYTES + 8);

    /* Per-lane buffers: large initial squeeze */
    static uint8_t buf[4][BUF_4X_SIZE];
    shake256x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3],
                             INIT_4X_NBLOCKS, &state4x);

    int16_t *r[4]        = {r0, r1, r2, r3};
    size_t pos[4], avail[4], coefcnt[4];

    /* Extract signs and initialize per-lane state */
    uint8_t signs[4][NGCC_SIGN_N / 8];
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

            int16_t z[NGCC_GAUSS_BATCH];
            sampler_sigma2(z, buf[j] + pos[j]);
            pos[j] += SIGMA2_RAND_BYTES;
            avail[j] -= SIGMA2_RAND_BYTES;

            for (int k = 0; k < NGCC_GAUSS_BATCH && coefcnt[j] < len[j]; k++) {
                int accepted = sample_gauss_sigma128(
                    &r[j][coefcnt[j]], buf[j] + pos[j], z[k],
                    signs[j], coefcnt[j]);
                pos[j] += GAUSS_RAND_BYTES;
                avail[j] -= GAUSS_RAND_BYTES;

                if (accepted) {
                    coefcnt[j]++;
                }
            }
        }
    }
}
