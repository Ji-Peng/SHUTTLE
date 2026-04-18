/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE (ref version).
 *
 * sigma = 128 = 2^7. Decomposition: k = 2^6, sigma_2 = 2.
 * Output: signed int16_t samples with |r| <= 11*128 = 1408.
 *
 * Unified buffer design: one SHAKE-256 stream per sample_gauss_N call
 * produces ALL randomness (signs, CDT, y, rejection). Processing in
 * mini-batches of GAUSS_BATCH attempts for low memory (~0.7 KB buffer).
 *
 * KAT-compatible with the AVX2 implementation: same SHAKE consumption order.
 */

#include "sampler.h"
#include "approx_exp.h"
#include <string.h>

/* ============================================================
 * RCDT Table: 93-bit cumulative probabilities, 3x31-bit limb format.
 *
 * Entry i = Pr[|x| >= i+1] * 2^93, for i = 0..21.
 * Layout: [0] = lowest 31 bits, [1] = middle 31 bits, [2] = upper 31 bits.
 * Identical to SM_Sign. Renyi Divergence: 1 + 2^{-93.85}
 * ============================================================ */
static const uint32_t RCDT_3x31[RCDT_ENTRIES][3] = {
    {0x0204FE5FU, 0x4E73D911U, 0x556D69B6U},
    {0x7624991FU, 0x1BCB49D6U, 0x2FDB7191U},
    {0x26A118E0U, 0x7DB923CBU, 0x16091DCDU},
    {0x70C18AE5U, 0x536A5628U, 0x0836DCFDU},
    {0x54119E8DU, 0x3D589AD0U, 0x0273E65EU},
    {0x2EB0249EU, 0x2EBE8BE2U, 0x00950CA4U},
    {0x0472F0ECU, 0x65A6660AU, 0x001BFA1DU},
    {0x70C91274U, 0x5945CB2DU, 0x000422EEU},
    {0x0517A411U, 0x4A81EF20U, 0x00007AFAU},
    {0x3A4CC14BU, 0x79C6FA82U, 0x00000B31U},
    {0x12B0F081U, 0x185517F0U, 0x000000CCU},
    {0x6073B7BEU, 0x2FB676CEU, 0x0000000BU},
    {0x1F24D098U, 0x3F53D508U, 0x00000000U},
    {0x5FBBC714U, 0x02267C01U, 0x00000000U},
    {0x0A24AB07U, 0x000E9522U, 0x00000000U},
    {0x61E4FCE8U, 0x00004D1EU, 0x00000000U},
    {0x7B65A4CBU, 0x0000013DU, 0x00000000U},
    {0x7EE574B9U, 0x00000003U, 0x00000000U},
    {0x04FF6C73U, 0x00000000U, 0x00000000U},
    {0x0009C080U, 0x00000000U, 0x00000000U},
    {0x00000ED2U, 0x00000000U, 0x00000000U},
    {0x00000011U, 0x00000000U, 0x00000000U}
};

/* ============================================================
 * Batched CDT sampler: sampler_sigma2 (pure function)
 *
 * Random data layout (AVX2-friendly, groups of 8 samples):
 *   Group k (samples 8k..8k+7):
 *     32 bytes: uint32_t[8] low  31-bit limbs
 *     32 bytes: uint32_t[8] mid  31-bit limbs
 *     32 bytes: uint32_t[8] high 31-bit limbs
 *
 * Internal computation in 32-bit. Output stored as int16_t.
 * ============================================================ */
int sampler_sigma2(int16_t *z_out, const uint8_t *rand) {
    for (int s = 0; s < GAUSS_BATCH; s++) {
        int group = s >> 3;          /* s / 8 */
        int lane  = s & 7;           /* s % 8 */
        int byte_base = group * 96 + lane * 4;

        uint32_t v0 = load_le32(rand + byte_base)      & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + byte_base + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + byte_base + 64) & 0x7FFFFFFFU;

        int32_t z = 0;
        for (int i = 0; i < RCDT_ENTRIES; i++) {
            uint32_t cc;
            cc = ct_lt_u32(v0, RCDT_3x31[i][0]);
            cc = ct_lt_u32(v1 - cc, RCDT_3x31[i][1]);
            cc = ct_lt_u32(v2 - cc, RCDT_3x31[i][2]);
            z += (int32_t)cc;
        }
        z_out[s] = (int16_t)z;
    }
    return GAUSS_BATCH;
}

/* ============================================================
 * Single Gaussian attempt: sample_gauss_sigma128
 *
 * CDT value x passed directly (from sampler_sigma2 output).
 * rand[0..8]: y (1 byte, 6 bits) + sign_r0+rejection (8 bytes).
 *
 * Returns 1 if accepted, 0 otherwise.
 * Sets *r to signed candidate (sign from pre-extracted sign bits).
 * ============================================================ */
static int sample_gauss_sigma128(int16_t *r, const uint8_t *rand,
                                 int16_t x, const uint8_t *signs,
                                 size_t idx) {
    /* y from rand[0], 6 bits */
    uint32_t y = rand[0] & 0x3FU;

    /* candidate = x * 64 + y (unsigned magnitude) */
    int32_t candidate = (int32_t)x * SHUTTLE_GAUSS_K + (int32_t)y;

    /* sign_r0 + rejection from rand[1..8] */
    uint64_t rand_tail = load_le64(rand + 1);
    uint8_t sign_r0 = rand_tail & 1;
    uint64_t rand_rej_63 = rand_tail >> 1;

    /* Rejection test:
     * exponent = y * (y + 2*k*x) / (2*sigma^2)
     * k = 64, sigma = 128, 2*sigma^2 = 32768 = 2^15
     * In Q60: exp_q60 = y * (y + 128*x) * 2^60 / 2^15 = y * (y + 128*x) << 45
     *
     * Max: 63*(63+128*22) = 63*2879 = 181377
     * 181377 << 45 = 6.38e18 < 2^63, fits in uint64_t.
     */
    uint32_t t = y + (uint32_t)(2 * SHUTTLE_GAUSS_K) * (uint32_t)x;
    uint64_t num = (uint64_t)y * (uint64_t)t;
    uint64_t exp_q60 = num << 45;  /* / 2^15 * 2^60 = << 45 */

    uint64_t exp_val_q63 = approx_exp(exp_q60);
    int accepted = (rand_rej_63 < exp_val_q63) ? 1 : 0;

    /* r = 0 case: accept with prob 1/2.
     * When candidate=0, exp(-0)=1, rejection always passes.
     * Must reject half to correct for folded distribution. */
    uint64_t r_is_zero = (candidate == 0) ? 1ULL : 0ULL;
    accepted &= (int)(1 - (r_is_zero & (1 - (uint64_t)sign_r0)));

    /* Apply sign from pre-extracted sign bits */
    uint8_t sign = (signs[idx / 8] >> (idx % 8)) & 1;
    int16_t cand16 = (int16_t)candidate;
    *r = sign ? (int16_t)(-cand16) : cand16;

    return accepted;
}

/* ============================================================
 * sample_gauss_N: unified SHAKE stream, mini-batch processing.
 *
 * Stream byte layout:
 *   [0 .. len/8-1]:    sign bytes (one bit per accepted sample)
 *   [len/8 .. ]:       repeating mini-batches:
 *     [SIGMA2_RAND_BYTES bytes]    CDT randomness (AVX2-friendly layout)
 *     [GAUSS_BATCH * 9 bytes] y + rejection per attempt
 *
 * Buffer management:
 *   ~0.7 KB buffer (for BATCH=16), refilled with SQUEEZE_NBLOCKS blocks.
 * ============================================================ */
void sample_gauss_N(int16_t *r,
                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                    uint64_t nonce, size_t len) {
    uint8_t buf[GAUSS_BUF_SIZE];
    stream256_state state;
    stream256_init(&state, seed, nonce);

    /* Initial squeeze: enough for signs + first mini-batch */
    size_t sign_bytes = len / 8;
    size_t init_needed = sign_bytes + MINIBATCH_RAND_BYTES;
    size_t init_nblocks = (init_needed + STREAM256_BLOCKBYTES - 1)
                          / STREAM256_BLOCKBYTES;
    stream256_squeezeblocks(buf, init_nblocks, &state);

    /* Extract sign bits */
    uint8_t signs[SHUTTLE_N / 8];  /* max len = N = 256, so 32 bytes */
    memcpy(signs, buf, sign_bytes);
    size_t pos = sign_bytes;
    size_t avail = init_nblocks * STREAM256_BLOCKBYTES - sign_bytes;

    /* Mini-batch loop */
    size_t coefcnt = 0;
    int16_t z[GAUSS_BATCH];

    while (coefcnt < len) {
        /* Refill buffer if not enough for a full mini-batch */
        if (avail < (size_t)MINIBATCH_RAND_BYTES) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            stream256_squeezeblocks(buf + avail, SQUEEZE_NBLOCKS, &state);
            avail += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }

        /* CDT batch: consume SIGMA2_RAND_BYTES */
        sampler_sigma2(z, buf + pos);
        pos += SIGMA2_RAND_BYTES;
        avail -= SIGMA2_RAND_BYTES;

        /* Process GAUSS_BATCH attempts */
        for (int j = 0; j < GAUSS_BATCH && coefcnt < len; j++) {
            int accepted = sample_gauss_sigma128(&r[coefcnt],
                                                  buf + pos, z[j],
                                                  signs, coefcnt);
            pos += GAUSS_RAND_BYTES;
            avail -= GAUSS_RAND_BYTES;

            if (accepted) {
                coefcnt++;
            }
        }
    }
}
