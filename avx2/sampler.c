/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE (AVX2 version).
 *
 * sampler_sigma2: AVX2 8-wide CDT comparison (2 groups x 8 = 16 samples).
 * All other functions identical to ref. KAT-compatible: same random byte
 * consumption order (signs -> CDT -> y+rej per mini-batch).
 */

#include "sampler.h"
#include "approx_exp.h"
#include <string.h>
#include <immintrin.h>

/* Compile-time check: AVX2 sampler_sigma2 is tuned for BATCH=16 (2x8) */
typedef char static_assert_batch16[(NGCC_GAUSS_BATCH == 16) ? 1 : -1];

/* ============================================================
 * RCDT Table: 8-lane AVX2 broadcast format.
 *
 * Each entry has 3 limbs (low, mid, high), each broadcast to 8 uint32_t
 * lanes for direct _mm256_load_si256.
 * Values identical to ref RCDT_3x31.
 * ============================================================ */
#define U32X8(x) {(x), (x), (x), (x), (x), (x), (x), (x)}

static const uint32_t RCDT_AVX2[RCDT_ENTRIES][3][8]
    __attribute__((aligned(32))) = {
    { U32X8(0x0204FE5FU), U32X8(0x4E73D911U), U32X8(0x556D69B6U) },
    { U32X8(0x7624991FU), U32X8(0x1BCB49D6U), U32X8(0x2FDB7191U) },
    { U32X8(0x26A118E0U), U32X8(0x7DB923CBU), U32X8(0x16091DCDU) },
    { U32X8(0x70C18AE5U), U32X8(0x536A5628U), U32X8(0x0836DCFDU) },
    { U32X8(0x54119E8DU), U32X8(0x3D589AD0U), U32X8(0x0273E65EU) },
    { U32X8(0x2EB0249EU), U32X8(0x2EBE8BE2U), U32X8(0x00950CA4U) },
    { U32X8(0x0472F0ECU), U32X8(0x65A6660AU), U32X8(0x001BFA1DU) },
    { U32X8(0x70C91274U), U32X8(0x5945CB2DU), U32X8(0x000422EEU) },
    { U32X8(0x0517A411U), U32X8(0x4A81EF20U), U32X8(0x00007AFAU) },
    { U32X8(0x3A4CC14BU), U32X8(0x79C6FA82U), U32X8(0x00000B31U) },
    { U32X8(0x12B0F081U), U32X8(0x185517F0U), U32X8(0x000000CCU) },
    { U32X8(0x6073B7BEU), U32X8(0x2FB676CEU), U32X8(0x0000000BU) },
    { U32X8(0x1F24D098U), U32X8(0x3F53D508U), U32X8(0x00000000U) },
    { U32X8(0x5FBBC714U), U32X8(0x02267C01U), U32X8(0x00000000U) },
    { U32X8(0x0A24AB07U), U32X8(0x000E9522U), U32X8(0x00000000U) },
    { U32X8(0x61E4FCE8U), U32X8(0x00004D1EU), U32X8(0x00000000U) },
    { U32X8(0x7B65A4CBU), U32X8(0x0000013DU), U32X8(0x00000000U) },
    { U32X8(0x7EE574B9U), U32X8(0x00000003U), U32X8(0x00000000U) },
    { U32X8(0x04FF6C73U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x0009C080U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x00000ED2U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x00000011U), U32X8(0x00000000U), U32X8(0x00000000U) }
};

/* ============================================================
 * AVX2 batched CDT sampler: sampler_sigma2
 *
 * 16 samples in 2 groups of 8 using AVX2 8-wide SIMD.
 * Internal 32-bit computation, output packed to int16_t.
 * ============================================================ */
int sampler_sigma2(int16_t *z_out, const uint8_t *rand) {
    const __m256i mask31 = _mm256_set1_epi32(0x7FFFFFFF);

    /* Load 192 bytes into 6 AVX2 registers (2 groups x 3 limbs) */
    __m256i rnd_low0 = _mm256_and_si256(
        _mm256_loadu_si256((const __m256i *)(rand +   0)), mask31);
    __m256i rnd_mid0 = _mm256_and_si256(
        _mm256_loadu_si256((const __m256i *)(rand +  32)), mask31);
    __m256i rnd_hig0 = _mm256_and_si256(
        _mm256_loadu_si256((const __m256i *)(rand +  64)), mask31);
    __m256i rnd_low1 = _mm256_and_si256(
        _mm256_loadu_si256((const __m256i *)(rand +  96)), mask31);
    __m256i rnd_mid1 = _mm256_and_si256(
        _mm256_loadu_si256((const __m256i *)(rand + 128)), mask31);
    __m256i rnd_hig1 = _mm256_and_si256(
        _mm256_loadu_si256((const __m256i *)(rand + 160)), mask31);

    __m256i z0 = _mm256_setzero_si256();
    __m256i z1 = _mm256_setzero_si256();

    for (int i = 0; i < RCDT_ENTRIES; i++) {
        __m256i tlo = _mm256_load_si256(
            (const __m256i *)RCDT_AVX2[i][0]);
        __m256i tmi = _mm256_load_si256(
            (const __m256i *)RCDT_AVX2[i][1]);
        __m256i thi = _mm256_load_si256(
            (const __m256i *)RCDT_AVX2[i][2]);

        __m256i cc0, cc1;

        /* Step 1: Low limb */
        cc0 = _mm256_cmpgt_epi32(tlo, rnd_low0);
        cc1 = _mm256_cmpgt_epi32(tlo, rnd_low1);

        /* Step 2: Mid limb with borrow */
        cc0 = _mm256_add_epi32(cc0, rnd_mid0);
        cc1 = _mm256_add_epi32(cc1, rnd_mid1);
        cc0 = _mm256_sub_epi32(cc0, tmi);
        cc1 = _mm256_sub_epi32(cc1, tmi);
        cc0 = _mm256_srli_epi32(cc0, 31);
        cc1 = _mm256_srli_epi32(cc1, 31);

        /* Step 3: High limb with borrow */
        cc0 = _mm256_sub_epi32(rnd_hig0, cc0);
        cc1 = _mm256_sub_epi32(rnd_hig1, cc1);
        cc0 = _mm256_sub_epi32(cc0, thi);
        cc1 = _mm256_sub_epi32(cc1, thi);
        cc0 = _mm256_srli_epi32(cc0, 31);
        cc1 = _mm256_srli_epi32(cc1, 31);

        z0 = _mm256_add_epi32(z0, cc0);
        z1 = _mm256_add_epi32(z1, cc1);
    }

    /* Pack 32-bit results to 16-bit.
     * _mm256_packs_epi32 packs with saturation (values [0,22] are fine).
     * Lane layout after pack: [z0_lo4, z1_lo4, z0_hi4, z1_hi4]
     * Permute to get [z0_0-7, z1_0-7] order. */
    __m256i packed = _mm256_packs_epi32(z0, z1);
    packed = _mm256_permute4x64_epi64(packed, 0xD8);  /* [0,2,1,3] */
    _mm256_storeu_si256((__m256i *)z_out, packed);

    return NGCC_GAUSS_BATCH;
}

/* ============================================================
 * Single Gaussian attempt (identical to ref)
 * ============================================================ */
static int sample_gauss_sigma128(int16_t *r, const uint8_t *rand,
                                 int16_t x, const uint8_t *signs,
                                 size_t idx) {
    uint32_t y = rand[0] & 0x3FU;
    int32_t candidate = (int32_t)x * SHUTTLE_GAUSS_K + (int32_t)y;

    uint64_t rand_tail = load_le64(rand + 1);
    uint8_t sign_r0 = rand_tail & 1;
    uint64_t rand_rej_63 = rand_tail >> 1;

    uint32_t t = y + (uint32_t)(2 * SHUTTLE_GAUSS_K) * (uint32_t)x;
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

/* ============================================================
 * sample_gauss_N (identical to ref)
 * ============================================================ */
void sample_gauss_N(int16_t *r,
                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                    uint64_t nonce, size_t len) {
    uint8_t buf[GAUSS_BUF_SIZE];
    stream256_state state;
    stream256_init(&state, seed, nonce);

    size_t sign_bytes = len / 8;
    size_t init_needed = sign_bytes + MINIBATCH_RAND_BYTES;
    size_t init_nblocks = (init_needed + STREAM256_BLOCKBYTES - 1)
                          / STREAM256_BLOCKBYTES;
    stream256_squeezeblocks(buf, init_nblocks, &state);

    uint8_t signs[SHUTTLE_N / 8];
    memcpy(signs, buf, sign_bytes);
    size_t pos = sign_bytes;
    size_t avail = init_nblocks * STREAM256_BLOCKBYTES - sign_bytes;

    size_t coefcnt = 0;
    int16_t z[NGCC_GAUSS_BATCH];

    while (coefcnt < len) {
        if (avail < (size_t)MINIBATCH_RAND_BYTES) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            stream256_squeezeblocks(buf + avail, SQUEEZE_NBLOCKS, &state);
            avail += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }

        sampler_sigma2(z, buf + pos);
        pos += SIGMA2_RAND_BYTES;
        avail -= SIGMA2_RAND_BYTES;

        for (int j = 0; j < NGCC_GAUSS_BATCH && coefcnt < len; j++) {
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
