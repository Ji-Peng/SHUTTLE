/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE (AVX2 version).
 *
 * Three modes selectable via SHUTTLE_SIGMA (see config.h / params.h):
 *
 *   sigma | small_sigma | RCDT entries | 2*sigma^2 | parameter set
 *   ------+-------------+--------------+-----------+---------------
 *    128  |   2.0000    |      22      |  32768    | legacy / reference
 *    101  |   1.5781    |      17      |  20402    | SHUTTLE-128  (default)
 *    149  |   2.3281    |      26      |  44402    | SHUTTLE-256
 *
 * sampler_sigma2: AVX2 8-wide CDT comparison (2 groups x 8 = 16 samples).
 * sampler_y:      scalar - see analysis comment above the function for
 *                 why an AVX2 variant is NOT used.
 * sample_gauss_attempt / sample_gauss_N: arithmetic identical to ref; the
 *                 same Q60 exponent is computed via a shift fast path for
 *                 sigma=128 and via mulh64 + precomputed Q64 reciprocal for
 *                 the other two modes.
 *
 * KAT-compatible with the ref implementation: same SHAKE consumption order
 * (signs -> CDT -> y-batch -> per-attempt rejection tails) and same math.
 */

#include "sampler.h"
#include "approx_exp.h"
#include <string.h>
#include <immintrin.h>

/* Compile-time check: AVX2 sampler_sigma2 is tuned for BATCH=16 (2x8) */
typedef char static_assert_batch16[(GAUSS_BATCH == 16) ? 1 : -1];

/* ============================================================
 * RCDT table: 8-lane AVX2 broadcast format.
 *
 * Each entry has 3 limbs (low, mid, high), each broadcast to 8 uint32_t
 * lanes for direct _mm256_load_si256. Values identical to ref RCDT_3x31.
 * ============================================================ */
#define U32X8(x) {(x), (x), (x), (x), (x), (x), (x), (x)}

#if SHUTTLE_SIGMA == 128
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
#elif SHUTTLE_SIGMA == 101
static const uint32_t RCDT_AVX2[RCDT_ENTRIES][3][8]
    __attribute__((aligned(32))) = {
    { U32X8(0x451749F4U), U32X8(0x18ACAE70U), U32X8(0x4C57A8B4U) },
    { U32X8(0x36957301U), U32X8(0x799737CEU), U32X8(0x2214D42AU) },
    { U32X8(0x78E3EE86U), U32X8(0x594E7B92U), U32X8(0x0AF10BD9U) },
    { U32X8(0x3F7F03BDU), U32X8(0x4C1BE93FU), U32X8(0x02763175U) },
    { U32X8(0x42E4A820U), U32X8(0x6029EBF2U), U32X8(0x0061BFAEU) },
    { U32X8(0x12E96D88U), U32X8(0x2C3F1462U), U32X8(0x000A584AU) },
    { U32X8(0x25E3DAE8U), U32X8(0x5DA4C386U), U32X8(0x0000BDF6U) },
    { U32X8(0x4FBF560BU), U32X8(0x16C867C8U), U32X8(0x00000932U) },
    { U32X8(0x6E585FDDU), U32X8(0x590E5C56U), U32X8(0x0000004CU) },
    { U32X8(0x7A9810AEU), U32X8(0x56D6F249U), U32X8(0x00000001U) },
    { U32X8(0x68B3CA36U), U32X8(0x0327875AU), U32X8(0x00000000U) },
    { U32X8(0x3BA50C27U), U32X8(0x0007F2CAU), U32X8(0x00000000U) },
    { U32X8(0x27184C49U), U32X8(0x00000D6BU), U32X8(0x00000000U) },
    { U32X8(0x1643B825U), U32X8(0x0000000FU), U32X8(0x00000000U) },
    { U32X8(0x05BEA231U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x0002E981U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x000000FCU), U32X8(0x00000000U), U32X8(0x00000000U) }
};
#elif SHUTTLE_SIGMA == 149
static const uint32_t RCDT_AVX2[RCDT_ENTRIES][3][8]
    __attribute__((aligned(32))) = {
    { U32X8(0x27D8432FU), U32X8(0x5CC00AEDU), U32X8(0x5A8CA8FBU) },
    { U32X8(0x6B59D50CU), U32X8(0x52D4B77EU), U32X8(0x38662F3FU) },
    { U32X8(0x0589FE19U), U32X8(0x52DE1A50U), U32X8(0x1E814196U) },
    { U32X8(0x19228C60U), U32X8(0x3DBC9E12U), U32X8(0x0E2DC05DU) },
    { U32X8(0x395F2C65U), U32X8(0x64BD9992U), U32X8(0x059E9128U) },
    { U32X8(0x21DD30A0U), U32X8(0x5D4446D4U), U32X8(0x01E35796U) },
    { U32X8(0x53D9CE7AU), U32X8(0x2423C1B7U), U32X8(0x0089147FU) },
    { U32X8(0x5BA66039U), U32X8(0x729CDEA1U), U32X8(0x0020B5B4U) },
    { U32X8(0x4EF21306U), U32X8(0x73B4FA72U), U32X8(0x00068CFFU) },
    { U32X8(0x26A028CEU), U32X8(0x7A947178U), U32X8(0x00011957U) },
    { U32X8(0x5D3622B3U), U32X8(0x1C22D94BU), U32X8(0x0000277BU) },
    { U32X8(0x6FED0352U), U32X8(0x06B4278AU), U32X8(0x000004A1U) },
    { U32X8(0x75A0969DU), U32X8(0x7E04A1C3U), U32X8(0x00000073U) },
    { U32X8(0x3CCE8069U), U32X8(0x3C0903B5U), U32X8(0x00000009U) },
    { U32X8(0x713308E2U), U32X8(0x527DF437U), U32X8(0x00000000U) },
    { U32X8(0x232E8D11U), U32X8(0x04ADABA1U), U32X8(0x00000000U) },
    { U32X8(0x4F10AD3CU), U32X8(0x003893F8U), U32X8(0x00000000U) },
    { U32X8(0x60487ED3U), U32X8(0x000239C1U), U32X8(0x00000000U) },
    { U32X8(0x462800E2U), U32X8(0x000012A8U), U32X8(0x00000000U) },
    { U32X8(0x18E3CA18U), U32X8(0x00000082U), U32X8(0x00000000U) },
    { U32X8(0x7A0221CEU), U32X8(0x00000002U), U32X8(0x00000000U) },
    { U32X8(0x07226B9EU), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x001CADF0U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x00005FE7U), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x0000010AU), U32X8(0x00000000U), U32X8(0x00000000U) },
    { U32X8(0x00000002U), U32X8(0x00000000U), U32X8(0x00000000U) }
};
#endif

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
     * _mm256_packs_epi32 packs with saturation (values [0, RCDT_ENTRIES] are fine).
     * Lane layout after pack: [z0_lo4, z1_lo4, z0_hi4, z1_hi4]
     * Permute to get [z0_0-7, z1_0-7] order. */
    __m256i packed = _mm256_packs_epi32(z0, z1);
    packed = _mm256_permute4x64_epi64(packed, 0xD8);  /* [0,2,1,3] */
    _mm256_storeu_si256((__m256i *)z_out, packed);

    return GAUSS_BATCH;
}

/* ============================================================
 * Batched uniform y sampler: sampler_y
 *
 * AVX2 NOTES (why scalar here):
 *   This function extracts 16 six-bit fields from 12 bytes = 96 bits.
 *   At 16 fields, the total scalar work is ~16 shifts + 16 masks + 2 loads,
 *   well under 20 cycles on modern cores. An AVX2 version would need:
 *     1) broadcasting / unpacking the 12 bytes into a 16-lane vector,
 *     2) a per-lane variable shift (_mm256_srlv_epi32 only comes in 32-bit
 *        lanes, so we'd need 2 ymm regs and recombination),
 *     3) a permute/blend to move the boundary-straddling field into place.
 *   The setup cost alone (~10-15 cy) makes SIMD at best break-even, while
 *   adding code complexity and a second KAT-critical code path. In the full
 *   sample_gauss_N loop this step is amortized over 16 approx_exp calls
 *   (~30 cy each) plus sampler_sigma2 (53 cy AVX2), so its contribution
 *   to total cycles is <1% either way.
 *
 * Byte-identical to the ref implementation - required for KAT parity.
 * ============================================================ */
static void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand) {
    const uint32_t y_mask = (1U << SHUTTLE_Y_BITS) - 1U;
    uint64_t lo = load_le64(rand);        /* bits [0, 64)  */
    uint32_t hi = load_le32(rand + 8);    /* bits [64, 96) */

    for (int i = 0; i < 10; i++)
        y_out[i] = (uint8_t)((lo >> (SHUTTLE_Y_BITS * i)) & y_mask);

    /* y[10] straddles the 64-bit boundary: low 4 bits from lo[60..63],
     * top 2 bits from hi[0..1]. */
    y_out[10] = (uint8_t)(((lo >> 60) | ((uint64_t)hi << 4)) & y_mask);

    for (int i = 11; i < 16; i++) {
        unsigned shift = (unsigned)(SHUTTLE_Y_BITS * i - 64);
        y_out[i] = (uint8_t)((hi >> shift) & y_mask);
    }
}

/* ============================================================
 * Single Gaussian attempt (identical math to ref).
 *
 * See the ref implementation for the full error analysis of the
 * non-shift path (sigma=101, 149). This file reproduces the same
 * comment block because the branch is security-critical and must be
 * auditable in both implementations.
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
    /* 2*sigma^2 = 2^15, so (num / 2^15) in Q60 = num << 45 (exact). */
    uint64_t exp_q60 = num << 45;
#else
    /* --------------------------------------------------------
     * Non-power-of-two divisor: use a Q64 reciprocal + mulh64.
     *
     * Goal: exp_q60 = round(num * 2^60 / (2*sigma^2)),   0 <= num <= 2^18.
     *
     * Let N = 2*sigma^2               (20402 for sigma=101;
     *                                  44402 for sigma=149)
     *     I = SHUTTLE_INV_2SIGMA2_Q64 = floor((2^64 - 1) / N)
     *
     * Because N is NOT a power of two, I equals floor(2^64 / N) exactly,
     * so 1/N - I / 2^64 = eps_I with 0 <= eps_I < 1 / 2^64.
     *
     * Let A = num << 44  (so 0 <= A <= 2^62, fits in uint64), and
     *     M = mulh64(A, I) = floor(A * I / 2^64).
     *
     * Error vs the ideal ratio A / N:
     *     A * I / 2^64     = A / N - A * eps_I
     *     |A * eps_I|      <= 2^62 * 2^-64 = 1/4
     *     |floor truncation| <= 1
     *   => |M - A/N| <= 1 + 1/4 < 2           (units of 1 in Q44)
     *
     * After the final left-shift by 16:
     *   exp_q60 = M << 16
     *   |exp_q60 - num * 2^60 / N| <= 2 * 2^16 = 2^17   (units in Q60)
     *
     * Relative error: exp_q60 reaches up to ~2^63, so 2^17 / 2^63 = 2^-46
     * worst case. approx_exp() has ~2^-56 intrinsic precision, so our
     * input noise does not degrade its output. For rejection sampling the
     * uniform 63-bit threshold masks errors below ~2^-40, so the
     * distribution is statistically indistinguishable from ideal for any
     * feasible sample count.
     * -------------------------------------------------------- */
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

/* ============================================================
 * sample_gauss_N: KAT-identical to ref; see ref for full commentary
 * on the byte schedule and buffer management.
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
    size_t pos   = sign_bytes;
    size_t avail = init_nblocks * STREAM256_BLOCKBYTES - sign_bytes;

    size_t coefcnt = 0;
    int16_t z[GAUSS_BATCH];
    uint8_t y[GAUSS_BATCH];

    while (coefcnt < len) {
        if (avail < (size_t)MINIBATCH_RAND_BYTES) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            stream256_squeezeblocks(buf + avail, SQUEEZE_NBLOCKS, &state);
            avail += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }

        sampler_sigma2(z, buf + pos);
        pos   += SIGMA2_RAND_BYTES;
        avail -= SIGMA2_RAND_BYTES;

        sampler_y(y, buf + pos);
        pos   += Y_RAND_BYTES;
        avail -= Y_RAND_BYTES;

        for (int j = 0; j < GAUSS_BATCH && coefcnt < len; j++) {
            int accepted = sample_gauss_attempt(&r[coefcnt],
                                                (uint32_t)z[j],
                                                (uint32_t)y[j],
                                                buf + pos,
                                                signs, coefcnt);
            pos   += GAUSS_RAND_BYTES;
            avail -= GAUSS_RAND_BYTES;

            if (accepted)
                coefcnt++;
        }
    }
}
