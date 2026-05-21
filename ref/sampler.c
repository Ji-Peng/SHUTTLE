/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE (ref version).
 *
 * Four modes selectable via SHUTTLE_SIGMA (see config.h / params.h):
 *
 *   sigma | small  |  k  | RCDT bits |  entries  | parameter set
 *   ------+--------+-----+-----------+-----------+---------------
 *    128  | 2.0000 |  64 |   93      |    22     | legacy / reference
 *    101  | 1.5781 |  64 |   93      |    17     | SHUTTLE-128
 *    149  | 2.3281 |  64 |   93      |    26     | SHUTTLE-256
 *    202  | 1.5781 | 128 |   96      |    18     | SHUTTLE-512
 *
 * SHUTTLE-512 doubles k from 64 to 128 (=> K_BITS=7, Y_BITS=7, TWO_K_BITS=8)
 * and switches the RCDT layout from 3x31-bit limbs to 3x32-bit limbs to
 * carry 96 bits of precision (matching tools/BaseSampler.ipynb).
 *
 * Output: signed int16_t samples with |r| <= 11*sigma. Max:
 *   sigma=101 -> 1111, sigma=149 -> 1639, sigma=202 -> 2222, sigma=128 -> 1408.
 *
 * Unified buffer design: one SHAKE-256 stream per sample_gauss_N call
 * produces ALL randomness, in repeating mini-batches.
 */

#include "sampler.h"
#include "approx_exp.h"
#include "approx_exp_constants.h"
#include <string.h>

/* ============================================================
 * RCDT table.
 *
 * Layout per entry: 3 limbs (low, mid, upper). Limb width is 31 bits
 * for SIGMA in {128,101,149} (high bit masked off in the comparison)
 * and 32 bits for SIGMA=202 (full uint32 used).
 *
 * Entry i encodes Pr[|X_small_sigma| >= i+1] * 2^RCDT_BITS for the
 * mode's small_sigma. RCDT_BITS = 93 for limb-31 modes and 96 for the
 * limb-32 mode (SHUTTLE-512).
 * ============================================================ */

#if SHUTTLE_SIGMA == 128
/* small_sigma = 2.0000, 22 entries, Renyi (1025) = 1 + 2^{-93.85}. */
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
#elif SHUTTLE_SIGMA == 101
/* small_sigma = 1.5781, 17 entries, Renyi (1025) = 1 + 2^{-93.40}. */
static const uint32_t RCDT_3x31[RCDT_ENTRIES][3] = {
    {0x451749F4U, 0x18ACAE70U, 0x4C57A8B4U},
    {0x36957301U, 0x799737CEU, 0x2214D42AU},
    {0x78E3EE86U, 0x594E7B92U, 0x0AF10BD9U},
    {0x3F7F03BDU, 0x4C1BE93FU, 0x02763175U},
    {0x42E4A820U, 0x6029EBF2U, 0x0061BFAEU},
    {0x12E96D88U, 0x2C3F1462U, 0x000A584AU},
    {0x25E3DAE8U, 0x5DA4C386U, 0x0000BDF6U},
    {0x4FBF560BU, 0x16C867C8U, 0x00000932U},
    {0x6E585FDDU, 0x590E5C56U, 0x0000004CU},
    {0x7A9810AEU, 0x56D6F249U, 0x00000001U},
    {0x68B3CA36U, 0x0327875AU, 0x00000000U},
    {0x3BA50C27U, 0x0007F2CAU, 0x00000000U},
    {0x27184C49U, 0x00000D6BU, 0x00000000U},
    {0x1643B825U, 0x0000000FU, 0x00000000U},
    {0x05BEA231U, 0x00000000U, 0x00000000U},
    {0x0002E981U, 0x00000000U, 0x00000000U},
    {0x000000FCU, 0x00000000U, 0x00000000U}
};
#elif SHUTTLE_SIGMA == 149
/* small_sigma = 2.3281, 26 entries, Renyi (1025) = 1 + 2^{-93.75}. */
static const uint32_t RCDT_3x31[RCDT_ENTRIES][3] = {
    {0x27D8432FU, 0x5CC00AEDU, 0x5A8CA8FBU},
    {0x6B59D50CU, 0x52D4B77EU, 0x38662F3FU},
    {0x0589FE19U, 0x52DE1A50U, 0x1E814196U},
    {0x19228C60U, 0x3DBC9E12U, 0x0E2DC05DU},
    {0x395F2C65U, 0x64BD9992U, 0x059E9128U},
    {0x21DD30A0U, 0x5D4446D4U, 0x01E35796U},
    {0x53D9CE7AU, 0x2423C1B7U, 0x0089147FU},
    {0x5BA66039U, 0x729CDEA1U, 0x0020B5B4U},
    {0x4EF21306U, 0x73B4FA72U, 0x00068CFFU},
    {0x26A028CEU, 0x7A947178U, 0x00011957U},
    {0x5D3622B3U, 0x1C22D94BU, 0x0000277BU},
    {0x6FED0352U, 0x06B4278AU, 0x000004A1U},
    {0x75A0969DU, 0x7E04A1C3U, 0x00000073U},
    {0x3CCE8069U, 0x3C0903B5U, 0x00000009U},
    {0x713308E2U, 0x527DF437U, 0x00000000U},
    {0x232E8D11U, 0x04ADABA1U, 0x00000000U},
    {0x4F10AD3CU, 0x003893F8U, 0x00000000U},
    {0x60487ED3U, 0x000239C1U, 0x00000000U},
    {0x462800E2U, 0x000012A8U, 0x00000000U},
    {0x18E3CA18U, 0x00000082U, 0x00000000U},
    {0x7A0221CEU, 0x00000002U, 0x00000000U},
    {0x07226B9EU, 0x00000000U, 0x00000000U},
    {0x001CADF0U, 0x00000000U, 0x00000000U},
    {0x00005FE7U, 0x00000000U, 0x00000000U},
    {0x0000010AU, 0x00000000U, 0x00000000U},
    {0x00000002U, 0x00000000U, 0x00000000U}
};
#elif SHUTTLE_SIGMA == 202
/* small_sigma = 1.5781, k=128, 18 entries, 96-bit precision (3x32 limbs),
 * Renyi (2049) = 1 + 2^{-96.31}. From tools/BaseSampler.ipynb. */
static const uint32_t RCDT_3x32[RCDT_ENTRIES][3] = {
    {0x28BA4FD7U, 0x62B2B9C2U, 0x98AF5168U},
    {0xB4AB983BU, 0xE65CDF39U, 0x4429A855U},
    {0xC71F7463U, 0x6539EE4BU, 0x15E217B3U},
    {0xFBF81E18U, 0x306FA4FDU, 0x04EC62EBU},
    {0x1725412BU, 0x80A7AFCAU, 0x00C37F5DU},
    {0x974B6C68U, 0xB0FC5188U, 0x0014B094U},
    {0x2F1ED768U, 0x76930E19U, 0x00017BEDU},
    {0x7DFAB07CU, 0x5B219F22U, 0x00001264U},
    {0x72C2FF0CU, 0x6439715BU, 0x00000099U},
    {0xD4C08593U, 0x5B5BC927U, 0x00000003U},
    {0x459E51D0U, 0x0C9E1D6BU, 0x00000000U},
    {0xDD286156U, 0x001FCB29U, 0x00000000U},
    {0x38C26260U, 0x000035ADU, 0x00000000U},
    {0xB21DC13BU, 0x0000003CU, 0x00000000U},
    {0x2DF51197U, 0x00000000U, 0x00000000U},
    {0x00174C15U, 0x00000000U, 0x00000000U},
    {0x000007E7U, 0x00000000U, 0x00000000U},
    {0x00000001U, 0x00000000U, 0x00000000U}
};
#endif

/* ============================================================
 * sampler_sigma2: batched small-sigma CDT sampler.
 *
 * Two implementations behind an #if -- the 3x31 path masks each limb's
 * high bit so the comparison stays in a strict 31-bit range (matches
 * the AVX2 layout). The 3x32 path uses the full 32-bit limb (no mask).
 * ============================================================ */
int sampler_sigma2(int16_t *z_out, const uint8_t *rand) {
    for (int s = 0; s < GAUSS_BATCH; s++) {
        int group = s >> 3;
        int lane  = s & 7;
        int byte_base = group * 96 + lane * 4;

#if SHUTTLE_RCDT_LIMB_BITS == 31
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
#elif SHUTTLE_RCDT_LIMB_BITS == 32
        uint32_t v0 = load_le32(rand + byte_base);
        uint32_t v1 = load_le32(rand + byte_base + 32);
        uint32_t v2 = load_le32(rand + byte_base + 64);

        int32_t z = 0;
        for (int i = 0; i < RCDT_ENTRIES; i++) {
            uint32_t cc;
            cc = ct_lt_u32(v0, RCDT_3x32[i][0]);
            cc = ct_lt_u32(v1 - cc, RCDT_3x32[i][1]);
            cc = ct_lt_u32(v2 - cc, RCDT_3x32[i][2]);
            z += (int32_t)cc;
        }
        z_out[s] = (int16_t)z;
#else
#  error "Unsupported SHUTTLE_RCDT_LIMB_BITS"
#endif
    }
    return GAUSS_BATCH;
}

/* ============================================================
 * sampler_y: batched uniform Y_BITS-bit sampler.
 *
 * Produces GAUSS_BATCH = 16 uniform Y_BITS-bit values from exactly
 *   Y_RAND_BYTES = (GAUSS_BATCH * Y_BITS + 7) / 8 = ceil(16*Y_BITS/8)
 * input bytes, packed little-endian bit-wise. Implementation is a
 * generic accumulator loop -- portable for any Y_BITS up to ~24.
 * ============================================================ */
static void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand) {
    const uint32_t y_mask = (1U << SHUTTLE_Y_BITS) - 1U;
    uint32_t acc = 0;
    unsigned int acc_bits = 0;
    unsigned int byte_idx = 0;
    for (int i = 0; i < GAUSS_BATCH; ++i) {
        while (acc_bits < (unsigned)SHUTTLE_Y_BITS) {
            acc |= (uint32_t)rand[byte_idx++] << acc_bits;
            acc_bits += 8;
        }
        y_out[i] = (uint8_t)(acc & y_mask);
        acc >>= SHUTTLE_Y_BITS;
        acc_bits -= SHUTTLE_Y_BITS;
    }
}

/* ============================================================
 * sample_gauss_attempt: one (x, y, sign_for_0, rej) attempt.
 *
 *   candidate = (x << K_BITS) | y                   in [0, 2^K_BITS * (RCDT_ENTRIES+1)]
 *   t         = y + (x << TWO_K_BITS)               = y + 2*k*x
 *   num       = y * t                                in uint32 (num < 2^20)
 *   a         = num / (2*sigma^2)
 *   accept iff rand_rej_63 < approx_exp(a_q60)       in Q63
 *
 * The input a_q60 = round(a * 2^60) to approx_exp is computed in
 * fixed-point with a Q80 reciprocal of 2*sigma^2 (see agent/ApproxExp/
 * ApproxExp.tex Section 4.1):
 *
 *   I       = floor(2^80 / (2*sigma^2))                 ~ 64..66 bits
 *   M       = num * I    (full 128-bit integer multiply)
 *   a_q60   = M >> 20
 *
 * Error chain (relative to the ideal a):
 *   I underflow:           |2^80/N - I| < 1                   (Q80 ULP)
 *   M underflow:           |M_ideal - M| <= num <= 2^18       (Q80 ULP)
 *   shift truncation:      <= 1                               (Q60 ULP)
 *   total                  |a_q60 - a*2^60| <= 5/4 + 1 < 2   (Q60 ULP)
 *   propagation to exp(-a): <= 2 / 2^60 * exp(-a) = 2^-59 * exp(-a)
 *
 * sigma=128 (legacy) is a power-of-two special case: 2*sigma^2 = 2^15
 * so a_q60 = num << 45 is exact, no Q80 reciprocal needed.
 *
 * For SHUTTLE-512 (k=128, RCDT_ENTRIES=18), max num <= 127 * (127 +
 * 256*18) = 601,345 < 2^20; the 128-bit M fits comfortably. For the
 * other modes num is smaller still (mode-128: 141057, mode-256:
 * 213633).
 * ============================================================ */
static int sample_gauss_attempt(int16_t *r,
                                uint32_t x,
                                uint32_t y,
                                const uint8_t *rej_rand,
                                const uint8_t *signs,
                                size_t idx) {
    int32_t candidate = ((int32_t)x << SHUTTLE_K_BITS) | (int32_t)y;

    uint64_t rand_tail = load_le64(rej_rand);
    uint8_t  sign_r0      = (uint8_t)(rand_tail & 1U);
    uint64_t rand_rej_63  = rand_tail >> 1;

    uint32_t t   = (uint32_t)y + ((uint32_t)x << SHUTTLE_TWO_K_BITS);
    uint64_t num = (uint64_t)y * (uint64_t)t;

#if SHUTTLE_SIGMA == 128
    /* 2*sigma^2 = 2^15, so (num / 2^15) in Q60 = num << 45 (exact). */
    uint64_t a_q60 = num << 45;
#else
    /* Q80 reciprocal path. Splitting I into (I_hi, I_lo) lets us cover
     * the 65--66 bit reciprocal of sigma=101/149 with two 64-bit lanes,
     * and degenerates cleanly when I < 2^64 (sigma=202: I_hi == 0).
     * The product is built as
     *
     *   M = num * I = num * I_lo  +  (num * I_hi) << 64
     *
     * `num` is bounded by 2^20 and `I_hi` is bounded by 2^2 (for the
     * largest case sigma=101), so `num * I_hi` always fits in a uint64
     * -- no nested 128-bit multiplication is required. */
    unsigned __int128 M_lo = (unsigned __int128)num * APPROX_EXP_I_LO;
    uint64_t          M_hi = (uint64_t)num * APPROX_EXP_I_HI;
    unsigned __int128 M = M_lo + ((unsigned __int128)M_hi << 64);
    uint64_t a_q60 = (uint64_t)(M >> 20);
#endif

    uint64_t exp_val_q63 = approx_exp(a_q60);
    int accepted = (rand_rej_63 < exp_val_q63) ? 1 : 0;

    /* r=0 case: accept with prob 1/2 (folded-distribution correction). */
    uint64_t r_is_zero = (candidate == 0) ? 1ULL : 0ULL;
    accepted &= (int)(1 - (r_is_zero & (1 - (uint64_t)sign_r0)));

    uint8_t sign = (signs[idx / 8] >> (idx % 8)) & 1;
    int16_t cand16 = (int16_t)candidate;
    *r = sign ? (int16_t)(-cand16) : cand16;

    return accepted;
}

/* ============================================================
 * sample_gauss_N: unified SHAKE stream, mini-batch processing.
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

    /* len <= SHUTTLE_N (== 1024 worst case for mode-512) so sign_bytes <= 128. */
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
