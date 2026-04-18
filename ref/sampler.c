/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE (ref version).
 *
 * Three modes selectable via SHUTTLE_SIGMA (see config.h / params.h):
 *
 *   sigma | small_sigma | RCDT entries | 2*sigma^2 | parameter set
 *   ------+-------------+--------------+-----------+---------------
 *    128  |   2.0000    |      22      |  32768    | legacy / reference
 *    101  |   1.5781    |      17      |  20402    | SHUTTLE-128  (default)
 *    149  |   2.3281    |      26      |  44402    | SHUTTLE-256
 *
 * Output: signed int16_t samples with |r| <= 11*sigma.
 *
 * Unified buffer design: one SHAKE-256 stream per sample_gauss_N call
 * produces ALL randomness (signs, CDT, y, rejection). Processing in
 * mini-batches of GAUSS_BATCH attempts.
 *
 * Stream byte order (per mini-batch):
 *   SIGMA2_RAND_BYTES (192) : CDT randomness (3x31-bit AVX2 layout)
 *   Y_RAND_BYTES      ( 12) : 16 six-bit y values packed into 96 bits
 *   GAUSS_RAND_BYTES*16=128 : per-attempt sign-for-zero + 63-bit rejection
 *
 * KAT-identical with the AVX2 implementation: same SHAKE consumption order
 * and same arithmetic path (shift fast-path for sigma=128; precomputed
 * Q64 reciprocal via mulh64 for the other two modes).
 */

#include "sampler.h"
#include "approx_exp.h"
#include <string.h>

/* ============================================================
 * RCDT table: 93-bit cumulative probabilities in 3x31-bit limb format.
 *
 * Entry i encodes Pr[|X_small_sigma| >= i+1] * 2^93, for i = 0 .. RCDT_ENTRIES-1.
 * Layout per entry: [0]=low 31 bits, [1]=mid 31 bits, [2]=upper 31 bits.
 * ============================================================ */

#if SHUTTLE_SIGMA == 128
/* small_sigma = 2.0000, 22 entries, Renyi (1025) = 1 + 2^{-93.85} */
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
/* small_sigma = 1.5781, 17 entries (truncated at 11*small_sigma),
 * Renyi (1025) = 1 + 2^{-93.40} */
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
/* small_sigma = 2.3281, 26 entries (truncated at 11*small_sigma),
 * Renyi (1025) = 1 + 2^{-93.75} */
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
#endif

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
 * Output range: [0, RCDT_ENTRIES], mode-dependent (see table above).
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
 * Batched uniform y sampler: sampler_y (pure function, static)
 *
 * Produces GAUSS_BATCH = 16 uniform six-bit values from exactly
 * Y_RAND_BYTES = 12 input bytes (= 16 * 6 bits / 8). Layout is little-endian
 * bit-packed: y[i] occupies bits [6*i, 6*i+6) of the 96-bit stream.
 *
 * Motivation: the previous design took 1 byte per y and discarded the
 * top 2 bits, wasting 2 out of every 8 bits. This function recovers that
 * waste while preserving a fixed, deterministic byte schedule for KAT.
 *
 * Implementation: two aligned loads (uint64 + uint32) cover all 96 bits;
 * each y is an (x >> shift) & mask of a 64-bit window.
 * ============================================================ */
static void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand) {
    const uint32_t y_mask = (1U << SHUTTLE_Y_BITS) - 1U;
    uint64_t lo = load_le64(rand);        /* bits [0, 64)  */
    uint32_t hi = load_le32(rand + 8);    /* bits [64, 96) */

    /* y[0..9]: six-bit fields fully inside the low 64 bits
     * (10 fields span 60 bits, within 64). */
    for (int i = 0; i < 10; i++)
        y_out[i] = (uint8_t)((lo >> (SHUTTLE_Y_BITS * i)) & y_mask);

    /* y[10] straddles the 64-bit boundary: bits [60, 66).
     *   low 4 bits come from lo[60..63], top 2 bits from hi[0..1]. */
    y_out[10] = (uint8_t)(((lo >> 60) | ((uint64_t)hi << 4)) & y_mask);

    /* y[11..15]: entirely within the high 32 bits, at offsets 2, 8, 14, 20, 26. */
    for (int i = 11; i < 16; i++) {
        unsigned shift = (unsigned)(SHUTTLE_Y_BITS * i - 64);
        y_out[i] = (uint8_t)((hi >> shift) & y_mask);
    }
}

/* ============================================================
 * Single Gaussian attempt: sample_gauss_attempt (static)
 *
 * Inputs:
 *   x         - CDT sample from sampler_sigma2, in [0, RCDT_ENTRIES].
 *   y         - 6-bit uniform value from sampler_y, in [0, 2^SHUTTLE_Y_BITS).
 *   rej_rand  - 8 bytes: low bit is the "sign when candidate==0" coin,
 *               upper 63 bits are the uniform rejection threshold.
 *   signs     - pre-extracted sign-bit buffer (one bit per accepted sample).
 *   idx       - index into signs[] to pull this attempt's sign.
 *
 * Output: *r receives the signed candidate (only meaningful if accepted).
 * Returns 1 if the attempt is accepted, 0 otherwise.
 *
 * Candidate & exponent are computed with shifts (no multiply by k):
 *   candidate = (x << K_BITS) | y
 *   t         = y + (x << TWO_K_BITS)    = y + 2*k*x
 *   num       = y * t
 *
 * Rejection probability = exp(-num / (2*sigma^2)).
 * ============================================================ */
static int sample_gauss_attempt(int16_t *r,
                                uint32_t x,
                                uint32_t y,
                                const uint8_t *rej_rand,
                                const uint8_t *signs,
                                size_t idx) {
    /* Candidate magnitude via pure shifts. y < 2^K_BITS is guaranteed by
     * sampler_y, so the OR is equivalent to an add. */
    int32_t candidate = ((int32_t)x << SHUTTLE_K_BITS) | (int32_t)y;

    /* Rejection-randomness split: 1 bit for the (r == 0) tie-breaker,
     * 63 bits compared against the accept probability in Q63. */
    uint64_t rand_tail = load_le64(rej_rand);
    uint8_t  sign_r0      = (uint8_t)(rand_tail & 1U);
    uint64_t rand_rej_63  = rand_tail >> 1;

    /* num = y * (y + 2*k*x). Max values per mode:
     *   sigma=128: 63*(63 + 128*22) = 181377   (~2^17.47)
     *   sigma=101: 63*(63 + 128*17) = 141057   (~2^17.11)
     *   sigma=149: 63*(63 + 128*26) = 213633   (~2^17.71)
     * All fit comfortably into uint32, and num into uint64. */
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

    /* r = 0 case: accept with prob 1/2 (folded-distribution correction).
     * When candidate = 0, exp(-0) = 1 and the rejection test always passes;
     * we must drop half of those to restore the correct two-sided density. */
    uint64_t r_is_zero = (candidate == 0) ? 1ULL : 0ULL;
    accepted &= (int)(1 - (r_is_zero & (1 - (uint64_t)sign_r0)));

    uint8_t sign = (signs[idx / 8] >> (idx % 8)) & 1;
    int16_t cand16 = (int16_t)candidate;
    *r = sign ? (int16_t)(-cand16) : cand16;

    return accepted;
}

/* ============================================================
 * sample_gauss_N: unified SHAKE stream, mini-batch processing.
 *
 * Stream byte layout (total per mini-batch = MINIBATCH_RAND_BYTES = 332):
 *   [0 .. len/8 - 1]   : sign bytes (one bit per accepted sample)
 *   [len/8 .. ]        : repeating mini-batches:
 *     SIGMA2_RAND_BYTES bytes  - CDT randomness   (x batch)
 *     Y_RAND_BYTES      bytes  - batched 6-bit y values (y batch)
 *     GAUSS_BATCH * GAUSS_RAND_BYTES bytes - per-attempt tails (sign0 + rej)
 *
 * Buffer: ~0.7 KB on-stack; refilled with SQUEEZE_NBLOCKS blocks at a time.
 * ============================================================ */
void sample_gauss_N(int16_t *r,
                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                    uint64_t nonce, size_t len) {
    uint8_t buf[GAUSS_BUF_SIZE];
    stream256_state state;
    stream256_init(&state, seed, nonce);

    /* Initial squeeze: enough for signs + first mini-batch. */
    size_t sign_bytes = len / 8;
    size_t init_needed = sign_bytes + MINIBATCH_RAND_BYTES;
    size_t init_nblocks = (init_needed + STREAM256_BLOCKBYTES - 1)
                          / STREAM256_BLOCKBYTES;
    stream256_squeezeblocks(buf, init_nblocks, &state);

    /* Extract sign bits. */
    uint8_t signs[SHUTTLE_N / 8];   /* max len = N = 256, so 32 bytes */
    memcpy(signs, buf, sign_bytes);
    size_t pos   = sign_bytes;
    size_t avail = init_nblocks * STREAM256_BLOCKBYTES - sign_bytes;

    size_t coefcnt = 0;
    int16_t z[GAUSS_BATCH];
    uint8_t y[GAUSS_BATCH];

    while (coefcnt < len) {
        /* Refill buffer if not enough for a full mini-batch. */
        if (avail < (size_t)MINIBATCH_RAND_BYTES) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            stream256_squeezeblocks(buf + avail, SQUEEZE_NBLOCKS, &state);
            avail += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }

        /* x batch: consume SIGMA2_RAND_BYTES. */
        sampler_sigma2(z, buf + pos);
        pos   += SIGMA2_RAND_BYTES;
        avail -= SIGMA2_RAND_BYTES;

        /* y batch: consume Y_RAND_BYTES. */
        sampler_y(y, buf + pos);
        pos   += Y_RAND_BYTES;
        avail -= Y_RAND_BYTES;

        /* GAUSS_BATCH attempts, each consumes GAUSS_RAND_BYTES. */
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
