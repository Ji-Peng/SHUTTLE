/*
 * polyvec.c -- vector-level XOF-driven sampling surface (AVX2 FORK).  P07
 * / M9.
 *
 * =========================================================================
 *  What this is
 * =========================================================================
 * This is the AVX2 fork of ref/polyvec.c.  It is BYTE-EXACT to the scalar
 * reference: every Expand / Sample routine consumes the IDENTICAL PRNG
 * stream (same nonces, same order, same amount) and produces the IDENTICAL
 * sampled coefficients, so make check-kat reproduces the EXACT ref
 * NGCC/SHA3 hashes. Only the INSTRUCTION SELECTION of two hot inner loops
 * changes:
 *
 *   1. ExpandA uniform-[0,q) rejection scan (uniform_reject_chunk):
 *      vectorized to a 16-wide AVX2 reject + BMI2 (pdep/pext) compaction,
 *      consuming the SAME contiguous 2-byte candidates in the SAME
 * ascending order with the SAME mask + unsigned-q test as the scalar
 *      us_next_candidate -> dst stream.  (ExpandA is ~25% of Sign and
 *      ~83% of Verify -- the headline lever.)
 *
 *   2. ExpandS / SampleY BaseSampler magnitude scan goes through the AVX2
 *      cdt_scan96 (the P05 sampler.c fork), which is bit-identical to the
 *      scalar borrow chain.  The wide-Gaussian finalize, the
 * zero-fold/sign logic and the ApproxExp Bernoulli
 * (approx_exp_accept_q64_x4) are kept scalar (cheap per candidate, already
 * bit-exact, low register pressure)
 *      -- they are called from this file gauss_stream_chunk /
 * noise_minibatch exactly as the scalar reference does.
 *
 * Everything else (ExpandSeeds/ExpandSigningSeeds/SampleC, the
 * gauss_stream orchestration, the per-attempt byte schedule, the
 * OUTPUT-indexed sign stream) is copied VERBATIM from ref/polyvec.c.  When
 * -DUSE_AVX2_SAMPLER is NOT defined the file is bit-for-bit the scalar
 * reference.
 *
 * Read polyvec.h FIRST for the binding NGCC-DRBG no-rate-cursor rule and
 * the 16-stream nonce layout.  The PINNED PRNG byte schedules are
 * documented in ref/polyvec.c and are UNCHANGED here.
 */
#include "polyvec.h"

#include <string.h>

#include "approx_exp.h" /* approx_exp_accept_q64_x4, approx_exp_accept_q64 */
#include "rcdt_tables.h" /* SHUTTLE_RCDT_Z, SHUTTLE_RCDT_NOISE_* */

#if defined(USE_AVX2_SAMPLER) && defined(__AVX2__)
#    include <immintrin.h>
#endif

/* ===================================================================== *
 *  Local little-endian byte helpers (data-independent schedule)         *
 * ===================================================================== */
static void put_le16(uint8_t *p, uint16_t x)
{
    p[0] = (uint8_t)(x & 0xFFu);
    p[1] = (uint8_t)((x >> 8) & 0xFFu);
}
static void put_le32(uint8_t *p, uint32_t x)
{
    p[0] = (uint8_t)(x & 0xFFu);
    p[1] = (uint8_t)((x >> 8) & 0xFFu);
    p[2] = (uint8_t)((x >> 16) & 0xFFu);
    p[3] = (uint8_t)((x >> 24) & 0xFFu);
}
static uint32_t get_le_masked(const uint8_t *p, unsigned nbytes,
                              uint32_t mask)
{
    uint32_t x = 0;
    unsigned i;
    for (i = 0; i < nbytes; i++)
        x |= (uint32_t)p[i] << (8 * i);
    return x & mask;
}

/* ===================================================================== *
 *  ExpandSeeds (DS 0x00, xof256, single one-shot)                       *
 * ===================================================================== */
void expand_seeds(uint8_t T[EXPAND_SEEDS_BYTES],
                  const uint8_t xi[SEEDBYTES], uint32_t kappa)
{
    /* absorb 0x00 || xi || IntToBytes(EM,2) || IntToBytes(kappa,4) */
    uint8_t in[1 + SEEDBYTES + 2 + 4];
    xof_ctx ctx;
    in[0] = DS_EXPAND_SEEDS;
    memcpy(in + 1, xi, SEEDBYTES);
    put_le16(in + 1 + SEEDBYTES, (uint16_t)EM);
    put_le32(in + 1 + SEEDBYTES + 2, kappa);
    xof256_init(&ctx, in, sizeof(in));
    xof256_squeeze(&ctx, T, EXPAND_SEEDS_BYTES);
}

/* ===================================================================== *
 *  ExpandSigningSeeds (DS 0x01, xof256, single one-shot, Sign-only)     *
 * ===================================================================== */
void expand_signing_seeds(uint8_t seedY[SEEDBYTES],
                          const uint8_t K[CHALLENGESEEDBYTES],
                          const uint8_t rnd[RNDBYTES],
                          const uint8_t mu[CHALLENGESEEDBYTES],
                          uint32_t kappa)
{
    /* absorb 0x01 || K || rnd || mu || IntToBytes(kappa,4) */
    uint8_t in[1 + CHALLENGESEEDBYTES + RNDBYTES + CHALLENGESEEDBYTES + 4];
    size_t off = 0;
    xof_ctx ctx;
    in[off++] = DS_EXPAND_SIGNING;
    memcpy(in + off, K, CHALLENGESEEDBYTES);
    off += CHALLENGESEEDBYTES;
    memcpy(in + off, rnd, RNDBYTES);
    off += RNDBYTES;
    memcpy(in + off, mu, CHALLENGESEEDBYTES);
    off += CHALLENGESEEDBYTES;
    put_le32(in + off, kappa);
    off += 4;
    xof256_init(&ctx, in, off);
    xof256_squeeze(&ctx, seedY, SEEDBYTES);
}

/* ===================================================================== *
 *  16-lane uniform-reject stream (ExpandA) -- one-squeeze-per-fill       *
 *                                                                       *
 *  Draws a fixed UNIFORM_BLOCK bytes for (lane, refill) in ONE xof128    *
 *  squeeze; on exhaustion re-inits a fresh ctx with refill+1.            *
 * ===================================================================== */
#define UNIFORM_BLOCK 4096

typedef struct {
    uint8_t
        nonce[1 + SEEDBYTES + 2 + 2]; /* tag||seed||LE16(lane)||LE16(rc) */
    size_t nonce_len;
    uint16_t lane, refill;
    size_t pos, avail;
    uint8_t buf[UNIFORM_BLOCK];
} uniform_stream;

static void us_fill(uniform_stream *us)
{
    xof_ctx ctx;
    put_le16(us->nonce + 1 + SEEDBYTES, us->lane);
    put_le16(us->nonce + 1 + SEEDBYTES + 2, us->refill);
    xof128_init(&ctx, us->nonce, us->nonce_len);
    xof128_squeeze(&ctx, us->buf, UNIFORM_BLOCK);
    us->pos = 0;
    us->avail = UNIFORM_BLOCK;
}

static void us_init(uniform_stream *us, const uint8_t *seedA,
                    unsigned lane)
{
    us->nonce[0] = DS_EXPAND_A;
    memcpy(us->nonce + 1, seedA, SEEDBYTES);
    us->nonce_len = 1 + SEEDBYTES + 2 + 2;
    us->lane = (uint16_t)lane;
    us->refill = 0;
    us_fill(us);
}

/* Pull BQ bytes (a uniform-Z_q candidate); refill across the block edge.
 */
static uint16_t us_next_candidate(uniform_stream *us)
{
    uint32_t mask = (1u << DQ_BITS) - 1u;
    uint16_t a;
    if (us->pos + BQ > us->avail) {
        us->refill++;
        us_fill(us);
    }
    a = (uint16_t)get_le_masked(us->buf + us->pos, BQ, mask);
    us->pos += BQ;
    return a;
}

#if defined(USE_AVX2_SAMPLER) && defined(__AVX2__)
/* compact8 -- BMI2 pdep/pext compaction of accepted uint16 candidates.
 * `m` is the 8-bit accept mask (bit k = candidate k of this 128-bit half).
 * Compacts the surviving uint16 to the front of `o`, returns popcount(m).
 * Identical to the Lithium reject_block_avx2 helper; stores a full 8 u16
 * (trailing lanes garbage the caller overwrites -- the bulk loop
 * guarantees
 * >= 8 slots of headroom).  Byte-exact to the scalar accept order: pext
 * preserves ascending lane order. */
static inline int compact8(uint16_t *o, __m128i f, unsigned m)
{
    uint64_t e = _pdep_u64(m, 0x0101010101010101ULL); /* byte k <- bit k */
    e = (e << 8) - e;                                 /* 0x00/0xFF/byte  */
    uint64_t off = _pext_u64(0x0E0C0A0806040200ULL, e); /* accepted offs */
    __m128i lo = _mm_cvtsi64_si128((long long)off);
    __m128i hi = _mm_add_epi8(lo, _mm_set1_epi8(1)); /* high byte of u16 */
    __m128i ctrl = _mm_unpacklo_epi8(lo, hi);        /* off0,off0+1,...  */
    _mm_storeu_si128((__m128i *)o, _mm_shuffle_epi8(f, ctrl));
    return _mm_popcnt_u32(m);
}

/* uniform_reject_block_avx2 -- 16-wide AVX2 rejection over a contiguous
 * run of 2-byte candidates already resident in the stream buffer,
 * byte-identical to a run of scalar us_next_candidate() calls.
 *
 * Consumes whole 16-candidate (32-byte) chunks from `src` while >= 16
 * input candidates and >= 16 output slots remain, masking `& (2^DQ-1)` and
 * accepting iff `< q` via the unsigned-compare identity  (f < Q) ==
 * (min_epu16(f,Q-1)==f). Accepted coeffs are compacted into dst[*cnt..) in
 * ascending input order (pext preserves order).  Returns the number of
 * INPUT candidates consumed (always a multiple of 16); the caller advances
 * the byte cursor by 2x that and finishes the sub-16 remainder / refill on
 * the scalar path so EVERY byte matches the scalar schedule.  `incap` =
 * candidates available, `want` = total output count for the chunk. */
static size_t uniform_reject_block_avx2(uint16_t *dst, size_t *cnt,
                                        size_t want, const uint8_t *src,
                                        size_t incap)
{
    const uint32_t qmask = (1u << DQ_BITS) - 1u;
    const __m256i maskv = _mm256_set1_epi16((short)qmask);
    const __m256i qm1 = _mm256_set1_epi16((short)((uint16_t)Q - 1u));
    size_t in = 0;
    while (in + 16 <= incap && *cnt + 16 <= want) {
        __m256i f = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(src + 2 * in)), maskv);
        /* unsigned (f < Q) == (min_epu16(f, Q-1) == f) */
        __m256i acc = _mm256_cmpeq_epi16(_mm256_min_epu16(f, qm1), f);
        uint32_t good = _pext_u32((uint32_t)_mm256_movemask_epi8(acc),
                                  0x55555555u); /* 1 bit / coeff */
        *cnt += (size_t)compact8(dst + *cnt, _mm256_castsi256_si128(f),
                                 good & 0xFFu);
        *cnt +=
            (size_t)compact8(dst + *cnt, _mm256_extracti128_si256(f, 1),
                             (good >> 8) & 0xFFu);
        in += 16;
    }
    return in;
}
#endif /* USE_AVX2_SAMPLER && __AVX2__ */

/* Sample `count` uniform coeffs in [0,q) (ascending index) into `dst`.
 *
 * Byte-exact to the scalar reference (a run of `count` accepted
 * us_next_candidate() draws).  Under AVX2 the bulk of the in-buffer
 * candidates is processed 16-at-a-time (uniform_reject_block_avx2); the
 * sub-16 remainder inside the current buffer and any refill boundary fall
 * back to the scalar us_next_candidate path -- so the consumed bytes, the
 * candidate order, the mask and the unsigned-q test are identical, and the
 * accepted stream matches the scalar oracle exactly (K10). */
static void uniform_reject_chunk(uniform_stream *us, uint16_t *dst,
                                 size_t count)
{
    size_t cnt = 0;
#if defined(USE_AVX2_SAMPLER) && defined(__AVX2__)
    while (cnt < count) {
        /* whole candidates currently resident in the buffer */
        size_t incap = (us->avail - us->pos) / (size_t)BQ;
        if (incap >= 16 && cnt + 16 <= count) {
            /* SIMD bulk: consume whole 16-candidate chunks (>=16 in, >=16
             * out slots).  Returns 0 only when no full chunk fit, in which
             * case the scalar step below makes progress. */
            size_t consumed = uniform_reject_block_avx2(
                dst, &cnt, count, us->buf + us->pos, incap);
            us->pos += consumed * (size_t)BQ;
            if (consumed != 0)
                continue; /* re-evaluate buffer/output capacity */
        }
        /* scalar step: handles the sub-16 input remainder, the sub-16
         * output remainder, and the refill boundary -- all byte-identical
         * to ref. */
        {
            uint16_t a;
            do {
                a = us_next_candidate(us);
            } while (a >= (uint16_t)Q);
            dst[cnt++] = a;
        }
    }
#else
    size_t k;
    for (k = 0; k < count; k++) {
        uint16_t a;
        do {
            a = us_next_candidate(us); /* REPEAT { 2B; mask } UNTIL a<q */
        } while (a >= (uint16_t)Q);
        dst[k] = a;
    }
    cnt = count;
#endif
    (void)cnt;
}

/* ===================================================================== *
 *  ExpandA (DS 0x02, xof128, 16 lanes; A_gen direct in NTT domain, K1)  *
 * ===================================================================== */
void expand_a(poly16 agen[EM], poly16 hAgen[EM * ELL],
              const uint8_t seedA[SEEDBYTES])
{
    /* Flat views over the two output arrays. */
    uint16_t *abar = (uint16_t *)agen;  /* EM*n coeffs   */
    uint16_t *hbar = (uint16_t *)hAgen; /* EM*ELL*n coeffs */
    const size_t na = (size_t)EM * N;
    const size_t nh = (size_t)EM * ELL * N;
    const size_t wa = na / XOF_STREAMS; /* agen chunk width  */
    const size_t wh = nh / XOF_STREAMS; /* hAgen chunk width */
    unsigned t;

    _Static_assert(((EM * N) % XOF_STREAMS) == 0,
                   "EM*n must split into 16 lanes");
    _Static_assert(((EM * ELL * N) % XOF_STREAMS) == 0,
                   "EM*ELL*n must split into 16 lanes");
    _Static_assert(sizeof(poly16) == 2 * N, "poly16 flat layout");

    for (t = 0; t < XOF_STREAMS; t++) {
        uniform_stream us;
        us_init(&us, seedA, t);
        uniform_reject_chunk(&us, abar + (size_t)t * wa, wa);
        uniform_reject_chunk(&us, hbar + (size_t)t * wh, wh);
    }
}

/* ===================================================================== *
 *  16-lane noise/gauss stream (ExpandS, SampleY) -- one-squeeze-per-fill *
 * ===================================================================== */
void gauss_stream_init(gauss_stream *gs, uint8_t tag, const uint8_t *seed,
                       unsigned lane)
{
    gs->tag = tag;
    gs->seed = seed;
    gs->lane = (uint16_t)lane;
    gs->refill = 0;
    gs->pos = 0;
    gs->avail = 0;
}

static void gs_fill(gauss_stream *gs)
{
    uint8_t nonce[1 + CHALLENGESEEDBYTES + 2 + 2];
    size_t seedlen =
        (gs->tag == DS_SAMPLE_Y) ? SEEDBYTES : CHALLENGESEEDBYTES;
    size_t nlen;
    xof_ctx ctx;
    nonce[0] = gs->tag;
    memcpy(nonce + 1, gs->seed, seedlen);
    put_le16(nonce + 1 + seedlen, gs->lane);
    put_le16(nonce + 1 + seedlen + 2, gs->refill);
    nlen = 1 + seedlen + 2 + 2;
    xof256_init(&ctx, nonce, nlen);
    xof256_squeeze(&ctx, gs->buf + gs->avail, GAUSS_STREAM_BLOCK);
    gs->avail += GAUSS_STREAM_BLOCK;
}

void gs_ensure(gauss_stream *gs, size_t need)
{
    if (gs->avail - gs->pos >= need)
        return;
    {
        size_t left = gs->avail - gs->pos;
        if (left)
            memmove(gs->buf, gs->buf + gs->pos, left);
        gs->pos = 0;
        gs->avail = left;
        gs->refill++;
        gs_fill(gs);
    }
}

/* ===================================================================== *
 *  ExpandS (DS 0x03, xof256, 16 lanes; BaseSampler + zero-fold + sign)  *
 * ===================================================================== */
/* One noise mini-batch: cdt_scan96 over `Z`/`entries` (NOISE_BATCH=32; the
 * AVX2 cdt_scan96 from the sampler.c fork, bit-identical to scalar), then
 * a 2-bit-per-candidate tail (bit0 sign, bit1 zero-fold), 4 cand/byte; the
 * WHOLE NOISE_MINIBATCH_RAND_BYTES (392) is consumed up front so the
 * cursor advances independently of the cnt==want early break (K6/K8).
 * Identical to ref. */
static void noise_minibatch(gauss_stream *gs, int32_t *dst, size_t *cnt,
                            size_t want, const uint32_t Z[][3],
                            int entries)
{
    int32_t mag[NOISE_BATCH];
    const uint8_t *tailp;
    int j;
    gs_ensure(gs, NOISE_MINIBATCH_RAND_BYTES);
    noise_magnitude_batch(mag, gs->buf + gs->pos, Z, entries);
    tailp = gs->buf + gs->pos + NOISE_CDT_BYTES; /* 2-bit fields */
    for (j = 0; j < NOISE_BATCH; j++) {
        uint32_t f = (uint32_t)(tailp[j >> 2] >> (2 * (j & 3))) & 3u;
        uint32_t reject = ct_is_zero_u32((uint32_t)mag[j]) & (f >> 1);
        int32_t r = ct_sel_i32(f & 1u, -mag[j], mag[j]);
        if (*cnt < want &&
            !reject) /* placement only; cursor already advanced */
            dst[(*cnt)++] = r;
    }
    gs->pos +=
        NOISE_MINIBATCH_RAND_BYTES; /* WHOLE tail consumed (K6/K8) */
}

void expand_s(poly s1s2[ELL + EM],
              const uint8_t seedsk[CHALLENGESEEDBYTES])
{
    int32_t *sbar = s1s2[0].coeffs;   /* ELL*n contiguous */
    int32_t *ebar = s1s2[ELL].coeffs; /* EM*n contiguous  */
    const size_t ns = (size_t)ELL * N;
    const size_t ne = (size_t)EM * N;
    const size_t ws = ns / XOF_STREAMS; /* s chunk width  */
    const size_t we = ne / XOF_STREAMS; /* e chunk width  */
    unsigned t;

    _Static_assert(((ELL * N) % XOF_STREAMS) == 0,
                   "ELL*n must split into 16 lanes");
    _Static_assert(((EM * N) % XOF_STREAMS) == 0,
                   "EM*n must split into 16 lanes");
    _Static_assert(sizeof(poly) == 4 * N, "poly flat layout");
    _Static_assert(ELL + EM >= 2,
                   "s1s2 holds ELL s-polys then EM e-polys");

    for (t = 0; t < XOF_STREAMS; t++) {
        gauss_stream gs;
        size_t cnt;
        gauss_stream_init(&gs, DS_EXPAND_S, seedsk, t);
        cnt = 0;
        while (cnt < ws)
            noise_minibatch(&gs, sbar + (size_t)t * ws, &cnt, ws,
                            RCDT_NOISE_S, RCDT_NOISE_S_ENTRIES);
        cnt = 0;
        while (cnt < we)
            noise_minibatch(&gs, ebar + (size_t)t * we, &cnt, we,
                            RCDT_NOISE_E, RCDT_NOISE_E_ENTRIES);
    }
}

/* ===================================================================== *
 *  SampleC (DS 0x07, xof256, single ctx; partial Fisher-Yates)          *
 * ===================================================================== */
void sample_c(poly *c, const uint8_t seedC[CHALLENGESEEDBYTES])
{
    uint8_t nonce[1 + CHALLENGESEEDBYTES + 2 + 2];
    uint8_t block[UNIFORM_BLOCK];
    uint16_t refill = 0;
    size_t pos = 0, avail = 0;
    uint32_t mask = (1u << DN_BITS) - 1u;
    int i;
    xof_ctx ctx;

    nonce[0] = DS_SAMPLE_C;
    memcpy(nonce + 1, seedC, CHALLENGESEEDBYTES);
    put_le16(nonce + 1 + CHALLENGESEEDBYTES, 0u); /* single stream idx 0 */

    for (i = 0; i < N; i++)
        c->coeffs[i] = 0;

    for (i = N - TAU; i < N; i++) {
        uint32_t j;
        do {
            if (pos + BN > avail) {
                put_le16(nonce + 1 + CHALLENGESEEDBYTES + 2, refill);
                refill++;
                xof256_init(&ctx, nonce, sizeof(nonce));
                xof256_squeeze(&ctx, block, UNIFORM_BLOCK);
                pos = 0;
                avail = UNIFORM_BLOCK;
            }
            j = get_le_masked(block + pos, BN, mask);
            pos += BN;
        } while (j > (uint32_t)i);
        c->coeffs[i] = c->coeffs[j];
        c->coeffs[j] = 1;
    }
}

/* ===================================================================== *
 *  SampleDGauss / SampleY (DS 0x08, xof256, 16 lanes; wide Gaussian r)  *
 * ===================================================================== */
/* gauss_stream_chunk: produce a flat run of `count` wide-Gaussian samples.
 *   1. up-front OUTPUT-indexed sign stream: (count+7)/8 bytes (K9).
 *   2. mini-batches of GAUSS_BATCH candidates until `count` accepted:
 *        - cdt_scan96(x[32], buf+pos) over RCDT_Z (AVX2 fork; advance 384)
 *        - y[32] = buf[pos..pos+32]   (advance 32; Y_BITS=8 byte copy)
 *        - p_hat[32] via approx_exp_accept_q64_x4 (8 groups of 4)
 *        - per-candidate gauss_finalize, tail at buf+pos+8*j (advance 256)
 *      The whole MINIBATCH_RAND_BYTES is consumed up front (K6/K8).
 * Identical to ref; only sampler_sigma2 (the cdt_scan96) is the AVX2 fork.
 */
void gauss_stream_chunk(gauss_stream *gs, int32_t *dst, size_t count)
{
    uint8_t signs[SIGN_BYTES_PER_CHUNK + SIGN_PAD_AVX512];
    size_t signbytes = (count + 7) / 8;
    size_t coefcnt = 0;

    memset(signs, 0,
           sizeof(signs)); /* pad zero for the AVX512 LE64 window */
    gs_ensure(gs, signbytes);
    memcpy(signs, gs->buf + gs->pos, signbytes);
    gs->pos += signbytes;

    while (coefcnt < count) {
        int32_t x[GAUSS_BATCH];
        int32_t yv[GAUSS_BATCH];
        uint64_t phat[GAUSS_BATCH];
        const uint8_t *yp, *tailp;
        int j;

        gs_ensure(gs, MINIBATCH_RAND_BYTES);
        sampler_sigma2(x, gs->buf + gs->pos);
        yp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES;
        for (j = 0; j < GAUSS_BATCH; j++)
            yv[j] = (int32_t)yp[j]; /* Y_BITS=8: plain byte copy */
        for (j = 0; j < GAUSS_BATCH; j += 4) {
            int xi[4], yi[4];
            uint64_t po[4];
            int g;
            for (g = 0; g < 4; g++) {
                xi[g] = (int)x[j + g];
                yi[g] = (int)yv[j + g];
            }
            approx_exp_accept_q64_x4(xi, yi, po);
            for (g = 0; g < 4; g++)
                phat[j + g] = po[g];
        }
        tailp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES + Y_RAND_BYTES;
        for (j = 0; j < GAUSS_BATCH; j++) {
            int32_t r;
            uint32_t sgn;
            size_t idx = coefcnt; /* OUTPUT index of the NEXT accept */
            sgn = (uint32_t)(signs[idx >> 3] >> (idx & 7)) & 1u;
            if (gauss_finalize(&r, x[j], yv[j], phat[j],
                               tailp + (size_t)j * GAUSS_RAND_BYTES,
                               sgn)) {
                if (coefcnt < count)
                    dst[coefcnt++] = r;
            }
        }
        gs->pos += MINIBATCH_RAND_BYTES; /* WHOLE tail consumed (K6/K8) */
    }
}

void sample_y(poly y[KVEC], const uint8_t seedY[SEEDBYTES])
{
    int32_t *ybar = y[0].coeffs; /* KVEC*n contiguous */
    const size_t ntot = (size_t)KVEC * N;
    const size_t wy = ntot / XOF_STREAMS; /* lane chunk width (coeffs) */
    unsigned t;

    _Static_assert(((KVEC * N) % XOF_STREAMS) == 0,
                   "KVEC*n must split into 16 lanes");
    _Static_assert(sizeof(poly) == 4 * N, "poly flat layout");

    for (t = 0; t < XOF_STREAMS; t++) {
        gauss_stream gs;
        gauss_stream_init(&gs, DS_SAMPLE_Y, seedY, t);
        gauss_stream_chunk(&gs, ybar + (size_t)t * wy, wy);
    }
}
