/*
 * polyvec.c -- vector-level XOF-driven sampling surface (AVX-512 FORK).
 * P07 / M9.
 *
 * =========================================================================
 *  What this is
 * =========================================================================
 * This is the AVX-512 fork of ref/polyvec.c.  It is BYTE-EXACT to the
 * scalar reference: every Expand / Sample routine consumes the IDENTICAL
 * PRNG stream (same nonces, same order, same amount) and produces the
 * IDENTICAL sampled coefficients, so make check-kat reproduces the EXACT
 * ref NGCC/SHA3 hashes.  Only the INSTRUCTION SELECTION of two hot inner
 * loops changes:
 *
 *   1. ExpandA uniform-[0,q) rejection scan (uniform_reject_chunk):
 *      vectorized to a 32-wide AVX-512 reject + vpcompressw compaction,
 *      consuming the SAME contiguous 2-byte candidates in the SAME
 * ascending order with the SAME mask + unsigned-q test as the scalar
 *      us_next_candidate -> dst stream.  (ExpandA is ~25% of Sign and
 *      ~83% of Verify -- the headline lever.)  This is the AVX-512 analog
 * of the avx2 BMI2 pdep/pext path: a native unsigned 16-bit compare
 * (_mm512_cmplt_epu16_mask) yields a 32-bit accept mask, and
 *      _mm512_mask_compressstoreu_epi16 (vpcompressw, requires VBMI2)
 * gathers the surviving uint16 to the front in ascending lane order.
 *
 *   2. ExpandS / SampleY BaseSampler magnitude scan goes through the
 * AVX-512 cdt_scan96 (the P05 sampler.c fork), which is bit-identical to
 * the scalar borrow chain.  The wide-Gaussian finalize, the
 * zero-fold/sign logic and the ApproxExp Bernoulli
 * (approx_exp_accept_q64_x4) are kept scalar (cheap per candidate, already
 * bit-exact, low register pressure)
 *      -- they are called from this file gauss_stream_chunk /
 * noise_minibatch exactly as the scalar reference does.
 *
 * Everything else (ExpandSeeds/ExpandSigningSeeds/SampleC, the
 * gauss_stream orchestration, the per-attempt byte schedule, the
 * OUTPUT-indexed sign stream) is copied VERBATIM from ref/polyvec.c.  When
 * -DUSE_AVX512_SAMPLER is NOT defined the file is bit-for-bit the scalar
 * reference.
 *
 * Read polyvec.h FIRST for the binding NGCC-DRBG no-rate-cursor rule and
 * the 16-stream nonce layout.  The PINNED PRNG byte schedules are
 * documented in ref/polyvec.c and are UNCHANGED here.
 */
#include "polyvec.h"

#include <stdlib.h> /* malloc/free for the N-way batched-fill scratch */
#include <string.h>

#include "approx_exp.h" /* approx_exp_accept_q64_x4, approx_exp_accept_q64 */
#include "rcdt_tables.h" /* SHUTTLE_RCDT_Z, SHUTTLE_RCDT_NOISE_* */
#include "test/prof.h" /* PT_G_SHAKE/BASESAMP/APPROXEXP/FINAL -- ((void)0) unless PROF_TIME */

#if defined(USE_AVX512_SAMPLER) && defined(__AVX512F__)
#    include <immintrin.h>
#endif

/* ===================================================================== *
 *  N-way batched XOF refill (M9, USE_AVX512_XOF_NWAY)                    *
 * ===================================================================== *
 *
 *  This mirrors the avx2 fork's xof_nway_fill16 / gs_batch_first_fill
 *  design at the AVX-512 lane width.  ExpandA / ExpandS / SampleY each
 *  partition their work into the fixed XOF_STREAMS (=16) logical streams
 *  (K6).  The scalar reference fills the 16 lanes one at a time with 16
 *  SEQUENTIAL single-stream xof128/256 squeezes -- and profiling shows that
 *  squeeze (the SM3/SHAKE XOF) is the dominant cost of ExpandA (~83% of
 *  Verify, ~25% of Sign) and a large part of SampleY.
 *
 *  This fork keeps the *bytes* identical but computes the 16 streams'
 *  initial blocks N-AT-A-TIME with the already-built, lane-equivalent N-way
 *  AVX-512 primitives (16-way SM3 under NGCC_MODE, 8-way SHAKE under
 *  SHA3_MODE), filling all 16 lane buffers in XOF_STREAMS/XOF_LANES_AVX512
 *  batched passes (a SINGLE pass of 16-way SM3, or 2 passes of 8-way SHAKE)
 *  instead of 16 sequential single-stream squeezes.
 *
 *  WHY THIS IS BYTE-EXACT (lane-equivalence -- P02-verified, K6):
 *  lane k of an N-way init+squeeze over the per-lane nonce
 *      tag || seed || LE16(stream_idx) || LE16(refill)
 *  produces the IDENTICAL byte stream as the scalar xof128/256 over the
 *  SAME absorbed nonce.  The 16-way SM3 / 8-way SHAKE kernels were proven
 *  byte-for-byte equal to the scalar reference (test_xof_nway, sm3-test, the
 *  existing Keccak gate).  We do NOT change which bytes a stream gets, only
 *  compute them N-at-a-time -- so the sampled coefficients and the KAT are
 *  unchanged.  The SM3 DRBG has no sub-call rate cursor (each lane draws its
 *  WHOLE block in one squeeze), which is exactly the one-squeeze-per-fill /
 *  no-rate-cursor discipline the scalar streams already use (see polyvec.h),
 *  so the 16-lane split is the natural batching.
 *
 *  Only the INITIAL block (refill==0) of every lane is batched here -- that
 *  is the block that is ALWAYS drawn (16 of them per Expand/Sample call).
 *  The statistically rare continuation refills (rc>=1, almost never hit;
 *  see UNIFORM_BLOCK / GAUSS_STREAM_BLOCK sizing) fall back to the scalar
 *  per-lane us_fill / gs_fill, byte-identical to ref.  This keeps the refill
 *  nonce/rc bookkeeping in exactly one place and matters not at all for
 *  throughput.
 */
#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
#    define SHUTTLE_XOF_DECLARE_AVX512 1
#    include "symmetric.h" /* xof_ctx_avx512, xof128/256_avx512_* */

/* Number of N-way passes to cover the 16 logical streams (1 for 16-way SM3,
 * 2 for 8-way SHAKE).  XOF_STREAMS is a multiple of XOF_LANES_AVX512 for both
 * MODEs (16 % 16 == 0, 16 % 8 == 0), so the passes tile the streams exactly
 * with no partial trailing pass. */
#    define XOF_NWAY_PASSES (XOF_STREAMS / XOF_LANES_AVX512)
_Static_assert((XOF_STREAMS % XOF_LANES_AVX512) == 0,
               "16 logical streams must tile the N-way lane count exactly");

/* Batch-fill the refill==0 block of all XOF_STREAMS lanes with the N-way
 * XOF.  `nonces[t]` is the fully-built absorbed nonce for stream t
 * (tag||seed||LE16(t)||LE16(0)); `nonce_len` is shared (the producer builds
 * them all the same length).  `dst[t]` receives `block_len` bytes for
 * stream t.  `use_xof128` selects the 128 (public, ExpandA) vs 256 family;
 * under NGCC_MODE the two collapse to the same SM3 DRBG (MS-C5).
 *
 * Lane k of pass p serves logical stream (p*XOF_LANES_AVX512 + k), matching
 * shuttle_xof_stream_of() -- so dst[stream_idx] gets lane k's bytes, which
 * the N-way primitive guarantees equals the scalar xof over nonces[stream]. */
static void xof_nway_fill16(uint8_t *const dst[XOF_STREAMS],
                            const uint8_t *const nonces[XOF_STREAMS],
                            size_t nonce_len, size_t block_len,
                            int use_xof128)
{
    unsigned pass, k;
    for (pass = 0; pass < (unsigned)XOF_NWAY_PASSES; pass++) {
        const uint8_t *seedp[XOF_LANES_AVX512];
        uint8_t *outp[XOF_LANES_AVX512];
        xof_ctx_avx512 ctx;
        for (k = 0; k < (unsigned)XOF_LANES_AVX512; k++) {
            unsigned s = pass * (unsigned)XOF_LANES_AVX512 + k;
            seedp[k] = nonces[s];
            outp[k] = dst[s];
        }
        if (use_xof128) {
            xof128_avx512_init(&ctx, seedp, nonce_len);
            xof128_avx512_squeeze(&ctx, outp, block_len);
        } else {
            xof256_avx512_init(&ctx, seedp, nonce_len);
            xof256_avx512_squeeze(&ctx, outp, block_len);
        }
    }
}
#endif /* USE_AVX512_XOF_NWAY && __AVX512F__ */

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

/* us_setup: build the per-lane nonce + state WITHOUT drawing the first
 * block (so the 16 lanes' initial fills can be batched N-way).  The nonce's
 * LE16(lane) / LE16(refill) fields are written here for refill==0 so the
 * batched N-way fill can absorb us->nonce directly. */
static void us_setup(uniform_stream *us, const uint8_t *seedA, unsigned lane)
{
    us->nonce[0] = DS_EXPAND_A;
    memcpy(us->nonce + 1, seedA, SEEDBYTES);
    us->nonce_len = 1 + SEEDBYTES + 2 + 2;
    us->lane = (uint16_t)lane;
    us->refill = 0;
    put_le16(us->nonce + 1 + SEEDBYTES, us->lane);       /* LE16(lane)   */
    put_le16(us->nonce + 1 + SEEDBYTES + 2, us->refill); /* LE16(rc=0)   */
    us->pos = 0;
    us->avail = 0; /* not yet filled */
}

static void us_init(uniform_stream *us, const uint8_t *seedA,
                    unsigned lane)
{
    us_setup(us, seedA, lane);
    us_fill(us); /* scalar single-stream fill (refill==0) */
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

#if defined(USE_AVX512_SAMPLER) && defined(__AVX512F__)
/* uniform_reject_block_avx512 -- 32-wide AVX-512 rejection over a contiguous
 * run of 2-byte candidates already resident in the stream buffer,
 * byte-identical to a run of scalar us_next_candidate() calls.
 *
 * Consumes whole 32-candidate (64-byte) chunks from `src` while >= 32
 * input candidates and >= 32 output slots remain, masking `& (2^DQ-1)` and
 * accepting iff `< q` via the NATIVE unsigned 16-bit compare
 * _mm512_cmplt_epu16_mask (the AVX-512 analog of the avx2
 * min_epu16==self trick).  The 32-bit accept mask drives
 * _mm512_mask_compressstoreu_epi16 (vpcompressw, VBMI2): the surviving
 * uint16 are written to dst[*cnt..) in ascending input order (vpcompressw
 * preserves lane order), exactly like the scalar accept stream.  Returns
 * the number of INPUT candidates consumed (always a multiple of 32); the
 * caller advances the byte cursor by 2x that and finishes the sub-32
 * remainder / refill on the scalar path so EVERY byte matches the scalar
 * schedule.  `incap` = candidates available, `want` = total output count
 * for the chunk. */
static size_t uniform_reject_block_avx512(uint16_t *dst, size_t *cnt,
                                          size_t want, const uint8_t *src,
                                          size_t incap)
{
    const uint16_t qmask = (uint16_t)((1u << DQ_BITS) - 1u);
    const __m512i maskv = _mm512_set1_epi16((short)qmask);
    const __m512i qv = _mm512_set1_epi16((short)(uint16_t)Q);
    size_t in = 0;
    while (in + 32 <= incap && *cnt + 32 <= want) {
        __m512i f = _mm512_and_si512(
            _mm512_loadu_si512((const __m512i *)(src + 2 * in)), maskv);
        /* native unsigned (f < Q) -> 32-bit accept mask */
        __mmask32 acc = _mm512_cmplt_epu16_mask(f, qv);
        /* compact surviving uint16 to the front of dst (ascending order) */
        _mm512_mask_compressstoreu_epi16(dst + *cnt, acc, f);
        *cnt += (size_t)_mm_popcnt_u32((unsigned)acc);
        in += 32;
    }
    return in;
}
#endif /* USE_AVX512_SAMPLER && __AVX512F__ */

/* Sample `count` uniform coeffs in [0,q) (ascending index) into `dst`.
 *
 * Byte-exact to the scalar reference (a run of `count` accepted
 * us_next_candidate() draws).  Under AVX-512 the bulk of the in-buffer
 * candidates is processed 32-at-a-time (uniform_reject_block_avx512); the
 * sub-32 remainder inside the current buffer and any refill boundary fall
 * back to the scalar us_next_candidate path -- so the consumed bytes, the
 * candidate order, the mask and the unsigned-q test are identical, and the
 * accepted stream matches the scalar oracle exactly (K10). */
static void uniform_reject_chunk(uniform_stream *us, uint16_t *dst,
                                 size_t count)
{
    size_t cnt = 0;
#if defined(USE_AVX512_SAMPLER) && defined(__AVX512F__)
    while (cnt < count) {
        /* whole candidates currently resident in the buffer */
        size_t incap = (us->avail - us->pos) / (size_t)BQ;
        if (incap >= 32 && cnt + 32 <= count) {
            /* SIMD bulk: consume whole 32-candidate chunks (>=32 in, >=32
             * out slots).  Returns 0 only when no full chunk fit, in which
             * case the scalar step below makes progress. */
            size_t consumed = uniform_reject_block_avx512(
                dst, &cnt, count, us->buf + us->pos, incap);
            us->pos += consumed * (size_t)BQ;
            if (consumed != 0)
                continue; /* re-evaluate buffer/output capacity */
        }
        /* scalar step: handles the sub-32 input remainder, the sub-32
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

#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
    /* N-WAY BATCHED INITIAL FILL: draw the refill==0 UNIFORM_BLOCK of all
     * XOF_STREAMS lanes N-at-a-time (16-way SM3 = 1 pass / 8-way SHAKE = 2
     * passes), byte-exact to 16 sequential scalar us_fill() (lane-equivalence,
     * K6).  ExpandA is the headline lever: this replaces 16 sequential xof128
     * squeezes -- the dominant cost of Verify (~83%) and a big chunk of Sign
     * -- with XOF_STREAMS/XOF_LANES_AVX512 batched squeezes.  Then the
     * (already vectorized) per-lane uniform_reject_chunk runs on each
     * pre-filled buffer, handling any rare continuation refill on the scalar
     * path. */
    {
        uniform_stream *us = (uniform_stream *)malloc(
            (size_t)XOF_STREAMS * sizeof(uniform_stream));
        const uint8_t *nonces[XOF_STREAMS];
        uint8_t *dst[XOF_STREAMS];
        if (us) {
            for (t = 0; t < XOF_STREAMS; t++) {
                us_setup(&us[t], seedA, t);
                nonces[t] = us[t].nonce;
                dst[t] = us[t].buf;
            }
            /* ExpandA uses xof128 (public material, tag 0x02). */
            xof_nway_fill16(dst, nonces, us[0].nonce_len, UNIFORM_BLOCK,
                            /*use_xof128=*/1);
            for (t = 0; t < XOF_STREAMS; t++) {
                us[t].pos = 0;
                us[t].avail = UNIFORM_BLOCK; /* batched fill complete */
                uniform_reject_chunk(&us[t], abar + (size_t)t * wa, wa);
                uniform_reject_chunk(&us[t], hbar + (size_t)t * wh, wh);
            }
            free(us);
            return;
        }
        /* malloc failure: fall through to the scalar per-lane path below. */
    }
#endif

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

/* Build the per-(lane,refill) absorbed nonce for a gauss_stream into
 * `nonce`; returns the nonce length.  Shared by gs_fill (scalar) and the
 * N-way batched first fill so the byte layout is defined in exactly one
 * place (tag || seed || LE16(lane) || LE16(refill)). */
static size_t gs_build_nonce(const gauss_stream *gs, uint8_t *nonce)
{
    size_t seedlen =
        (gs->tag == DS_SAMPLE_Y) ? SEEDBYTES : CHALLENGESEEDBYTES;
    nonce[0] = gs->tag;
    memcpy(nonce + 1, gs->seed, seedlen);
    put_le16(nonce + 1 + seedlen, gs->lane);
    put_le16(nonce + 1 + seedlen + 2, gs->refill);
    return 1 + seedlen + 2 + 2;
}

static void gs_fill(gauss_stream *gs)
{
    uint8_t nonce[1 + CHALLENGESEEDBYTES + 2 + 2];
    size_t nlen;
    xof_ctx ctx;
    nlen = gs_build_nonce(gs, nonce);
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

#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
/* gs_batch_first_fill: draw the FIRST GAUSS_STREAM_BLOCK of all XOF_STREAMS
 * gauss_streams N-at-a-time (xof256: 16-way SM3 = 1 pass / 8-way SHAKE = 2
 * passes).
 *
 * BYTE-EXACT to 16 sequential scalar first fills.  Note the scalar discipline
 * (gauss_stream_init sets refill=0; the first gs_ensure does refill++ THEN
 * gs_fill): the first block is therefore drawn with refill==1 in the absorbed
 * nonce.  We replicate that EXACTLY -- set each lane's refill to 1, build the
 * nonce with LE16(refill=1), batch-fill buf[0..GAUSS_STREAM_BLOCK), and set
 * pos=0 / avail=GAUSS_STREAM_BLOCK.  Subsequent (rare) continuation refills
 * then run the scalar gs_ensure path with refill = 2,3,... -- byte-identical
 * to ref.  Lane k of the N-way squeeze over nonce(stream s, rc=1) equals the
 * scalar xof256 over the same nonce (lane-equivalence, K6).
 *
 * `gss[t]` must already be gauss_stream_init()'d (tag/seed/lane set). */
static void gs_batch_first_fill(gauss_stream gss[XOF_STREAMS])
{
    uint8_t nonces[XOF_STREAMS][1 + CHALLENGESEEDBYTES + 2 + 2];
    const uint8_t *noncep[XOF_STREAMS];
    uint8_t *dst[XOF_STREAMS];
    size_t nlen = 0;
    unsigned t;
    for (t = 0; t < XOF_STREAMS; t++) {
        gss[t].refill = 1; /* matches the scalar refill++ on the 1st fill */
        nlen = gs_build_nonce(&gss[t], nonces[t]);
        noncep[t] = nonces[t];
        dst[t] = gss[t].buf;
    }
    /* All gauss streams use xof256 (ExpandS tag 0x03 / SampleY tag 0x08). */
    xof_nway_fill16(dst, noncep, nlen, GAUSS_STREAM_BLOCK, /*use_xof128=*/0);
    for (t = 0; t < XOF_STREAMS; t++) {
        gss[t].pos = 0;
        gss[t].avail = GAUSS_STREAM_BLOCK;
    }
}
#endif /* USE_AVX512_XOF_NWAY && __AVX512F__ */

/* ===================================================================== *
 *  ExpandS (DS 0x03, xof256, 16 lanes; BaseSampler + zero-fold + sign)  *
 * ===================================================================== */
/* One noise mini-batch: cdt_scan96 over `Z`/`entries` (NOISE_BATCH=32; the
 * AVX-512 cdt_scan96 from the sampler.c fork, bit-identical to scalar),
 * then a 2-bit-per-candidate tail (bit0 sign, bit1 zero-fold), 4 cand/byte;
 * the WHOLE NOISE_MINIBATCH_RAND_BYTES (392) is consumed up front so the
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

#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
    /* N-WAY BATCHED INITIAL FILL: prime all 16 gauss_streams, then draw
     * their first GAUSS_STREAM_BLOCK N-at-a-time (xof256: 16-way SM3 = 1 pass
     * / 8-way SHAKE = 2 passes) instead of 16 sequential single-stream
     * squeezes.  Byte-exact to ref (lane-equivalence, K6); the per-lane
     * BaseSampler scan + zero-fold/sign logic is then run on each pre-filled
     * buffer exactly as scalar. */
    {
        gauss_stream *gss = (gauss_stream *)malloc(
            (size_t)XOF_STREAMS * sizeof(gauss_stream));
        if (gss) {
            for (t = 0; t < XOF_STREAMS; t++)
                gauss_stream_init(&gss[t], DS_EXPAND_S, seedsk, t);
            gs_batch_first_fill(gss);
            for (t = 0; t < XOF_STREAMS; t++) {
                size_t cnt = 0;
                while (cnt < ws)
                    noise_minibatch(&gss[t], sbar + (size_t)t * ws, &cnt,
                                    ws, RCDT_NOISE_S, RCDT_NOISE_S_ENTRIES);
                cnt = 0;
                while (cnt < we)
                    noise_minibatch(&gss[t], ebar + (size_t)t * we, &cnt,
                                    we, RCDT_NOISE_E, RCDT_NOISE_E_ENTRIES);
            }
            free(gss);
            return;
        }
        /* malloc failure: fall through to the scalar per-lane path. */
    }
#endif

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
 *        - cdt_scan96(x[32], buf+pos) over RCDT_Z (AVX-512 fork; advance 384)
 *        - y[32] = buf[pos..pos+32]   (advance 32; Y_BITS=8 byte copy)
 *        - p_hat[32] via approx_exp_accept_q64_x4 (8 groups of 4)
 *        - per-candidate gauss_finalize, tail at buf+pos+8*j (advance 256)
 *      The whole MINIBATCH_RAND_BYTES is consumed up front (K6/K8).
 * Identical to ref; only sampler_sigma2 (the cdt_scan96) is the AVX-512 fork.
 */
void gauss_stream_chunk(gauss_stream *gs, int32_t *dst, size_t count)
{
    uint8_t signs[SIGN_BYTES_PER_CHUNK + SIGN_PAD_AVX512];
    size_t signbytes = (count + 7) / 8;
    size_t coefcnt = 0;

    memset(signs, 0,
           sizeof(signs)); /* pad zero for the AVX512 LE64 window */
    {
        PROF_START(t_sh0);
        gs_ensure(gs, signbytes);
        PROF_STOP(PT_G_SHAKE, t_sh0);
    }
    memcpy(signs, gs->buf + gs->pos, signbytes);
    gs->pos += signbytes;

    while (coefcnt < count) {
        int32_t x[GAUSS_BATCH];
        int32_t yv[GAUSS_BATCH];
        uint64_t phat[GAUSS_BATCH];
        const uint8_t *yp, *tailp;
        int j;

        {
            PROF_START(t_sh);
            gs_ensure(gs, MINIBATCH_RAND_BYTES);
            PROF_STOP(PT_G_SHAKE, t_sh);
        }
        {
            PROF_START(t_bs);
            sampler_sigma2(x, gs->buf + gs->pos);
            PROF_STOP(PT_G_BASESAMP, t_bs);
        }
        yp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES;
        for (j = 0; j < GAUSS_BATCH; j++)
            yv[j] = (int32_t)yp[j]; /* Y_BITS=8: plain byte copy */
        {
            PROF_START(t_ae);
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
            PROF_STOP(PT_G_APPROXEXP, t_ae);
        }
        tailp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES + Y_RAND_BYTES;
        {
            PROF_START(t_fin);
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
            PROF_STOP(PT_G_FINAL, t_fin);
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

#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
    /* N-WAY BATCHED INITIAL FILL: prime all 16 wide-Gaussian streams, draw
     * their first GAUSS_STREAM_BLOCK N-at-a-time (xof256: 16-way SM3 = 1 pass
     * / 8-way SHAKE = 2 passes), then run the per-lane wide-Gaussian chunk on
     * each pre-filled buffer.  Byte-exact to ref (lane-equivalence, K6);
     * SampleY is a large fraction of Sign, and its XOF squeeze shrinks here. */
    {
        gauss_stream *gss = (gauss_stream *)malloc(
            (size_t)XOF_STREAMS * sizeof(gauss_stream));
        if (gss) {
            for (t = 0; t < XOF_STREAMS; t++)
                gauss_stream_init(&gss[t], DS_SAMPLE_Y, seedY, t);
            {
                PROF_START(t_nway);
                gs_batch_first_fill(gss); /* N-way XOF first fill */
                PROF_STOP(PT_G_SHAKE, t_nway);
            }
            for (t = 0; t < XOF_STREAMS; t++)
                gauss_stream_chunk(&gss[t], ybar + (size_t)t * wy, wy);
            free(gss);
            return;
        }
        /* malloc failure: fall through to the scalar per-lane path. */
    }
#endif

    for (t = 0; t < XOF_STREAMS; t++) {
        gauss_stream gs;
        gauss_stream_init(&gs, DS_SAMPLE_Y, seedY, t);
        gauss_stream_chunk(&gs, ybar + (size_t)t * wy, wy);
    }
}
