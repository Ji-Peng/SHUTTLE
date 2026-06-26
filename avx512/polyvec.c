/*
 * polyvec.c -- vector-level XOF-driven sampling surface (AVX-512 FORK).
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
 * AVX-512 cdt_scan96 (the sampler.c fork), which is bit-identical to
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
 *  N-way batched XOF refill (USE_AVX512_XOF_NWAY)                        *
 * ===================================================================== *
 *
 *  This mirrors the avx2 fork's xof_nway_fill16 / gs_batch_first_fill
 *  design at the AVX-512 lane width.  ExpandA / ExpandS / SampleY each
 *  partition their work into the fixed XOF_STREAMS (=16) logical streams.
 *  The scalar reference fills the 16 lanes one at a time with 16
 *  SEQUENTIAL single-stream xof128/256 squeezes -- and profiling shows
 * that squeeze (the SM3/SHAKE XOF) is the dominant cost of ExpandA (~83%
 * of Verify, ~25% of Sign) and a large part of SampleY.
 *
 *  This fork keeps the *bytes* identical but computes the 16 streams'
 *  initial blocks N-AT-A-TIME with the already-built, lane-equivalent
 * N-way AVX-512 primitives (16-way SM3 under NGCC_MODE, 8-way SHAKE under
 *  SHA3_MODE), filling all 16 lane buffers in XOF_STREAMS/XOF_LANES_AVX512
 *  batched passes (a SINGLE pass of 16-way SM3, or 2 passes of 8-way
 * SHAKE) instead of 16 sequential single-stream squeezes.
 *
 *  WHY THIS IS BYTE-EXACT (lane-equivalence, test-verified):
 *  lane k of an N-way init+squeeze over the per-lane nonce
 *      tag || seed || LE16(stream_idx) || LE16(refill)
 *  produces the IDENTICAL byte stream as the scalar xof128/256 over the
 *  SAME absorbed nonce.  The 16-way SM3 / 8-way SHAKE kernels were proven
 *  byte-for-byte equal to the scalar reference (test_xof_nway, sm3-test,
 * the existing Keccak gate).  We do NOT change which bytes a stream gets,
 * only compute them N-at-a-time -- so the sampled coefficients and the KAT
 * are unchanged.  The SM3 DRBG has no sub-call rate cursor (each lane
 * draws its WHOLE block in one squeeze), which is exactly the
 * one-squeeze-per-fill / no-rate-cursor discipline the scalar streams
 * already use (see polyvec.h), so the 16-lane split is the natural
 * batching.
 *
 *  Only the INITIAL block (refill==0) of every lane is batched here --
 * that is the block that is ALWAYS drawn (16 of them per Expand/Sample
 * call). The statistically rare continuation refills (rc>=1, almost never
 * hit; see UNIFORM_BLOCK / GAUSS_STREAM_BLOCK sizing) fall back to the
 * scalar per-lane us_fill / gs_fill, byte-identical to ref.  This keeps
 * the refill nonce/rc bookkeeping in exactly one place and matters not at
 * all for throughput.
 */
#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
#    define SHUTTLE_XOF_DECLARE_AVX512 1
#    include "symmetric.h" /* xof_ctx_avx512, xof128/256_avx512_* */

/* At most XOF_STREAMS/XOF_LANES_AVX512 N-way passes cover the 16 logical
 * streams (1 for 16-way SM3, 2 for 8-way SHAKE).  XOF_STREAMS is a
 * multiple of XOF_LANES_AVX512 for both MODEs (16 % 16 == 0, 16 % 8 == 0),
 * so the passes tile the streams exactly with no partial trailing pass; a
 * caller that needs only nfill streams drives ceil(nfill / lanes) of them.
 */
_Static_assert(
    (XOF_STREAMS % XOF_LANES_AVX512) == 0,
    "16 logical streams must tile the N-way lane count exactly");

/* Batch-fill the refill==0 block of the leading `nfill` lanes (slots
 * 0..nfill-1) with the N-way XOF.  `nonces[t]` is the fully-built absorbed
 * nonce for stream t (tag||seed||LE16(t)||LE16(0)); `nonce_len` is shared
 * (the producer builds them all the same length).  `dst[t]` receives
 * `block_len` bytes for stream t.  `use_xof128` selects the 128 (public,
 * ExpandA) vs 256 family; under NGCC_MODE the two collapse to the same SM3
 * DRBG (the 128/256 XOF collapse).
 *
 * Callers compact their active lanes into the leading `nfill` slots, so we
 * issue only ceil(nfill / lanes) passes: a pass whose lanes are all >=
 * nfill would fill nothing the caller reads.  When the lane width covers
 * all 16 streams in one pass (16-way SM3) this is always a single pass;
 * the win is in the 8-way SHAKE path, where a caller needing <= 8 lanes
 * skips the second pass entirely.
 *
 * Lane k of pass p serves logical stream (p*XOF_LANES_AVX512 + k),
 * matching shuttle_xof_stream_of() -- so dst[stream_idx] gets lane k's
 * bytes, which the N-way primitive guarantees equals the scalar xof over
 * nonces[stream]. */
static void xof_nway_fill16(uint8_t *const dst[XOF_STREAMS],
                            const uint8_t *const nonces[XOF_STREAMS],
                            size_t nonce_len, size_t block_len,
                            int use_xof128, unsigned nfill)
{
    const unsigned passes = (nfill + (unsigned)XOF_LANES_AVX512 - 1) /
                            (unsigned)XOF_LANES_AVX512;
    unsigned pass, k;
    for (pass = 0; pass < passes; pass++) {
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
 *  16-lane uniform-reject stream (ExpandA)                              *
 *                                                                       *
 *  A per-(lane,refill) XOF instance supplies UNIFORM_BLOCK logical      *
 *  bytes; the sampler crosses to the next instance (refill+1) only after *
 *  all of them are consumed.  We draw the block LAZILY in UNIFORM_DRAW   *
 *  slices from the same persisted ctx -- squeezing N then M bytes yields *
 *  the same bytes as one N+M squeeze (SHAKE/SM3-DRBG are streams), so    *
 *  the consumed byte sequence and the refill boundary are byte-identical *
 *  to one big squeeze; only the wasted tail is never squeezed.          *
 * ===================================================================== */
#define UNIFORM_BLOCK 4096

/* Right-sized slice: ceil-to-granularity of ~2x the mean per-lane byte
 * need (per-lane coeffs * BQ * 2, mean acceptance ~q/2^DQ).  Covers the
 * largest mode in one slice with margin; lazy continuation almost never
 * runs.  The refill boundary is governed by the UNIFORM_BLOCK cap. */
#define UNIFORM_LANE_COEFFS ((size_t)EM * N * (1u + ELL) / XOF_STREAMS)
#define UNIFORM_DRAW_RAW (UNIFORM_LANE_COEFFS * (size_t)BQ * 2u)
#define UNIFORM_DRAW                                                    \
    (((UNIFORM_DRAW_RAW + (size_t)XOF_SQUEEZE_GRANULARITY_BYTES - 1u) / \
      (size_t)XOF_SQUEEZE_GRANULARITY_BYTES) *                          \
     (size_t)XOF_SQUEEZE_GRANULARITY_BYTES)

/* Right-sized first/continuation slice for SampleC's partial Fisher-Yates.
 * The mean candidate need is sum_{i=n-tau}^{n-1} 2^DN/(i+1) BN-byte draws
 * (well under 256 B for every mode); TAU*BN*4 ceil-to-granularity gives a
 * comfortable margin so the first slice covers the whole challenge in
 * essentially every call.  The instance still caps at SAMPLEC_BLOCK bytes,
 * preserving the per-refill XOF boundary byte-for-byte. */
#define SAMPLEC_BLOCK UNIFORM_BLOCK /* logical per-refill instance cap */
#define SAMPLEC_DRAW_RAW ((size_t)TAU * (size_t)BN * 2u)
#define SAMPLEC_DRAW                                                    \
    (((SAMPLEC_DRAW_RAW + (size_t)XOF_SQUEEZE_GRANULARITY_BYTES - 1u) / \
      (size_t)XOF_SQUEEZE_GRANULARITY_BYTES) *                          \
     (size_t)XOF_SQUEEZE_GRANULARITY_BYTES)

typedef struct {
    uint8_t
        nonce[1 + SEEDBYTES + 2 + 2]; /* tag||seed||LE16(lane)||LE16(rc) */
    size_t nonce_len;
    uint16_t lane, refill;
    size_t pos, avail; /* cursor / valid bytes within buf            */
    size_t drawn;      /* bytes squeezed from the current XOF instance */
    int ctx_ready;     /* persisted ctx initialised for this instance  */
    xof_ctx ctx;       /* live XOF instance for incremental slices      */
    uint8_t buf[UNIFORM_DRAW];
} uniform_stream;

/* Squeeze the next slice of the CURRENT XOF instance into buf.  If the
 * persisted ctx is not yet live (the initial slice was filled by the N-way
 * batched path), init it from the nonce and fast-forward past the `drawn`
 * bytes already consumed -- identical bytes either way (stream). */
static void us_draw_slice(uniform_stream *us)
{
    size_t want = UNIFORM_BLOCK - us->drawn;
    if (want > UNIFORM_DRAW)
        want = UNIFORM_DRAW;
    if (!us->ctx_ready) {
        put_le16(us->nonce + 1 + SEEDBYTES, us->lane);
        put_le16(us->nonce + 1 + SEEDBYTES + 2, us->refill);
        xof128_init(&us->ctx, us->nonce, us->nonce_len);
        if (us->drawn) {
            uint8_t skip[UNIFORM_DRAW];
            size_t left = us->drawn;
            while (left) {
                size_t s = left < UNIFORM_DRAW ? left : UNIFORM_DRAW;
                xof128_squeeze(&us->ctx, skip, s);
                left -= s;
            }
        }
        us->ctx_ready = 1;
    }
    xof128_squeeze(&us->ctx, us->buf, want);
    us->pos = 0;
    us->avail = want;
    us->drawn += want;
}

static void us_fill(uniform_stream *us)
{
    us->ctx_ready = 0;
    us->drawn = 0;
    us_draw_slice(us);
}

/* us_setup: build the per-lane nonce + state WITHOUT drawing the first
 * slice (so the 16 lanes' initial fills can be batched N-way).  The
 * nonce's LE16(lane) / LE16(refill) fields are written here for refill==0
 * so the batched N-way fill can absorb us->nonce directly. */
static void us_setup(uniform_stream *us, const uint8_t *seedA,
                     unsigned lane)
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
    us->drawn = 0;
    us->ctx_ready = 0;
}

static void us_init(uniform_stream *us, const uint8_t *seedA,
                    unsigned lane)
{
    us_setup(us, seedA, lane);
    us_fill(us); /* scalar single-stream fill (refill==0) */
}

/* Pull BQ bytes (a uniform-Z_q candidate).  When the current slice runs
 * out, draw the next slice of the same instance; only after the whole
 * UNIFORM_BLOCK is consumed do we advance to the next instance (refill+1).
 */
static uint16_t us_next_candidate(uniform_stream *us)
{
    uint32_t mask = (1u << DQ_BITS) - 1u;
    uint16_t a;
    if (us->pos + BQ > us->avail) {
        if (us->drawn < UNIFORM_BLOCK) {
            us_draw_slice(us); /* same instance, next slice */
        } else {
            us->refill++; /* instance exhausted: next instance */
            us_fill(us);
        }
    }
    a = (uint16_t)get_le_masked(us->buf + us->pos, BQ, mask);
    us->pos += BQ;
    return a;
}

#if defined(USE_AVX512_SAMPLER) && defined(__AVX512F__)
/* uniform_reject_block_avx512 -- 32-wide AVX-512 rejection over a
 * contiguous run of 2-byte candidates already resident in the stream
 * buffer, byte-identical to a run of scalar us_next_candidate() calls.
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
        /* compact surviving uint16 to the front of dst (ascending order)
         */
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
 * accepted stream matches the scalar oracle exactly (constant-time). */
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
 *  ExpandA (DS 0x02, xof128, 16 lanes; A_gen direct in NTT domain)     *
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
    /* N-WAY BATCHED INITIAL FILL: draw the refill==0 FIRST SLICE
     * (UNIFORM_DRAW bytes) of all XOF_STREAMS lanes N-at-a-time (16-way
     * SM3 = 1 pass / 8-way SHAKE = 2 passes), byte-exact to 16 sequential
     * scalar us_fill() slices (lane-equivalence).  ExpandA is the headline
     * lever: this replaces 16 sequential xof128 squeezes -- the dominant
     * cost of Verify (~83%) and a big chunk of Sign -- with
     * XOF_STREAMS/XOF_LANES_AVX512 batched squeezes, now right-sized so
     * the wasted tail is never squeezed.  Then the (already vectorized)
     * per-lane uniform_reject_chunk runs on each pre-filled buffer,
     * drawing any further slice on the scalar path.  UNIFORM_DRAW is a
     * whole multiple of the per-MODE squeeze granularity, so the N-way
     * squeeze is rate-aligned and byte-exact. */
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
            xof_nway_fill16(dst, nonces, us[0].nonce_len, UNIFORM_DRAW,
                            /*use_xof128=*/1, /*nfill=*/XOF_STREAMS);
            for (t = 0; t < XOF_STREAMS; t++) {
                /* First slice resident; ctx not yet live -- a continuation
                 * slice (rare) lazily inits the scalar ctx and skips
                 * drawn. */
                us[t].pos = 0;
                us[t].avail = UNIFORM_DRAW;
                us[t].drawn = UNIFORM_DRAW;
                us[t].ctx_ready = 0;
                uniform_reject_chunk(&us[t], abar + (size_t)t * wa, wa);
                uniform_reject_chunk(&us[t], hbar + (size_t)t * wh, wh);
            }
            free(us);
            return;
        }
        /* malloc failure: fall through to the scalar per-lane path below.
         */
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
    size_t blk = gauss_block_bytes(gs->tag);
    size_t nlen;
    xof_ctx ctx;
    nlen = gs_build_nonce(gs, nonce);
    xof256_init(&ctx, nonce, nlen);
    xof256_squeeze(&ctx, gs->buf + gs->avail, blk);
    gs->avail += blk;
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
/* gs_batch_first_fill: draw the FIRST GAUSS_STREAM_BLOCK of all
 * XOF_STREAMS gauss_streams N-at-a-time (xof256: 16-way SM3 = 1 pass /
 * 8-way SHAKE = 2 passes).
 *
 * BYTE-EXACT to 16 sequential scalar first fills.  Note the scalar
 * discipline (gauss_stream_init sets refill=0; the first gs_ensure does
 * refill++ THEN gs_fill): the first block is therefore drawn with
 * refill==1 in the absorbed nonce.  We replicate that EXACTLY -- set each
 * lane's refill to 1, build the nonce with LE16(refill=1), batch-fill
 * buf[0..GAUSS_STREAM_BLOCK), and set pos=0 / avail=GAUSS_STREAM_BLOCK.
 * Subsequent (rare) continuation refills then run the scalar gs_ensure
 * path with refill = 2,3,... -- byte-identical to ref.  Lane k of the
 * N-way squeeze over nonce(stream s, rc=1) equals the scalar xof256 over
 * the same nonce (lane-equivalence).
 *
 * `gss[t]` must already be gauss_stream_init()'d (tag/seed/lane set). */
static void gs_batch_first_fill(gauss_stream gss[XOF_STREAMS])
{
    uint8_t nonces[XOF_STREAMS][1 + CHALLENGESEEDBYTES + 2 + 2];
    const uint8_t *noncep[XOF_STREAMS];
    uint8_t *dst[XOF_STREAMS];
    /* All lanes in one ExpandS/SampleY call share the tag, hence the same
     * tag-tuned per-refill block size (no over-squeezed tail). */
    size_t blk = gauss_block_bytes(gss[0].tag);
    size_t nlen = 0;
    unsigned t;
    for (t = 0; t < XOF_STREAMS; t++) {
        gss[t].refill =
            1; /* matches the scalar refill++ on the 1st fill */
        nlen = gs_build_nonce(&gss[t], nonces[t]);
        noncep[t] = nonces[t];
        dst[t] = gss[t].buf;
    }
    /* All gauss streams use xof256 (ExpandS tag 0x03 / SampleY tag 0x08).
     */
    xof_nway_fill16(dst, noncep, nlen, blk,
                    /*use_xof128=*/0, /*nfill=*/XOF_STREAMS);
    for (t = 0; t < XOF_STREAMS; t++) {
        gss[t].pos = 0;
        gss[t].avail = blk;
    }
}
#endif /* USE_AVX512_XOF_NWAY && __AVX512F__ */

/* ===================================================================== *
 *  ExpandS (DS 0x03, xof256, 16 lanes; BaseSampler + zero-fold + sign)  *
 * ===================================================================== */
/* noise_consume_minibatch: process ONE NOISE_BATCH mini-batch from the
 * lane buffer (cdt_scan96 over `Z`/`entries`, NOISE_BATCH=32; the AVX-512
 * cdt_scan96 from the sampler.c fork, bit-identical to scalar), then a
 * 2-bit-per-candidate tail (bit0 sign, bit1 zero-fold), 4 cand/byte.  PURE
 * buffer operation -- the caller MUST have ensured >=
 * NOISE_MINIBATCH_RAND_BYTES are resident at gs->buf+gs->pos.  The WHOLE
 * NOISE_MINIBATCH_RAND_BYTES (392) is consumed up front so the cursor
 * advances independently of the *cnt==want early break.  Shared by the
 * scalar (A) and N-way round-robin (B) ExpandS paths so they cannot drift.
 */
static void noise_consume_minibatch(gauss_stream *gs, int32_t *dst,
                                    size_t *cnt, size_t want,
                                    const uint32_t Z[][3], int entries)
{
    int32_t mag[NOISE_BATCH];
    const uint8_t *tailp;
    int j;
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
    gs->pos += NOISE_MINIBATCH_RAND_BYTES; /* WHOLE tail consumed */
}

/* noise_minibatch: scalar single-lane mini-batch (A) -- ensure one
 * mini-batch is resident (one scalar xof256 squeeze on refill) then
 * consume it.  Byte-identical to the shared pure-buffer body above; kept
 * for the non-USE_AVX512_XOF_NWAY build and the malloc-failure fallback.
 */
static void noise_minibatch(gauss_stream *gs, int32_t *dst, size_t *cnt,
                            size_t want, const uint32_t Z[][3],
                            int entries)
{
    gs_ensure(gs, NOISE_MINIBATCH_RAND_BYTES);
    noise_consume_minibatch(gs, dst, cnt, want, Z, entries);
}

#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
/* gs_lane_refill_prep (defined with the SampleY driver below) replicates
 * the scalar gs_ensure(need) bookkeeping for one lane WITHOUT issuing the
 * scalar squeeze, so the ExpandS round-robin driver can batch the squeeze
 * N-way. */
static size_t gs_lane_refill_prep(gauss_stream *gs, uint8_t *nonce,
                                  size_t *fill_at);

/* ------------------------------------------------------------------- *
 *  N-way round-robin ExpandS driver (B).                              *
 * ------------------------------------------------------------------- *
 *
 *  ExpandS's per-lane work is a fixed, PUBLIC sequence of NOISE_BATCH
 *  mini-batches: ceil-driven runs that fill `ws` s-coeffs then `we`
 *  e-coeffs (the run lengths are bounded by the PUBLIC widths -- each
 *  mini-batch's whole NOISE_MINIBATCH_RAND_BYTES tail is consumed
 *  up front, so the number of mini-batches per lane is data-independent).
 *  Like the SampleY driver, every refill across the 16 lanes is issued
 *  N-at-a-time (16-way SM3 / 8-way SHAKE) instead of one scalar squeeze
 * per lane.  Byte-exact to ref (consumed-bytes invariance: each lane's
 *  gauss_stream is the same deterministic nonce schedule, consumed by the
 *  same noise_consume_minibatch body; only WHERE the bytes are produced
 *  moves to the N-way batch). */
typedef struct {
    int32_t *sdst; /* this lane's s-output slice            */
    int32_t *edst; /* this lane's e-output slice            */
    size_t ws;     /* s-coeffs to produce (lane width)      */
    size_t we;     /* e-coeffs to produce (lane width)      */
    size_t scnt;   /* s-coeffs produced so far              */
    size_t ecnt;   /* e-coeffs produced so far              */
    int phase;     /* 0 = filling s, 1 = filling e, 2 done  */
} noise_lane_state;

/* noise_lane_advance_phase: roll the lane to its next non-finished phase.
 */
static void noise_lane_advance_phase(noise_lane_state *st)
{
    if (st->phase == 0 && st->scnt >= st->ws)
        st->phase = (st->we > 0) ? 1 : 2;
    if (st->phase == 1 && st->ecnt >= st->we)
        st->phase = 2;
}

/* expands_roundrobin: drive all XOF_STREAMS ExpandS lanes in lockstep with
 * N-way batched refills.  `gss` are init'd (refill=0, first block already
 * filled by gs_batch_first_fill); `lst` carry the per-lane s/e slices.
 * Loops in rounds:
 *   (1) for each not-done lane whose buffer is short for its next
 *       mini-batch, prep its refill (memmove + refill++ + nonce) and queue
 *       an N-way fill;
 *   (2) one N-way pass (xof_nway_fill16) tops up every queued lane;
 *   (3) each not-done lane consumes ONE mini-batch from its now-resident
 *       buffer (s-table or e-table per its phase).
 * Round count is bounded by the PUBLIC per-lane mini-batch count -- no
 * secret-dependent loop bound. */
static void expands_roundrobin(gauss_stream gss[XOF_STREAMS],
                               noise_lane_state lst[XOF_STREAMS])
{
    int all_done = 0;
    /* ExpandS: every lane shares tag DS_EXPAND_S -> the same per-refill
     * block size (NOISE_MINIBATCH_RAND_BYTES, no over-squeezed tail). */
    size_t blk = gauss_block_bytes(gss[0].tag);
    unsigned t;

    for (t = 0; t < XOF_STREAMS; t++)
        noise_lane_advance_phase(&lst[t]);

    while (!all_done) {
        uint8_t nonces[XOF_STREAMS][1 + CHALLENGESEEDBYTES + 2 + 2];
        const uint8_t *noncep[XOF_STREAMS];
        uint8_t *dstp[XOF_STREAMS];
        /* Per-round throwaway destinations for the inactive N-way slots,
         * so no two slots ever alias a live buffer (each padding slot gets
         * its OWN region; its bytes are discarded). */
        uint8_t scratch[XOF_STREAMS][GAUSS_STREAM_BLOCK];
        size_t nlen = 0;
        unsigned nfill = 0;
        unsigned fill_lane[XOF_STREAMS];

        /* (1) collect the lanes that need a refill for their next
         * mini-batch. */
        for (t = 0; t < XOF_STREAMS; t++) {
            if (lst[t].phase == 2)
                continue;
            if (gss[t].avail - gss[t].pos <
                (size_t)NOISE_MINIBATCH_RAND_BYTES) {
                size_t fill_at;
                nlen =
                    gs_lane_refill_prep(&gss[t], nonces[nfill], &fill_at);
                noncep[nfill] = nonces[nfill];
                dstp[nfill] = gss[t].buf + fill_at;
                fill_lane[nfill] = t;
                nfill++;
            }
        }

        /* (2) the queued lanes get their next block in ceil(nfill / lanes)
         * N-way passes.  The `nfill` lanes needing a fill occupy the
         * leading slots; the trailing slots of the LAST (partial) pass are
         * padded to DISTINCT throwaway buffers (never consumed) so each
         * issued pass is a whole N-way squeeze.  Passes that would serve
         * only slots >= nfill are skipped inside xof_nway_fill16 --
         * byte-neutral, since the caller never reads those slots.
         * Correctness depends only on the leading `nfill` slots, whose
         * bytes land in the real lane buffers and match the scalar
         * single-stream squeeze over the SAME nonce (lane-equivalence). */
        if (nfill) {
            const unsigned padto =
                ((nfill + (unsigned)XOF_LANES_AVX512 - 1) /
                 (unsigned)XOF_LANES_AVX512) *
                (unsigned)XOF_LANES_AVX512;
            PROF_START(t_sh);
            for (t = nfill; t < padto; t++) {
                noncep[t] =
                    noncep[0]; /* any valid nonce; output discarded */
                dstp[t] = scratch[t];
            }
            xof_nway_fill16(dstp, noncep, nlen, blk,
                            /*use_xof128=*/0, nfill);
            for (t = 0; t < nfill; t++)
                gss[fill_lane[t]].avail += blk;
            PROF_STOP(PT_G_SHAKE, t_sh);
        }

        /* (3) each not-done lane consumes ONE mini-batch from the buffer.
         */
        all_done = 1;
        for (t = 0; t < XOF_STREAMS; t++) {
            noise_lane_state *st = &lst[t];
            if (st->phase == 2)
                continue;
            if (st->phase == 0)
                noise_consume_minibatch(&gss[t], st->sdst, &st->scnt,
                                        st->ws, RCDT_NOISE_S,
                                        RCDT_NOISE_S_ENTRIES);
            else
                noise_consume_minibatch(&gss[t], st->edst, &st->ecnt,
                                        st->we, RCDT_NOISE_E,
                                        RCDT_NOISE_E_ENTRIES);
            noise_lane_advance_phase(st);
            if (st->phase != 2)
                all_done = 0;
        }
    }
}
#endif /* USE_AVX512_XOF_NWAY && __AVX512F__ */

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
    /* N-WAY ROUND-ROBIN (B): keep ALL of ExpandS's XOF on the N-way path.
     * Prime all 16 gauss_streams, batch their first block N-at-a-time
     * (xof256: 16-way SM3 = 1 pass / 8-way SHAKE = 2 passes), then drive
     * the 16 lanes in lockstep so EVERY subsequent refill -- not just
     * block 0 -- is computed N-way instead of one scalar squeeze per lane.
     * Byte-exact to ref (consumed-bytes invariance, see
     * expands_roundrobin). */
    {
        gauss_stream *gss = (gauss_stream *)malloc((size_t)XOF_STREAMS *
                                                   sizeof(gauss_stream));
        noise_lane_state *lst = (noise_lane_state *)malloc(
            (size_t)XOF_STREAMS * sizeof(noise_lane_state));
        if (gss && lst) {
            for (t = 0; t < XOF_STREAMS; t++) {
                gauss_stream_init(&gss[t], DS_EXPAND_S, seedsk, t);
                lst[t].sdst = sbar + (size_t)t * ws;
                lst[t].edst = ebar + (size_t)t * we;
                lst[t].ws = ws;
                lst[t].we = we;
                lst[t].scnt = 0;
                lst[t].ecnt = 0;
                lst[t].phase = 0;
            }
            gs_batch_first_fill(gss);
            expands_roundrobin(gss, lst);
            free(gss);
            free(lst);
            return;
        }
        /* malloc failure: fall through to the scalar per-lane path. */
        free(gss);
        free(lst);
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
    /* Single-stream rejection sampler over the public seedC.  Each refill
     * is a distinct XOF instance (tag||seedC||LE16(0)||LE16(refill)) that
     * logically supplies SAMPLEC_BLOCK bytes; we draw those bytes LAZILY
     * in SAMPLEC_DRAW-byte slices from the SAME persisted ctx and only
     * roll to the next instance (refill+1) once all SAMPLEC_BLOCK bytes of
     * the current one are consumed.  Squeezing N then M bytes from one
     * stream yields the same bytes as one N+M squeeze, so the consumed-
     * byte sequence and the refill boundary are byte-identical to the old
     * single 4096 B squeeze -- only the wasted tail is never squeezed.
     * Variable-time loop is fine here (seedC is public). */
    uint8_t nonce[1 + CHALLENGESEEDBYTES + 2 + 2];
    uint8_t block[SAMPLEC_DRAW];
    uint16_t refill = 0;
    /* cursor / valid bytes within block */
    size_t pos = 0, avail = 0;
    /* bytes squeezed from the current instance */
    size_t drawn = 0;
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
                size_t want;
                if (drawn >= SAMPLEC_BLOCK) {
                    /* current instance exhausted: advance to refill+1 */
                    refill++;
                    drawn = 0;
                }
                if (drawn == 0) {
                    put_le16(nonce + 1 + CHALLENGESEEDBYTES + 2, refill);
                    xof256_init(&ctx, nonce, sizeof(nonce));
                }
                want = SAMPLEC_BLOCK - drawn;
                if (want > SAMPLEC_DRAW)
                    want = SAMPLEC_DRAW;
                xof256_squeeze(&ctx, block, want);
                pos = 0;
                avail = want;
                drawn += want;
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
 * ===================================================================== *
 *
 *  Two byte-IDENTICAL code paths share the SAME per-buffer consumption
 *  kernels (gs_load_sign_stream / gs_consume_minibatch) so they cannot
 *  drift:
 *
 *   (A) scalar per-lane: gauss_stream_chunk() drives ONE lane to
 * completion, refilling its own buffer with the scalar single-stream
 * gs_ensure -> gs_fill (one xof256 squeeze per refill).  This is the
 * scalar reference path and the malloc-failure fallback.
 *
 *   (B) N-way round-robin (the headline lever): sample_y() drives ALL
 * 16 lanes in LOCKSTEP -- every round, the lanes that still need bytes are
 *       refilled N-AT-A-TIME (16-way SM3 = 1 pass / 8-way SHAKE = 2
 * passes) via xof_nway_fill16, then each not-yet-done lane consumes from
 * its freshly topped-up buffer.  Lanes finish at different rounds
 * (per-lane rejection makes R_k vary); finished lanes simply stop
 * consuming.
 *
 *  WHY (B) IS BYTE-EXACT TO (A) (consumed-bytes invariance):
 *  lane k's gauss_stream is the deterministic byte schedule
 *      xof256(0x08 || seedY || LE16(k) || LE16(refill)),  refill = 1,2,...
 *  consumed in a fixed sign-stream prefix + GAUSS_BATCH mini-batches with
 * a fully-consumed-up-front tail (the cursor advances regardless of the
 * early coefcnt==count break).  Lane k consumes refills 1..R_k.  Both
 * paths use the IDENTICAL nonce schedule (gs_build_nonce) and the
 * IDENTICAL per-buffer consumption kernels; only WHERE the bytes are
 * computed moves -- (A) makes them one scalar squeeze at a time, (B) makes
 * the not-done lanes' next blocks N-at-a-time.  The N-way SM3/SHAKE
 * kernels were proven byte-for-byte equal to the scalar xof256 over the
 * same nonce (lane-equivalence, test_xof_nway).  So the bytes lane k
 * CONSUMES -- and therefore every sampled coefficient and the KAT hash --
 * are unchanged.  Extra N-way bytes produced for already-finished lanes in
 * a round are simply UNUSED (harmless: the KAT depends only on the bytes
 * each lane CONSUMES). */

/* gs_load_sign_stream: consume the up-front OUTPUT-indexed sign stream
 * from the lane buffer into `signs` (signbytes = (count+7)/8).  PURE
 * buffer operation -- the caller MUST have ensured >= signbytes are
 * resident.  The AVX512 +SIGN_PAD_AVX512 over-read window stays zero
 * (caller zero-inits signs). */
static void gs_load_sign_stream(gauss_stream *gs, uint8_t *signs,
                                size_t signbytes)
{
    memcpy(signs, gs->buf + gs->pos, signbytes);
    gs->pos += signbytes;
}

/* gs_consume_minibatch: process ONE GAUSS_BATCH mini-batch from the lane
 * buffer, appending accepted coeffs to dst[*coefcnt..count).  PURE buffer
 * operation -- the caller MUST have ensured >= MINIBATCH_RAND_BYTES are
 * resident.  The whole MINIBATCH_RAND_BYTES is consumed up front (cursor
 * advances independent of the coefcnt==count early break).  This is
 * the EXACT scalar body extracted verbatim from gauss_stream_chunk,
 * shared by the scalar (A) and N-way (B) paths so they cannot diverge. */
static void gs_consume_minibatch(gauss_stream *gs, int32_t *dst,
                                 size_t *coefcnt, size_t count,
                                 const uint8_t *signs)
{
    int32_t x[GAUSS_BATCH];
    int32_t yv[GAUSS_BATCH];
    uint64_t phat[GAUSS_BATCH];
    const uint8_t *yp, *tailp;
    int j;

    {
        PROF_START(t_bs);
        sampler_sigma2(x, gs->buf + gs->pos);
        PROF_STOP(PT_G_BASESAMP, t_bs);
    }
    yp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES;
    for (j = 0; j < GAUSS_BATCH; j++)
        yv[j] = (int32_t)yp[j]; /* Y_BITS=8: plain byte copy */
    tailp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES + Y_RAND_BYTES;

    /* Two-path body.  When more than a full GAUSS_BATCH accepts are still
     * outstanding (count - *coefcnt > GAUSS_BATCH) every one of the 32
     * candidates is needed -- coefcnt cannot reach count this mini-batch
     * -- so we run the full-width SIMD path, byte-identical to before and
     * with the common-case throughput unchanged.  Otherwise this MIGHT be
     * the last mini-batch: we fuse approx_exp + finalize candidate-major
     * in SIMD groups of 8 (the AVX-512 gauss_finalize_batch width) and
     * BREAK the moment `count` accepts have landed, eliding the approx_exp
     * / finalize for the unused tail of the GAUSS_BATCH window.  Both
     * paths leave any candidate past `count` DISCARDED (the dst write was
     * already
     * `*coefcnt < count`-gated), and the byte cursor always advances the
     * WHOLE MINIBATCH_RAND_BYTES, so the per-lane byte schedule is byte-
     * identical to processing all 32.  The break is on the PUBLIC accept-
     * count (same data-dependence as the mini-batch count). */
    if (count - *coefcnt > (size_t)GAUSS_BATCH) {
        /* ----- common path: every candidate needed, full-width SIMD -----
         */
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
        {
            PROF_START(t_fin);
#if defined(USE_AVX512_SAMPLER) && defined(__AVX512F__) && \
    !defined(GAUSS_FINALIZE_SCALAR)
            int32_t cand[GAUSS_BATCH], negcand[GAUSS_BATCH];
            int32_t accept[GAUSS_BATCH], z0[GAUSS_BATCH];
            gauss_finalize_batch(cand, negcand, accept, z0, x, yv, phat,
                                 tailp, GAUSS_BATCH);
            for (j = 0; j < GAUSS_BATCH; j++) {
                size_t idx =
                    *coefcnt; /* OUTPUT index of the NEXT accept */
                uint32_t sgn =
                    (uint32_t)(signs[idx >> 3] >> (idx & 7)) & 1u;
                uint32_t keep =
                    (uint32_t)accept[j] & (1u ^ ((uint32_t)z0[j] & sgn));
                if (keep)
                    dst[(*coefcnt)++] = sgn ? negcand[j] : cand[j];
            }
#else
            for (j = 0; j < GAUSS_BATCH; j++) {
                int32_t r;
                uint32_t sgn;
                size_t idx =
                    *coefcnt; /* OUTPUT index of the NEXT accept */
                sgn = (uint32_t)(signs[idx >> 3] >> (idx & 7)) & 1u;
                if (gauss_finalize(&r, x[j], yv[j], phat[j],
                                   tailp + (size_t)j * GAUSS_RAND_BYTES,
                                   sgn))
                    dst[(*coefcnt)++] = r;
            }
#endif
            PROF_STOP(PT_G_FINAL, t_fin);
        }
    } else {
        /* ----- last-batch path: cover the candidate window [base,32) in
         * at most two WIDE SIMD passes, breaking the compaction once
         * `count` accepts land.  The first pass spans the candidates we
         * expect to need (remaining accepts rounded up to the 8-candidate
         * SIMD group); the second pass (taken only if a rare in-window
         * rejection left us short) covers the rest, so we NEVER
         * under-process and force an extra mini-batch.  Eliding the tail
         * of approx_exp / finalize past the window is what saves work;
         * keeping ONE wide finalize_batch per pass keeps the per-candidate
         * SIMD efficiency of the common path. Byte-exact: candidates are
         * processed in order with the same per- candidate math and the
         * same coefcnt progression as the full path. ----- */
        size_t remain = count - *coefcnt;     /* >0, <= GAUSS_BATCH */
        int win = (int)((remain + 7u) & ~7u); /* round up to SIMD group */
        int base = 0;
        if (win > GAUSS_BATCH)
            win = GAUSS_BATCH;
        while (base < GAUSS_BATCH && *coefcnt < count) {
            int n = win - base;
            int e = win;
            {
                PROF_START(t_ae);
                for (j = base; j < e; j += 4) {
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
            {
                PROF_START(t_fin);
#if defined(USE_AVX512_SAMPLER) && defined(__AVX512F__) && \
    !defined(GAUSS_FINALIZE_SCALAR)
                int32_t cand[GAUSS_BATCH], negcand[GAUSS_BATCH];
                int32_t accept[GAUSS_BATCH], z0[GAUSS_BATCH];
                gauss_finalize_batch(
                    cand + base, negcand + base, accept + base, z0 + base,
                    x + base, yv + base, phat + base,
                    tailp + (size_t)base * GAUSS_RAND_BYTES, n);
                for (j = base; j < e && *coefcnt < count; j++) {
                    size_t idx =
                        *coefcnt; /* OUTPUT index of NEXT accept */
                    uint32_t sgn =
                        (uint32_t)(signs[idx >> 3] >> (idx & 7)) & 1u;
                    uint32_t keep = (uint32_t)accept[j] &
                                    (1u ^ ((uint32_t)z0[j] & sgn));
                    if (keep)
                        dst[(*coefcnt)++] = sgn ? negcand[j] : cand[j];
                }
#else
                for (j = base; j < e && *coefcnt < count; j++) {
                    int32_t r;
                    uint32_t sgn;
                    size_t idx =
                        *coefcnt; /* OUTPUT index of NEXT accept */
                    sgn = (uint32_t)(signs[idx >> 3] >> (idx & 7)) & 1u;
                    if (gauss_finalize(
                            &r, x[j], yv[j], phat[j],
                            tailp + (size_t)j * GAUSS_RAND_BYTES, sgn))
                        dst[(*coefcnt)++] = r;
                }
#endif
                PROF_STOP(PT_G_FINAL, t_fin);
            }
            base = win;
            win =
                GAUSS_BATCH; /* second pass (if needed) covers the rest */
        }
    }
    gs->pos += MINIBATCH_RAND_BYTES; /* WHOLE tail consumed */
}

/* gauss_stream_chunk: scalar per-lane path (A) -- produce a flat run of
 * `count` wide-Gaussian samples, refilling this lane's buffer one scalar
 * xof256 squeeze at a time.  Byte-identical to the scalar reference;
 * kept for the non-USE_AVX512_XOF_NWAY build and the malloc-failure
 * fallback. */
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
    gs_load_sign_stream(gs, signs, signbytes);

    while (coefcnt < count) {
        {
            PROF_START(t_sh);
            gs_ensure(gs, MINIBATCH_RAND_BYTES);
            PROF_STOP(PT_G_SHAKE, t_sh);
        }
        gs_consume_minibatch(gs, dst, &coefcnt, count, signs);
    }
}

#if defined(USE_AVX512_XOF_NWAY) && defined(__AVX512F__)
/* ------------------------------------------------------------------- *
 *  N-way round-robin SampleY driver (B).                              *
 * ------------------------------------------------------------------- *
 *
 *  Per-lane resumable consumption state.  Each lane runs the SAME logical
 *  sequence as gauss_stream_chunk (load signs, then mini-batches until
 *  `count` accepted), but yields back to the driver whenever its buffer is
 *  too short for the next step, so the driver can batch that refill N-way
 *  across all lanes that need one.  `signs` is per-lane (OUTPUT-indexed
 * sign stream). */
typedef struct {
    int32_t *dst;     /* this lane's output slice              */
    size_t count;     /* coeffs to produce for this lane (wy)  */
    size_t coefcnt;   /* coeffs produced so far                */
    size_t signbytes; /* (count+7)/8                           */
    int signs_loaded; /* 0 until the up-front sign stream read */
    int done;         /* coefcnt == count                      */
    uint8_t signs[SIGN_BYTES_PER_CHUNK + SIGN_PAD_AVX512];
} gauss_lane_state;

/* gs_lane_bytes_needed: how many contiguous bytes lane `t` needs resident
 * in its buffer to make its NEXT step (sign-stream load, then each
 * mini-batch). This mirrors the gs_ensure(need) request in
 * gauss_stream_chunk EXACTLY, so the refill decision -- and therefore the
 * refill counter and consumed bytes
 * -- match the scalar path byte-for-byte. */
static size_t gs_lane_bytes_needed(const gauss_lane_state *st)
{
    return st->signs_loaded ? (size_t)MINIBATCH_RAND_BYTES : st->signbytes;
}

/* gs_lane_refill_prep: replicate the gs_ensure(need) bookkeeping for one
 * lane WITHOUT issuing the scalar single-stream squeeze.  Memmoves the
 * unconsumed leftover to the front of the buffer, bumps the refill
 * counter, builds the continuation nonce, and reports the
 * destination/length so the driver can fill
 * gs->buf[avail..avail+GAUSS_STREAM_BLOCK) N-way.  Byte-identical to the
 * scalar gs_ensure(need) path (same memmove, same refill++, same nonce),
 * only the squeeze itself is deferred to the N-way batch.  Returns the
 * nonce length; `*fill_at` receives the buffer offset to fill. */
static size_t gs_lane_refill_prep(gauss_stream *gs, uint8_t *nonce,
                                  size_t *fill_at)
{
    size_t left = gs->avail - gs->pos;
    if (left)
        memmove(gs->buf, gs->buf + gs->pos, left);
    gs->pos = 0;
    gs->avail = left;
    gs->refill++; /* matches scalar gs_ensure: refill++ THEN fill        */
    *fill_at = gs->avail;
    return gs_build_nonce(gs, nonce); /* continuation nonce (rc=refill)  */
}

/* gauss_roundrobin: drive all XOF_STREAMS lanes' SampleY in lockstep with
 * N-way batched refills.  `gss` are init'd (refill=0); `lst` carry the
 * per-lane dst/count.  Loops in rounds:
 *   (1) for each not-done lane whose buffer is short for its next step,
 *       prep its refill (memmove+refill++ +nonce) and queue an N-way fill;
 *   (2) one N-way pass (xof_nway_fill16) tops up every queued lane;
 *   (3) each not-done lane consumes its next step (sign load or one
 *       mini-batch) from the now-resident buffer.
 * Terminates when every lane is done.  Round count is bounded by the
 * PUBLIC per-lane work (sum of mini-batches) -- no secret-dependent loop
 * bound. */
static void gauss_roundrobin(gauss_stream gss[XOF_STREAMS],
                             gauss_lane_state lst[XOF_STREAMS])
{
    int all_done = 0;
    /* SampleY: every lane shares the tag -> the same tag-tuned per-refill
     * block size (right-sized, no over-squeezed tail). */
    size_t blk = gauss_block_bytes(gss[0].tag);
    while (!all_done) {
        uint8_t nonces[XOF_STREAMS][1 + CHALLENGESEEDBYTES + 2 + 2];
        const uint8_t *noncep[XOF_STREAMS];
        uint8_t *dstp[XOF_STREAMS];
        /* Per-round throwaway destinations for the inactive N-way slots,
         * so no two slots ever alias a live buffer (each padding slot gets
         * its OWN region; its bytes are discarded). */
        uint8_t scratch[XOF_STREAMS][GAUSS_STREAM_BLOCK];
        size_t nlen = 0;
        unsigned t, nfill = 0;
        unsigned fill_lane[XOF_STREAMS];

        /* (1) collect the lanes that need a refill for their next step. */
        for (t = 0; t < XOF_STREAMS; t++) {
            size_t need;
            if (lst[t].done)
                continue;
            need = gs_lane_bytes_needed(&lst[t]);
            if (gss[t].avail - gss[t].pos < need) {
                size_t fill_at;
                nlen =
                    gs_lane_refill_prep(&gss[t], nonces[nfill], &fill_at);
                noncep[nfill] = nonces[nfill];
                dstp[nfill] = gss[t].buf + fill_at;
                fill_lane[nfill] = t;
                nfill++;
            }
        }

        /* (2) the queued lanes get their next block in ceil(nfill / lanes)
         * N-way passes.  xof_nway_fill16 wants XOF_STREAMS-length
         * nonce/dst arrays indexed by physical N-way slot; the `nfill`
         * lanes needing a fill occupy the leading slots and the trailing
         * slots of the LAST (partial) pass are padded to DISTINCT
         * throwaway buffers (never consumed) so each issued pass is a
         * whole N-way squeeze with the SIMD fully utilised.  Passes that
         * would serve only slots >= nfill are skipped inside
         * xof_nway_fill16 -- byte-neutral, since the caller never reads
         * those slots.  Each padding slot reuses an active nonce (any
         * valid absorb) but its own scratch dst, so there is no aliasing
         * -- correctness depends only on the leading `nfill` slots, whose
         * bytes land in the real lane buffers and match the scalar
         * single-stream squeeze over the SAME nonce (lane-equivalence). */
        if (nfill) {
            const unsigned padto =
                ((nfill + (unsigned)XOF_LANES_AVX512 - 1) /
                 (unsigned)XOF_LANES_AVX512) *
                (unsigned)XOF_LANES_AVX512;
            PROF_START(t_sh);
            for (t = nfill; t < padto; t++) {
                noncep[t] =
                    noncep[0]; /* any valid nonce; output discarded */
                dstp[t] = scratch[t];
            }
            xof_nway_fill16(dstp, noncep, nlen, blk,
                            /*use_xof128=*/0, nfill);
            for (t = 0; t < nfill; t++)
                gss[fill_lane[t]].avail += blk;
            PROF_STOP(PT_G_SHAKE, t_sh);
        }

        /* (3) each not-done lane consumes its next step from the buffer.
         */
        all_done = 1;
        for (t = 0; t < XOF_STREAMS; t++) {
            gauss_lane_state *st = &lst[t];
            if (st->done)
                continue;
            if (!st->signs_loaded) {
                gs_load_sign_stream(&gss[t], st->signs, st->signbytes);
                st->signs_loaded = 1;
            } else {
                gs_consume_minibatch(&gss[t], st->dst, &st->coefcnt,
                                     st->count, st->signs);
            }
            if (st->coefcnt >= st->count)
                st->done = 1;
            else
                all_done = 0;
        }
    }
}
#endif /* USE_AVX512_XOF_NWAY && __AVX512F__ */

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
    /* N-WAY ROUND-ROBIN (B): keep ALL of SampleY's XOF on the N-way path.
     * Every refill across the 16 lanes is computed N-at-a-time (16-way SM3
     * = 1 pass / 8-way SHAKE = 2 passes), not one scalar squeeze per lane
     * -- so the BULK ~74 KB/sig of gauss XOF (PT_G_SHAKE, ~93% of SampleY)
     * runs at N-way throughput.  Byte-exact to ref (consumed-bytes
     * invariance, see the section comment). */
    {
        gauss_stream *gss = (gauss_stream *)malloc((size_t)XOF_STREAMS *
                                                   sizeof(gauss_stream));
        gauss_lane_state *lst = (gauss_lane_state *)malloc(
            (size_t)XOF_STREAMS * sizeof(gauss_lane_state));
        if (gss && lst) {
            for (t = 0; t < XOF_STREAMS; t++) {
                gauss_stream_init(&gss[t], DS_SAMPLE_Y, seedY, t);
                lst[t].dst = ybar + (size_t)t * wy;
                lst[t].count = wy;
                lst[t].coefcnt = 0;
                lst[t].signbytes = (wy + 7) / 8;
                lst[t].signs_loaded = 0;
                lst[t].done = 0;
                /* zero the AVX512 LE64 over-read window in the sign array
                 */
                memset(lst[t].signs, 0, sizeof(lst[t].signs));
            }
            gauss_roundrobin(gss, lst);
            free(gss);
            free(lst);
            return;
        }
        /* malloc failure: fall through to the scalar per-lane path. */
        free(gss);
        free(lst);
    }
#endif

    for (t = 0; t < XOF_STREAMS; t++) {
        gauss_stream gs;
        gauss_stream_init(&gs, DS_SAMPLE_Y, seedY, t);
        gauss_stream_chunk(&gs, ybar + (size_t)t * wy, wy);
    }
}
