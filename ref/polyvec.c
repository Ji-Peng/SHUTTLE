/*
 * polyvec.c -- vector-level XOF-driven sampling surface (PURE SCALAR
 *              REFERENCE; the KAT oracle).  P07.
 *
 * Implements ExpandSeeds / ExpandSigningSeeds / ExpandA / ExpandS /
 * SampleC / SampleDGauss / SampleY.  The AVX2/AVX512 forks (M6) must be
 * byte-identical to this scalar oracle within each MODE.
 *
 * Read polyvec.h FIRST for the binding NGCC-DRBG no-rate-cursor rule and
 * the 16-stream nonce layout: every logical stream draws its whole byte
 * budget in ONE squeeze; refills re-init a fresh ctx with a 2-byte refill
 * counter appended to the nonce.  This is what keeps ref==avx2==avx512
 * byte-exact (K6) and what makes the cursor advance independent of any
 * early rejection break (K6/K8).
 *
 * ===================================================================== *
 *  PINNED PRNG byte schedules (KAT-NORMATIVE)                            *
 * ===================================================================== *
 * ExpandSeeds        : absorb 0x00 || xi || LE16(EM) || LE32(kappa);
 *                      squeeze 5*LAMBDA/8 bytes. * ExpandSigningSeeds :
 * absorb 0x01 || K || rnd || mu || LE32(kappa);     * squeeze SEEDBYTES
 * bytes.                           * ExpandA (per lane) : nonce 0x02 ||
 * seedA || LE16(t) || LE16(rc);       * per coeff REPEAT{ 2 bytes; a = LE
 * & (2^DQ-1) }     * UNTIL a < q.  agen chunk THEN hAgen chunk.         *
 * ExpandS (per lane) : nonce 0x03 || seedsk || LE16(t) || LE16(rc);      *
 *                      per coeff via a NOISE-style mini-batch: 384-byte *
 *                      grouped rho_u (cdt_scan96) + 8-byte 2-bit tail *
 *                      (bit0 sign, bit1 zero-fold), whole tail up front. *
 *                      s chunk (RCDT_NOISE_S) THEN e chunk
 * (RCDT_NOISE_E).* SampleC            : nonce 0x07 || seedC || LE16(0) ||
 * LE16(rc);       * per draw REPEAT{ BN bytes; j = LE & (2^DN-1) }     *
 *                      UNTIL j <= i. * SampleY (per lane) : nonce 0x08 ||
 * seedY || LE16(t) || LE16(rc);       * per poly: SIGN_BYTES_PER_POLY
 * up-front sign bits,  * then mini-batches of GAUSS_BATCH candidates: *
 *                      384-byte rho_u (cdt_scan96 RCDT_Z) + 32-byte y + *
 *                      256-byte (32x8) Bernoulli tails, whole tail up *
 *                      front; accept via gauss_finalize. *
 */
#include "polyvec.h"

#include <string.h>

#include "approx_exp.h" /* approx_exp_accept_q64_x4, approx_exp_accept_q64 */
#include "rcdt_tables.h" /* SHUTTLE_RCDT_Z, SHUTTLE_RCDT_NOISE_* */

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
/* One block must comfortably exceed the expected per-coeff cost.  Each
 * coeff needs ~q/2^DQ acceptance => on average <2 candidates of BQ=2
 * bytes. The largest per-lane chunk is hAgen's w_h = EM*ELL*n/16 =
 * 144/192/384 coeffs.  4096 bytes / 2 = 2048 candidates covers >>5x the
 * worst chunk; refills are essentially never exercised but are implemented
 * + tested. */
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

/* Sample `count` uniform coeffs in [0,q) (ascending index) into `dst`. */
static void uniform_reject_chunk(uniform_stream *us, uint16_t *dst,
                                 size_t count)
{
    size_t k;
    for (k = 0; k < count; k++) {
        uint16_t a;
        do {
            a = us_next_candidate(us); /* REPEAT { 2B; mask } UNTIL a<q */
        } while (a >= (uint16_t)Q);
        dst[k] = a;
    }
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

    /* Each lane fills its agen chunk THEN its hAgen chunk, from the SAME
     * per-lane stream (the contiguous flattening per the spec).  hAgen is
     * written ascending-index = canonical NTT order (ref needs no
     * nttunpack; the AVX forks poly_ntt_import the canonical poly). */
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
/* gauss_stream services ExpandS (noise mini-batches) and SampleY (wide
 * mini-batches).  Both consume their per-lane stream as a sequence of
 * fixed-size mini-batches whose whole tail is drawn up front (K6/K8).  The
 * struct is declared in polyvec.h.  GAUSS_STREAM_BLOCK is sized to hold a
 * SampleY poly's sign stream + one wide mini-batch, the larger of the two
 * uses; ExpandS mini-batches (NOISE_MINIBATCH_RAND_BYTES = 392) fit too.
 */

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

/* Build the per-(lane,refill) nonce and draw GAUSS_STREAM_BLOCK bytes in
 * ONE xof256 squeeze, appending after any leftover already memmoved to the
 * front of buf.  seed length differs by tag: SampleY uses SEEDBYTES,
 * ExpandS uses CHALLENGESEEDBYTES.  We carry the length via the tag. */
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
    /* memmove leftover to the front, bump refill, draw a fresh block. */
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
/* One noise mini-batch: cdt_scan96 over `Z`/`entries` (NOISE_BATCH=32),
 * then a 2-bit-per-candidate tail (bit0 sign, bit1 zero-fold), 4
 * cand/byte; the WHOLE NOISE_MINIBATCH_RAND_BYTES (392) is consumed up
 * front so the cursor advances independently of the cnt==want early break
 * (K6/K8). Appends accepted signed coeffs to dst[*cnt..]; stops at `want`.
 */
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
        /* f = (tailp[j>>2] >> (2*(j&3))) & 3 : bit0 = sign, bit1 = fold */
        uint32_t f = (uint32_t)(tailp[j >> 2] >> (2 * (j & 3))) & 3u;
        uint32_t reject = ct_is_zero_u32((uint32_t)mag[j]) & (f >> 1);
        int32_t r = ct_sel_i32(f & 1u, -mag[j], mag[j]);
        if (*cnt < want && !reject) /* placement only; cursor already
                                       advanced */
            dst[(*cnt)++] = r;
    }
    gs->pos +=
        NOISE_MINIBATCH_RAND_BYTES; /* WHOLE tail consumed (K6/K8) */
}

void expand_s(poly s1s2[ELL + EM],
              const uint8_t seedsk[CHALLENGESEEDBYTES])
{
    /* Flat views over the s (ELL*n) and e (EM*n) blocks of s1s2. */
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
        /* s chunk uses RCDT_NOISE_S; e chunk uses RCDT_NOISE_E.  Both ride
         * the SAME per-lane stream (s first, then e), so the byte cursor
         * threads s->e exactly like the spec's flattened lane chunk. */
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
    /* Single-stream rejection sampler over the public seedC.  Re-uses the
     * uniform-stream one-squeeze-per-fill discipline but with seedC's
     * CHALLENGESEEDBYTES length and the SampleC tag.  Variable-time loop
     * is fine here (seedC is public). */
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

    /* for i = n-tau .. n-1: REPEAT{ BN bytes; j = LE & (2^DN-1) } UNTIL
     * j<=i; then c[i]=c[j]; c[j]=1. */
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
/* gauss_stream_chunk: produce a flat run of `count` wide-Gaussian samples
 * (the lane chunk KVEC*n/16).
 *   1. up-front OUTPUT-indexed sign stream: (count+7)/8 bytes (K9) -- ONE
 *      bit per output, indexed by the running accept counter.
 *   2. mini-batches of GAUSS_BATCH candidates until `count` accepted:
 *        - cdt_scan96(x[32], buf+pos) over RCDT_Z (advance 384)
 *        - y[32] = buf[pos..pos+32]   (advance 32; Y_BITS=8 byte copy)
 *        - p_hat[32] via approx_exp_accept_q64_x4 (8 groups of 4)
 *        - per-candidate gauss_finalize, tail at buf+pos+8*j (advance 256)
 *      The whole MINIBATCH_RAND_BYTES is consumed up front (K6/K8). */
void gauss_stream_chunk(gauss_stream *gs, int32_t *dst, size_t count)
{
    uint8_t signs[SIGN_BYTES_PER_CHUNK + SIGN_PAD_AVX512];
    size_t signbytes = (count + 7) / 8;
    size_t coefcnt = 0;

    /* (1) up-front sign bits (one per OUTPUT sample in this chunk). */
    memset(signs, 0,
           sizeof(signs)); /* pad zero for the AVX512 LE64 window */
    gs_ensure(gs, signbytes);
    memcpy(signs, gs->buf + gs->pos, signbytes);
    gs->pos += signbytes;

    /* (2) wide mini-batches. */
    while (coefcnt < count) {
        int32_t x[GAUSS_BATCH];
        int32_t yv[GAUSS_BATCH];
        uint64_t phat[GAUSS_BATCH];
        const uint8_t *yp, *tailp;
        int j;

        gs_ensure(gs, MINIBATCH_RAND_BYTES);
        /* magnitude x via RCDT_Z (grouped 96B/8-sample layout) */
        sampler_sigma2(x, gs->buf + gs->pos);
        yp = gs->buf + gs->pos + SIGMA_S_RAND_BYTES;
        for (j = 0; j < GAUSS_BATCH; j++)
            yv[j] = (int32_t)yp[j]; /* Y_BITS=8: plain byte copy */
        /* p_hat = approx_exp accept threshold, 8 groups of 4 (x4 batch).
         */
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
    /* y flattened to ybar[0..KVEC*n-1]; chunk width w_y = KVEC*n/16.  Lane
     * t fills the FLAT chunk ybar[t*w_y .. (t+1)*w_y).  The flattening is
     * by COEFFICIENT index (the spec's "for each flat index k in chunk t:
     * ybar[k] <- SampleDGauss"), so a lane chunk is a contiguous run that
     * may cross poly boundaries -- which is fine, the output poly array is
     * one contiguous int32 buffer. */
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
