/*
 * shuttle_rans.c - Byte-wise rANS encoder/decoder for SHUTTLE.
 *
 * Two contexts, two public API pairs (see shuttle_rans.h for the full
 * design notes; agent/rANS/SHUTTLE_rANS.tex §2.2-2.3 for the algebra):
 *
 *   z-hi  : HighBits_{alpha_0'}(z^(0)) ⨁ HighBits_{alpha_r}(z^(1..lenS))
 *   hint  : MakeHint output
 *
 * Two distinct frequency tables per mode (ZHI and HINT) — the earlier
 * "mode-128 shares one unified table" optimization was removed once we
 * discovered that the hint distribution is NOT a discrete Gaussian
 * (SHUTTLE_rANS.tex §4.5 errata). Even though mode-128's nominal
 * sigma_h coincides with sigma_zhi, the actual hint PMF (bucket-crossing,
 * wider than D_{Z, sigma_h}) differs in shape from the z-hi PMF, so a
 * shared table would be sub-optimal.
 *
 * State is uint32_t; bytes are emitted to the TAIL of the output buffer
 * during encode (LIFO), then compacted with memmove so the final layout
 * is a prefix of length *out_len. The 4-byte flush at the end is written
 * little-endian, matching ryg_rans / HAETAE. Decoding reads the leading
 * 4 little-endian bytes to prime the state, then walks the stream.
 *
 * Final-state check on decode: after the last symbol is consumed, x must
 * equal L_ren (= 2^23). Any byte-level corruption that did not blow up
 * the renorm loop will, with overwhelming probability, land here.
 */

#include <assert.h>
#include <string.h>

#include "params.h"
#include "shuttle_rans.h"
#include "rans_tables.h"

/* ============================================================
 * Per-mode table bindings.
 *
 * All three modes have separate ZHI_ and HINT_ tables in rans_tables.h.
 * ============================================================ */

#if SHUTTLE_MODE == 128
#  define SRANS_ZHI_NUM_SYMS    SHUTTLE128_RANS_ZHI_NUM_SYMS
#  define SRANS_ZHI_SYM_MIN     SHUTTLE128_RANS_ZHI_SYM_MIN
#  define SRANS_ZHI_SYM_MAX     SHUTTLE128_RANS_ZHI_SYM_MAX
#  define srans_zhi_syms        shuttle128_rans_zhi_syms
#  define srans_zhi_freqs       shuttle128_rans_zhi_freqs

#  define SRANS_HINT_NUM_SYMS   SHUTTLE128_RANS_HINT_NUM_SYMS
#  define SRANS_HINT_SYM_MIN    SHUTTLE128_RANS_HINT_SYM_MIN
#  define SRANS_HINT_SYM_MAX    SHUTTLE128_RANS_HINT_SYM_MAX
#  define srans_hint_syms       shuttle128_rans_hint_syms
#  define srans_hint_freqs      shuttle128_rans_hint_freqs
#elif SHUTTLE_MODE == 256
#  define SRANS_ZHI_NUM_SYMS    SHUTTLE256_RANS_ZHI_NUM_SYMS
#  define SRANS_ZHI_SYM_MIN     SHUTTLE256_RANS_ZHI_SYM_MIN
#  define SRANS_ZHI_SYM_MAX     SHUTTLE256_RANS_ZHI_SYM_MAX
#  define srans_zhi_syms        shuttle256_rans_zhi_syms
#  define srans_zhi_freqs       shuttle256_rans_zhi_freqs

#  define SRANS_HINT_NUM_SYMS   SHUTTLE256_RANS_HINT_NUM_SYMS
#  define SRANS_HINT_SYM_MIN    SHUTTLE256_RANS_HINT_SYM_MIN
#  define SRANS_HINT_SYM_MAX    SHUTTLE256_RANS_HINT_SYM_MAX
#  define srans_hint_syms       shuttle256_rans_hint_syms
#  define srans_hint_freqs      shuttle256_rans_hint_freqs
#elif SHUTTLE_MODE == 512
#  define SRANS_ZHI_NUM_SYMS    SHUTTLE512_RANS_ZHI_NUM_SYMS
#  define SRANS_ZHI_SYM_MIN     SHUTTLE512_RANS_ZHI_SYM_MIN
#  define SRANS_ZHI_SYM_MAX     SHUTTLE512_RANS_ZHI_SYM_MAX
#  define srans_zhi_syms        shuttle512_rans_zhi_syms
#  define srans_zhi_freqs       shuttle512_rans_zhi_freqs

#  define SRANS_HINT_NUM_SYMS   SHUTTLE512_RANS_HINT_NUM_SYMS
#  define SRANS_HINT_SYM_MIN    SHUTTLE512_RANS_HINT_SYM_MIN
#  define SRANS_HINT_SYM_MAX    SHUTTLE512_RANS_HINT_SYM_MAX
#  define srans_hint_syms       shuttle512_rans_hint_syms
#  define srans_hint_freqs      shuttle512_rans_hint_freqs
#else
#  error "Unsupported SHUTTLE_MODE for rANS tables"
#endif

/* ============================================================
 * Per-table context + lazy init.
 *
 * For mode-128 both contexts point to the same underlying tables, but
 * each owns a *separate* cdf[] / sym_lookup[] scratch slot. That is
 * harmless (they would derive identical contents anyway) and keeps the
 * code path identical to mode-256/512.
 * ============================================================ */

typedef struct {
    const int16_t  *syms;       /* length num_syms, contiguous SYM_MIN..SYM_MAX */
    const uint16_t *freqs;      /* length num_syms, sum == 2^PROB_BITS */
    uint16_t       *cdf;        /* length num_syms+1, lazy-built */
    uint16_t       *sym_lookup; /* length 2^PROB_BITS, lazy-built */
    int16_t         sym_min;
    int16_t         sym_max;
    uint16_t        num_syms;
    uint8_t         initialized;
} rans_ctx_t;

/* z-hi scratch storage. */
static uint16_t g_zhi_cdf[SRANS_ZHI_NUM_SYMS + 1];
static uint16_t g_zhi_sym_lookup[SHUTTLE_RANS_PROB_TOTAL];
static rans_ctx_t g_ctx_zhi = {
    .syms       = srans_zhi_syms,
    .freqs      = srans_zhi_freqs,
    .cdf        = g_zhi_cdf,
    .sym_lookup = g_zhi_sym_lookup,
    .sym_min    = SRANS_ZHI_SYM_MIN,
    .sym_max    = SRANS_ZHI_SYM_MAX,
    .num_syms   = SRANS_ZHI_NUM_SYMS,
    .initialized = 0,
};

/* hint scratch storage. */
static uint16_t g_hint_cdf[SRANS_HINT_NUM_SYMS + 1];
static uint16_t g_hint_sym_lookup[SHUTTLE_RANS_PROB_TOTAL];
static rans_ctx_t g_ctx_hint = {
    .syms       = srans_hint_syms,
    .freqs      = srans_hint_freqs,
    .cdf        = g_hint_cdf,
    .sym_lookup = g_hint_sym_lookup,
    .sym_min    = SRANS_HINT_SYM_MIN,
    .sym_max    = SRANS_HINT_SYM_MAX,
    .num_syms   = SRANS_HINT_NUM_SYMS,
    .initialized = 0,
};

static void rans_ctx_init(rans_ctx_t *ctx) {
    uint32_t running;
    unsigned i, k;

    if (ctx->initialized) return;

    running = 0;
    ctx->cdf[0] = 0;
    for (i = 0; i < ctx->num_syms; ++i) {
        running += ctx->freqs[i];
        ctx->cdf[i + 1] = (uint16_t)running;
    }
    /* running == 2^PROB_BITS by construction of the static table. */

    for (i = 0; i < ctx->num_syms; ++i) {
        uint32_t start = ctx->cdf[i];
        uint32_t end   = ctx->cdf[i + 1];
        for (k = start; k < end; ++k)
            ctx->sym_lookup[k] = (uint16_t)i;
    }

    ctx->initialized = 1;
}

/* Map a signed integer symbol to its index in ctx->syms.
 *
 * The theoretical vocabulary covers [SYM_MIN, SYM_MAX] for every legal
 * signature (SHUTTLE_rANS.tex §3 tight bound from the 11*sigma SampleY
 * truncation), so OOV cannot occur for a faithful signer. We still
 * assert it in debug builds so a generation-side bug (e.g. parameter
 * drift or stale rans_tables.h) shows up immediately. */
static inline int sym_to_index(const rans_ctx_t *ctx, int32_t sym) {
    int32_t idx = sym - ctx->sym_min;
    assert(idx >= 0 && idx < (int32_t)ctx->num_syms);
    assert(ctx->syms[idx] == sym);
    return (int)idx;
}

/* ============================================================
 * Core encoder / decoder (table-parameterized)
 * ============================================================ */

static int encode_core(rans_ctx_t *ctx,
                       uint8_t *out, size_t *out_len, size_t max_bytes,
                       const int32_t *syms, size_t n)
{
    uint32_t x = SHUTTLE_RANS_L;
    size_t idx_i;
    size_t write_pos;

    rans_ctx_init(ctx);

    write_pos = max_bytes;

    for (idx_i = n; idx_i > 0; --idx_i) {
        int si = sym_to_index(ctx, syms[idx_i - 1]);

        uint32_t freq  = ctx->freqs[si];
        uint32_t start = ctx->cdf[si];

        uint32_t x_max = ((SHUTTLE_RANS_L >> SHUTTLE_RANS_PROB_BITS)
                          << SHUTTLE_RANS_BYTE_BITS) * freq;
        while (x >= x_max) {
            if (write_pos == 0) return -2;
            out[--write_pos] = (uint8_t)(x & 0xFFu);
            x >>= SHUTTLE_RANS_BYTE_BITS;
        }

        x = ((x / freq) << SHUTTLE_RANS_PROB_BITS) + start + (x % freq);
    }

    /* Final flush: write x as 4 little-endian bytes. */
    if (write_pos < 4) return -2;
    write_pos -= 4;
    out[write_pos + 0] = (uint8_t)(x >>  0);
    out[write_pos + 1] = (uint8_t)(x >>  8);
    out[write_pos + 2] = (uint8_t)(x >> 16);
    out[write_pos + 3] = (uint8_t)(x >> 24);

    {
        size_t bytes_used = max_bytes - write_pos;
        memmove(out, out + write_pos, bytes_used);
        *out_len = bytes_used;
    }
    return 0;
}

static int decode_core(rans_ctx_t *ctx,
                       int32_t *syms, size_t n,
                       const uint8_t *in, size_t in_len)
{
    rans_ctx_init(ctx);

    if (in_len < 4) return -1;

    /* Final flush was little-endian; read it back in matching order. */
    uint32_t x = ((uint32_t)in[0] <<  0)
               | ((uint32_t)in[1] <<  8)
               | ((uint32_t)in[2] << 16)
               | ((uint32_t)in[3] << 24);
    size_t pos = 4;

    const uint32_t prob_mask = SHUTTLE_RANS_PROB_TOTAL - 1u;

    for (size_t k = 0; k < n; ++k) {
        uint32_t c = x & prob_mask;
        unsigned si = ctx->sym_lookup[c];

        uint32_t freq  = ctx->freqs[si];
        uint32_t start = ctx->cdf[si];

        syms[k] = ctx->syms[si];

        x = freq * (x >> SHUTTLE_RANS_PROB_BITS) + (c - start);

        while (x < SHUTTLE_RANS_L) {
            if (pos >= in_len) return -1;
            x = (x << SHUTTLE_RANS_BYTE_BITS) | in[pos];
            ++pos;
        }
    }

    /* Final-state verification: a faithful round-trip leaves x == L_ren.
     * Catches truncation, byte-level corruption and length-prefix
     * tampering that survived the renorm loop. */
    if (x != SHUTTLE_RANS_L) return -1;

    return 0;
}

/* ============================================================
 * Public API thin wrappers
 * ============================================================ */

int shuttle_rans_encode_zhi(uint8_t *out, size_t *out_len, size_t max_bytes,
                            const int32_t *syms, size_t n)
{
    return encode_core(&g_ctx_zhi, out, out_len, max_bytes, syms, n);
}

int shuttle_rans_decode_zhi(int32_t *syms, size_t n,
                            const uint8_t *in, size_t in_len)
{
    return decode_core(&g_ctx_zhi, syms, n, in, in_len);
}

int shuttle_rans_encode_hint(uint8_t *out, size_t *out_len, size_t max_bytes,
                             const int32_t *syms, size_t n)
{
    return encode_core(&g_ctx_hint, out, out_len, max_bytes, syms, n);
}

int shuttle_rans_decode_hint(int32_t *syms, size_t n,
                             const uint8_t *in, size_t in_len)
{
    return decode_core(&g_ctx_hint, syms, n, in, in_len);
}
