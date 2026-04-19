/*
 * shuttle_rans.c - Byte-wise rANS encoder/decoder for SHUTTLE.
 *
 * Supports two frequency tables at compile time (SHUTTLE_MODE selects the
 * active mode; within a mode we expose both a HINT table and a Z1 table):
 *
 *   - hint table : distribution of the mod-2(q-1) signed hint h (Alg 2).
 *   - z1 table   : distribution of HighBits(z[1..L]) = round(z / alpha_h)
 *                   used by Phase 6c to compress z[1..L].
 *
 * Both tables share the same encoder/decoder core, differing only in the
 * CDF, symbol array, and symbol-range metadata. A ctx struct captures the
 * per-table data; lazy initialization builds both CDFs and both reverse
 * lookups on first use.
 *
 * State x is uint32_t; bytes emitted to TAIL of the output buffer during
 * encode (LIFO), then compacted to the front via memmove so the final
 * layout is a prefix of length *out_len. Decoder primes from the leading
 * 4 bytes (the last 4 written by the encoder).
 */

#include <string.h>

#include "params.h"
#include "shuttle_rans.h"
#include "rans_tables.h"

/* ============================================================
 * Per-mode table bindings. We bind both hint and z1 tables.
 * ============================================================ */

#if SHUTTLE_MODE == 128
#  define SRANS_HINT_PROB_BITS  SHUTTLE128_RANS_PROB_BITS
#  define SRANS_HINT_NUM_SYMS   SHUTTLE128_RANS_NUM_SYMS
#  define SRANS_HINT_SYM_MIN    SHUTTLE128_RANS_SYM_MIN
#  define SRANS_HINT_SYM_MAX    SHUTTLE128_RANS_SYM_MAX
#  define srans_hint_syms       shuttle128_rans_syms
#  define srans_hint_freqs      shuttle128_rans_freqs

#  define SRANS_Z1_PROB_BITS    SHUTTLE128_RANS_Z1_PROB_BITS
#  define SRANS_Z1_NUM_SYMS     SHUTTLE128_RANS_Z1_NUM_SYMS
#  define SRANS_Z1_SYM_MIN      SHUTTLE128_RANS_Z1_SYM_MIN
#  define SRANS_Z1_SYM_MAX      SHUTTLE128_RANS_Z1_SYM_MAX
#  define srans_z1_syms         shuttle128_rans_z1_syms
#  define srans_z1_freqs        shuttle128_rans_z1_freqs
#elif SHUTTLE_MODE == 256
#  define SRANS_HINT_PROB_BITS  SHUTTLE256_RANS_PROB_BITS
#  define SRANS_HINT_NUM_SYMS   SHUTTLE256_RANS_NUM_SYMS
#  define SRANS_HINT_SYM_MIN    SHUTTLE256_RANS_SYM_MIN
#  define SRANS_HINT_SYM_MAX    SHUTTLE256_RANS_SYM_MAX
#  define srans_hint_syms       shuttle256_rans_syms
#  define srans_hint_freqs      shuttle256_rans_freqs

#  define SRANS_Z1_PROB_BITS    SHUTTLE256_RANS_Z1_PROB_BITS
#  define SRANS_Z1_NUM_SYMS     SHUTTLE256_RANS_Z1_NUM_SYMS
#  define SRANS_Z1_SYM_MIN      SHUTTLE256_RANS_Z1_SYM_MIN
#  define SRANS_Z1_SYM_MAX      SHUTTLE256_RANS_Z1_SYM_MAX
#  define srans_z1_syms         shuttle256_rans_z1_syms
#  define srans_z1_freqs        shuttle256_rans_z1_freqs
#else
#  error "Unsupported SHUTTLE_MODE for rANS tables"
#endif

/* Both tables use the same prob_bits = 12 by construction. We still keep
 * per-table macros so a future asymmetric layout works without edits. */
#define SRANS_HINT_PROB_TOTAL   (1u << SRANS_HINT_PROB_BITS)
#define SRANS_Z1_PROB_TOTAL     (1u << SRANS_Z1_PROB_BITS)

/* ============================================================
 * Per-table context + lazy init.
 * ============================================================ */

typedef struct {
    const int16_t  *syms;       /* length num_syms, contiguous SYM_MIN..SYM_MAX */
    const uint16_t *freqs;      /* length num_syms, sum == prob_total */
    const uint16_t *cdf;        /* length num_syms+1, lazy-built */
    const uint16_t *sym_lookup; /* length prob_total, lazy-built */
    int16_t         sym_min;
    int16_t         sym_max;
    uint16_t        num_syms;
    uint8_t         prob_bits;
    uint8_t         initialized;
    uint32_t        prob_total;
} rans_ctx_t;

/* Hint-table derived storage. */
static uint16_t g_hint_cdf[SRANS_HINT_NUM_SYMS + 1];
static uint16_t g_hint_sym_lookup[SRANS_HINT_PROB_TOTAL];
static rans_ctx_t g_ctx_hint = {
    .syms       = srans_hint_syms,
    .freqs      = srans_hint_freqs,
    .cdf        = g_hint_cdf,
    .sym_lookup = g_hint_sym_lookup,
    .sym_min    = SRANS_HINT_SYM_MIN,
    .sym_max    = SRANS_HINT_SYM_MAX,
    .num_syms   = SRANS_HINT_NUM_SYMS,
    .prob_bits  = SRANS_HINT_PROB_BITS,
    .initialized = 0,
    .prob_total = SRANS_HINT_PROB_TOTAL,
};

/* z1-table derived storage. */
static uint16_t g_z1_cdf[SRANS_Z1_NUM_SYMS + 1];
static uint16_t g_z1_sym_lookup[SRANS_Z1_PROB_TOTAL];
static rans_ctx_t g_ctx_z1 = {
    .syms       = srans_z1_syms,
    .freqs      = srans_z1_freqs,
    .cdf        = g_z1_cdf,
    .sym_lookup = g_z1_sym_lookup,
    .sym_min    = SRANS_Z1_SYM_MIN,
    .sym_max    = SRANS_Z1_SYM_MAX,
    .num_syms   = SRANS_Z1_NUM_SYMS,
    .prob_bits  = SRANS_Z1_PROB_BITS,
    .initialized = 0,
    .prob_total = SRANS_Z1_PROB_TOTAL,
};

static void rans_ctx_init(rans_ctx_t *ctx) {
    uint16_t *cdf;
    uint16_t *lookup;
    uint32_t running;
    unsigned i, k;

    if (ctx->initialized) return;

    /* Safe cast: the storage slots are writable, even though the public
     * struct view exposes them as const. */
    cdf    = (uint16_t *)(uintptr_t)ctx->cdf;
    lookup = (uint16_t *)(uintptr_t)ctx->sym_lookup;

    running = 0;
    cdf[0] = 0;
    for (i = 0; i < ctx->num_syms; ++i) {
        running += ctx->freqs[i];
        cdf[i + 1] = (uint16_t)running;
    }
    /* running == ctx->prob_total by construction of the table. */

    for (i = 0; i < ctx->num_syms; ++i) {
        uint32_t start = cdf[i];
        uint32_t end   = cdf[i + 1];
        for (k = start; k < end; ++k)
            lookup[k] = (uint16_t)i;
    }

    ctx->initialized = 1;
}

/* Map a signed integer symbol to its index in ctx->syms.
 * Returns -1 if the symbol is outside [SYM_MIN, SYM_MAX]. */
static int sym_to_index(const rans_ctx_t *ctx, int32_t sym) {
    int32_t idx = sym - ctx->sym_min;
    if (idx < 0 || idx >= (int32_t)ctx->num_syms)
        return -1;
    if (ctx->syms[idx] != sym)
        return -1;
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
        int32_t sym = syms[idx_i - 1];
        int si = sym_to_index(ctx, sym);
        if (si < 0) return -1;

        uint32_t freq  = ctx->freqs[si];
        uint32_t start = ctx->cdf[si];

        uint32_t x_max = ((SHUTTLE_RANS_L >> ctx->prob_bits)
                          << SHUTTLE_RANS_BYTE_BITS) * freq;
        while (x >= x_max) {
            if (write_pos == 0) return -2;
            out[--write_pos] = (uint8_t)(x & 0xFF);
            x >>= SHUTTLE_RANS_BYTE_BITS;
        }

        x = ((x / freq) << ctx->prob_bits) + start + (x % freq);
    }

    if (write_pos < 4) return -2;
    write_pos -= 4;
    out[write_pos + 0] = (uint8_t)(x >> 24);
    out[write_pos + 1] = (uint8_t)(x >> 16);
    out[write_pos + 2] = (uint8_t)(x >>  8);
    out[write_pos + 3] = (uint8_t)(x >>  0);

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

    uint32_t x = ((uint32_t)in[0] << 24)
               | ((uint32_t)in[1] << 16)
               | ((uint32_t)in[2] <<  8)
               | ((uint32_t)in[3] <<  0);
    size_t pos = 4;

    const uint32_t prob_mask = ctx->prob_total - 1;

    for (size_t k = 0; k < n; ++k) {
        uint32_t c = x & prob_mask;
        unsigned si = ctx->sym_lookup[c];

        uint32_t freq  = ctx->freqs[si];
        uint32_t start = ctx->cdf[si];

        syms[k] = ctx->syms[si];

        x = freq * (x >> ctx->prob_bits) + (c - start);

        while (x < SHUTTLE_RANS_L) {
            if (pos >= in_len) return -1;
            x = (x << SHUTTLE_RANS_BYTE_BITS) | in[pos];
            ++pos;
        }
    }

    return 0;
}

/* ============================================================
 * Public API thin wrappers
 * ============================================================ */

int shuttle_rans_encode(uint8_t *out, size_t *out_len, size_t max_bytes,
                        const int32_t *syms, size_t n)
{
    return encode_core(&g_ctx_hint, out, out_len, max_bytes, syms, n);
}

int shuttle_rans_decode(int32_t *syms, size_t n,
                        const uint8_t *in, size_t in_len)
{
    return decode_core(&g_ctx_hint, syms, n, in, in_len);
}

int shuttle_rans_encode_z1(uint8_t *out, size_t *out_len, size_t max_bytes,
                           const int32_t *syms, size_t n)
{
    return encode_core(&g_ctx_z1, out, out_len, max_bytes, syms, n);
}

int shuttle_rans_decode_z1(int32_t *syms, size_t n,
                           const uint8_t *in, size_t in_len)
{
    return decode_core(&g_ctx_z1, syms, n, in, in_len);
}
