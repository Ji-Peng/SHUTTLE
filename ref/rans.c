/*
 * rans.c -- static byte-renormalized rANS codec for SHUTTLE.
 *
 * Byte-exact to the Python golden model tools/rans.py.  32-bit state,
 * L = 2^23, 8-bit renorm, prob_bits = 10, RANS_N = 2 interleaved streams.
 * The flat symbol sequence is Q0 then Qs then hint (logical order); symbol
 * t is coded on its power-of-two interleave state s = t & MASK.  Encode
 * pushes in REVERSE t (so decode pulls forward t), all states sharing one
 * backward-written byte stream.  N=2 hides the ~5-cyc reciprocal-multiply
 * latency at the cost of only 4 extra flush bytes (negligible vs
 * RANS_RESERVED_BYTES); see rans.h for the engine rationale.
 *
 * CONSTANT-TIME: this operates on the PUBLIC, post-signing signature
 * data, so the renorm `while` loops and `pos` bounds are not a CT
 * violation.
 */
#include "rans.h"

#include <stdint.h>
#include <string.h>

#define PB RANS_PROB_BITS /* 10 */
#define PSCALE (1u << PB)
#define RN RANS_INTERLEAVED_STREAMS
#define RMASK RANS_INTERLEAVED_MASK

static inline unsigned rans_stream_index(size_t t)
{
    return (unsigned)(t & RMASK);
}

static inline int rans_state_in_range(uint32_t x)
{
    return x >= RANS_L && x < (RANS_L << 8);
}

/* one encode step (push): updates *x, writes renorm bytes to out[--*pos].
 * The quotient/remainder is replaced by a multiply-by-reciprocal (ryg rANS
 * / Granlund-Montgomery; rcp/rsh/bias are generated + verified in
 * tools/gen_rans_tables.py), bit-exact for the encoder's x range, so no
 * hardware divide touches the (public) signature data and the byte stream
 * -- hence the KAT -- is unchanged. */
static inline void enc_put(uint32_t *x, uint8_t *out, size_t *pos,
                           uint32_t freq, uint32_t rcp, uint32_t rsh,
                           uint32_t bias)
{
    uint32_t x_max = ((RANS_L >> PB) << 8) * freq;
    while (*x >= x_max) {
        out[--(*pos)] = (uint8_t)(*x & 0xff);
        *x >>= 8;
    }
    uint32_t q = (uint32_t)(((uint64_t)*x * rcp) >> 32) >> rsh;
    *x = *x + bias + q * (PSCALE - freq);
}

/* one encode step keyed by a contiguous-alphabet table; returns -2 if the
 * symbol is out of the modelled support (Sign restarts). */
#define ENC_ONE(TBL, X, OUT, POS, SYM)                                  \
    do {                                                                \
        int slot_ = (int)((SYM)-TBL##_LO);                              \
        if (slot_ < 0 || slot_ >= TBL##_N)                              \
            return -2;                                                  \
        enc_put((X), (OUT), (POS), TBL##_FREQ[slot_], TBL##_RCP[slot_], \
                TBL##_RSH[slot_], TBL##_BIAS[slot_]);                   \
    } while (0)

int shuttle_rans_encode(uint8_t *out, size_t *out_len, size_t cap,
                        const int32_t *q0, const int32_t *qs,
                        const int32_t *h, size_t nq0, size_t nqs,
                        size_t nh)
{
    uint32_t x[RN];
    for (int s = 0; s < RN; s++)
        x[s] = RANS_L;
    size_t pos = cap; /* write backwards */
    size_t ntot = nq0 + nqs + nh;

    /* flat S[t]: q0[t] (t<nq0, Q0) then qs[t-nq0] (Qs) then h[...] (HINT).
     * Push in REVERSE t; decode mirrors forward t with the same schedule.
     */
    for (size_t t = ntot; t-- > 0;) {
        unsigned s = rans_stream_index(t);
        if (pos < (size_t)(4 * RN + 4))
            return -2;
        if (t < nq0) {
            ENC_ONE(RANS_Q0, &x[s], out, &pos, q0[t]);
        } else if (t < nq0 + nqs) {
            ENC_ONE(RANS_QS, &x[s], out, &pos, qs[t - nq0]);
        } else {
            ENC_ONE(RANS_HINT, &x[s], out, &pos, h[t - nq0 - nqs]);
        }
    }
    /* flush all RN states, state 0 first (lands last in the stream, so
     * decode reads state RN-1 first). */
    for (int s = 0; s < RN; s++)
        for (int b = 0; b < 4; b++) {
            if (pos == 0)
                return -2;
            out[--pos] = (uint8_t)(x[s] & 0xff);
            x[s] >>= 8;
        }
    size_t used = cap - pos;
    memmove(out, out + pos, used);
    *out_len = used;
    return 0;
}

/* one decode step keyed by a table: emit the symbol value, update *x. */
#define DEC_ONE(TBL, X, DST)                                            \
    do {                                                                \
        uint32_t val_ = (*(X)) & (PSCALE - 1);                          \
        unsigned slot_ = TBL##_SLOT[val_];                              \
        (DST) = TBL##_LO + (int)slot_;                                  \
        *(X) =                                                          \
            TBL##_FREQ[slot_] * (*(X) >> PB) + val_ - TBL##_CDF[slot_]; \
    } while (0)

int shuttle_rans_decode(int32_t *q0, int32_t *qs, int32_t *h, size_t nq0,
                        size_t nqs, size_t nh, const uint8_t *in,
                        size_t in_len)
{
    size_t ntot = nq0 + nqs + nh;
    if (in_len < (size_t)(4 * RN))
        return -1;
    uint32_t x[RN];
    size_t bp = 0;
    /* init states in reverse flush order: stream front is state RN-1, ...,
     * then state 0 (mirror of the encoder's flush).  Range-check each (the
     * canonical initial-state check). */
    for (int s = RN - 1; s >= 0; s--) {
        uint32_t v = 0;
        for (int b = 0; b < 4; b++)
            v = (v << 8) | in[bp++];
        x[s] = v;
        if (!rans_state_in_range(x[s]))
            return -1;
    }

    for (size_t t = 0; t < ntot; t++) {
        unsigned s = rans_stream_index(t);
        /* val in [0, PSCALE) -> slot via O(1) direct lookup; CDF[0]=0,
         * CDF[N]=PSCALE so every val maps to a valid slot (no CDF hole).
         * The decoded SYMBOL VALUE is range-checked by the caller (the
         * per-block support + hint range-check live in unpack_sig). */
        if (t < nq0) {
            DEC_ONE(RANS_Q0, &x[s], q0[t]);
        } else if (t < nq0 + nqs) {
            DEC_ONE(RANS_QS, &x[s], qs[t - nq0]);
        } else {
            DEC_ONE(RANS_HINT, &x[s], h[t - nq0 - nqs]);
        }
        while (x[s] < RANS_L) {
            if (bp >= in_len)
                return -1;
            x[s] = (x[s] << 8) | in[bp++];
        }
    }
    /* Canonical stream check: all bytes consumed and all interleaved
     * states rewind to the encoder's initial state L. */
    if (bp != in_len)
        return -1;
    for (int s = 0; s < RN; s++)
        if (x[s] != RANS_L)
            return -1;
    return 0;
}
