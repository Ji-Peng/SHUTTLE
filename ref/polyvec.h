/*
 * polyvec.h -- vector-level XOF-driven sampling surface for SHUTTLE (P07).
 *
 * This header owns the deterministic seed expanders and the 16-lane
 * matrix/vector samplers that sit between the unified XOF layer (P02) and
 * the top-level KeyGen/Sign orchestration (P11):
 *
 *     expand_seeds          (DS 0x00, xof256, single one-shot)
 *     expand_signing_seeds  (DS 0x01, xof256, single one-shot)
 *     expand_a              (DS 0x02, xof128, 16 lanes; A direct in NTT)
 *     expand_s              (DS 0x03, xof256, 16 lanes; BaseSampler+fold)
 *     sample_c              (DS 0x07, xof256, single ctx; Fisher-Yates)
 *     sample_y              (DS 0x08, xof256, 16 lanes; wide Gaussian r)
 *
 * It pins the exact per-attempt PRNG byte schedule, the 16-lane stream
 * partition, the OUTPUT-indexed sign stream (K9), the cursor-advance
 * determinism rule (K6/K8), and the batch-orchestration plumbing
 * (gauss_stream) so that ref == avx2 == avx512 is byte-exact within each
 * MODE.  It is the KAT-fragile heart of the scheme (K2/K6/K7/K8/K9).
 *
 * ===================================================================== *
 *  THE NGCC DRBG NO-RATE-CURSOR RULE (binding; K6) -- read this.         *
 * ===================================================================== *
 * Under NGCC_MODE the SM3 Hash-DRBG behind xof256_squeeze has NO sub-call
 * rate cursor: every squeeze is a FRESH `SM3_DRNG_Generate` that produces
 * bytes from the current internal state V and then advances V by exactly
 * ONE step, regardless of how many bytes were requested.  Two consecutive
 * xof256_squeeze calls therefore DO NOT concatenate like SHAKE's
 * rate-buffered squeeze: the first byte of the second squeeze restarts
 * from a fresh generate, not from where the first left off.  (Verified
 * empirically: a 64-byte squeeze's prefix equals a 32-byte squeeze, but
 * its tail differs from a following 32-byte squeeze.)
 *
 * CONSEQUENCE: every producer here draws each logical stream's WHOLE byte
 * budget in ONE xof_squeeze call.  When a rejection-sampling stream needs
 * more than its initial buffer (statistically rare), it does NOT chain a
 * second squeeze on the same ctx; instead it RE-INITS a fresh ctx whose
 * absorbed nonce carries a 2-byte little-endian REFILL COUNTER appended
 * after the stream index, and squeezes the next block one-shot.  The same
 * refill discipline is applied under SHA3_MODE so both MODEs are
 * internally self-consistent (they are two DISTINCT KAT sets; within each
 * MODE, ref/avx2/avx512 must be byte-exact).
 *
 * The absorbed nonce layout for the 16-lane / single-ctx producers is thus
 *
 *     tag || seed || IntegerToBytes(stream_idx, 2) || IntegerToBytes(rc,
 * 2)
 *
 * with rc = 0 for the first (almost always only) block; the single-context
 * one-shots (ExpandSeeds/ExpandSigningSeeds) need no rc because their
 * output length is fixed and small.  SampleC is single-stream but
 * rejection-driven, so it uses stream_idx = 0 and a refill counter.
 */
#ifndef SHUTTLE_POLYVEC_H
#define SHUTTLE_POLYVEC_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "sampler.h"
#include "xof.h"

/* RNDBYTES: the per-signature randomness length fed to ExpandSigningSeeds.
 * For deterministic/KAT signing this is a fixed value owned by P11 (Sign);
 * P07 only fixes the byte LENGTH (= SEEDBYTES, the Dilithium rnd width).
 */
#ifndef RNDBYTES
#    define RNDBYTES SEEDBYTES
#endif

/* ===================================================================== *
 *  Deterministic seed expanders (single-context one-shot, xof256)       *
 * ===================================================================== */

/* T (len = 5*LAMBDA/8) = ExpandSeeds(0x00 || xi || IntToBytes(EM,2) ||
 *                                    IntToBytes(kappa,4)).
 * Slice T into seedA(SEEDBYTES) | seedsk(CHALLENGESEEDBYTES) |
 * masterSeed K(CHALLENGESEEDBYTES).  KeyGen passes kappa>=1. */
#define EXPAND_SEEDS_BYTES (5 * LAMBDA / 8) /* 80 / 160 / 320 */
void expand_seeds(uint8_t T[EXPAND_SEEDS_BYTES],
                  const uint8_t xi[SEEDBYTES], uint32_t kappa);

/* seedY (len = SEEDBYTES) = ExpandSigningSeeds(0x01 || K || rnd || mu ||
 *                                              IntToBytes(kappa,4)).
 * K = masterSeed (CHALLENGESEEDBYTES), rnd (RNDBYTES), mu
 * (CHALLENGESEEDBYTES).  Sign passes the CURRENT kappa (first iter 0). */
void expand_signing_seeds(uint8_t seedY[SEEDBYTES],
                          const uint8_t K[CHALLENGESEEDBYTES],
                          const uint8_t rnd[RNDBYTES],
                          const uint8_t mu[CHALLENGESEEDBYTES],
                          uint32_t kappa);

/* ===================================================================== *
 *  ExpandA (xof128, tag 0x02, 16 lanes; A_gen direct in NTT domain)     *
 * ===================================================================== */
/* agen : coeff-domain, EM polys, every coeff in [0,q).
 * hAgen: NTT-domain, EM*ELL polys, every coeff in [0,q).  hAgen is written
 * in canonical (ref ascending-index == bit-reversed) NTT order (K1); the
 * scalar ref needs no nttunpack (the AVX forks import via
 * poly_ntt_import).
 */
void expand_a(poly16 agen[EM], poly16 hAgen[EM * ELL],
              const uint8_t seedA[SEEDBYTES]);

/* ===================================================================== *
 *  ExpandS (xof256, tag 0x03, 16 lanes; BaseSampler + zero-fold + sign) *
 * ===================================================================== */
/* s1s2: ONE contiguous poly[ELL+EM]; the caller aliases
 *   poly *s = s1s2;  poly *e = s1s2 + ELL;
 * |s|<=L_sigma1, |e|<=L_sigma2.  Coefficient domain (signed int32). */
void expand_s(poly s1s2[ELL + EM],
              const uint8_t seedsk[CHALLENGESEEDBYTES]);

/* ===================================================================== *
 *  SampleC (xof256, tag 0x07, single ctx; partial Fisher-Yates)         *
 * ===================================================================== */
/* Binary challenge c in {0,1}^n, Hamming weight TAU, all set bits +1. */
void sample_c(poly *c, const uint8_t seedC[CHALLENGESEEDBYTES]);

/* ===================================================================== *
 *  SampleY / SampleDGauss (xof256, tag 0x08, 16 lanes; D_{Z,r=825})     *
 * ===================================================================== */
/* y: KVEC (= 1+ELL+EM) polys, flattened across 16 lanes; each coeff is a
 * wide discrete-Gaussian sample (scheme/coefficient domain, signed int32).
 */
void sample_y(poly y[KVEC], const uint8_t seedY[SEEDBYTES]);

/* ===================================================================== *
 *  gauss_stream -- persistent per-lane XOF buffer for the wide sampler   *
 *                                                                       *
 *  The scalar reference services one lane at a time (the degenerate      *
 *  LANES=1 fold of the fixed 16-stream flow).  gauss_stream wraps the    *
 *  ONE-SQUEEZE-PER-FILL discipline: gs_fill draws GAUSS_STREAM_BLOCK     *
 *  bytes for the current (lane, refill-counter) in a single xof256       *
 *  squeeze; gs_ensure refills (bumping the refill counter and re-initing *
 *  a fresh ctx) when the cursor would run past `avail`.  The whole       *
 *  mini-batch tail is requested at once so the cursor advances           *
 *  identically across backends (K6/K8). */
/* The wide-sampler lane chunk is KVEC*n/16 coeffs; its OUTPUT-indexed sign
 * stream is therefore (KVEC*n/16 + 7)/8 bytes (+ AVX512 LE64 read pad).  A
 * stream block must hold that sign stream plus one full mini-batch so the
 * "draw the whole tail up front" rule (K6/K8) fits in one squeeze. */
#define SIGN_BYTES_PER_CHUNK (((KVEC * N / XOF_STREAMS) + 7) / 8)
#define GAUSS_STREAM_BLOCK \
    (MINIBATCH_RAND_BYTES + SIGN_BYTES_PER_CHUNK + SIGN_PAD_AVX512 + 64)
typedef struct {
    xof_ctx ctx;         /* current backend xof256 state               */
    uint8_t tag;         /* domain-separation tag (DS_SAMPLE_Y)        */
    const uint8_t *seed; /* lane seed (SEEDBYTES)                      */
    uint16_t lane;       /* logical stream index 0..15                 */
    uint16_t refill;     /* refill counter (continuation nonce)        */
    size_t pos, avail;   /* cursor into buf                            */
    uint8_t buf[2 * GAUSS_STREAM_BLOCK]; /* one block + headroom */
} gauss_stream;

void gauss_stream_init(gauss_stream *gs, uint8_t tag, const uint8_t *seed,
                       unsigned lane);
/* Ensure at least `need` contiguous bytes are available at
 * gs->buf+gs->pos; memmoves the leftover to the front and draws a fresh
 * one-shot block. */
void gs_ensure(gauss_stream *gs, size_t need);
/* Produce a flat run of `count` wide-Gaussian samples into dst[0..count).
 * The OUTPUT-indexed sign stream covers exactly these `count` outputs and
 * is squeezed up front (K9); the wide mini-batches then fill them. `count`
 * is the lane chunk width KVEC*n/16 (<= GAUSS_STREAM_BLOCK-worth of
 * signs).
 */
void gauss_stream_chunk(gauss_stream *gs, int32_t *dst, size_t count);

#endif /* SHUTTLE_POLYVEC_H */
