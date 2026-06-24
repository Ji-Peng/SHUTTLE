/*
 * symmetric_avx2.c -- AVX2 lane-batched XOF variant bodies (P02-T6).
 *
 * Implements xof128/256_avx2_init/squeeze on top of the already-built
 * N-way primitives.  The signature shape mirrors the existing N-way
 * contract EXACTLY: lengths are SHARED across lanes; the pointer arrays
 * are PER-LANE; NOTHING is validated.  These wrappers do not defeat the
 * N-way perf: they reuse one ctx across refills and apply the bytes->bits
 * *8 shim in the NGCC squeeze.
 *
 *   NGCC_MODE: lanes = XOF_LANES_AVX2 = 8 (SM3 DRBG).  Body =
 *       init_random_number_avx2 / get_random_number_avx2(.., out_len*8).
 *       xof128_avx2_* alias xof256_avx2_* (the 128/256 collapse).
 *
 *   SHA3_MODE: lanes = XOF_LANES_AVX2 = 4 (Keccak).  Body =
 *       shake{128,256}x4_absorb_once over the shared seed, then full-rate
 *       shake{128,256}x4_squeezeblocks with a buffered tail for the final
 *       partial block.  xof128 = SHAKE128 (rate 168), xof256 = SHAKE256
 * (136).
 *
 * Compiled only in the avx2 build (USE_AVX2_SHAKE4X selects the SHA3 path;
 * the NGCC path is the default).  Lane-equivalence (lane k == scalar xof
 * for stream k) is exercised by test_xof_nway (the SM3 N-way path is
 * already proven byte-exact to scalar by the existing sm3-test sweep; this
 * wrapper only adds the pointer marshaling + the *8 shim, which
 * test_xof_nway re-checks).
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Pull in the xof_ctx_avx2 typedef + the four xof*_avx2_* prototypes
 * (declared in symmetric.h under this guard); our definitions below must
 * match them. */
#define SHUTTLE_XOF_DECLARE_AVX2 1
#include "symmetric.h"

#if defined(SHA3_MODE)

/* 4-lane SHAKE squeeze with a buffered final partial block.  out[k]
 * receives out_len bytes for lane k.  We squeeze whole rate blocks
 * straight into the destination, then one extra block into a scratch row
 * to copy the tail. */
static void shakex4_squeeze_buffered(uint8_t *const out[4], size_t out_len,
                                     keccakx4_state *st, unsigned int rate,
                                     void (*sqblk)(uint8_t *, uint8_t *,
                                                   uint8_t *, uint8_t *,
                                                   size_t,
                                                   keccakx4_state *))
{
    size_t nblocks = out_len / rate;
    size_t off = nblocks * rate;
    size_t tail = out_len - off;
    uint8_t t0[SHAKE128_RATE], t1[SHAKE128_RATE];
    uint8_t t2[SHAKE128_RATE], t3[SHAKE128_RATE];

    if (nblocks) {
        sqblk(out[0], out[1], out[2], out[3], nblocks, st);
    }
    if (tail) {
        sqblk(t0, t1, t2, t3, 1, st);
        memcpy(out[0] + off, t0, tail);
        memcpy(out[1] + off, t1, tail);
        memcpy(out[2] + off, t2, tail);
        memcpy(out[3] + off, t3, tail);
    }
}

void xof128_avx2_init(xof_ctx_avx2 *ctx,
                      const uint8_t *const seed[XOF_LANES_AVX2],
                      size_t seed_len)
{
    shake128x4_absorb_once(ctx, seed[0], seed[1], seed[2], seed[3],
                           seed_len);
}

void xof128_avx2_squeeze(xof_ctx_avx2 *ctx,
                         uint8_t *const out[XOF_LANES_AVX2],
                         size_t out_len)
{
    shakex4_squeeze_buffered(out, out_len, ctx, SHAKE128_RATE,
                             shake128x4_squeezeblocks);
}

void xof256_avx2_init(xof_ctx_avx2 *ctx,
                      const uint8_t *const seed[XOF_LANES_AVX2],
                      size_t seed_len)
{
    shake256x4_absorb_once(ctx, seed[0], seed[1], seed[2], seed[3],
                           seed_len);
}

void xof256_avx2_squeeze(xof_ctx_avx2 *ctx,
                         uint8_t *const out[XOF_LANES_AVX2],
                         size_t out_len)
{
    shakex4_squeeze_buffered(out, out_len, ctx, SHAKE256_RATE,
                             shake256x4_squeezeblocks);
}

#else /* NGCC_MODE: 8-way SM3 DRBG ====================================== \
       */

#    include "drng_avx2.h"

void xof256_avx2_init(xof_ctx_avx2 *ctx,
                      const uint8_t *const seed[XOF_LANES_AVX2],
                      size_t seed_len)
{
    (void)init_random_number_avx2(ctx, seed, (unsigned long long)seed_len);
}

void xof256_avx2_squeeze(xof_ctx_avx2 *ctx,
                         uint8_t *const out[XOF_LANES_AVX2],
                         size_t out_len)
{
    /* bytes -> bits *8 shim, identical to the scalar path; per-lane out[].
     */
    (void)get_random_number_avx2(ctx, out,
                                 (unsigned long long)out_len * 8u);
}

void xof128_avx2_init(xof_ctx_avx2 *ctx,
                      const uint8_t *const seed[XOF_LANES_AVX2],
                      size_t seed_len)
{
    xof256_avx2_init(ctx, seed, seed_len);
}

void xof128_avx2_squeeze(xof_ctx_avx2 *ctx,
                         uint8_t *const out[XOF_LANES_AVX2],
                         size_t out_len)
{
    xof256_avx2_squeeze(ctx, out, out_len);
}

#endif /* SHA3_MODE / NGCC_MODE */
