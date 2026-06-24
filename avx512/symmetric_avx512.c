/*
 * symmetric_avx512.c -- AVX512 lane-batched XOF variant bodies (P02-T6).
 *
 * Implements xof128/256_avx512_init/squeeze on top of the N-way
 * primitives. Same contract as symmetric_avx2.c (lengths shared, pointer
 * arrays per-lane, nothing validated), at the AVX512 lane width.
 *
 *   NGCC_MODE: lanes = XOF_LANES_AVX512 = 16 (SM3 DRBG).  Body =
 *       init_random_number_avx512 / get_random_number_avx512(..,
 * out_len*8). In the fixed 16-stream flow this is a SINGLE pass.
 *
 *   SHA3_MODE: lanes = XOF_LANES_AVX512 = 8 (Keccak).  Body =
 *       shake{128,256}x8_absorb_once + buffered
 * shake{128,256}x8_squeezeblocks. In the 16-stream flow this is TWO passes
 * (8 lanes x 2).
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define SHUTTLE_XOF_DECLARE_AVX512 1
#include "symmetric.h"

#if defined(SHA3_MODE)

#    include "fips202x8.h"

static void shakex8_squeeze_buffered(
    uint8_t *const out[8], size_t out_len, keccakx8_state *st,
    unsigned int rate,
    void (*sqblk)(uint8_t *, uint8_t *, uint8_t *, uint8_t *, uint8_t *,
                  uint8_t *, uint8_t *, uint8_t *, size_t,
                  keccakx8_state *))
{
    size_t nblocks = out_len / rate;
    size_t off = nblocks * rate;
    size_t tail = out_len - off;
    uint8_t t[8][SHAKE128_RATE];

    if (nblocks) {
        sqblk(out[0], out[1], out[2], out[3], out[4], out[5], out[6],
              out[7], nblocks, st);
    }
    if (tail) {
        size_t k;
        sqblk(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], 1, st);
        for (k = 0; k < 8; k++) {
            memcpy(out[k] + off, t[k], tail);
        }
    }
}

void xof128_avx512_init(xof_ctx_avx512 *ctx,
                        const uint8_t *const seed[XOF_LANES_AVX512],
                        size_t seed_len)
{
    shake128x8_absorb_once(ctx, seed[0], seed[1], seed[2], seed[3],
                           seed[4], seed[5], seed[6], seed[7], seed_len);
}

void xof128_avx512_squeeze(xof_ctx_avx512 *ctx,
                           uint8_t *const out[XOF_LANES_AVX512],
                           size_t out_len)
{
    shakex8_squeeze_buffered(out, out_len, ctx, SHAKE128_RATE,
                             shake128x8_squeezeblocks);
}

void xof256_avx512_init(xof_ctx_avx512 *ctx,
                        const uint8_t *const seed[XOF_LANES_AVX512],
                        size_t seed_len)
{
    shake256x8_absorb_once(ctx, seed[0], seed[1], seed[2], seed[3],
                           seed[4], seed[5], seed[6], seed[7], seed_len);
}

void xof256_avx512_squeeze(xof_ctx_avx512 *ctx,
                           uint8_t *const out[XOF_LANES_AVX512],
                           size_t out_len)
{
    shakex8_squeeze_buffered(out, out_len, ctx, SHAKE256_RATE,
                             shake256x8_squeezeblocks);
}

#else /* NGCC_MODE: 16-way SM3 DRBG ==================================== \
       */

#    include "drng_avx512.h"

void xof256_avx512_init(xof_ctx_avx512 *ctx,
                        const uint8_t *const seed[XOF_LANES_AVX512],
                        size_t seed_len)
{
    (void)init_random_number_avx512(ctx, seed,
                                    (unsigned long long)seed_len);
}

void xof256_avx512_squeeze(xof_ctx_avx512 *ctx,
                           uint8_t *const out[XOF_LANES_AVX512],
                           size_t out_len)
{
    (void)get_random_number_avx512(ctx, out,
                                   (unsigned long long)out_len * 8u);
}

void xof128_avx512_init(xof_ctx_avx512 *ctx,
                        const uint8_t *const seed[XOF_LANES_AVX512],
                        size_t seed_len)
{
    xof256_avx512_init(ctx, seed, seed_len);
}

void xof128_avx512_squeeze(xof_ctx_avx512 *ctx,
                           uint8_t *const out[XOF_LANES_AVX512],
                           size_t out_len)
{
    xof256_avx512_squeeze(ctx, out, out_len);
}

#endif /* SHA3_MODE / NGCC_MODE */
