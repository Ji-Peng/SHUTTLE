/* sm3.c */
// author: https://github.com/8891689
#include "sm3.h"
#include <string.h>

/* Optimized 32-bit left rotation */
#define ROTL32(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

/* Boolean functions */
#define FF0(a, b, c) ((a) ^ (b) ^ (c))
#define FF1(a, b, c) (((a) & (b)) | ((a) & (c)) | ((b) & (c)))
#define GG0(x, y, z) ((x) ^ (y) ^ (z))
#define GG1(x, y, z) (((x) & (y)) | ((~(x)) & (z)))

/* Permutation functions */
#define P0(x) ((x) ^ ROTL32((x), 9) ^ ROTL32((x), 17))
#define P1(x) ((x) ^ ROTL32((x), 15) ^ ROTL32((x), 23))

/* Initialization vector */
static const uint32_t IV[8] = {
    0x7380166F, 0x4914B2B9,
    0x172442D7, 0xDA8A0600,
    0xA96F30BC, 0x163138AA,
    0xE38DEE4D, 0xB0FB0E4E
};

/* Big-endian load */
static inline uint32_t load_be32(const uint8_t in[4]) {
    return ((uint32_t)in[0] << 24) | 
           ((uint32_t)in[1] << 16) | 
           ((uint32_t)in[2] << 8) | 
           (uint32_t)in[3];
}

/* Big-endian store */
static inline void store_be32(uint8_t out[4], uint32_t w) {
    out[0] = (w >> 24) & 0xFF;
    out[1] = (w >> 16) & 0xFF;
    out[2] = (w >> 8) & 0xFF;
    out[3] = w & 0xFF;
}

/* Compression function */
static void sm3_compress(uint32_t state[8], const uint8_t block[64]) {
    uint32_t W[68];
    uint32_t A = state[0], B = state[1], C = state[2], D = state[3];
    uint32_t E = state[4], F = state[5], G = state[6], H = state[7];
    
    // 1. Message expansion
    for (int j = 0; j < 16; j++) {
        W[j] = load_be32(block + j * 4);
    }

    // Message expansion loop
    for (int j = 16; j < 68; j++) {
        uint32_t tmp = W[j-16] ^ W[j-9] ^ ROTL32(W[j-3], 15);
        W[j] = P1(tmp) ^ ROTL32(W[j-13], 7) ^ W[j-6];
    }
    
    // 2. Message expansion (W')
    uint32_t WW[64];
    for (int j = 0; j < 64; j++) {
        WW[j] = W[j] ^ W[j+4];
    }
    
    // 3. Main compression loop
    for (int j = 0; j < 64; j++) {
        // Compute T_j constant (direct calculation)
        uint32_t T_j = (j < 16) ? 0x79CC4519 : 0x7A879D8A;
        T_j = ROTL32(T_j, j);

        // Precompute A rotation results
        uint32_t A_rot12 = ROTL32(A, 12);
        uint32_t tmp1 = A_rot12 + E + T_j;
        uint32_t SS1 = ROTL32(tmp1, 7);
        uint32_t SS2 = SS1 ^ A_rot12;
        
        uint32_t TT1, TT2;
        if (j < 16) {
            TT1 = FF0(A, B, C) + D + SS2 + WW[j];
            TT2 = GG0(E, F, G) + H + SS1 + W[j];
        } else {
            TT1 = FF1(A, B, C) + D + SS2 + WW[j];
            TT2 = GG1(E, F, G) + H + SS1 + W[j];
        }
        
        // State update
        D = C;
        C = ROTL32(B, 9);
        B = A;
        A = TT1;
        H = G;
        G = ROTL32(F, 19);
        F = E;
        E = P0(TT2);
    }
    
    // 4. Update state
    state[0] ^= A;
    state[1] ^= B;
    state[2] ^= C;
    state[3] ^= D;
    state[4] ^= E;
    state[5] ^= F;
    state[6] ^= G;
    state[7] ^= H;
}

/* Initialization */
void sm3_init(sm3_ctx_t *ctx) {
    memcpy(ctx->state, IV, sizeof(IV));
    ctx->total_len = 0;
    ctx->buf_len = 0;
}

/* Optimized batch data processing */
void sm3_update(sm3_ctx_t *ctx, const uint8_t *data, size_t len) {
    const size_t block_size = 64;
    ctx->total_len += len;
    
    size_t buf_len = ctx->buf_len;
    
    // Fill buffer
    if (buf_len > 0) {
        size_t space = block_size - buf_len;
        if (len < space) {
            memcpy(ctx->buffer + buf_len, data, len);
            ctx->buf_len = buf_len + len;
            return;
        }
        
        memcpy(ctx->buffer + buf_len, data, space);
        sm3_compress(ctx->state, ctx->buffer);
        data += space;
        len -= space;
        buf_len = 0;
    }
    
    // Process full blocks
    while (len >= block_size) {
        sm3_compress(ctx->state, data);
        data += block_size;
        len -= block_size;
    }
    
    // Save remaining data
    if (len > 0) {
        memcpy(ctx->buffer, data, len);
        buf_len = len;
    }
    
    ctx->buf_len = buf_len;
}

/* Finalization */
void sm3_final(sm3_ctx_t *ctx, uint8_t digest[32]) {
    // Compute total message length in bits
    uint64_t bit_len = (uint64_t)ctx->total_len * 8;
    size_t buf_len = ctx->buf_len;

    // Compute padding length
    size_t pad_len = (buf_len < 56) ? 56 - buf_len : 120 - buf_len;

    // Build padding block directly
    uint8_t padding[128] = {0};
    padding[0] = 0x80;

    // Append length information
    store_be32(padding + pad_len - 8, (uint32_t)(bit_len >> 32));
    store_be32(padding + pad_len - 4, (uint32_t)bit_len);

    // Process padding
    sm3_update(ctx, padding, pad_len);

    // Ensure final block is processed
    if (ctx->buf_len > 0) {
        sm3_compress(ctx->state, ctx->buffer);
        ctx->buf_len = 0;
    }

    // Output digest
    for (int i = 0; i < 8; i++) {
        store_be32(digest + i * 4, ctx->state[i]);
    }
}

/* High-performance one-shot interface */
void sm3(const uint8_t *data, size_t len, uint8_t digest[32]) {
    // Create local context to avoid struct overhead
    uint32_t state[8];
    memcpy(state, IV, sizeof(IV));
    size_t total_blocks = len / 64;
    size_t tail_len = len % 64;

    // Process full blocks
    for (size_t i = 0; i < total_blocks; i++) {
        sm3_compress(state, data + i * 64);
    }

    // Prepare tail data
    uint8_t tail_block[128] = {0};
    size_t pad_pos = tail_len;

    // Copy tail data
    if (tail_len > 0) {
        memcpy(tail_block, data + len - tail_len, tail_len);
    }

    // Append padding bit
    tail_block[pad_pos++] = 0x80;

    // Compute padding length
    size_t pad_len;
    if (tail_len < 56) {
        pad_len = 56 - tail_len;
    } else {
        pad_len = 120 - tail_len;
    }

    // Append zero bytes for padding
    if (pad_len > 1) {
        memset(tail_block + pad_pos, 0, pad_len - 1);
        pad_pos += pad_len - 1;
    }

    // Append bit length
    uint64_t bit_len = (uint64_t)len * 8;
    store_be32(tail_block + pad_pos, (uint32_t)(bit_len >> 32));
    store_be32(tail_block + pad_pos + 4, (uint32_t)bit_len);

    // Process padded blocks
    sm3_compress(state, tail_block);

    // If there is a second padding block
    if (tail_len + pad_len + 8 > 64) {
        sm3_compress(state, tail_block + 64);
    }

    // Output digest
    for (int i = 0; i < 8; i++) {
        store_be32(digest + i * 4, state[i]);
    }
}
