// author: https://github.com/8891689
#ifndef SM3_AVX_H
#define SM3_AVX_H

#include <stdint.h>
#include <stddef.h>
#include <immintrin.h>

// --- Single-channel SM3 context and functions ---
// SM3 single-channel context structure
typedef struct {
    uint32_t state[8];
} sm3_context;

// Single-channel SM3 function declarations
void sm3_starts(sm3_context *ctx);
void sm3_compress(uint32_t state[8], const unsigned char block[64]);
void sm3_single(const unsigned char *input, size_t ilen, unsigned char *output);


// --- 8-channel AVX2 SM3 context and functions ---
// SM3 8-channel AVX2 context structure
typedef struct {
    __m256i state[8];
    __m256i active_mask;
} sm3_8x_context;

// 8-channel SM3 function declarations (public interface)
void sm3_8x_starts(sm3_8x_context *ctx);
void sm3_8x_compress(sm3_8x_context *ctx, const unsigned char blocks[8][64]);
void sm3_8x_final(sm3_8x_context *ctx, unsigned char outputs[8][32]);
void sm3_8x(const unsigned char *inputs[8], size_t ilens[8], unsigned char outputs[8][32]);

#endif // SM3_AVX_H
