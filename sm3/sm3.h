/* sm3.h */
// author: https://github.com/8891689
#ifndef SM3_H
#define SM3_H

#include <stddef.h>
#include <stdint.h>

/* SM3 computation context */
typedef struct {
    uint32_t total_len;      /* Total bytes processed (lower 32 bits) */
    uint32_t state[8];       /* Intermediate hash state */
    uint8_t  buffer[64];     /* Unprocessed partial block data */
    size_t   buf_len;        /* Number of bytes currently in buffer */
} sm3_ctx_t;

/**
 * @brief Initialize SM3 context (set initial IV)
 */
void sm3_init(sm3_ctx_t *ctx);

/**
 * @brief Feed data into SM3 context
 * @param ctx  SM3 context
 * @param data Input data
 * @param len  Data length in bytes
 */
void sm3_update(sm3_ctx_t *ctx, const uint8_t *data, size_t len);

/**
 * @brief Finalize hash computation and output 32-byte digest
 * @param ctx    SM3 context with all message data fed
 * @param digest Output buffer (at least 32 bytes)
 */
void sm3_final(sm3_ctx_t *ctx, uint8_t digest[32]);

/**
 * @brief One-shot SM3 hash computation
 * @param data   Input data
 * @param len    Input length
 * @param digest Output buffer (32 bytes)
 */
void sm3(const uint8_t *data, size_t len, uint8_t digest[32]);

#endif /* SM3_H */
