/*
 * rng.h - Deterministic SHAKE-256 DRBG used for test-vector generation.
 *
 * Replaces the default /dev/urandom-backed randombytes() from
 * ../randombytes.c. Output bytes depend only on the 48-byte seed,
 * so test vectors are byte-for-byte reproducible.
 */
#ifndef SHUTTLE_VECTOR_RNG_H
#define SHUTTLE_VECTOR_RNG_H

#include <stddef.h>
#include <stdint.h>

/* Re-seed the DRBG. Call before any randombytes() usage in a phase. */
void randombytes_seed(const uint8_t seed[48]);

/* randombytes() with the same signature as ../randombytes.h. */
void randombytes(uint8_t *out, size_t outlen);

#endif
