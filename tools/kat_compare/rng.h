/*
 * rng.h - Deterministic SHAKE-256 DRBG for KAT comparison.
 *
 * Replaces the default /dev/urandom-backed randombytes() from
 * ../../ref/randombytes.c with a seeded PRF so identical inputs
 * deterministically produce identical outputs in both ref/ and
 * avx2/ builds. Output order and byte-for-byte content depend only
 * on the 48-byte seed, so the two implementations can be compared
 * on the same randomness schedule.
 */
#ifndef KAT_RNG_H
#define KAT_RNG_H

#include <stddef.h>
#include <stdint.h>

/* Re-seed the DRBG. Call before any randombytes() usage; the 48-byte
 * size is inherited from the NIST PQCgenKAT_sign convention. */
void randombytes_seed(const uint8_t seed[48]);

/* randombytes() with the same signature as ref/randombytes.h so it
 * link-replaces the default implementation. */
void randombytes(uint8_t *out, size_t outlen);

#endif
