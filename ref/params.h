/*
 * params.h - Parameters for SHUTTLE signature scheme.
 *
 * Ring: R_q = Z_q[x]/(x^n + 1), q = 15361, n = 256.
 * Module: secret key [1, s, e] with s in R^l, e in R^m.
 * Gaussian sampler: sigma selectable among {101, 128, 149} via
 * SHUTTLE_SIGMA (default 101 = SHUTTLE-128 parameter set).
 */

#ifndef SHUTTLE_PARAMS_H
#define SHUTTLE_PARAMS_H

#include <stdint.h>
#include "config.h"

/* ============================================================
 * Ring parameters
 * ============================================================ */
#define SHUTTLE_N              256    /* Polynomial ring degree */
#define SHUTTLE_Q              15361  /* Modulus (prime) */
#define SHUTTLE_ROOT_OF_UNITY  5301   /* Primitive 512-th root of unity mod q */

/* ============================================================
 * Module structure
 * ============================================================ */
#define SHUTTLE_L       3   /* Number of polynomial components in s */
#define SHUTTLE_M       2   /* Number of polynomial components in e */
#define SHUTTLE_VECLEN  (1 + SHUTTLE_L + SHUTTLE_M)  /* 6: full sk vector [1,s,e] */

/* ============================================================
 * Secret key distribution
 * ============================================================ */
#define SHUTTLE_ETA     1   /* CBD parameter for s and e */

/* ============================================================
 * Signing parameters
 * ============================================================ */
#define SHUTTLE_TAU       30    /* Hamming weight of challenge polynomial */
#define SHUTTLE_ALPHA_H   128   /* Compression parameter for hint */
#define SHUTTLE_ALPHA_1   1     /* Compression parameter for y[0] (1 = no compression) */

/* ============================================================
 * Norm bounds
 * ============================================================ */
#define SHUTTLE_BK       26         /* Secret key norm upper bound */
#define SHUTTLE_BK_LOW   0          /* Secret key norm lower bound */
#define SHUTTLE_BS_SQ    17923948   /* floor(4233.70^2) - signing bound squared */
#define SHUTTLE_BV_SQ    25000000   /* 5000^2 - verification bound squared */

/* ============================================================
 * Hash and seed sizes (bytes)
 * ============================================================ */
#define SHUTTLE_SEEDBYTES    32
#define SHUTTLE_CRHBYTES     64
#define SHUTTLE_TRBYTES      64
#define SHUTTLE_RNDBYTES     32
#define SHUTTLE_CTILDEBYTES  32

/* ============================================================
 * Packing sizes (bytes)
 * ============================================================ */
/* eta=1: coefficients in {-1,0,1}, encode as 2 bits/coeff */
#define SHUTTLE_POLYETA_PACKEDBYTES   64   /* 256 * 2 / 8 */

/* Public key b: coefficients in [0, q-1], 14 bits/coeff */
#define SHUTTLE_POLYPK_PACKEDBYTES    448  /* 256 * 14 / 8 */

/* Signature z1: Gaussian bounded coeff, ~13 bits/coeff + sign */
#define SHUTTLE_POLYZ_PACKEDBYTES     448  /* 256 * 14 / 8 (signed 14-bit) */

/* High bits of commitment for hashing */
#define SHUTTLE_POLYW1_PACKEDBYTES    192  /* 256 * 6 / 8 (6-bit high) */

/* Key sizes */
#define SHUTTLE_PUBLICKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYPK_PACKEDBYTES)
    /* 32 + 2*448 = 928 */

#define SHUTTLE_SECRETKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_TRBYTES \
                                   + SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_L * SHUTTLE_POLYETA_PACKEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYETA_PACKEDBYTES)
    /* 32 + 64 + 32 + 3*64 + 2*64 = 448 */

/* Hint encoding: per-coefficient 1 bit in each of m polys.
 * Reserved for future compressed signature format. */
#define SHUTTLE_POLYVECH_PACKEDBYTES  (SHUTTLE_M * (SHUTTLE_N / 8 + 1))
    /* 2 * 33 = 66 */

/* IRS sign bits: ceil(TAU/8) bytes encoding the per-monomial sign choices
 * from iterative rejection sampling. The verifier needs these to reconstruct
 * the effective challenge c_eff = sum(sigma_i * c_i * x^{j_i}). */
#define SHUTTLE_IRS_SIGNBYTES  ((SHUTTLE_TAU + 7) / 8)
    /* ceil(30/8) = 4 */

/* Signature format (Version 1 - full response):
 * c_tilde (32B) || irs_signs (4B) || polyz_pack(z[0..VECLEN-1]) (6*448B)
 * = 2724B
 *
 * This publishes the full VECLEN=6 response vector z plus the IRS sign bits.
 * The verifier reconstructs the effective challenge from (seed_c, irs_signs),
 * then checks w = B*z - c_eff*b has HighBits matching the commitment.
 * A future version will compress z2 for smaller signatures. */
#define SHUTTLE_BYTES  (SHUTTLE_CTILDEBYTES \
                          + SHUTTLE_IRS_SIGNBYTES \
                          + SHUTTLE_VECLEN * SHUTTLE_POLYZ_PACKEDBYTES)
    /* 32 + 4 + 6*448 = 2724 */

/* ============================================================
 * Gaussian sampler parameters (mode-dependent).
 *
 * SHUTTLE_SIGMA selects the discrete Gaussian standard deviation and is
 * defined in config.h (default 101). Three modes are supported:
 *
 *   sigma | small sigma | RCDT | 2*sigma^2 |   Renyi (1025)  | parameter set
 *   ------+------------+------+-----------+-----------------+---------------
 *    128  |    2.0000  |  22  |  32768    | 1 + 2^-93.85    | legacy / ref
 *    101  |    1.5781  |  17  |  20402    | 1 + 2^-93.40    | SHUTTLE-128
 *    149  |    2.3281  |  26  |  44402    | 1 + 2^-93.75    | SHUTTLE-256
 *
 * Sampling (BLISS-style convolution):
 *   1) x ~ D_{small_sigma} via 93-bit RCDT (3x31-bit limbs)
 *   2) y ~ Uniform{0, ..., 2^SHUTTLE_Y_BITS - 1}
 *   3) candidate r = (x << SHUTTLE_K_BITS) | y
 *   4) accept with probability exp(-y*(y + 2*k*x) / (2*sigma^2))
 *   5) random sign
 *
 * All three modes share k = 64 (so K_BITS = 6) and truncation = 11*sigma.
 * ============================================================ */

#if SHUTTLE_SIGMA == 128
#  define RCDT_ENTRIES 22
#elif SHUTTLE_SIGMA == 101
#  define RCDT_ENTRIES 17
#elif SHUTTLE_SIGMA == 149
#  define RCDT_ENTRIES 26
#else
#  error "Unsupported SHUTTLE_SIGMA (expected 128, 101, or 149)"
#endif

/* All three modes share the same convolution factor k = 64. Defined
 * inside this block so all sampler-related constants stay mode-scoped
 * and audit-friendly. */
#define SHUTTLE_GAUSS_K      64
#define SHUTTLE_K_BITS       6    /* log2(k) */
#define SHUTTLE_Y_BITS       6    /* # of random bits per y; y in [0, 2^6) */
#define SHUTTLE_TWO_K_BITS   7    /* log2(2k); used for 2kx in the rejection exponent */
#define SHUTTLE_TRUNC        11
#define SHUTTLE_BOUND        (SHUTTLE_TRUNC * SHUTTLE_SIGMA)

/* Reciprocal of 2*sigma^2 in Q64. Used by the rejection-sampling step
 * when 2*sigma^2 is NOT a power of two (true for sigma=101 and 149,
 * where N = 20402, 44402 respectively). For sigma=128 this macro is
 * still defined but unused, since 2*sigma^2 = 2^15 admits an exact
 * left-shift fast path.
 *
 * Derivation: (uint64_t)(-1ULL / N) = floor((2^64 - 1) / N), which
 * equals floor(2^64 / N) whenever N does NOT divide 2^64. Both 20402
 * and 44402 are non-powers-of-two, so the equality holds. */
#define SHUTTLE_INV_2SIGMA2_Q64 \
    ((uint64_t)(-1ULL / (2ULL * (uint64_t)SHUTTLE_SIGMA * (uint64_t)SHUTTLE_SIGMA)))

#endif /* SHUTTLE_PARAMS_H */
