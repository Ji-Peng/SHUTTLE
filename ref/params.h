/*
 * params.h - Parameters for SHUTTLE signature scheme.
 *
 * Ring: R_q = Z_q[x]/(x^n + 1), q = 15361, n = 256.
 * Module: secret key [1, s, e] with s in R^l, e in R^m.
 * Gaussian sampler: sigma = 128 = 2^7.
 */

#ifndef SHUTTLE_PARAMS_H
#define SHUTTLE_PARAMS_H

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
 * Gaussian sampler parameters
 * ============================================================ */
#define SHUTTLE_SIGMA     128
#define SHUTTLE_GAUSS_K   64     /* 2^6, convolution decomposition factor */
#define SHUTTLE_K_BITS    6      /* log2(GAUSS_K) */
#define SHUTTLE_TRUNC     11     /* truncation factor */
#define SHUTTLE_BOUND     (SHUTTLE_TRUNC * SHUTTLE_SIGMA)  /* 1408 */

#endif /* SHUTTLE_PARAMS_H */
