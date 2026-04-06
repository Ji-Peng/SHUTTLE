/*
 * params.h - Parameters for NGCC_SIGN signature scheme.
 *
 * Ring: R_q = Z_q[x]/(x^n + 1), q = 15361, n = 256.
 * Module: secret key [1, s, e] with s in R^l, e in R^m.
 * Gaussian sampler: sigma = 128 = 2^7.
 */

#ifndef NGCC_SIGN_PARAMS_H
#define NGCC_SIGN_PARAMS_H

#include "config.h"

/* ============================================================
 * Ring parameters
 * ============================================================ */
#define NGCC_SIGN_N              256    /* Polynomial ring degree */
#define NGCC_SIGN_Q              15361  /* Modulus (prime) */
#define NGCC_SIGN_ROOT_OF_UNITY  5301   /* Primitive 512-th root of unity mod q */

/* ============================================================
 * Module structure
 * ============================================================ */
#define NGCC_SIGN_L       3   /* Number of polynomial components in s */
#define NGCC_SIGN_M       2   /* Number of polynomial components in e */
#define NGCC_SIGN_VECLEN  (1 + NGCC_SIGN_L + NGCC_SIGN_M)  /* 6: full sk vector [1,s,e] */

/* ============================================================
 * Secret key distribution
 * ============================================================ */
#define NGCC_SIGN_ETA     1   /* CBD parameter for s and e */

/* ============================================================
 * Signing parameters
 * ============================================================ */
#define NGCC_SIGN_TAU       30    /* Hamming weight of challenge polynomial */
#define NGCC_SIGN_ALPHA_H   128   /* Compression parameter for hint */
#define NGCC_SIGN_ALPHA_1   1     /* Compression parameter for y[0] (1 = no compression) */

/* ============================================================
 * Norm bounds
 * ============================================================ */
#define NGCC_SIGN_BK       26         /* Secret key norm upper bound */
#define NGCC_SIGN_BK_LOW   0          /* Secret key norm lower bound */
#define NGCC_SIGN_BS_SQ    17923948   /* floor(4233.70^2) - signing bound squared */
#define NGCC_SIGN_BV_SQ    25000000   /* 5000^2 - verification bound squared */

/* ============================================================
 * Hash and seed sizes (bytes)
 * ============================================================ */
#define NGCC_SIGN_SEEDBYTES    32
#define NGCC_SIGN_CRHBYTES     64
#define NGCC_SIGN_TRBYTES      64
#define NGCC_SIGN_RNDBYTES     32
#define NGCC_SIGN_CTILDEBYTES  32

/* ============================================================
 * Packing sizes (bytes)
 * ============================================================ */
/* eta=1: coefficients in {-1,0,1}, encode as 2 bits/coeff */
#define NGCC_SIGN_POLYETA_PACKEDBYTES   64   /* 256 * 2 / 8 */

/* Public key b: coefficients in [0, q-1], 14 bits/coeff */
#define NGCC_SIGN_POLYPK_PACKEDBYTES    448  /* 256 * 14 / 8 */

/* Signature z1: Gaussian bounded coeff, ~13 bits/coeff + sign */
#define NGCC_SIGN_POLYZ_PACKEDBYTES     448  /* 256 * 14 / 8 (signed 14-bit) */

/* High bits of commitment for hashing */
#define NGCC_SIGN_POLYW1_PACKEDBYTES    192  /* 256 * 6 / 8 (6-bit high) */

/* Key sizes */
#define NGCC_SIGN_PUBLICKEYBYTES  (NGCC_SIGN_SEEDBYTES \
                                   + NGCC_SIGN_M * NGCC_SIGN_POLYPK_PACKEDBYTES)
    /* 32 + 2*448 = 928 */

#define NGCC_SIGN_SECRETKEYBYTES  (NGCC_SIGN_SEEDBYTES \
                                   + NGCC_SIGN_TRBYTES \
                                   + NGCC_SIGN_SEEDBYTES \
                                   + NGCC_SIGN_L * NGCC_SIGN_POLYETA_PACKEDBYTES \
                                   + NGCC_SIGN_M * NGCC_SIGN_POLYETA_PACKEDBYTES)
    /* 32 + 64 + 32 + 3*64 + 2*64 = 448 */

/* Hint encoding: per-coefficient 1 bit in each of m polys.
 * Reserved for future compressed signature format. */
#define NGCC_SIGN_POLYVECH_PACKEDBYTES  (NGCC_SIGN_M * (NGCC_SIGN_N / 8 + 1))
    /* 2 * 33 = 66 */

/* IRS sign bits: ceil(TAU/8) bytes encoding the per-monomial sign choices
 * from iterative rejection sampling. The verifier needs these to reconstruct
 * the effective challenge c_eff = sum(sigma_i * c_i * x^{j_i}). */
#define NGCC_SIGN_IRS_SIGNBYTES  ((NGCC_SIGN_TAU + 7) / 8)
    /* ceil(30/8) = 4 */

/* Signature format (Version 1 - full response):
 * c_tilde (32B) || irs_signs (4B) || polyz_pack(z[0..VECLEN-1]) (6*448B)
 * = 2724B
 *
 * This publishes the full VECLEN=6 response vector z plus the IRS sign bits.
 * The verifier reconstructs the effective challenge from (seed_c, irs_signs),
 * then checks w = B*z - c_eff*b has HighBits matching the commitment.
 * A future version will compress z2 for smaller signatures. */
#define NGCC_SIGN_BYTES  (NGCC_SIGN_CTILDEBYTES \
                          + NGCC_SIGN_IRS_SIGNBYTES \
                          + NGCC_SIGN_VECLEN * NGCC_SIGN_POLYZ_PACKEDBYTES)
    /* 32 + 4 + 6*448 = 2724 */

/* ============================================================
 * Gaussian sampler parameters
 * ============================================================ */
#define NGCC_SIGN_SIGMA     128
#define NGCC_SIGN_GAUSS_K   64     /* 2^6, convolution decomposition factor */
#define NGCC_SIGN_K_BITS    6      /* log2(GAUSS_K) */
#define NGCC_SIGN_TRUNC     11     /* truncation factor */
#define NGCC_SIGN_BOUND     (NGCC_SIGN_TRUNC * NGCC_SIGN_SIGMA)  /* 1408 */

#endif /* NGCC_SIGN_PARAMS_H */
