/*
 * params.h - Parameters for NGCC_SIGN discrete Gaussian sampler.
 *
 * Gaussian sampler: sigma = 128 = 2^7.
 * Decomposition: sigma = k * sigma_2, k = 2^6, sigma_2 = 2.
 * Truncation: 11 * sigma = 1408.
 */

#ifndef NGCC_SIGN_PARAMS_H
#define NGCC_SIGN_PARAMS_H

#define NGCC_SIGN_N 256  /* Polynomial ring degree */

/* Seed bytes for SHAKE-256 stream */
#define NGCC_SIGN_SEEDBYTES 32

/* Gaussian parameters */
#define NGCC_SIGN_SIGMA     128
#define NGCC_SIGN_K         64     /* 2^6 */
#define NGCC_SIGN_K_BITS    6      /* log2(k) */
#define NGCC_SIGN_TRUNC     11     /* truncation factor */
#define NGCC_SIGN_BOUND     (NGCC_SIGN_TRUNC * NGCC_SIGN_SIGMA)  /* 1408 */

#endif /* NGCC_SIGN_PARAMS_H */
