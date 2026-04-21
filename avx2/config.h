#ifndef SHUTTLE_CONFIG_H
#define SHUTTLE_CONFIG_H

/* ============================================================
 * SHUTTLE parameter-set selector (AVX2 implementation).
 *
 * Mirrors ref/config.h but uses the _avx2 namespace suffix so both
 * implementations can co-exist in the same process.
 *
 * Valid values:
 *   128 - SHUTTLE-128 (n=256, sigma=101)
 *   256 - SHUTTLE-256 (n=512, sigma=149)
 *   512 - SHUTTLE-512 (parameters TBD; currently #error in params.h)
 *
 * Override at build time with: -DSHUTTLE_MODE=128|256|512
 * ============================================================ */
#ifndef SHUTTLE_MODE
#define SHUTTLE_MODE 128
#endif

#if SHUTTLE_MODE == 128
#  define CRYPTO_ALGNAME       "SHUTTLE-128"
#  define SHUTTLE_NAMESPACETOP shuttle128_avx2
#  define SHUTTLE_NAMESPACE(s) shuttle128_avx2_##s
#elif SHUTTLE_MODE == 256
#  define CRYPTO_ALGNAME       "SHUTTLE-256"
#  define SHUTTLE_NAMESPACETOP shuttle256_avx2
#  define SHUTTLE_NAMESPACE(s) shuttle256_avx2_##s
#elif SHUTTLE_MODE == 512
#  define CRYPTO_ALGNAME       "SHUTTLE-512"
#  define SHUTTLE_NAMESPACETOP shuttle512_avx2
#  define SHUTTLE_NAMESPACE(s) shuttle512_avx2_##s
#else
#  error "Unsupported SHUTTLE_MODE (expected 128, 256, or 512)"
#endif

/* ============================================================
 * Discrete Gaussian sampler standard deviation.
 *
 * SHUTTLE_SIGMA is normally derived from SHUTTLE_MODE. It can still be
 * overridden on the command line (e.g. -DSHUTTLE_SIGMA=128) to exercise
 * the legacy sigma=128 RCDT table kept in sampler.c for regression/audit.
 * ============================================================ */
#ifndef SHUTTLE_SIGMA
#  if SHUTTLE_MODE == 128
#    define SHUTTLE_SIGMA 101
#  elif SHUTTLE_MODE == 256
#    define SHUTTLE_SIGMA 149
#  elif SHUTTLE_MODE == 512
#    error "SHUTTLE-512 sigma not yet specified (see NGCC-Signature Table 2)"
#  endif
#endif

#endif
