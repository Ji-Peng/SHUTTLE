#ifndef SHUTTLE_CONFIG_H
#define SHUTTLE_CONFIG_H

/* ============================================================
 * SHUTTLE parameter-set selector.
 *
 * Valid values:
 *   128 - SHUTTLE-128 (n=256,  q=13313, sigma=101)
 *   256 - SHUTTLE-256 (n=512,  q=32257, sigma=149)
 *   512 - SHUTTLE-512 (n=1024, q=64513, sigma=202)
 *
 * Override at build time with: -DSHUTTLE_MODE=128|256|512
 * ============================================================ */
#ifndef SHUTTLE_MODE
#define SHUTTLE_MODE 128
#endif

#if SHUTTLE_MODE == 128
#  define CRYPTO_ALGNAME       "SHUTTLE-128"
#  define SHUTTLE_NAMESPACETOP shuttle128_ref
#  define SHUTTLE_NAMESPACE(s) shuttle128_ref_##s
#elif SHUTTLE_MODE == 256
#  define CRYPTO_ALGNAME       "SHUTTLE-256"
#  define SHUTTLE_NAMESPACETOP shuttle256_ref
#  define SHUTTLE_NAMESPACE(s) shuttle256_ref_##s
#elif SHUTTLE_MODE == 512
#  define CRYPTO_ALGNAME       "SHUTTLE-512"
#  define SHUTTLE_NAMESPACETOP shuttle512_ref
#  define SHUTTLE_NAMESPACE(s) shuttle512_ref_##s
#else
#  error "Unsupported SHUTTLE_MODE (expected 128, 256, or 512)"
#endif

/* ============================================================
 * Discrete Gaussian sampler standard deviation.
 *
 * SHUTTLE_SIGMA is derived from SHUTTLE_MODE. Override on the command
 * line (e.g. -DSHUTTLE_SIGMA=128) to exercise the legacy sigma=128
 * RCDT table kept in sampler.c for regression / audit.
 * ============================================================ */
#ifndef SHUTTLE_SIGMA
#  if SHUTTLE_MODE == 128
#    define SHUTTLE_SIGMA 101
#  elif SHUTTLE_MODE == 256
#    define SHUTTLE_SIGMA 149
#  elif SHUTTLE_MODE == 512
#    define SHUTTLE_SIGMA 202
#  endif
#endif

#endif
