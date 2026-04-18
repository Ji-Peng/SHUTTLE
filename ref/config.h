#ifndef SHUTTLE_CONFIG_H
#define SHUTTLE_CONFIG_H

#define SHUTTLE_NAMESPACE(s) shuttle_ref_##s

/* ============================================================
 * Discrete Gaussian sampler mode selector.
 *
 * Valid values:
 *   101 - SHUTTLE-128 parameter set (default)
 *   149 - SHUTTLE-256 parameter set
 *   128 - legacy / reference variant (sigma = 2^7, kept for regression)
 *
 * Override at build time with: -DSHUTTLE_SIGMA=128|101|149
 * ============================================================ */
#ifndef SHUTTLE_SIGMA
#define SHUTTLE_SIGMA 101
#endif

#endif
