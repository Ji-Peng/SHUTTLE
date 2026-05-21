#ifndef SHUTTLE_NTT_H
#define SHUTTLE_NTT_H

#include <stdint.h>
#include "params.h"

#define ntt SHUTTLE_NAMESPACE(ntt)
void ntt(int32_t a[SHUTTLE_N]);

#define invntt_tomont SHUTTLE_NAMESPACE(invntt_tomont)
void invntt_tomont(int32_t a[SHUTTLE_N]);

#if SHUTTLE_BASE_DEG == 2
/* Kyber-style basemul over all N/2 deg-2 basecases. */
#define poly_basemul_montgomery_native SHUTTLE_NAMESPACE(poly_basemul_montgomery_native)
void poly_basemul_montgomery_native(int32_t r[SHUTTLE_N],
                                    const int32_t a[SHUTTLE_N],
                                    const int32_t b[SHUTTLE_N]);
#endif

#endif
