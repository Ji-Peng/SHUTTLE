#ifndef SHUTTLE_NTT_H
#define SHUTTLE_NTT_H

#include <stdint.h>
#include "params.h"

#define ntt SHUTTLE_NAMESPACE(ntt)
void ntt(int32_t a[SHUTTLE_N]);

#define invntt_tomont SHUTTLE_NAMESPACE(invntt_tomont)
void invntt_tomont(int32_t a[SHUTTLE_N]);

#endif
