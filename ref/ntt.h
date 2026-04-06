#ifndef NGCC_SIGN_NTT_H
#define NGCC_SIGN_NTT_H

#include <stdint.h>
#include "params.h"

#define ntt NGCC_SIGN_NAMESPACE(ntt)
void ntt(int32_t a[NGCC_SIGN_N]);

#define invntt_tomont NGCC_SIGN_NAMESPACE(invntt_tomont)
void invntt_tomont(int32_t a[NGCC_SIGN_N]);

#endif
