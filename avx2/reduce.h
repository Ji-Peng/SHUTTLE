#ifndef NGCC_SIGN_REDUCE_H
#define NGCC_SIGN_REDUCE_H

#include <stdint.h>
#include "params.h"

#define NGCC_SIGN_MONT  974          /* 2^32 mod q */
#define NGCC_SIGN_QINV  1309656065   /* q^(-1) mod 2^32 */

#define montgomery_reduce NGCC_SIGN_NAMESPACE(montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define reduce32 NGCC_SIGN_NAMESPACE(reduce32)
int32_t reduce32(int32_t a);

#define caddq NGCC_SIGN_NAMESPACE(caddq)
int32_t caddq(int32_t a);

#define freeze NGCC_SIGN_NAMESPACE(freeze)
int32_t freeze(int32_t a);

#endif
