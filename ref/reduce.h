#ifndef SHUTTLE_REDUCE_H
#define SHUTTLE_REDUCE_H

#include <stdint.h>
#include "params.h"

#define SHUTTLE_MONT  974          /* 2^32 mod q */
#define SHUTTLE_QINV  1309656065   /* q^(-1) mod 2^32 */

#define montgomery_reduce SHUTTLE_NAMESPACE(montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define reduce32 SHUTTLE_NAMESPACE(reduce32)
int32_t reduce32(int32_t a);

#define caddq SHUTTLE_NAMESPACE(caddq)
int32_t caddq(int32_t a);

#define freeze SHUTTLE_NAMESPACE(freeze)
int32_t freeze(int32_t a);

#endif
