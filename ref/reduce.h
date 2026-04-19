#ifndef SHUTTLE_REDUCE_H
#define SHUTTLE_REDUCE_H

#include <stdint.h>
#include "params.h"

/* Montgomery constants for q = 13313, base R = 2^32.
 *   MONT  = 2^32 mod q
 *   QINV  = q^{-1} mod 2^32  (used as uint32_t to ensure exact arithmetic)
 */
#define SHUTTLE_MONT  7114U
#define SHUTTLE_QINV  3398421505U

/* Barrett precomputation: V = round(2^32 / q) = 322585 (for q = 13313).
 * Used by reduce32 to compute a quotient estimate good to +/- 1. */
#define SHUTTLE_BARRETT_V  322585

#define montgomery_reduce SHUTTLE_NAMESPACE(montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define reduce32 SHUTTLE_NAMESPACE(reduce32)
int32_t reduce32(int32_t a);

#define caddq SHUTTLE_NAMESPACE(caddq)
int32_t caddq(int32_t a);

#define freeze SHUTTLE_NAMESPACE(freeze)
int32_t freeze(int32_t a);

#endif
