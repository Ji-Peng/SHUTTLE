#ifndef SHUTTLE_REDUCE_H
#define SHUTTLE_REDUCE_H

#include <stdint.h>
#include "params.h"

/* Montgomery / Barrett constants are now per-mode. They are provided by
 * params.h via ntt_constants.h (auto-generated). The base R = 2^32 is
 * uniform across all modes.
 *
 *   SHUTTLE_MONT      = 2^32 mod q
 *   SHUTTLE_QINV      = q^{-1} mod 2^32   (uint32_t)
 *   SHUTTLE_BARRETT_V = round(2^32 / q)   (for AVX2 ports)
 */

#define montgomery_reduce SHUTTLE_NAMESPACE(montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define reduce32 SHUTTLE_NAMESPACE(reduce32)
int32_t reduce32(int32_t a);

#define caddq SHUTTLE_NAMESPACE(caddq)
int32_t caddq(int32_t a);

#define freeze SHUTTLE_NAMESPACE(freeze)
int32_t freeze(int32_t a);

/* ============================================================
 * mod 2q helpers
 * ============================================================ */

#define caddq2 SHUTTLE_NAMESPACE(caddq2)
int32_t caddq2(int32_t a);

#define reduce_mod_2q SHUTTLE_NAMESPACE(reduce_mod_2q)
int32_t reduce_mod_2q(int32_t a);

#endif
