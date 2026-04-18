#ifndef SHUTTLE_ROUNDING_H
#define SHUTTLE_ROUNDING_H

#include <stdint.h>
#include "params.h"

#define decompose SHUTTLE_NAMESPACE(decompose)
int32_t decompose(int32_t *a0, int32_t a);

#define make_hint SHUTTLE_NAMESPACE(make_hint)
unsigned int make_hint(int32_t z2, int32_t comY);

#define use_hint SHUTTLE_NAMESPACE(use_hint)
int32_t use_hint(int32_t a, unsigned int hint);

#endif
