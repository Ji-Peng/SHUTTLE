#ifndef PRINT_SPEED_H
#define PRINT_SPEED_H

#include <stddef.h>
#include <stdint.h>

uint64_t cycles_median(uint64_t *t, size_t tlen);
void print_results(const char *s, uint64_t *t, size_t tlen);

#endif
