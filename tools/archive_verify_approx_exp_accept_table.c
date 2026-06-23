#include <math.h>
#include <quadmath.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "approx_exp_accept_table.h"

#define BLISS_R 825
#define BLISS_K 256
#define TARGET_BITS 53.0Q

static __float128 reference_exp(int x, int y)
{
    uint64_t num = (uint64_t)y * (uint64_t)(y + 2 * BLISS_K * x);
    uint64_t den = 2ULL * BLISS_R * BLISS_R;
    return expq(-((__float128)num) / ((__float128)den));
}

static __float128 table_probability(int x, int y)
{
    if (y == 0) {
        return 1.0Q;
    }
    uint64_t q = shuttle_exp_accept_q64(x, y);
    return ((__float128)q) / ldexpq(1.0Q, 64);
}

int main(void)
{
    __float128 max_rel = 0.0Q;
    int worst_x = 0;
    int worst_y = 0;
    for (int x = 0; x <= SHUTTLE_EXP_ACCEPT_X_MAX; x++) {
        for (int y = 0; y <= SHUTTLE_EXP_ACCEPT_Y_MAX; y++) {
            __float128 ref = reference_exp(x, y);
            __float128 got = table_probability(x, y);
            __float128 rel = fabsq(got - ref) / ref;
            if (rel > max_rel) {
                max_rel = rel;
                worst_x = x;
                worst_y = y;
            }
        }
    }
    __float128 bits = -logq(max_rel) / logq(2.0Q);
    char err_s[128];
    char bits_s[128];
    quadmath_snprintf(err_s, sizeof(err_s), "%.40Qe", max_rel);
    quadmath_snprintf(bits_s, sizeof(bits_s), "%.20Qf", bits);
    printf("max accept relative error = %s\n", err_s);
    printf("accept precision          = %s bits\n", bits_s);
    printf("worst point               = x=%d y=%d\n", worst_x, worst_y);
    printf("table size                = %u bytes\n",
           (unsigned)((SHUTTLE_EXP_ACCEPT_X_MAX + 1) *
                      (SHUTTLE_EXP_ACCEPT_Y_MAX + 1) * sizeof(uint64_t)));
    return bits >= TARGET_BITS ? EXIT_SUCCESS : EXIT_FAILURE;
}
