#include <float.h>
#include <math.h>
#include <quadmath.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define BLISS_R 825
#define BLISS_K 256
#define X_MAX 36
#define Y_MAX 255
#define TAYLOR_CUTOFF 0x1p-5Q

static const int64_t kChebCoeff[17] = {
    3012773964870532179LL,
    1317854137785297206LL,
    238888013567701743LL,
    36678453356540424LL,
    4861731890245871LL,
    565337493569314LL,
    58449486838614LL,
    5433446436407LL,
    458445434344LL,
    35391153625LL,
    2516952771LL,
    165883110LL,
    10183987LL,
    585046LL,
    31576LL,
    1607LL,
    77LL,
};

static int64_t high64_s64(int64_t a, int64_t b)
{
    return (int64_t)(((__int128)a * (__int128)b) >> 64);
}

static int64_t mul_q63_hi(int64_t a, int64_t b)
{
    return (int64_t)((__int128)high64_s64(a, b) * 2);
}

static int64_t mul_q61_to_q59_hi(int64_t a, int64_t b)
{
    return (int64_t)((__int128)high64_s64(a, b) * 2);
}

static __int128 mul_q59_q64_to_q64_hi(int64_t a, int64_t b)
{
    return ((__int128)high64_s64(a, b) * 32);
}

static uint64_t div_round_u64_scaled(uint64_t num, uint64_t den, int bits)
{
    uint64_t rem = num % den;
    __uint128_t q = (__uint128_t)(num / den) << bits;
    for (int i = bits - 1; i >= 0; i--) {
        rem <<= 1;
        if (rem >= den) {
            rem -= den;
            q |= ((__uint128_t)1) << i;
        }
    }
    rem <<= 1;
    if (rem >= den) {
        q++;
    }
    return (uint64_t)q;
}

static int64_t exp_arg_q61(int x, int y)
{
    uint64_t num = (uint64_t)y * (uint64_t)(y + 2 * BLISS_K * x);
    uint64_t den = 2ULL * BLISS_R * BLISS_R;
    return -(int64_t)div_round_u64_scaled(num, den, 61);
}

static int64_t z_q63_from_xy(int x, int y)
{
    uint64_t num = (uint64_t)y * (uint64_t)(y + 2 * BLISS_K * x);
    uint64_t den = (uint64_t)BLISS_R * BLISS_R * 351ULL;
    uint64_t delta = div_round_u64_scaled(num * 100ULL, den, 63);
    return (int64_t)((1ULL << 63) - delta);
}

static int64_t eval_r_q63(int64_t z)
{
    int64_t b1 = 0;
    int64_t b2 = 0;
    for (int i = 17; i >= 1; i--) {
        int64_t two_z_b1 = (int64_t)((__int128)mul_q63_hi(z, b1) * 2);
        int64_t b0 = kChebCoeff[i] + two_z_b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return kChebCoeff[0] + mul_q63_hi(z, b1) - b2;
}

static __float128 taylor_small_expq(__float128 p)
{
    __float128 term = 1.0Q;
    __float128 sum = 1.0Q;
    for (int i = 1; i <= 10; i++) {
        term *= p / (__float128)i;
        sum += term;
    }
    return sum;
}

static __float128 approx_exp_bliss(int x, int y)
{
    if (y == 0) {
        return 1.0Q;
    }
    uint64_t num = (uint64_t)y * (uint64_t)(y + 2 * BLISS_K * x);
    uint64_t den = 2ULL * BLISS_R * BLISS_R;
    __float128 p = -((__float128)num) / ((__float128)den);
    if (-p < TAYLOR_CUTOFF) {
        return taylor_small_expq(p);
    }
    int64_t p_q61 = exp_arg_q61(x, y);
    int64_t z_q63 = z_q63_from_xy(x, y);
    int64_t r_q63 = eval_r_q63(z_q63);
    int64_t p2_q59 = mul_q61_to_q59_hi(p_q61, p_q61);
    int64_t p2r_q63 = mul_q59_q64_to_q64_hi(p2_q59, r_q63);
    __int128 acc =
        ((__int128)1 << 63) + (((__int128)p_q61) << 2) + p2r_q63;
    return ((__float128)acc) / ((__float128)((uint64_t)1 << 63));
}

int main(void)
{
    __float128 max_rel = 0.0Q;
    __float128 max_rej = 0.0Q;
    int wr_x = 0, wr_y = 0, wj_x = 0, wj_y = 0;
    for (int x = 0; x <= X_MAX; x++) {
        for (int y = 0; y <= Y_MAX; y++) {
            uint64_t num = (uint64_t)y * (uint64_t)(y + 2 * BLISS_K * x);
            uint64_t den = 2ULL * BLISS_R * BLISS_R;
            __float128 p = -((__float128)num) / ((__float128)den);
            __float128 ref = expq(p);
            __float128 got = approx_exp_bliss(x, y);
            __float128 err = fabsq(got - ref);
            __float128 rel = err / ref;
            if (rel > max_rel) {
                max_rel = rel;
                wr_x = x;
                wr_y = y;
            }
            __float128 rej_den = fabsq(1.0Q - ref);
            if (rej_den > 0.0Q) {
                __float128 rej = err / rej_den;
                if (rej > max_rej) {
                    max_rej = rej;
                    wj_x = x;
                    wj_y = y;
                }
            }
        }
    }
    __float128 rel_bits = -logq(max_rel) / logq(2.0Q);
    __float128 rej_bits = -logq(max_rej) / logq(2.0Q);
    char rel_s[128], rej_s[128], relb_s[128], rejb_s[128];
    quadmath_snprintf(rel_s, sizeof(rel_s), "%.40Qe", max_rel);
    quadmath_snprintf(rej_s, sizeof(rej_s), "%.40Qe", max_rej);
    quadmath_snprintf(relb_s, sizeof(relb_s), "%.20Qf", rel_bits);
    quadmath_snprintf(rejb_s, sizeof(rejb_s), "%.20Qf", rej_bits);
    printf("max relative error     = %s\n", rel_s);
    printf("relative precision     = %s bits at x=%d y=%d\n", relb_s, wr_x,
           wr_y);
    printf("max rejection error    = %s\n", rej_s);
    printf("rejection precision    = %s bits at x=%d y=%d\n", rejb_s, wj_x,
           wj_y);
    char comb_s[128];
    quadmath_snprintf(comb_s, sizeof(comb_s), "%.20Qf",
                      rel_bits < rej_bits ? rel_bits : rej_bits);
    printf("combined precision     = %s bits\n", comb_s);
    return (rel_bits >= 52.0Q && rej_bits >= 52.0Q) ? EXIT_SUCCESS
                                                    : EXIT_FAILURE;
}
