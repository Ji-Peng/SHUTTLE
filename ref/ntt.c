#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"

static const int32_t zetas[SHUTTLE_N] = {
         0,    -6140,     5618,     3413,     6932,    -5375,    -3513,    -7157,
      2405,     3859,     7532,    -5530,     -484,     -387,    -6766,     3540,
     -1285,      972,     5205,    -7105,     5017,     -400,     3168,     5326,
      4823,    -2142,    -2083,    -1126,      607,    -3101,    -1861,     4193,
     -6109,     -854,    -6754,     5073,      744,     2880,     1768,    -4553,
      5213,    -6083,     5781,     5035,     1774,     3894,    -5034,    -5612,
     -1955,     -135,    -2003,    -6267,     6557,    -3358,     -440,     5234,
      2317,    -7383,      716,     -697,     -511,        4,    -1875,    -5276,
      1350,    -4189,     1226,    -4669,     2858,     4126,    -6257,    -4400,
     -2975,    -7552,     6970,     7160,      -40,    -5110,     6677,    -3389,
      -614,     6047,    -4267,    -3634,     7143,     2379,    -5324,    -4257,
      1922,     7440,     -553,     2319,    -4647,    -6096,      354,     6821,
      7493,    -6672,    -6144,    -1485,     3870,    -4840,    -4678,    -6216,
      7667,    -7485,    -5235,    -4408,    -3408,    -5264,    -5621,       44,
      6059,     2147,    -4101,    -5469,      288,     6070,     4153,    -3249,
      4000,     4087,    -7177,      958,    -5755,     5967,     2511,    -5641,
      1878,     1819,    -3961,    -2945,     3020,     1780,    -4881,     2413,
      -725,    -4293,     3893,    -5742,     -397,     6887,     1369,    -5602,
     -6862,     6637,     3349,     1567,     5226,     -582,     3995,     -392,
      6019,    -2963,     2576,     6503,     7258,    -2131,    -3399,     -274,
     -2821,     4441,     3537,    -5138,    -3833,    -1954,     1958,    -3322,
      -326,    -3244,     -114,    -6883,     3042,    -3090,    -3177,     5045,
      5220,     6332,    -3452,     4476,    -3286,     2641,     2432,     3468,
     -6383,     2545,     1349,     7204,    -5275,     5843,     -808,     4305,
     -1876,     6117,     1323,    -3798,     4312,    -2138,    -3958,    -6402,
      5272,    -2386,     4765,    -1871,     3014,    -6667,     3033,     7281,
      1718,    -3260,     7386,    -1140,      178,     -302,    -4367,    -1048,
      4179,    -7608,     2498,     4219,     5297,     4648,     2512,    -1673,
     -3153,    -7250,    -4024,    -7153,    -7426,    -3970,    -5424,    -1671,
     -2439,     -522,     6592,    -2727,    -1272,     6473,     3419,     2829,
     -1092,    -1254,    -3586,    -4962,     5949,    -4225,     2740,    -3268,
      5820,     6177,     3920,    -6133,     -309,     2768,    -7176,     4926
};

/*************************************************
* Name:        ntt
*
* Description: Forward NTT, in-place. No modular reduction is performed after
*              additions or subtractions. Output vector is in bitreversed order.
*
* Arguments:   - int32_t a[SHUTTLE_N]: input/output coefficient array
**************************************************/
void ntt(int32_t a[SHUTTLE_N]) {
  unsigned int len, start, j, k;
  int32_t zeta, t;

  k = 0;
  for(len = 128; len > 0; len >>= 1) {
    for(start = 0; start < SHUTTLE_N; start = j + len) {
      zeta = zetas[++k];
      for(j = start; j < start + len; ++j) {
        t = montgomery_reduce((int64_t)zeta * a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}

/*************************************************
* Name:        invntt_tomont
*
* Description: Inverse NTT and multiplication by Montgomery factor 2^32.
*              In-place. No modular reductions after additions or
*              subtractions; input coefficients need to be smaller than
*              Q in absolute value. Output coefficients are smaller than Q in
*              absolute value.
*
* Arguments:   - int32_t a[SHUTTLE_N]: input/output coefficient array
**************************************************/
void invntt_tomont(int32_t a[SHUTTLE_N]) {
  unsigned int start, len, j, k;
  int32_t t, zeta;
  const int32_t f = 7306; /* mont^2/256 mod q */

  k = 256;
  for(len = 1; len < SHUTTLE_N; len <<= 1) {
    for(start = 0; start < SHUTTLE_N; start = j + len) {
      zeta = -zetas[--k];
      for(j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
      }
    }
  }

  for(j = 0; j < SHUTTLE_N; ++j) {
    a[j] = montgomery_reduce((int64_t)f * a[j]);
  }
}
