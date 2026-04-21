# SHUTTLE: Implementation and Benchmark Report

This repository contains the reference (`ref/`) and AVX2-optimized (`avx2/`)
implementations of the **SHUTTLE** lattice-based post-quantum digital
signature scheme. Two parameter sets are supported:

- **SHUTTLE-128** (NIST category I target): `n=256, q=13313, sigma=101, tau=30`
- **SHUTTLE-256** (NIST category III target): `n=512, q=13313, sigma=149, tau=58`

Both modes share a single algorithmic source tree; mode selection is a
compile-time switch driven by `-DSHUTTLE_MODE={128,256}` (see `config.h`).
The AVX2 tree keeps only the SIMD-specific files as real sources
(`fips202x4.c`, `f1600x4.S`, `sampler.c`, `sampler_4x.c`); everything else is
a symlink to `ref/` so the two trees cannot drift.

## Test Environment

All results were collected on April 21, 2026 on the following host:

- **OS**: Linux 6.17.0-1012-gcp x86_64 GNU/Linux
- **CPU**: AMD EPYC 9B45 (Turin-class Zen 5)
- **Topology**: 1 socket, 8 cores, 16 hardware threads (KVM guest)
- **Cache**: 384 KiB L1d, 256 KiB L1i, 8 MiB L2, 32 MiB L3
- **Compiler**: gcc 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)
- **ISA flags observed**: `avx2`, `bmi2`, `popcnt`, `aes`, `fma`, `avx512f/bw/vl`

Build flags:

- **ref**: `-Wall -Wextra -Wpedantic -Wmissing-prototypes -O3 -std=c99 -fPIC`
- **avx2**: additionally `-mavx2 -mbmi2 -mpopcnt -march=native -mtune=native`

Timing is based on `rdtsc`, reported in CPU cycles/ticks. Because the machine
is virtualized, medians (1000 runs for sign, 10000 for sampler micro-benchmarks)
are more reliable than averages.

---

## 1. Scheme Parameters

| Parameter | SHUTTLE-128 | SHUTTLE-256 | Description |
|-----------|------------:|------------:|-------------|
| `n`       | 256   | 512   | Polynomial ring degree |
| `q`       | 13313 | 13313 | Modulus (prime, 14-bit; `q = 13 * 1024 + 1`) |
| `l`       | 3     | 3     | Number of polynomial components in **s** |
| `m`       | 2     | 2     | Number of polynomial components in **e** |
| `eta`     | 1     | 1     | Centered binomial distribution parameter |
| `tau`     | 30    | 58    | Hamming weight of challenge polynomial |
| `sigma`   | 101   | 149   | Discrete Gaussian standard deviation |
| `alpha_h` | 128   | 256   | Compression parameter for the hint |
| `alpha_1` | 8     | 16    | Compression parameter for response slot 0 |
| `B_k`     | 26    | 37    | Secret key norm upper bound |
| `B_s`     | 3893.66 | 9090.00 | Signing norm bound |
| `B_v`     | 4663.00 | 11202.00 | Verification norm bound |

## 2. Key and Signature Sizes

| Component   | SHUTTLE-128 (bytes) | SHUTTLE-256 (bytes) |
|-------------|--------------------:|--------------------:|
| Public key  | 928  | 1824 |
| Secret key  | 448  | 768  |
| Signature   | 1314 | 2578 |

The signature is a fixed-size container built from `seedC (32 B)`,
`irs_signs (ceil(tau/8) B)`, and three rANS-compressed streams
(`Z_0`, HighBits(`z[1..L]`), and the hint `h`), each preceded by a 2-byte
length prefix and zero-padded to a per-mode reserved budget:

| Stream  | SHUTTLE-128 reserved (B) | SHUTTLE-256 reserved (B) |
|---------|------------------------:|------------------------:|
| `Z_0`    | 202 | 360 |
| `z1_hi`  | 198 | 306 |
| `hint`   | 200 | 330 |

The per-block total rANS failure probability is targeted at `2^-20`
(overflow `2^-21` plus OOV `2^-21`), giving a combined per-signature
failure of about `3 * 2^-20 ~= 2^-18.4`, negligible compared to the IRS
rejection rate (~0.5).

---

## 3. Build and Run

Both trees use `make`, which emits per-mode binaries into
`test/out/` with a numeric suffix (`128` or `256`).

### Reference implementation

```bash
cd ref
make clean && make
# Correctness
./test/out/test_sign128
./test/out/test_sign256
# KeyGen/Sign/Verify throughput
./test/out/speed_sign128
./test/out/speed_sign256
# Per-component sign profiling
./test/out/profile_sign128
./test/out/profile_sign256
# Optional: sampler micro-benchmarks
./test/out/speed_sampler128
./test/out/speed_sampler256
```

The convenience target `make check` runs the full functional suite for both
modes (`test_sampler`, `test_ntt`, `test_mod2q`, `test_rans`,
`test_poly_z1`, `test_sig`, `test_sign`).

### AVX2 implementation

```bash
cd avx2
make clean && make
./test/out/test_sign128 && ./test/out/test_sign256
./test/out/speed_sign128 && ./test/out/speed_sign256
./test/out/profile_sign128 && ./test/out/profile_sign256
```

---

## 4. Correctness Tests

Both `ref` and `avx2` pass all five test suites for both modes (shared test
sources):

| Test | ref-128 | avx2-128 | ref-256 | avx2-256 |
|------|:-------:|:--------:|:-------:|:--------:|
| Key generation (norm bounds)             | PASS | PASS | PASS | PASS |
| Sign + Verify (single round)             | PASS | PASS | PASS | PASS |
| Forgery detection (bit flips, msg tamper) | PASS | PASS | PASS | PASS |
| Stress test (100 rounds, random msg len)  | PASS | PASS | PASS | PASS |
| Combined `crypto_sign` / `crypto_sign_open` | PASS | PASS | PASS | PASS |

The single-round Test 2 also prints the resulting signed-message length,
confirming `SHUTTLE_BYTES` is honored end-to-end: 1330 B for mode-128
(1314 sig + 16 B message) and 2594 B for mode-256 (2578 sig + 16 B).

---

## 5. KeyGen / Sign / Verify Benchmarks

Median / average of 1000 runs.

### 5.1 SHUTTLE-128

| Operation | ref (cycles) | avx2 (cycles) | Speedup |
|-----------|-------------:|--------------:|--------:|
| `crypto_sign_keypair`   | 103,031   | 94,904    | **1.09x** |
| `crypto_sign_signature` | 657,854   | 461,267   | **1.43x** |
| `crypto_sign_verify`    | 138,698   | 128,708   | **1.08x** |

Raw output - ref (mode 128):

```
crypto_sign_keypair:
  median:  103031 cycles/ticks
  average: 104984 cycles/ticks

crypto_sign_signature:
  median:  657854 cycles/ticks
  average: 659248 cycles/ticks

crypto_sign_verify:
  median:  138698 cycles/ticks
  average: 139462 cycles/ticks
```

Raw output - avx2 (mode 128):

```
crypto_sign_keypair:
  median:   94904 cycles/ticks
  average:  96541 cycles/ticks

crypto_sign_signature:
  median:  461267 cycles/ticks
  average: 462434 cycles/ticks

crypto_sign_verify:
  median:  128708 cycles/ticks
  average: 129454 cycles/ticks
```

### 5.2 SHUTTLE-256

| Operation | ref (cycles) | avx2 (cycles) | Speedup |
|-----------|-------------:|--------------:|--------:|
| `crypto_sign_keypair`   | 203,336   | 189,323   | **1.07x** |
| `crypto_sign_signature` | 1,646,108 | 1,014,983 | **1.62x** |
| `crypto_sign_verify`    | 275,453   | 262,196   | **1.05x** |

Raw output - ref (mode 256):

```
crypto_sign_keypair:
  median:   203336 cycles/ticks
  average:  205567 cycles/ticks

crypto_sign_signature:
  median:  1646108 cycles/ticks
  average: 1645813 cycles/ticks

crypto_sign_verify:
  median:   275453 cycles/ticks
  average:  277228 cycles/ticks
```

Raw output - avx2 (mode 256):

```
crypto_sign_keypair:
  median:   189323 cycles/ticks
  average:  190444 cycles/ticks

crypto_sign_signature:
  median:  1014983 cycles/ticks
  average: 1009500 cycles/ticks

crypto_sign_verify:
  median:   262196 cycles/ticks
  average:  263520 cycles/ticks
```

---

## 6. Sign Profiling: Component Breakdown

`profile_sign` instruments the signing flow at eight stages (precomputation,
Gaussian sampling, mod-2q commitment, challenge derivation, IRS, finalize,
pack plumbing, rANS encode) and reports cycle counts for each. Only the
**first successful** signing attempt is measured per run; IRS / norm / rANS
rejections trigger a retry whose overhead is excluded to isolate per-stage
cost.

### 6.1 SHUTTLE-128

Median of 1000 runs:

| # | Component | ref cycles | ref % | avx2 cycles | avx2 % | Speedup |
|:-:|-----------|-----------:|------:|------------:|-------:|--------:|
| 1 | Precomputation (unpack, ExpandA, NTT, mu, rhoprime, `b_hat`) | 97,605  | 15.1% | 91,638  | 20.0% | 1.07x |
| 2 | Gaussian sampling (`sample_gauss_N_4x` x2)                    | 314,982 | 48.8% | 214,406 | 46.9% | **1.47x** |
| 3 | Commitment (mod-2q NTT matrix-vector multiply)                | 46,332  |  7.2% | 45,441  |  9.9% | 1.02x |
| 4 | Challenge derivation (HighBits, SHAKE, SampleC)               | 8,127   |  1.3% | 7,290   |  1.6% | 1.11x |
| 5 | IRS (tau steps: ApproxLn, IntervalParity, inner product)      | 150,741 | 23.3% | 72,387  | 15.8% | **2.08x** |
| 6 | Finalize (CompressY z[0], norm check, MakeHint)               | 7,398   |  1.1% | 5,670   |  1.2% | 1.31x |
| 7 | Pack plumbing (headers, z1 lo split+pack, length, pad)        | 1,457   |  0.2% | 1,242   |  0.3% | 1.17x |
| 8 | rANS encode (Z_0 + z1_hi + hint)                              | 18,954  |  2.9% | 19,521  |  4.3% | 0.97x |
|   | **Sum of components**                                         | **645,596** | **100%** | **457,595** | **100%** | 1.41x |
|   | **Total (end-to-end)**                                        | **647,838** |      | **460,323** |      | 1.41x |

rANS sub-breakdown (median, avx2 shown; ref numbers are within noise):

| Stream | cycles | fraction of rANS total |
|--------|-------:|-----------------------:|
| Z_0    | 3,348  | 17.2% |
| z1_hi  | 9,693  | 49.7% |
| hint   | 6,480  | 33.2% |

Raw output - ref (mode 128):

```
Component                             median    average pct(med)
-------------------------------------------------------------------
1. Precomputation                      97605      98429    15.1%
2. Gaussian sampling (y)              314982     317213    48.8%
3. Commitment (mod-2q)                 46332      46789     7.2%
4. Challenge derivation                 8127       8199     1.3%
5. IRS (rejection sampling)           150741     152666    23.3%
6. Finalize (norm + MakeHint)           7398       7536     1.1%
7. Pack plumbing                        1457       1463     0.2%
8. rANS encode (total)                 18954      19035     2.9%
     Z_0  (shuttle_rans_encode_z0)      3267       3288   (17.2% of rANS)
     z1_hi(shuttle_rans_encode_z1)      9396       9446   (49.6% of rANS)
     hint (shuttle_rans_encode)         6291       6301   (33.2% of rANS)
-------------------------------------------------------------------
Sum of components                     645596     651330
Total (end-to-end)                    647838     651334
```

Raw output - avx2 (mode 128):

```
Component                             median    average pct(med)
-------------------------------------------------------------------
1. Precomputation                      91638      96718    20.0%
2. Gaussian sampling (y)              214406     219310    46.9%
3. Commitment (mod-2q)                 45441      47827     9.9%
4. Challenge derivation                 7290       7625     1.6%
5. IRS (rejection sampling)            72387      78635    15.8%
6. Finalize (norm + MakeHint)           5670       5938     1.2%
7. Pack plumbing                        1242       1298     0.3%
8. rANS encode (total)                 19521      19770     4.3%
     Z_0  (shuttle_rans_encode_z0)      3348       3382   (17.2% of rANS)
     z1_hi(shuttle_rans_encode_z1)      9693       9808   (49.7% of rANS)
     hint (shuttle_rans_encode)         6480       6579   (33.2% of rANS)
-------------------------------------------------------------------
Sum of components                     457595     477121
Total (end-to-end)                    460323     477126
```

### 6.2 SHUTTLE-256

Median of 1000 runs:

| # | Component | ref cycles | ref % | avx2 cycles | avx2 % | Speedup |
|:-:|-----------|-----------:|------:|------------:|-------:|--------:|
| 1 | Precomputation                                                | 200,151   | 12.2% | 189,702   | 17.4% | 1.06x |
| 2 | Gaussian sampling (`sample_gauss_N_4x` x2)                    | 604,476   | 36.8% | 380,484   | 34.9% | **1.59x** |
| 3 | Commitment (mod-2q)                                           | 94,851    |  5.8% | 97,848    |  9.0% | 0.97x |
| 4 | Challenge derivation                                          | 14,526    |  0.9% | 13,230    |  1.2% | 1.10x |
| 5 | IRS (tau=58 steps)                                            | 674,406   | 41.1% | 359,073   | 32.9% | **1.88x** |
| 6 | Finalize                                                      | 14,391    |  0.9% | 11,421    |  1.0% | 1.26x |
| 7 | Pack plumbing                                                 | 1,620     |  0.1% | 1,026     |  0.1% | 1.58x |
| 8 | rANS encode (Z_0 + z1_hi + hint)                              | 37,530    |  2.3% | 37,773    |  3.5% | 0.99x |
|   | **Sum of components**                                         | **1,641,951** | **100%** | **1,090,557** | **100%** | 1.51x |
|   | **Total (end-to-end)**                                        | **1,645,297** |      | **1,089,963** |      | 1.51x |

Mode-256 uses `tau=58` (vs 30 for mode-128), which is why IRS climbs from
23% to 41% of total cycles in the reference build. The AVX2 version
absorbs most of that growth and brings IRS back below Gaussian sampling.

Raw output - ref (mode 256):

```
Component                             median    average pct(med)
-------------------------------------------------------------------
1. Precomputation                     200151     201156    12.2%
2. Gaussian sampling (y)              604476     605650    36.8%
3. Commitment (mod-2q)                 94851      95348     5.8%
4. Challenge derivation                14526      14651     0.9%
5. IRS (rejection sampling)           674406     673884    41.1%
6. Finalize (norm + MakeHint)          14391      14486     0.9%
7. Pack plumbing                        1620       1648     0.1%
8. rANS encode (total)                 37530      37719     2.3%
     Z_0  (shuttle_rans_encode_z0)      6453       6512   (17.2% of rANS)
     z1_hi(shuttle_rans_encode_z1)     18414      18541   (49.1% of rANS)
     hint (shuttle_rans_encode)        12663      12666   (33.7% of rANS)
-------------------------------------------------------------------
Sum of components                    1641951    1644542
Total (end-to-end)                   1645297    1644546
```

Raw output - avx2 (mode 256):

```
Component                             median    average pct(med)
-------------------------------------------------------------------
1. Precomputation                     189702     192625    17.4%
2. Gaussian sampling (y)              380484     383161    34.9%
3. Commitment (mod-2q)                 97848      99379     9.0%
4. Challenge derivation                13230      13471     1.2%
5. IRS (rejection sampling)           359073     328371    32.9%
6. Finalize (norm + MakeHint)          11421      11591     1.0%
7. Pack plumbing                        1026       1036     0.1%
8. rANS encode (total)                 37773      37973     3.5%
     Z_0  (shuttle_rans_encode_z0)      6453       6497   (17.1% of rANS)
     z1_hi(shuttle_rans_encode_z1)     18603      18675   (49.2% of rANS)
     hint (shuttle_rans_encode)        12744      12800   (33.7% of rANS)
-------------------------------------------------------------------
Sum of components                    1090557    1067607
Total (end-to-end)                   1089963    1067610
```

---

## 7. Analysis and Optimization Opportunities

### 7.1 Gaussian sampling is the #1 bottleneck in mode-128 (~47%)

In both `ref` and `avx2`, the discrete Gaussian sampler
(`sample_gauss_N_4x` called twice for 6 polynomials of `n` coefficients
each) accounts for roughly half of total signing cost in mode-128 and
about 35% in mode-256. The AVX2 version achieves 1.47x / 1.59x speedup
here through 4-lane parallel SHAKE-256 (`fips202x4` with AVX2 Keccak
assembly) plus a vectorized CDT comparison loop, but this remains the
single largest bucket for mode-128.

Optimization paths:

- AVX2 vectorization of `sampler_sigma2` (already done; 2.5x-2.7x speedup
  on that kernel alone — see section 8).
- Reducing SHAKE squeeze overhead through larger initial buffer
  allocation.
- Alternative sampler designs (e.g., Knuth-Yao, lazy sampling).

### 7.2 IRS dominates mode-256 signing (~33%-41%)

The iterative rejection sampling executes `tau` steps (30 for mode-128, 58
for mode-256), each involving:

- Cyclic shift of the secret key vector (`VECLEN * n` coefficients).
- Inner product computation (~`VECLEN * n` multiply-accumulate, i.e. 1536
  MACs for mode-128 and 3072 for mode-256).
- `approx_neg_ln` call (~30 cycles).
- IntervalParity constant-time check (100 addition+comparison steps).
- Branchless sign decision and vector update.

The 2.08x (mode-128) / 1.88x (mode-256) AVX2 speedup in IRS comes mainly
from faster SHAKE squeeze for random byte generation and general
`-march=native` compiler benefits; the inner-product kernel itself is
still scalar.

Optimization paths:

- AVX2 vectorization of the inner product (8-wide int32 MAC with
  `_mm256_mullo_epi32`).
- AVX2 vectorization of the cyclic shift.
- Precomputing all `tau` cyclic shifts before the IRS loop.

### 7.3 Commitment and challenge are already cheap

The NTT-based mod-2q matrix-vector multiply costs only 6%-10% of total
signing. Since `q = 13313` is 14-bit (vs Dilithium's 23-bit
`q = 8380417`), each coefficient fits comfortably in 16 bits and
arithmetic has no overflow issues even without frequent reductions.
Challenge derivation adds another 1%-2%; both are far below the sampling
and IRS buckets and warrant no dedicated AVX2 effort today.

Potential future path: a Dilithium-style AVX2 assembly NTT for `q =
13313` (requires regenerating the twiddle-factor tables in Dilithium's
AVX2 register layout).

### 7.4 rANS encode is near-free (<5%)

The three rANS streams (Z_0, HighBits z[1..L], hint) together cost only
~3% of total signing in both modes, and AVX2 is within noise of ref
because the encoder is already a tight scalar loop over a small table.
No optimization effort is warranted here.

### 7.5 Other components are negligible

Finalize (~1%) and pack plumbing (<0.3%) together account for less than
2% of total cost. No optimization effort is warranted.

---

## 8. Sampler Benchmarks (Micro-level)

Collected from `./test/out/speed_sampler{128,256}` (10000 iterations each).

### 8.1 SHUTTLE-128 (`sigma = 101`, `N = 256`)

| Benchmark | ref | avx2 | Speedup |
|-----------|----:|-----:|--------:|
| `sampler_sigma2` (BATCH=16)   | 134 cycles | 53 cycles | **2.53x** |
| `approx_exp` (single call)    | 26 cycles  | 26 cycles | 1.00x |
| SHAKE256 throughput (1 lane)  | 5.00 cy/byte | 4.16 cy/byte | 1.20x |
| SHAKE256x4 throughput (AVX2)  | -          | 1.93 cy/byte | **2.60x vs 1-lane ref** |
| `sample_gauss_N` (N=256)      | 53,351 cycles | 43,955 cycles | **1.21x** |
| `sample_gauss_N_4x` (4 x 256) | 209,465 cycles | 122,282 cycles | **1.71x** |
| Per-sample (from 4x)          | ~52,366 cy | ~30,571 cy | 1.71x |

### 8.2 SHUTTLE-256 (`sigma = 149`, `N = 512`)

| Benchmark | ref | avx2 | Speedup |
|-----------|----:|-----:|--------:|
| `sampler_sigma2` (BATCH=16)   | 215 cycles | 80 cycles | **2.69x** |
| `approx_exp` (single call)    | 26 cycles  | 26 cycles | 1.00x |
| SHAKE256 throughput (1 lane)  | 5.00 cy/byte | 4.16 cy/byte | 1.20x |
| SHAKE256x4 throughput (AVX2)  | -          | 1.93 cy/byte | **2.60x vs 1-lane ref** |
| `sample_gauss_N` (N=512)      | 99,008 cycles | 80,729 cycles | **1.23x** |
| `sample_gauss_N_4x` (4 x 512) | 394,604 cycles | 218,537 cycles | **1.81x** |
| Per-sample (from 4x)          | ~98,651 cy | ~54,634 cy | 1.81x |

Mode-256 doubles `N` (from 256 to 512) and uses a larger `sigma` (149 vs
101), so per-poly costs roughly double; the AVX2 speedup ratio is nearly
identical across the two modes because both are driven by the same
`fips202x4` and CDT kernels.

---

## 9. Architecture Overview

### Signing flow (per `sign.c::crypto_sign_signature`)

```
Sign(sk, msg):
  1. Precomputation: decode sk, build sk_full, compute mu, rhoprime;
                      ExpandA and NTT(A); compute b = A*s + e, NTT(b).
  2. [Loop until IRS accepts, norm < B_s, rANS OOV not hit]
     a. Gaussian sampling: y <- D_{sigma}^{(1+L+M) * N} via sample_gauss_N_4x
                            (two 4-way calls; last two lanes idle).
     b. Commitment:        Y = y with y[0] stretched by alpha_1;
                            w = mat-vec multiply over modulus 2q (NTT domain).
     c. Challenge:         w_h = HighBits_{2q}(w), w_0 = LSB bitmap,
                            seed_c = SHAKE256(w_h || w_0 || mu || rho);
                            c = SampleC(seed_c) with weight tau.
     d. IRS:               z = y; run tau log-domain rejection steps over
                            sk_full, updating z with +-c_i * x^{j_i} * sk.
     e. Finalize:          CompressY z[0], polyvec_sq_norm <? B_s^2,
                            h = MakeHint_{2q}(w, z_2).
     f. Pack:              sig = seed_c || irs_signs || rANS(Z_0)
                            || polyz1_lo(z[1..L]) || rANS(z1_hi)
                            || rANS(h), each stream length-prefixed and
                            zero-padded to the reserved budget.
```

### Key design choices

- **IRS method**: simplified log-domain rejection (0 exp, 1 ln per step),
  constant-time 100-step IntervalParity, branchless sign update.
- **NTT**: reference Cooley-Tukey for both ref and avx2 (`q = 13313` is
  small enough that assembly NTT is not on the critical path).
- **rANS-compressed signature**: fixed-size container with per-stream
  reserved budgets calibrated by `SHUTTLE/tools/calibrate_rans.py`
  (5000 / 2500 trials for mode-128 / mode-256). Per-block failure target
  `p_block ~= 2^-20`.
- **AVX2 advantage**: primarily from 4-lane SHAKE-256 (`f1600x4.S` Keccak
  assembly, `fips202x4.c`) and vectorized CDT sampler; NTT and packing
  code are shared with `ref/` via symlinks.

---

## 10. Reproducibility

```bash
# Reference
cd ref
make clean && make -j
./test/out/test_sign128      && ./test/out/test_sign256
./test/out/speed_sign128     && ./test/out/speed_sign256
./test/out/profile_sign128   && ./test/out/profile_sign256
./test/out/speed_sampler128  && ./test/out/speed_sampler256

# AVX2
cd ../avx2
make clean && make -j
./test/out/test_sign128      && ./test/out/test_sign256
./test/out/speed_sign128     && ./test/out/speed_sign256
./test/out/profile_sign128   && ./test/out/profile_sign256
./test/out/speed_sampler128  && ./test/out/speed_sampler256
```

All raw data in this report was produced by the commands above on
2026-04-21 with the environment listed in the "Test Environment" section.
