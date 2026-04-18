# SHUTTLE: Implementation and Benchmark Report

This repository contains the reference (`ref/`) and AVX2-optimized (`avx2/`) implementations of the **SHUTTLE** lattice-based post-quantum digital signature scheme.

## Test Environment

All results were collected on April 6, 2026 on the following host:

- **OS**: Linux 6.17.0-1010-gcp x86_64 GNU/Linux
- **CPU**: AMD EPYC 9B45
- **Topology**: 1 socket, 8 cores, 16 hardware threads (KVM guest)
- **Cache**: 384 KiB L1d, 256 KiB L1i, 8 MiB L2, 32 MiB L3
- **Compiler**: gcc 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)
- **ISA flags**: `avx2`, `bmi2`, `popcnt`, `aes`, `fma`, `avx512f/bw/vl`

Build flags:
- **ref**: `-Wall -Wextra -Wpedantic -O3 -std=c99`
- **avx2**: additionally `-mavx2 -mbmi2 -mpopcnt -march=native -mtune=native`

Timing is based on `rdtsc`, reported in CPU cycles/ticks. Because the machine is virtualized, medians (1000 runs) are more reliable than averages.

---

## 1. Scheme Parameters

| Parameter | Value | Description |
|-----------|------:|-------------|
| `n` | 256 | Polynomial ring degree |
| `q` | 15361 | Modulus (prime, 14-bit) |
| `l` | 3 | Number of polynomial components in **s** |
| `m` | 2 | Number of polynomial components in **e** |
| `η` | 1 | Centered binomial distribution parameter |
| `τ` | 30 | Hamming weight of challenge polynomial |
| `σ` | 128 | Discrete Gaussian standard deviation |
| `α_h` | 128 | Compression parameter for hint |
| `B_k` | 26 | Secret key norm upper bound |
| `B_s` | 4233.70 | Signing norm bound |
| `B_v` | 5000 | Verification norm bound |

## 2. Key and Signature Sizes

| Component | Size (bytes) |
|-----------|------------:|
| Public key | 928 |
| Secret key | 448 |
| Signature | 2724 |

---

## 3. Build and Run

### Reference implementation

```bash
cd ref
# Build all targets
gcc -O3 -std=c99 -o test/test_sign test/test_sign.c sign.c packing.c \
    rejsample.c approx_log.c poly.c polyvec.c ntt.c reduce.c rounding.c \
    sampler.c sampler_4x.c approx_exp.c fips202.c symmetric-shake.c randombytes.c

gcc -O3 -std=c99 -o test/speed_sign test/speed_sign.c test/cpucycles.c \
    test/speed_print.c sign.c packing.c rejsample.c approx_log.c poly.c \
    polyvec.c ntt.c reduce.c rounding.c sampler.c sampler_4x.c approx_exp.c \
    fips202.c symmetric-shake.c randombytes.c

gcc -O3 -std=c99 -o test/profile_sign test/profile_sign.c test/cpucycles.c \
    test/speed_print.c sign.c packing.c rejsample.c approx_log.c poly.c \
    polyvec.c ntt.c reduce.c rounding.c sampler.c sampler_4x.c approx_exp.c \
    fips202.c symmetric-shake.c randombytes.c

# Run
./test/test_sign       # Correctness tests
./test/speed_sign      # Overall KeyGen/Sign/Verify benchmarks
./test/profile_sign    # Sign component profiling
```

### AVX2 implementation

```bash
cd avx2
make clean && make
./test/out/test_sign       # Correctness tests
./test/out/speed_sign      # Overall benchmarks
./test/out/profile_sign    # Sign component profiling
```

---

## 4. Correctness Tests

Both `ref` and `avx2` pass all five test suites (identical test code):

| Test | ref | avx2 |
|------|:---:|:----:|
| Key generation (norm bounds) | PASS | PASS |
| Sign + Verify (single round) | PASS | PASS |
| Forgery detection (bit flips) | PASS | PASS |
| Stress test (100 rounds) | PASS | PASS |
| Combined `crypto_sign`/`crypto_sign_open` | PASS | PASS |

---

## 5. KeyGen / Sign / Verify Benchmarks

Median of 1000 runs:

| Operation | ref (cycles) | avx2 (cycles) | Speedup |
|-----------|------------:|-------------:|--------:|
| `crypto_sign_keypair` | 116,072 | 104,840 | **1.11×** |
| `crypto_sign_signature` | 563,894 | 382,103 | **1.48×** |
| `crypto_sign_verify` | 141,776 | 130,085 | **1.09×** |

### Raw output — ref

```
crypto_sign_keypair:
  median: 116072 cycles/ticks
  average: 118117 cycles/ticks

crypto_sign_signature:
  median: 563894 cycles/ticks
  average: 567231 cycles/ticks

crypto_sign_verify:
  median: 141776 cycles/ticks
  average: 142775 cycles/ticks
```

### Raw output — avx2

```
crypto_sign_keypair:
  median: 104840 cycles/ticks
  average: 106811 cycles/ticks

crypto_sign_signature:
  median: 382103 cycles/ticks
  average: 383959 cycles/ticks

crypto_sign_verify:
  median: 130085 cycles/ticks
  average: 130969 cycles/ticks
```

---

## 6. Sign Profiling: Component Breakdown

The profiler instruments the signing flow at six stages and reports cycle counts for each component. Only the first **successful** signing attempt is measured per run (IRS rejections trigger a retry; the retry overhead is excluded to isolate per-component cost).

### 6.1 ref — Component Breakdown

Median of 1000 runs:

| # | Component | Cycles | % |
|:-:|-----------|-------:|--:|
| 1 | **Precomputation** (sk decode, ExpandA, NTT, mu, rhoprime) | 61,857 | 10.9% |
| 2 | **Gaussian sampling** (`sample_gauss_N_4x` ×2, 6 polys) | 293,841 | **51.9%** |
| 3 | **Commitment** (NTT, matrix-vector multiply, INTT) | 41,877 | 7.4% |
| 4 | **Challenge derivation** (decompose, hash, SampleC) | 6,507 | 1.1% |
| 5 | **IRS** (30 steps: ApproxLn, IntervalParity, inner product) | 158,841 | **28.1%** |
| 6 | **Norm check + pack sig** | 3,132 | 0.6% |
| | **Sum of components** | **566,055** | **100%** |
| | **Total (end-to-end)** | **567,918** | |

```
Component                        median    average pct(med)
--------------------------------------------------------------
1. Precomputation                 61857      62427    10.9%
   (sk decode, ExpandA,
    NTT, mu, rhoprime)
2. Gaussian sampling (y)         293841     295865    51.9%
   (sample_gauss_N_4x x2)
3. Commitment (B*y)               41877      42212     7.4%
   (NTT, mat-vec mul, INTT)
4. Challenge derivation            6507       6575     1.1%
   (decompose, hash, SampleC)
5. IRS (rejection sampling)      158841     159456    28.1%
   (30 steps: ApproxLn,
    IntervalParity, inner prod)
6. Norm check + pack sig           3132       3166     0.6%
--------------------------------------------------------------
Sum of components                566055     569701
Total (end-to-end)               567918     569743
```

### 6.2 avx2 — Component Breakdown

Median of 1000 runs:

| # | Component | Cycles | % |
|:-:|-----------|-------:|--:|
| 1 | **Precomputation** | 53,433 | 13.9% |
| 2 | **Gaussian sampling** | 202,905 | **52.9%** |
| 3 | **Commitment** | 40,905 | 10.7% |
| 4 | **Challenge derivation** | 6,021 | 1.6% |
| 5 | **IRS** | 78,408 | **20.4%** |
| 6 | **Norm check + pack sig** | 2,241 | 0.6% |
| | **Sum of components** | **383,913** | **100%** |
| | **Total (end-to-end)** | **383,454** | |

```
Component                        median    average pct(med)
--------------------------------------------------------------
1. Precomputation                 53433      53749    13.9%
   (sk decode, ExpandA,
    NTT, mu, rhoprime)
2. Gaussian sampling (y)         202905     204286    52.9%
   (sample_gauss_N_4x x2)
3. Commitment (B*y)               40905      41196    10.7%
   (NTT, mat-vec mul, INTT)
4. Challenge derivation            6021       6056     1.6%
   (decompose, hash, SampleC)
5. IRS (rejection sampling)       78408      77526    20.4%
   (30 steps: ApproxLn,
    IntervalParity, inner prod)
6. Norm check + pack sig           2241       2257     0.6%
--------------------------------------------------------------
Sum of components                383913     385070
Total (end-to-end)               383454     385102
```

### 6.3 Component-Level Comparison (ref vs avx2)

| Component | ref (cycles) | avx2 (cycles) | Speedup | Notes |
|-----------|------------:|-------------:|--------:|-------|
| Precomputation | 61,857 | 53,433 | 1.16× | ExpandA dominates (8 poly_uniform calls) |
| Gaussian sampling | 293,841 | 202,905 | **1.45×** | SHAKE4x parallelism in sampler_4x |
| Commitment | 41,877 | 40,905 | 1.02× | Same ref NTT used in both |
| Challenge derivation | 6,507 | 6,021 | 1.08× | Mostly SHAKE hash |
| IRS | 158,841 | 78,408 | **2.03×** | Inner product + SHAKE squeeze benefit |
| Norm check + pack | 3,132 | 2,241 | 1.40× | Marginal |
| **Total** | **567,918** | **383,454** | **1.48×** | |

---

## 7. Analysis and Optimization Opportunities

### 7.1 Signing is dominated by Gaussian sampling (~52%)

In both ref and avx2, the discrete Gaussian sampler (`sample_gauss_N_4x` called twice for 6 polynomials × 256 coefficients) accounts for over half of the total signing cost. The AVX2 version achieves a 1.45× speedup here through 4-lane parallel SHAKE-256 (`fips202x4` with AVX2 Keccak assembly), but this remains the primary bottleneck.

**Optimization paths**:
- AVX2 vectorization of the CDT comparison loop (already implemented, providing 3.55× on `sampler_sigma2`)
- Reducing SHAKE squeeze overhead through larger initial buffer allocation
- Alternative sampler designs (e.g., Knuth-Yao, lazy sampling)

### 7.2 IRS is the second bottleneck (~20–28%)

The iterative rejection sampling executes 30 steps (one per challenge monomial), each involving:
- Cyclic shift of the secret key vector (6 × 256 coefficients)
- Inner product computation (6 × 256 multiply-accumulate → 1536 MACs)
- `approx_neg_ln` call (~30 cycles)
- IntervalParity (100 addition+comparison steps)
- Branchless sign decision and vector update

The 2.03× AVX2 speedup in IRS comes mainly from faster SHAKE squeeze for random byte generation and general compilation benefits (`-march=native`).

**Optimization paths**:
- AVX2 vectorization of the inner product (8-wide int32 MAC with `_mm256_mullo_epi32`)
- AVX2 vectorization of the cyclic shift
- Precomputing all 30 cyclic shifts before the IRS loop

### 7.3 Commitment computation is already fast (~7–11%)

The NTT-based matrix-vector multiplication uses a reference Cooley-Tukey NTT for `q = 15361`. Since the modulus is only 14 bits (vs Dilithium's 23-bit `q = 8380417`), each coefficient fits comfortably in 16 bits and the arithmetic has no overflow issues even without frequent reductions.

**Optimization path**:
- AVX2 assembly NTT adapted for `q = 15361` (requires regenerating the twiddle factor tables in Dilithium's AVX2 register layout)

### 7.4 Other components are negligible

Challenge derivation (1.1–1.6%) and signature packing (0.6%) together account for less than 3% of total cost. No optimization effort is warranted here.

---

## 8. Sampler Benchmarks (Micro-level)

These were collected separately from the sampler-only benchmark suite.

| Benchmark | ref | avx2 | Speedup |
|-----------|----:|-----:|--------:|
| `sampler_sigma2` (BATCH=16) | 188 cycles | 53 cycles | **3.55×** |
| `approx_exp` (single call) | 26 cycles | 26 cycles | 1.00× |
| SHAKE throughput | 5.02 cy/byte | 1.93 cy/byte | **2.60×** |
| `sample_gauss_N` (N=256) | 49,976 cycles | 41,174 cycles | **1.21×** |
| `sample_gauss_N_4x` (4×N=256) | 197,207 cycles | 113,939 cycles | **1.73×** |
| Per-sample (from 4x) | 49,302 cy/sample | 28,485 cy/sample | **1.73×** |

---

## 9. Architecture Overview

### Signing flow

```
Sign(sk, msg):
  1. Precomputation: decode sk, ExpandA, NTT(A), compute mu, derive rhoprime
  2. [Loop]
     a. Gaussian sampling: y ← D_{σ=128}^{6×256} via sample_gauss_N_4x
     b. Commitment: w = B·y (NTT domain matrix-vector multiply)
     c. Challenge: seed_c = H(HighBits(w) ‖ mu), c = SampleC(seed_c)
     d. IRS: z = y ± Σ(c_i · x^{j_i} · sk)  (30 log-domain rejection steps)
     e. Norm check: ‖z₁‖² < Bs²
     f. Pack: σ = (seed_c, irs_signs, z)
```

### Key design choices

- **IRS method**: Simplified log-domain (0 exp, 1 ln per step), constant-time 100-step IntervalParity
- **NTT**: Reference Cooley-Tukey for both ref and avx2 (q = 15361 is small enough that assembly NTT is not the bottleneck)
- **AVX2 advantage**: Primarily from 4-lane SHAKE-256 (`f1600x4.S`) and vectorized CDT sampler

---

## 10. Reproducibility

```bash
# Run all benchmarks
cd ref
./test/test_sign          # Correctness
./test/speed_sign         # KeyGen/Sign/Verify
./test/profile_sign       # Sign profiling

cd ../avx2
make clean && make
./test/out/test_sign
./test/out/speed_sign
./test/out/profile_sign
```

Raw profiling outputs are saved in:
- `agent/data/shuttle_sign_profile_ref.txt`
- `agent/data/shuttle_sign_profile_avx2.txt`
- `agent/data/shuttle_ref_benchmark_20260406.txt`
- `agent/data/shuttle_avx2_benchmark_20260406.txt`
