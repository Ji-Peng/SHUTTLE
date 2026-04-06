# NGCC-Signature-Code Benchmark Report

This README records a local benchmark run for the NGCC sampler implementation in both [ref](./ref) and [avx2](./avx2). The goal is to compare the reference path against the AVX2-optimized path on the current machine and to preserve the exact commands and raw outputs used for the measurement.

## Test Environment

All results below were collected on April 6, 2026 on the following host:

- OS: Linux 6.17.0-1010-gcp x86_64 GNU/Linux
- CPU: AMD EPYC 9B45
- Topology: 1 socket, 8 cores, 16 hardware threads visible to the VM
- Cache summary: 384 KiB L1d, 256 KiB L1i, 8 MiB L2, 32 MiB L3
- Hypervisor: KVM
- Compiler: gcc 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)
- Relevant ISA flags: `avx2`, `bmi2`, `popcnt`, `aes`, `fma`, `avx512f`, `avx512bw`, `avx512vl`

Notes:

- The AVX2 build uses `-mavx2 -mbmi2 -mpopcnt -march=native -mtune=native`.
- Measurements are based on the built-in `cpucycles` timing harness, so the numbers are reported in CPU cycles/ticks.
- Because the machine is virtualized, the median is more reliable than the average when short-term scheduling noise is present.

## Benchmark Targets

The benchmark executable in each directory is `test/out/speed_sampler`. It covers the following components:

1. `sampler_sigma2`
2. `approx_exp`
3. SHAKE-256 throughput
4. `sample_gauss_N` for `N = 256`
5. `sample_gauss_N_4x` for four parallel samples of length 256

## Build And Run Commands

Reference implementation:

```bash
cd ref
make clean
make
./test/out/speed_sampler
```

AVX2 implementation:

```bash
cd avx2
make clean
make
./test/out/speed_sampler
```

The raw outputs from this benchmark session were also saved in the repository for traceability:

- `../agent/data/ngcc_ref_benchmark_20260406.txt`
- `../agent/data/ngcc_avx2_benchmark_20260406.txt`

## Raw Results

### Reference (`ref`)

```text
=== NGCC_SIGN Component Benchmarks (sigma=128) ===

sampler_sigma2 (BATCH=16):
median: 188 cycles/ticks
average: 187 cycles/ticks

approx_exp (V2, single call):
median: 26 cycles/ticks
average: 29 cycles/ticks

stream256_squeezeblocks (16 blocks = 2176 bytes):
median: 10988 cycles/ticks
average: 11051 cycles/ticks

stream256 avg per-byte: 5.02 cycles/byte (2176 bytes/squeeze)

sample_gauss_N (N=256):
median: 49976 cycles/ticks
average: 49477 cycles/ticks

sample_gauss_N_4x (4 x N=256):
median: 197207 cycles/ticks
average: 197349 cycles/ticks
```

### AVX2 (`avx2`)

```text
=== NGCC_SIGN Component Benchmarks - AVX2 (sigma=128) ===

sampler_sigma2 (AVX2, BATCH=16):
median: 53 cycles/ticks
average: 62 cycles/ticks

approx_exp (V2, single call):
median: 26 cycles/ticks
average: 28 cycles/ticks

shake256x4_squeezeblocks (16 blk x 4 lanes = 8704 bytes):
median: 16766 cycles/ticks
average: 16853 cycles/ticks

shake256x4 avg per-byte: 1.93 cycles/byte (8704 total bytes/squeeze)

stream256_squeezeblocks (16 blocks = 2176 bytes, single-lane ref):
median: 9125 cycles/ticks
average: 9164 cycles/ticks

stream256 (single-lane) avg per-byte: 4.18 cycles/byte (2176 bytes/squeeze)

sample_gauss_N (N=256):
median: 41174 cycles/ticks
average: 40720 cycles/ticks

sample_gauss_N_4x (4 x N=256):
median: 113939 cycles/ticks
average: 114416 cycles/ticks
```

## Summary Table

The table below compares the main medians from the two builds. For the SHAKE path, the per-byte cost is the most meaningful normalized metric.

| Benchmark | Reference | AVX2 | Speedup |
| --- | ---: | ---: | ---: |
| `sampler_sigma2` | 188 cycles | 53 cycles | 3.55x |
| `approx_exp` | 26 cycles | 26 cycles | 1.00x |
| SHAKE throughput | 5.02 cycles/byte | 1.93 cycles/byte | 2.60x |
| `sample_gauss_N` | 49,976 cycles | 41,174 cycles | 1.21x |
| `sample_gauss_N_4x` total cost | 197,207 cycles | 113,939 cycles | 1.73x |
| `sample_gauss_N_4x` normalized per sample | 49,302 cycles/sample | 28,485 cycles/sample | 1.73x |

## Interpretation

The main observations from this machine are:

1. The strongest improvement appears in `sampler_sigma2`, where the AVX2 path reduces the median cost from 188 cycles to 53 cycles. This is a 3.55x gain, which confirms that the batched Gaussian core benefits substantially from vectorization.
2. `approx_exp` is effectively unchanged. That is expected because the approximation routine remains scalar in both builds, so the optimized implementation does not materially change this component.
3. The SHAKE-based randomness generation shows a clear benefit from the four-lane Keccak implementation. The measured cost drops from 5.02 cycles/byte in the reference path to 1.93 cycles/byte in the AVX2 x4 path, for a 2.60x throughput improvement.
4. For a complete `sample_gauss_N` over 256 coefficients, the end-to-end speedup is more modest at 1.21x. This indicates that only part of the full sampler is accelerated by AVX2, while the remaining scalar work and control overhead still matter.
5. The four-way sampler is where the optimized design pays off most at the full-kernel level. `sample_gauss_N_4x` drops from 197,207 cycles to 113,939 cycles, which corresponds to 1.73x better total throughput, or about 28.5k cycles per sample instead of 49.3k cycles per sample.

## Reproducibility Notes

- Timing is based on a single local benchmark session on a virtualized cloud host.
- No explicit CPU pinning, turbo control, or frequency locking was applied.
- Both builds completed successfully, but `gcc` reported a warning in `randombytes.c` about the `syscall` declaration. This warning did not prevent the benchmark binaries from running.

For higher-confidence microbenchmarking, rerun the same commands on a pinned core and repeat the experiment multiple times before using the numbers in a paper.
