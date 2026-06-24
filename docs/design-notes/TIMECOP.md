# TIMECOP / Valgrind Variable-Latency Notes (SHUTTLE-NGCC)

This note records the TIMECOP setup for SHUTTLE and its methodology. Patched Valgrind is **dynamic** confirmation of the KyberSlash variable-latency property; it complements, and does NOT replace, the **static** object-code constant-time scan (`ct_scan_matrix.sh` -> `tools/ct_scan.py`) and the statistical `dudect-sign` harness. SHUTTLE port of `Lithium-Code/docs/design-notes/TIMECOP.md`.

The `ref` and `avx2` backends run under Valgrind; the `avx512` backend cannot (Valgrind 3.23 SIGILLs on EVEX/AVX-512/VBMI2) and is covered statically instead.

## Tools (`tools/timecop/`)

- One-shot driver (build-if-needed + positive control + smoke per runnable backend x mode): `tools/timecop/run_all.sh`
- Build patched Valgrind: `tools/timecop/build_valgrind_varlat.sh`
- Positive-control secret-division detector: `tools/timecop/prove_varlat.sh --valgrind <vg>`
- Per-backend signing smoke (`--backend ref|avx2|avx512 --mode 128|256|512`): `tools/timecop/run_sign_smoke.sh`
- Constant-time cmov suppression: `tools/timecop/shuttle-ct.supp`
- Harness: `tools/dudect/timecop_smoke.c` (backend-agnostic; built against each backend's signing sources).

From the SHUTTLE root, `sh tools/timecop/run_all.sh` is the only command needed; it builds the patched Valgrind on first use and runs every probe. `TIMECOP_BACKENDS`, `TIMECOP_MODES`, and `TIMECOP_TRY_AVX512=1` tune what runs. Opt-in via `TIMECOP=1 ./run_tests.sh`.

## Methodology

The patched Memcheck (KyberSlash `valgrind-varlat` patch) adds `--variable-latency-errors=yes`, which flags when a secret/undefined operand reaches a **variable-latency instruction** (`div`/`mod`/...). This is the exact class KyberSlash exploits.

- **Scope of the smoke = the variable-latency property only.** The gate is: *fail iff Memcheck reports `Variable-latency instruction operand ...`.* Secret-dependent **control flow** is owned separately by the object-code CT scan (`tools/ct_scan.py` rejects `idiv|div|sdiv|udiv|vdiv` and secret gather/scatter in disassembly) and by `dudect-sign`.
- **`--undef-value-errors` must stay at its default (`yes`).** The variable-latency check is implemented on top of Memcheck's undefined-value instrumentation; `--undef-value-errors=no` silently disables it (tests nothing).
- **`Memcheck:Cond` is suppressed; `Memcheck:Value` is not.** SHUTTLE's sampler / ApproxExp / ApproxLog / IRS / packing use branchless constant-time idioms (the 96-bit `cdt_scan96` full-table scan, masked selects, `ct_lt_u32` / `ct_sel_i32` / `ct_is_zero_u32`) that GCC `-O2`/`-O3` lowers to `cmov`/`setcc`. Memcheck reports each secret-conditioned `cmov` as `Conditional jump or move depends on uninitialised value` — a known false positive on constant-time code — so those are suppressed via `shuttle-ct.supp`. The patched Valgrind records a variable-latency hit as an `Err_Value` (same kind as a plain use-of-uninitialised value), so a blanket `Memcheck:Value` suppression would blind the gate; we keep `Memcheck:Value` un-suppressed and let the message-text grep isolate the variable-latency signal. Each smoke run saves the **unsuppressed** stderr as an audit artifact.

## Smoke Taint Scope (SHUTTLE re-key)

`tools/dudect/timecop_smoke.c` marks the signing-secret region of the packed `sk` and the fresh signing randomness as undefined/secret; sampler and IRS inputs are tainted by propagation through the normal signing data flow. SHUTTLE's `skEncode` layout is

```
seedA(SEEDBYTES) + EM*POLYPK_PACKEDBYTES (b body)
  + masterSeed(CHALLENGESEEDBYTES) + tr(CHALLENGESEEDBYTES)
  + ELL*(secret s) + EM*(secret e')
```

so the harness marks `masterSeed` and the trailing `ELL*(secret s) + EM*(secret e')` region UNDEFINED, and keeps `seedA`, the `b` body, and `tr` DEFINED. Under NGCC_MODE the signing randomness is drawn from `drng_algorithm` (SM3 Hash-DRBG); the harness taints the `rnd` bytes that feed `rhoprime` and the IRS `0x09 || seed_y` stream (Overview 4.2, K2). Public `seedA`, cached `bn`, message, verifier inputs, and the produced signature (after signing) are kept defined/declassified.

## avx512 = SKIP (exit 3), not FAIL

Valgrind 3.23's AVX-512 support is incomplete; it cannot decode the EVEX/VBMI2 instructions the backend emits (and `-march=native` even auto-vectorizes scalar Keccak/SM3 to AVX-512). The native binary runs fine, but under Valgrind the guest aborts with `SIGILL — Illegal opcode` before signing starts. This is a Valgrind limitation, not a code finding, and is not fixed by a newer Valgrind (AVX512IFMA/VBMI2 remain unsupported) — so `run_sign_smoke.sh --backend avx512` reports **SKIP** (exit 3). The avx512 no-variable-latency-division property is enforced statically: `tools/ct_scan.py` disassembles the avx512 objects across `-O3`/`-Os` and `gcc`/`clang`, and `make -C avx512 dudect-sign` (when wired) measures the avx512 signing path statistically with real `rdtsc`. SHUTTLE-512's avx512 superblock (`q59393n1024`) makes the SIGILL even more certain.

## Patch Sources

- KyberSlash artifacts page: `https://kyberslash.cr.yp.to/papers.html`
- Preferred patch: `valgrind-varlat-patch-20240808.txt` (sha256 `b4fe3b9f1badcbec078e2dda08718ccc1571d3cf4f79c4b80ef6adb157d3999e`; native target Valgrind 3.23.0, sha256 `c5c34a3380457b9b75606df890102e7df2c702b9420c2ebef9540f8b5d56264d`)
- Fallback patch: `valgrind-varlat-patch-20250805.txt`

The patch is cut from the Valgrind git tree, so its first hunk edits `.gitignore`, which the released source tarballs do not ship. `build_valgrind_varlat.sh` filters out hunks that modify files absent from the extracted tree (keeping new-file creations), applies the remainder with fuzz, and verifies success by capability (`--variable-latency-errors` present) rather than by `patch`'s raw exit code.

## Current Decision

TIMECOP is an **optional, integrated** audit workflow:

- It is **not** part of the default hard gate (`./run_tests.sh`, `./security_audit_matrix.sh`), because it requires a one-time local patched-Valgrind build (download + ~3 min compile).
- It is wired in as an **opt-in**: `TIMECOP=1 ./run_tests.sh` appends `tools/timecop/run_all.sh` (ref + avx2 across all modes) after the standard gates.
- The static object-code CT scan remains the mandatory owner of the no-variable-latency-division and no-secret-gather properties for **all** backends; the TIMECOP smoke is the dynamic cross-check for the backends Valgrind can run (`ref`, `avx2`).

## Float-exp variant note (P13-T12)

SHUTTLE is integer-only by default (Overview 4.7). There is **no opt-in float / IFMA exp path**: no `approx_exp_avx2.c` / `approx_exp_avx512ifma.c`, no `-DAVX2_EXP` / `-DAVX512IFMA_EXP`. Therefore `tools/ct_scan.py --include-optional` is a documented no-op (adds no variant), and the `check-kat-avx2exp` / `check-kat-avx512ifma` gates do not exist for SHUTTLE.

## Environment

Recorded on the machine where this harness was authored (2026-06-24).

| Component | Value |
|---|---|
| OS | Ubuntu 22.04 (WSL2) |
| Kernel | Linux 5.15.153.1-microsoft-standard-WSL2 |
| Compilers | gcc 11.4.0, clang 14 |
| ISA | AVX2 + AVX-512 present |

Build dependencies for the patched Valgrind (`apt-get install`): `build-essential`, `automake`, `autoconf`, `libtool`, `libc6-dbg`, `pkg-config`.

## Result Report

(Populated by `tools/timecop/{prove_varlat.sh,run_sign_smoke.sh}` on each run; the patched-Valgrind build is one-time and local. Until built, the variable-latency property is owned by `ct_scan_matrix.sh` + `run_dudect_all.sh`.)
