# SHUTTLE build notes

This file records the build-system conventions and submission mapping for the SHUTTLE reference + AVX2 + AVX-512 implementation. It is the durable home for the build-system decisions; the final packaging step assembles the submission tree.

## Three backends, one source of truth

`ref/` is pure C99 scalar and is the single source of truth (the KAT oracle). `avx2/` and `avx512/` symlink every shared scalar file from `ref/` (`tools/mklinks.sh` materializes the links) and fork only the files that need vectorization. The fork set per backend is: `config.h`, `polyvec.c`, `sampler.c`, `irs.c`, `rounding.c`, `sign.c`, plus the AVX-only XOF/NTT files and the existing SM3 N-way files. `params.h`, `namespace.h`, `api.h` are backend-independent and are relative symlinks (`../ref/<f>`); `config.h` is a real fork carrying the `_ref` / `_avx2` / `_avx512` namespace tail so all 3 sets times 3 backends (9 libraries) co-link without symbol clash.

## MODE axes

`SHUTTLE_MODE` (128/256/512) is `-D`-injected per built binary (`out/<test>_<mode>`); there is no runtime mode dispatch. Each binary is one parameter set. Orthogonally, the XOF MODE is a second build-time axis: `make MODE=SHA3 ...` injects `-DSHA3_MODE` (SHAKE backend); the default is the NGCC SM3 Hash-DRBG backend (no extra flag). The two axes multiply into six KAT sets; each must be byte-exact across ref/avx2/avx512.

## NGCC submission folder mapping

- `ref/` maps to `Implementations/Reference_Implementation/SHUTTLE-<set>/`
- `avx2/` maps to `Implementations/Optimized_Implementation/SHUTTLE-<set>/`
- `avx512/` maps to `Implementations/Additional_Implementation/SHUTTLE-<set>/`

The AVX-512 build steps outside the literal NGCC perf-flag string (it adds `-mavx512*`), so it ships as "Additional". The NativeUbuntu bench box has no AVX-512, so AVX-512 performance is local-only. One `Test_Vectors/KAT_SIG_SHUTTLE-<set>.txt` per set.

## Exact compile flags per backend

- Reference (`ref/`): `gcc -std=c99 -Wpedantic -Wall -Wextra -O2` — pure ISO C99, no asm/intrinsics/extensions, MUST be warning-clean.
- Optimized (`avx2/`): `gcc -O3 -march=x86-64 -mavx2 -mtune=native -flto -fomit-frame-pointer -std=c99 -Wpedantic -Wall -Wextra -mbmi2 -mpopcnt -mno-avx512f`. The `-mno-avx512f` forces a pure-AVX2 path (it transitively disables every AVX-512 subset while keeping AVX2/BMI2/POPCNT) and is enforced by `make check-no-avx512`.
- Additional (`avx512/`): the Optimized flags minus `-mno-avx512f`, plus `-mavx512f -mavx512bw -mavx512dq -mavx512vl -mavx512vbmi2` (VBMI2 is required by the existing SM3 N-way code).

## The -Wpedantic-vs-intrinsics decision

The mandated perf flags carry `-std=c99 -Wpedantic`, which warns on `__m256i`, inline asm, and `__attribute__`. Decision: tolerate the pedantic warnings on the avx2/avx512 perf builds (they are unavoidable for SIMD and the NGCC flag set is fixed), but keep the reference build strictly warning-clean — that is the hard "no warnings" expectation, enforced by `make check-warn` in `ref/` (adds `-Werror`). Where a single intrinsic TU emits a noisy warning, localize with `#pragma GCC diagnostic push/ignored "-Wpedantic"` around the intrinsic include only. Do NOT add `-Wno-pedantic` globally. If the NGCC submission rules require a clean perf build too, the fallback is per-TU pragma silencing; flag this at packaging time.

## Reproducible constants

Every runtime constant lives in an `@@AUTOGEN:<name>@@ BEGIN/END` region inside a hand-written header, emitted by a `tools/gen_*.py` generator with an audit log in `tools/log/`. `make tables` regenerates all regions; `make check-consts` re-runs the generators and `cmp`s against the committed bytes (drift gate). `gen_params.py` (derived NTT/mask/size constants) and `gen_bounds.py` (integer squared-norm thresholds) both write `params.h` and `tools/log/{params,bounds}_derivation.txt`. No magic numbers.

## Open placeholders (pinned once the corresponding component lands)

- `CRYPTO_SECRETKEYBYTES` is a skEncode-layout estimate (confirmed against the KAT once keygen/packing land).
- `CRYPTO_BYTES` is a generous placeholder over the spec-table sig-size target (the rANS encoder pins `RANS_RESERVED_BYTES`).
- The KAT hashes in the avx2/avx512 Makefiles are `TBD` until the KAT recording step fills them in.
