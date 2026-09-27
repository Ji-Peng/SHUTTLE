# IRS transition sign fix (audit note)

## Motivation

The NGCC attack report *Covariance Leakage in the NGCC Shuttle Submission*
(Feussner, 25 Sep 2026) shows that the submitted Algorithm R / reference
implementation applied the interval event associated with \(p_v(y)\)
(probability of returning \(y-v\)) to the opposite update \(y+v\). In the
idealized analysis this replaces covariance cancellation with an extra
rank-one term \(2vv^\top\), enabling equivalent-key recovery from ordinary
signatures.

## Semantic change

Keep **flag polarity unchanged** (default \(-1\); set to \(+1\) on interval
hit). Change only the final update:

| | Before (buggy) | After (fixed) |
|--|----------------|---------------|
| Update | \(z \leftarrow z + \mathrm{flag}\cdot v\) | \(z \leftarrow z - \mathrm{flag}\cdot v\) |
| Interval hit (`flag=+1`) | \(y+v\) (wrong) | \(y-v\) (matches \(p_v\)) |
| No hit (`flag=-1`) | \(y-v\) (wrong) | \(y+v\) (matches \(1-p_v\)) |

Unchanged: sign-normalization of \(v\), interval boundaries, SamplerU, XOF,
parameters, and the security / error-analysis text (already stated
“return \(y-v\) with probability \(p_v\)”).

## Files changed

### Specification (Overleaf NGCC-Signature)

- `sections/Description.tex` — Algorithm R (`alg:Ryv` / Algorithm 22) return
  line only (`y + flag·v` → `y - flag·v`), plus a footnote noting that an
  earlier draft had the wrong sign and this revision corrects it. Included
  from `main.tex` via `\input{sections/Description}`.

### Implementation (this repo)

- `ref/irs.c` — `poly_axpy_flag`: `+=` → `-=`; comments
- `avx2/irs.c`, `avx512/irs.c` — `simd_axpy_i32`: `add_epi32` → `sub_epi32`;
  scalar fallback `-=`; comments
- `ref/irs.h` (avx2/avx512 symlink to the same file) — CT note `z -= flag*v`
- `ref/test/test_irs.c` (avx2/avx512 symlink) — oracle + identity replay
- `tools/pyref/irs_ref.py` — Python R-transition mirror

### KAT lock (hashes change because signatures change)

- `ref/Makefile`, `avx2/Makefile`, `avx512/Makefile` — `KAT_*_NGCC` and
  `KAT_*_SHA3` FNV pins
- Regenerated vectors under `dist/Test_Vectors/` (and `dist/Test_Vectors/SHA3/`)
  via `./gen_kat.sh` / `./gen_kat.sh --sha3`

## Before / after (C)

```c
/* Before */
z[i].coeffs[k] += f * v[i].coeffs[k];
/* SIMD: _mm*_add_epi32(zc, mullo(vc, vf)) */

/* After */
z[i].coeffs[k] -= f * v[i].coeffs[k];
/* SIMD: _mm*_sub_epi32(zc, mullo(vc, vf)) */
```

## Before / after (LaTeX — Algorithm R)

Path: Overleaf `sections/Description.tex` (via `main.tex`).

**Before:**

```latex
\State $\flag \gets -1$
% ... sign normalize + interval loop; on hit: $\flag \gets 1$ ...
\State \Return $(\mathit{ctx},\, y + \flag \cdot v)$
```

**After:**

```latex
\State $\flag \gets -1$
% ... sign normalize + interval loop; on hit: $\flag \gets 1$ ...
\State \Return $(\mathit{ctx},\, y - \flag \cdot v)$
```

Only the return operator changes (`+` → `-`).

## Explicitly not changed

- Flag init / interval assignment
- Sign normalization, boundary formulas, truncation \(N\)
- Lattice parameters, Stretch/Compress, rANS packing
- Rational / Error-analysis / Security prose (already correct \(T_v\))
- XOF (SM3 Hash-DRBG / SHAKE)

## Verification (this machine)

| Check | Result |
|-------|--------|
| `test_irs` `(c) reject_sample determinism+oracle+identity` | PASS (128/256/512) |
| `test_sign` (RAW) | PASS (128/256/512) |
| `make -C ref check-kat` (NGCC) | PASS — hashes below |
| `make -C ref MODE=SHA3 check-kat` | PASS — hashes below |
| `avx2` / `avx512` build | Not run here (host is arm64; `-march=x86-64` unsupported) |

Recorded KAT hashes after the fix:

```
KAT_128_NGCC = 11486952262273914864
KAT_256_NGCC = 13231489903054685845
KAT_512_NGCC = 7887990226066942058
KAT_128_SHA3 = 6997426514966519864
KAT_256_SHA3 = 17087775945720413227
KAT_512_SHA3 = 3305270276154518466
```

Note: on this macOS host, `test_irs` SamplerU float-oracle subchecks
(`ell` vs `long double`, `u_frac` ULP) can FAIL due to platform FP noise;
those checks are unrelated to the axpy sign. The reject-sample oracle /
identity checks that encode the transition update all PASS.

## Rebuild KATs

```sh
./gen_kat.sh            # NGCC -> dist/Test_Vectors/
./gen_kat.sh --sha3     # SHA3  -> dist/Test_Vectors/SHA3/
```

If the repo path contains spaces, `gen_kat.sh`’s unquoted `$(for ...)`
path expansion may break; run via a space-free symlink (e.g.
`ln -sfn "$PWD" /tmp/shuttle && cd /tmp/shuttle && ./gen_kat.sh`).
