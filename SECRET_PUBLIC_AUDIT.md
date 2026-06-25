# Secret/Public Audit Baseline (SHUTTLE-NGCC)

This note is the secret/public classification baseline for SHUTTLE. Any future division, modulus, lookup, branch, gather allowlist, rejection-timing decision, or declassification in secret-handling code must cite one of the entries below, or extend this note in the same style, before the code change lands. It is the human-auditable source of truth that every `tools/ct_scan.py` allowlisted forbidden-instruction exemption points at.

SHUTTLE's parameter sets are `128/256/512` (module dims `1+ELL+EM = KVEC = 7/6/6`), the default XOF is the NGCC SM3 Hash-DRBG (a `XOF / DRNG` section, replacing a SHAKE-only path), and the core rejection sampler is the rejection-FREE Iterated Rejection Sampler (IRS), so the per-message iteration count is NOT a secret-dependent loop count (see the IRS section).

## Classification Labels

- `secret`: Derived from signing secret material, signing randomness, Gaussian/noise samples, IRS randomness, or unreleased intermediate values whose distribution or value must not control variable-latency operations.
- `public`: Known to the verifier, attacker-controlled, deterministic from public inputs, or already serialized in public artifacts.
- `public length`: A byte count, mode constant, loop bound, fixed retry budget, or one-shot buffer size that may drive branches or divisions because it is independent of secret values.
- `declassified`: A value derived from secret data but intentionally released by the scheme through the public key, signature, or accept/reject output after the documented check.

## Keygen

- The keygen seed `xi`, expanded sub-seeds, the master/key seed, the secret factors `f0`, `f`, `s0`, `e`, the assembled secret vector `s` (the `ELL` secret polys and `EM` secret error polys), the keygen 2-norm `nrm`, and the norm-window / invertibility retry state are `secret` until `pk` and `sk` are packed.
- `seedA`, the expanded public matrix `A` (NTT domain), loop bounds, module dims (`ELL`, `EM`, `KVEC`), polynomial dimension `N`, and pack sizes are `public` or `public length`.
- `b` and the canonical cached `bn` are `declassified` when packed into `pk` (`seedA || EM * POLYPK_PACKEDBYTES`) and copied into `sk`; after that point `bn` may be treated as public cached key material.
- Keygen rejection timing (the `B_k' <= ||.||_2 <= B_k` norm-window gate, ~34/47/38 % accept per set) is deployment-classified as a **provisioning-time** retry, not a per-message side channel, but the arithmetic inside each attempt still follows `secret` rules (constant-time per attempt, accepted or not).

## Sign

- The unpacked secret key, the secret vector `s`, signing randomness `rnd`, `rhoprime`, `seed_y`, the Gaussian `y`, the pre-release response `z` / `z1` / `z2'`, the IRS RNG state, the IRS state `V` / `t` / `ell` / `flag`, the 2-norm value, all rejection decisions, and the rANS overflow-retry state are `secret`.
- Message bytes and message length are attacker-controlled `public`; `tr`, `mu`, and the challenge-hash absorption length are public once bound to the public key and message transcript.
- `seedA`, cached public `bn`, expanded `A`, mode constants, `CRYPTO_BYTES`, and rANS reserve lengths are `public` or `public length`.
- The challenge `c`, the hint `h`, `z1`, the rANS `z`-head/low symbols, the rANS stream length, and the zero padding are `declassified` ONLY through `pack_sig`; before successful packing they remain secret-derived and must not drive variable-latency operations except through documented fixed-format writers.

## Verify

- Public key, message, signature bytes, signature length, the unpacked `c` / `h` / `z1`, the reconstructed `w0`-bits / `z2'` / 2-norm, and the verifier-side narrower `A1`-hat (the `1+ELL`-column matrix WITHOUT the `2*I_m` block) are `public`, because verification operates only on adversary-provided or public-key-derived data.
- Verify accept/reject codes and early exits (length mismatch `-2`, invalid-sig `-1`, accept `0`) are public parser behavior; they must be deterministic and covered by the negative-parser corpus, but they are NOT secret-dependent branches.
- The reconstructed `z2'` and the verify-side 2-norm are `public` because they are functions of public signature and public-key data.

## Sampler

- Gaussian and noise sample values, sign bits, CDT decisions, ApproxExp acceptance values, refill state that depends on accepted samples, and the generated secret-key / signing vectors are `secret`.
- XOF output block counts derived from fixed sampler batch sizes, initial refill sizes, lane counts, and mode constants are `public length`.
- The 96-bit RCDT tables, the ApproxExp / ApproxLog polynomial tables, and the fixed distribution tables are `public`; indexes into them are allowed only through full scans (the branchless 96-bit `cdt_scan96` full-table scan + borrow-fold equality mask), fixed-lane SIMD selection, or documented constant-time table selection. SamplerU's `__builtin_clzll` (MSB-first bit extraction) lowers to `lzcnt`/`bsr` (NOT in the forbidden set) and consumes no float conversion (SHUTTLE is integer-only).

## IRS (see the "ISOCHRONY / LEAKAGE" notes in irs.h)

- `seed_y`, the stretched secret `sk_tilde`, the IRS accumulator `t`, the state `V`, `ell`, and `flag`, the random bits from the IRS RNG, and the adjusted `z` values are `secret` until signature packing declassifies the accepted output.
- The challenge weight `tau` (the support size), the polynomial dimension `N`, and the IRS loop trip count are `public length`.
- **The IRS is rejection-FREE: RejectSample applies EXACTLY `tau` R-transitions per call and never aborts, so there is NO secret-dependent loop count and the SamplerU byte consumption is a fixed `18*tau` bytes (deterministic).** The single documented public-branch exception is the ascending-`j` `if (c[j])` test inside RejectSample: `c` is recomputed by the verifier from `seedC`, so the challenge support pattern is `public`; gating on `c[j]` reveals only the public challenge support, not secret data. The only timing observable in the whole signing path is the **outer** Sign-loop iteration count driven by the PUBLIC norm test `||(z1,z2')||_2 <= B_v` — that count is a function of public/declassified quantities. Any future shortcut in IRS rejection or coefficient adjustment must cite this section and prove it does not expose secret-dependent variable latency.

## Challenge

- Sign-side `buckets`, `w0`-bits, `mu`, `seedC`, and `c` are secret-derived until the final signature declassifies `c`; hash absorption lengths and mode constants are `public length`.
- Verify-side challenge recomputation inputs are `public`; a challenge mismatch is public reject behavior.
- Challenge packing is fixed-size and public after serialization; malformed challenge bytes are parser-negative corpus inputs.

## rANS

- Sign-side `z`-head and hint symbols, the rANS encoder states, the stream bytes before publication, the rANS used length, and the overflow-retry decisions are `secret` until the signature is successfully emitted. The rANS renormalization is a division acting DIRECTLY on the secret symbol stream — it MUST be implemented as precomputed-reciprocal multiplication (the single most dangerous division site, KyberSlash class); `tools/ct_scan.py` enforces that `rans.c` contains ZERO division mnemonics (any div here is a real bug, NOT a public-length exception — `rans.c` carries no allowlist entry).
- rANS frequency / CDF / reciprocal tables, the interleave count `RANS_N`, the reserve size, and the final fixed signature size are `public` or `public length`.
- Verify-side rANS bytes, decoded symbols, the non-canonical-state checks (initial state out of range, incomplete consumption, terminal state != L, nonzero padding, CDF-hole), the reserve-overflow check, and the padding checks are `public` parser behavior and must stay in the parser-negative corpus.

## Packing

- Key packing treats the secret-key seed and the secret polynomials (`s`, `e'`) as `secret`; public-key bytes and the cached `bn` are `declassified` public-key material. The SHUTTLE `skEncode` layout is `seedA (SEEDBYTES) + EM*POLYPK_PACKEDBYTES (b body) + masterSeed (CHALLENGESEEDBYTES) + tr (CHALLENGESEEDBYTES) + ELL*(secret s) + EM*(secret e')`, summing to `CRYPTO_SECRETKEYBYTES = 2288 / 3680 / 7104`.
- Sign-side signature packing handles secret-derived `c`, `h`, and `z1` until the complete signature is emitted; zero padding is `declassified` and authenticated by verify.
- Verify-side unpacking treats all signature bytes as attacker-controlled `public`; branches on malformed lengths, padding, rANS decode failure, the non-power-of-2 hint range (a decoded hint `>= H_h = 30/120/58` must REJECT, never `mod H_h`-wrap), and signature size are public parser behavior.

## NTT Boundary

- Secret-key and signing vectors entering NTT or inverse-NTT routines are `secret`; the public matrix `A`, the public-key `bn`, and verifier-side imported NTT data are `public`.
- Backend-native NTT slot order is an INTERNAL representation detail, NOT a declassification event; only `pack_pk` / `pack_sk` / `pack_sig` define serialized declassification boundaries.
- SIMD gathers or shuffles may be allowlisted only when their addresses depend on fixed permutations, public lane layout, public SHAKE lane pointers, or mode constants; secret-indexed gather remains forbidden. SHUTTLE-512's avx512 superblock (`q59393n1024/ntt_avx512.S`) uses only fixed permutes; the avx2 8-block NTT must contain no `zmm` (enforced by `make -C avx2 check-no-avx512`).

## XOF / DRNG

- SM3 Hash-DRBG state and SHAKE state are `secret` when seeded from secret material: the signing randomness `rnd` is drawn from `drng_algorithm` (SM3 Hash-DRBG) under NGCC_MODE and feeds `rhoprime` and the IRS `0x09 || seed_y` stream; that DRBG state is `secret`.
- Block-count divisions on PUBLIC output byte counts (how many SM3 / SHAKE blocks to squeeze for a fixed-length request) are `public length` and MAY drive divisions — but they are allowlisted in `tools/ct_scan.py` only inside the named DRBG/XOF block-count functions (`drng_*`, `get_random_number`, `*xof*`, `sm3_*`, `*_shake128/256`), never blanket. Preferred outcome: division-free (reciprocal-multiply); any live allowlisted division here cites this section.
- Public `seedA`, cached `bn`, message, and verifier inputs stay `public`; the produced signature is `declassified`.

## Current Allowlist Citations

- `tools/ct_scan.py` `ALLOW_GATHER_SOURCES` (avx2/`fips202x4.c`, avx512/`fips202x8.c`) cite the **NTT Boundary** + **Sampler** entries: SHAKE lane pointer/offset selection is public and fixed by the backend schedule.
- `tools/ct_scan.py` `ALLOW_PUBLIC_DIVS` for `fips202.c` `*_shake128/256`, `fips202x4.c` `*x4`, `fips202x8.c` `*x8` cite the **XOF / DRNG** public-length entry (one-shot / N-way SHAKE output-length block counts).
- `tools/ct_scan.py` `ALLOW_PUBLIC_DIVS` for `drng.c` / `auxfunc.c` (`drng_*`, `get_random_number`, `sm3_*`, `*xof*`) cite the **XOF / DRNG** public-length entry (SM3 Hash-DRBG / pseudoXOF output-length block counts). Default expectation: none emitted.
- `tools/ct_scan.py` `ALLOW_PUBLIC_DIVS` for `sampler.c` (`sample_gauss_*`, `noise_magnitude_batch`, `sampler_sigma2`) cite the **Sampler** public-length entry (fixed Gaussian batch block counts).
- `rans.c` carries NO allowlist: the rANS reciprocal-multiply renorm must be division-free (cite **rANS**).
- Verify-side parser early exits cite the **Verify**, **rANS**, and **Packing** public-parser-behavior entries (negative-corpus covered).

## Implementation Risk -> Covering Gate Traceability

| Risk | Covering gate(s) |
|---|---|
| NTT slot order divergence | trace-KAT (`security_audit_matrix.sh`) + NTT-Boundary classification + `check-no-avx512` |
| rnd/seed derivation order | KAT (deterministic SM3 DRBG) + diff-fuzz |
| SamplerU MSB-first CLZ | ct_scan (no `cvtsi2sd`; `lzcnt`/`bsr` not forbidden) + KAT |
| keygen kappa timing | provisioning-time classification (Keygen) + dudect-sign |
| XOF refill cursor determinism | `check-refill` + KAT |
| sampler bit-consumption | trace-KAT + KAT |
| sign-bit fold | KAT + diff-fuzz |
| SIMD == scalar byte-exact + division-free | KAT cross-backend equality + ct_scan (all 3 backends) + trace-KAT + diff-fuzz |
| 96-bit `cdt_scan96` branchless | ct_scan (no div/branch on secret index) + dudect-components |
| base-2 exact exponent term | check-consts (ApproxLog table) + KAT |
| mod-2q parity | diff-fuzz + trace-KAT |
| non-pow-2 hint range | parser-negative (decoded hint `>= H_h` must REJECT) |
| rANS canonical decode | parser-negative (>= 6 rANS mutations) + KAT |
