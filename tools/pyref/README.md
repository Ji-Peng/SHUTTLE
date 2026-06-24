# SHUTTLE Python reference (P12)

Independent, pure-Python (stdlib only) mirror of the deterministic,
KAT-load-bearing layers of the SHUTTLE signature scheme, plus the byte-exact
cross-checks against the C reference (`SHUTTLE/ref/`). The C ref is the
oracle; the Python ref leads it for correctness review and is asserted
byte-for-byte equal. A divergence is a bug.

## Modules (the mirrors)

| file | mirrors | what it proves |
|---|---|---|
| `params.py` | `ref/params.h`, `ref/reduce.h` | per-set constant table; `pk/sk` sizes, `ZETA^n=-1`, `n*NINV=1`, `H_h` integer |
| `reduce_ref.py` | `ref/reduce.{h,c}` | `montgomery_reduce16`/`fqmul16`/`addm16`/`subm16` + `reduce32`/`freeze`/`reduce_mod_2q` |
| `ntt_ref.py` | `ref/ntt/<qset>/ntt_ref.c` | complete CT/GS NTT; `ntt->pointwise->invntt_tomont == schoolbook` |
| `packing_ref.py` | `ref/packing.c` | `integer/poly_to_bytes` (LSB-first), `bytes_to_bits` (MSB-first, K3), `pack_pk/sk/com`, `encode_com_vec` (MS-A7) |
| `drng_ref.py` | `ref/drng.c` | SM3 Hash-DRBG (`init/get_random_number`); SM3("abc") KAT |
| `rans_ref.py` | `ref/rans.c` (via `tools/rans.py`) | SHUTTLE block-plan rANS encode/decode + re-encode equality (K15) |
| `sampler_ref.py` | `ref/sampler.c`, `ref/sampler_u.c` | `cdt_scan96` (96-bit 3x32 limb, borrow-fold), `sampler_u` (MSB-first 80-bit CLZ + 57-bit mantissa + Q62 ApproxLog) |

`sampler_ref.py` parses the committed `tools/rcdt_tables.h` (RCDT tables) and
`tools/approx_log_poly.h` (Q62 log2 coeffs) -- no re-derivation, so it reads
the SAME integers the C does.

## Cross-checks (Python == C, byte-exact)

| file | C oracle | coverage |
|---|---|---|
| `xcheck.py` + `xcheck_dump.c` | reduce/poly/poly_ntt/packing/drng | DRNG draws, `reduce32/freeze/reduce_mod_2q`, NTT fwd+inv, `pk/sk/com` |
| `xcheck_rans.py` + `rans_oracle.c` | `shuttle_rans_encode` | `rans_ref == committed golden == C` on real response vectors + round-trip + re-encode |
| `xcheck_sampler.py` + `sampler_oracle.c` | `sampler_sigma2`/`sampler_u` | `cdt_scan96` byte-exact; `sampler_u` (a,frac) byte-exact from raw 18-byte stream |

The C dumpers seed a deterministic xorshift64 mirrored exactly in Python, so
both sides see identical inputs; the Python side recomputes every dumped value
and asserts equality.

## M4 empirical rANS validation

`m4_validate.py` (+ `m4_dump.c`) signs N real messages with the C rANS path,
parses the realized rANS com length and the Q0/Qs/h symbols, and checks the
`RANS_RESERVED_BYTES` reserve holds (no overflow) + the empirical laws match
the model PMFs (`tools/rans_model.py`).

## Run

    python3 run_pyref.py --selftest    # every module self-test
    python3 run_pyref.py --xcheck      # build+run the 3 C byte-exact xchecks
    python3 run_pyref.py --m4 10000    # M4 empirical reserve check (N sigs)
    python3 run_pyref.py --all         # selftest + xcheck

C dumpers (`*_dump_*`, `*_oracle_*`) are build artifacts written next to the
sources; they are gitignored.
