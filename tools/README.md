# SHUTTLE Python Tooling (Phase 6b-4)

Offline helpers for Phase 6b-4 of the SHUTTLE implementation: empirical
hint sampling, rANS frequency-table generation, and a reference rANS
encoder/decoder used for cross-validation of the C-side implementation in
Phase 6b-5.

## Files

| file | role |
|---|---|
| `rans.py` | Pure-Python rANS encoder/decoder (reference; for cross-check with C). |
| `sample_hints.py` | Simulate SHUTTLE signing to collect a hint histogram per mode. |
| `gen_rans_tables.py` | Turn a histogram JSON into a quantized frequency table and C header. |
| `test_rans_endtoend.py` | Validate the full round-trip: sample -> table -> encode -> decode. |

## Dependencies

- Python 3.8+
- `numpy` (for polynomial arithmetic in the simulator)

## Quick start

```bash
cd SHUTTLE/tools

# 1. Sample hint distributions (250-500 trials each; takes a few minutes).
python3 sample_hints.py --mode 128 --trials 500 --out hints_mode128.json
python3 sample_hints.py --mode 256 --trials 250 --out hints_mode256.json

# 2. Generate rANS tables. prob-bits=12 means frequencies sum to 4096.
python3 gen_rans_tables.py hints_mode128.json hints_mode256.json \
    --prob-bits 12 --format c --out ../ref/rans_tables.h

# 3. End-to-end validation (Python only).
python3 test_rans_endtoend.py \
    --hist128 hints_mode128.json --hist256 hints_mode256.json \
    --fresh-trials 50
```

Expected end-to-end result:

```
=== SHUTTLE-128 end-to-end ===
  round-trip: PASS  (25600 symbols, 8812 bytes)
  actual rate:    2.754 bits/symbol
  sizes per signing round (512 coef = M*N):
    naive 1-byte:   512.0 B
    rANS:           176.2 B  (saves 336 B)

=== SHUTTLE-256 end-to-end ===
  round-trip: PASS  (51200 symbols, 14983 bytes)
  actual rate:    2.341 bits/symbol
  sizes per signing round (1024 coef = M*N):
    naive 1-byte:   1024.0 B
    rANS:           299.7 B  (saves 724 B)
```

## How the simulator works

`sample_hints.py` replicates the SHUTTLE signing pipeline in Python, using
schoolbook polynomial arithmetic (no NTT) for clarity. Per trial:

1. Sample `a_gen` and `A_gen` uniformly mod q.
2. Sample `s, e` from CBD(η).
3. Compute `b = a_gen + A_gen · s + e (mod q)`.
4. Sample `y` with discrete Gaussian (σ, truncated at 11σ).
5. CompressY: `Y_0 = round(y[0] / α_1)`.
6. Compute `comY = (a_gen − b)·Y_0 + A_gen·y[1..L] + y[L+1..L+M]` and lift
   to mod 2q via `2·U + q·(Y_0 & 1)`.
7. Sample challenge `c` (τ ones, rest zero).
8. Simulate IRS output `z = y + c·sk_stretched` (idealized; no rejection).
9. Split `z_2 = z[L+1..L+M]`.
10. Compute `h = MakeHint(z_2, comY)` per spec Alg 9, with the raw
    difference centered to `(-HINT_MAX/2, HINT_MAX/2]`.

The simulator does NOT do CompressY on z[0] here because the hint's
distribution only depends on z_2. Realistic IRS rejection is also skipped
(it only affects the overall signature-rate probability, not the hint's
conditional distribution given successful signing).

**Caveat**: Because IRS rejection is idealized away, the distribution seen
here is an approximation of the real one. For production tables, a
follow-up pass should sample from the actual C signing code after Phase
6b-5 is complete.

## Frequency table generation

`gen_rans_tables.py` enforces symmetry (`p(v) ← (p(v) + p(−v))/2`),
quantizes real frequencies to integers summing to `2^prob_bits`, and
outputs either a text summary or a C header with:

```c
#define SHUTTLE128_RANS_PROB_BITS 12
#define SHUTTLE128_RANS_NUM_SYMS  15
static const int16_t  shuttle128_rans_syms[15];   /* ascending symbol values */
static const uint16_t shuttle128_rans_freqs[15];  /* integer freqs, sum = 2^12 */
```

The C-side rANS engine (Phase 6b-5) will build the `cdf[]` and `sym_lookup[]`
decode tables at load time from these constants.

## rANS reference implementation

`rans.py` follows Fabian Giesen's `rANS_byte` pattern (which HAETAE also
uses): 32-bit state `x`, byte-wise renormalization, encode pushes symbols
in reverse order and the finalize step flushes `x` as 4 bytes. Decode
mirrors: pulls the 4 state bytes first, then pops symbols left-to-right.

```python
from rans import RansTable, encode, decode
table = RansTable(syms=[...], freqs=[...], prob_bits=12)
stream = encode([0, -1, 2, 0, 1, ...], table)
recovered = decode(stream, table, len(msg))
```

Use `python3 rans.py` to run the built-in self-test.

## Known limitations

- **Idealized IRS**: rejection sampling is skipped in the simulator; real
  hint distribution conditional on successful signing may differ slightly.
- **Fixed-range symbol set**: the quantized table includes only symbols
  observed during training. Production code must handle tail events (a
  rare symbol not in the table) either by enlarging the table with
  symmetric padding or by a fallback path (rejection + resample).
- **No IND-CCA2 considerations**: the table is purely a statistical
  compressor; constant-time decoding on the C side is a Phase 6b-5
  concern.

## Known sizes (Phase 6b-4 sampling)

After running `sample_hints.py` with sufficient trials, the hint
distribution should converge to:

| Mode | σ | 2σ/α_h | empirical std | entropy | rANS rate |
|---:|---:|---:|---:|---:|---:|
| 128 | 101 | 1.58 | 1.63 | 2.75 bits | 2.76 bits |
| 256 | 149 | 1.16 | 1.23 | 2.35 bits | 2.35 bits |

These match the theoretical Gaussian approximation in the design memo
(`docs/NGCC_Sign/SHUTTLE_draft.md` §14) to within 5%.
