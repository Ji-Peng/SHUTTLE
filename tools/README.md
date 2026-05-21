# SHUTTLE Python Tooling

Offline helpers for the SHUTTLE rANS layer: theoretical frequency-table
generation, a pure-Python reference encoder/decoder used to cross-check the
C implementation, and a KAT comparator for byte-identical regression.

The design rationale and all numeric targets live in
`agent/rANS/SHUTTLE_rANS.tex`.

## Files

| file | role |
|---|---|
| `rans.py` | Pure-Python rANS encoder/decoder (reference; cross-checks the C side). |
| `gen_rans_tables.py` | Emit `ref/rans_tables.h` from the theoretical discrete Gaussian PMF. |
| `kat_compare/` | KAT regression helpers. |

## Dependencies

Python 3.8+ (no external libraries).

## Quick start

```bash
cd SHUTTLE/tools

# Regenerate the C frequency tables (one-shot, deterministic).
python3 gen_rans_tables.py --out ../ref/rans_tables.h

# Inspect the generated tables in a human-readable format.
python3 gen_rans_tables.py --format text

# Run the Python reference self-test.
python3 rans.py
```

## Where the frequency tables come from

We do *not* sample any empirical histogram. The discrete Gaussian PMF
is evaluated analytically:

  z-hi context:  sigma = r / alpha_r,  M_voc = ceil((11*r + tau*eta)/alpha_r)
  hint context:  sigma = 2r / alpha_h, M_voc = floor(2*(11*r + tau*eta)/alpha_h) + 1

SampleY hard-truncates `|y_k| <= 11*sigma`, so the M_voc bound holds with
probability 1 — the vocabulary covers every possible coefficient, and
out-of-vocabulary failure is mathematically impossible (`p_OOV ≡ 0`).
mode-128's `alpha_h = 2*alpha_r` makes `sigma_zhi == sigma_hint`, so a
single shared table doubles as both contexts; mode-256/512 keep two
separate tables.

Quantization is the standard "round each `p(s) * 2^t` to an integer,
nudge the largest bucket to absorb the slack" recipe, with `t = 10`.

## Reservation budgets

The per-stream byte reservations in `ref/params.h` come from
`tools/SigSize.py::compute_two_stream_reservation` (Gauss CLT for z-hi,
Poisson model for narrow-sigma hint streams in mode-256/512). Re-run it
when changing any spec parameter — the numbers feed straight into
`SHUTTLE_ZHI_RESERVED_BYTES` / `SHUTTLE_HINT_RESERVED_BYTES`.

## rANS reference

`rans.py` is a literal transliteration of the C engine. State is 32-bit,
renormalization is byte-wise, the final-state flush is little-endian
(matching ryg_rans / HAETAE), and the decoder verifies `x == RANS_L`
after the last symbol — any byte-level corruption that survives the
renorm loop will overwhelmingly land here.

```python
from rans import RansTable, encode, decode
table = RansTable(syms=[...], freqs=[...], prob_bits=10)
stream = encode([0, -1, 2, 0, 1, ...], table)
recovered = decode(stream, table, len(msg))   # raises on final-state mismatch
```
