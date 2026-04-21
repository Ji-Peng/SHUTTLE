#!/usr/bin/env python3
"""End-to-end validation of Phase 6b-4 Python rANS tooling.

For each mode (128, 256):
  1. Sample hints via sample_hints.py (or reuse a saved JSON).
  2. Generate frequency table via gen_rans_tables.
  3. Draw a fresh sample sequence (independent of training data).
  4. Encode via rans.py; decode; verify round-trip.
  5. Compare compressed size to Shannon entropy lower bound.

Also ensures all observed hint values are representable in the frequency
table (no "out-of-vocabulary" samples); flags if the training distribution
was too narrow.
"""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from rans import RansTable, encode, decode                          # noqa: E402
from sample_hints import MODE_PARAMS, simulate_one_round             # noqa: E402
from gen_rans_tables import load_histogram, symmetrize, quantize     # noqa: E402

import numpy as np


def sample_round_vec(rng: np.random.Generator, mode: int):
    """Return a flat list of hint coefficients from one signing round."""
    h = simulate_one_round(rng, MODE_PARAMS[mode])
    return [int(v) for v in h.ravel()]


def build_table(histogram_path: str, prob_bits: int = 12):
    _, _, raw, _ = load_histogram(histogram_path)
    sym = symmetrize(raw)
    syms, freqs = quantize(sym, prob_bits)
    return RansTable(syms=syms, freqs=freqs, prob_bits=prob_bits), syms, freqs


def endtoend(mode: int, histogram_path: str,
             fresh_trials: int = 20, prob_bits: int = 12):
    print(f"\n=== SHUTTLE-{mode} end-to-end ===")
    table, syms, freqs = build_table(histogram_path, prob_bits)
    print(f"Table: {len(syms)} symbols in [{syms[0]}, {syms[-1]}], prob_total = 2^{prob_bits}")

    params = MODE_PARAMS[mode]
    rng = np.random.default_rng(0xAABBCCDD ^ mode)
    n_coef_per_round = params['M'] * params['n']

    # Fresh sample (not from training data).
    msg = []
    for _ in range(fresh_trials):
        msg.extend(sample_round_vec(rng, mode))

    # Check vocabulary coverage.
    unseen = set(msg) - set(syms)
    if unseen:
        print(f"  WARN: {len(unseen)} unseen symbols in fresh sample: {sorted(unseen)[:20]}")
        print(f"        These are tail events that happened after the training window.")
        print(f"        A production build would either extend the table (symmetric padding)")
        print(f"        or apply a fallback (rejection / plain byte) path.")
        return False

    # Encode / decode.
    enc_bytes = encode(msg, table)
    recovered = decode(enc_bytes, table, len(msg))
    if recovered != msg:
        print("  FAIL: round-trip mismatch")
        # Find first mismatch.
        for i, (a, b) in enumerate(zip(msg, recovered)):
            if a != b:
                print(f"    first diff at {i}: msg={a} dec={b}")
                break
        return False

    # Stats.
    n = len(msg)
    n_bits = 8 * len(enc_bytes)
    bits_per_sym = n_bits / n
    prob_total = 1 << prob_bits
    fit_entropy = -sum(f/prob_total * math.log2(f/prob_total) for f in freqs if f > 0)
    emp_entropy = _empirical_entropy(msg)

    size_naive_byte = n                      # current `polyveck_hint_pack_basic`
    size_naive_4bit = (n * 4 + 7) // 8       # if we used 4 bits/coef
    size_rans       = len(enc_bytes)

    size_per_round_rans  = size_rans / fresh_trials
    size_per_round_naive = size_naive_byte / fresh_trials

    print(f"  round-trip: PASS  ({n} symbols, {len(enc_bytes)} bytes)")
    print(f"  actual rate:         {bits_per_sym:6.3f} bits/symbol")
    print(f"  fitted rate:         {fit_entropy:6.3f} bits/symbol (floor)")
    print(f"  empirical entropy:   {emp_entropy:6.3f} bits/symbol")
    print(f"  sizes per signing round ({n_coef_per_round} coef = M*N):")
    print(f"    naive 1-byte:      {size_per_round_naive:6.1f} B")
    print(f"    naive 4-bit:       {size_naive_4bit/fresh_trials:6.1f} B")
    print(f"    rANS:              {size_per_round_rans:6.1f} B  "
          f"(vs current {size_per_round_naive:.0f} B -> saves {size_per_round_naive - size_per_round_rans:.0f} B)")
    return True


def _empirical_entropy(msg):
    from collections import Counter
    cnt = Counter(msg)
    tot = len(msg)
    return sum(-c/tot * math.log2(c/tot) for c in cnt.values())


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist128", required=True)
    ap.add_argument("--hist256", required=True)
    ap.add_argument("--fresh-trials", type=int, default=20,
                    help="number of fresh signing rounds to encode")
    args = ap.parse_args()

    ok = True
    ok &= endtoend(128, args.hist128, fresh_trials=args.fresh_trials)
    ok &= endtoend(256, args.hist256, fresh_trials=args.fresh_trials)
    print(f"\n=== end-to-end {'PASSED' if ok else 'FAILED'} ===")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
