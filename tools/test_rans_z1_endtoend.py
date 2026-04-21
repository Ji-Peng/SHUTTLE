#!/usr/bin/env python3
"""End-to-end validation of the Phase 6c z[1..L] HighBits rANS tooling.

Workflow mirrors ``test_rans_endtoend.py`` but for the z1 high distribution
produced by ``sample_z1.py``:

  1. Load a z1 histogram JSON.
  2. Build the quantized frequency table.
  3. Draw fresh z1 rounds (independent of training data).
  4. Encode via ``rans.py``; decode; verify round-trip.
  5. Report the per-round rANS size and the resulting signature-size
     contribution of z[1..L] (= L * n * log2(alpha_h) bits of low-part
     + rANS bits of high-part).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from rans import RansTable, encode, decode                          # noqa: E402
from sample_hints import MODE_PARAMS                                # noqa: E402
from sample_z1 import sample_z1_one_round                           # noqa: E402
from gen_rans_tables import load_histogram, symmetrize, quantize    # noqa: E402

import numpy as np


def build_table(histogram_path: str, prob_bits: int = 12):
    mode, kind, raw, _ = load_histogram(histogram_path)
    if kind != "z1_high":
        raise ValueError(f"expected kind=z1_high, got {kind!r} in {histogram_path}")
    sym = symmetrize(raw)
    syms, freqs = quantize(sym, prob_bits)
    return mode, RansTable(syms=syms, freqs=freqs, prob_bits=prob_bits), syms, freqs


def endtoend(mode: int, histogram_path: str,
             fresh_trials: int = 20, prob_bits: int = 12):
    print(f"\n=== SHUTTLE-{mode} z[1..L] rANS end-to-end ===")
    ret_mode, table, syms, freqs = build_table(histogram_path, prob_bits)
    assert ret_mode == mode, f"histogram mode mismatch {ret_mode} vs {mode}"
    print(f"Table: {len(syms)} symbols in [{syms[0]}, {syms[-1]}], prob_total = 2^{prob_bits}")

    params = MODE_PARAMS[mode]
    n_per_round = params['L'] * params['n']
    low_bits = int(math.log2(params['alpha_h']))
    assert (1 << low_bits) == params['alpha_h'], "alpha_h must be power of 2"

    rng = np.random.default_rng(0xBADC0FFE ^ mode)
    msg = []
    for _ in range(fresh_trials):
        z_hi = sample_z1_one_round(rng, params)
        msg.extend(int(v) for v in z_hi.ravel())

    unseen = set(msg) - set(syms)
    if unseen:
        print(f"  WARN: {len(unseen)} unseen symbols in fresh sample: {sorted(unseen)[:20]}")
        print(f"        These are tail events; training window was too small.")
        return False

    enc_bytes = encode(msg, table)
    recovered = decode(enc_bytes, table, len(msg))
    if recovered != msg:
        print("  FAIL: round-trip mismatch")
        for i, (a, b) in enumerate(zip(msg, recovered)):
            if a != b:
                print(f"    first diff at {i}: msg={a} dec={b}")
                break
        return False

    n = len(msg)
    bits_per_sym = 8 * len(enc_bytes) / n
    prob_total = 1 << prob_bits
    fit_entropy = -sum(f / prob_total * math.log2(f / prob_total) for f in freqs if f > 0)
    # z1 size contribution per signing round:
    size_rans_high = len(enc_bytes) / fresh_trials          # rANS-compressed high part
    size_low       = n_per_round * low_bits / 8              # uniform low part
    size_z1_total  = size_rans_high + size_low
    # Uncompressed baseline (current v2 path): L*n*14 bits.
    size_baseline  = n_per_round * 14 / 8

    print(f"  round-trip: PASS  ({n} symbols, {len(enc_bytes)} bytes)")
    print(f"  actual rate:         {bits_per_sym:6.3f} bits/symbol")
    print(f"  fitted rate:         {fit_entropy:6.3f} bits/symbol (rANS floor)")
    print(f"  per-round z[1..L] size (L={params['L']}, n={params['n']}, alpha_h={params['alpha_h']}):")
    print(f"    v2 baseline (14 bit/coef): {size_baseline:7.1f} B")
    print(f"    v3 low part ({low_bits} bit/coef uniform): {size_low:7.1f} B")
    print(f"    v3 high part (rANS):       {size_rans_high:7.1f} B")
    print(f"    v3 total z[1..L]:          {size_z1_total:7.1f} B  "
          f"(-{size_baseline - size_z1_total:.1f} B vs v2)")
    return True


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
    print(f"\n=== z1 end-to-end {'PASSED' if ok else 'FAILED'} ===")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
