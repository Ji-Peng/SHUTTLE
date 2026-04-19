#!/usr/bin/env python3
"""Generate static rANS frequency tables from empirical hint histograms.

Reads a hint histogram JSON produced by sample_hints.py, quantizes to
integer frequencies summing to 2^prob_bits (typically 2^12), and emits
a C header snippet embedding the table.

Quantization procedure (safe against zero-frequency symbols):

  1. Enforce symmetry: replace p(v) with (p(v) + p(-v)) / 2.
  2. Compute exact real frequencies  f(v) = p(v) * PROB_TOTAL.
  3. Round each f(v) to nearest integer, minimum 1 for observed symbols.
  4. Adjust the single largest bucket's frequency up or down so the sum
     equals PROB_TOTAL exactly.

Output format is controlled by --format:
  --format text : human-readable table + stats
  --format c    : emits C header text on stdout

Usage:
  python3 gen_rans_tables.py hints_mode128.json --prob-bits 12 --format c

Combined-mode usage (emit a single C header for both modes):
  python3 gen_rans_tables.py hints_mode128.json hints_mode256.json \
      --prob-bits 12 --format c --out rans_tables.h
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Dict, List, Tuple


def load_histogram(path: str) -> Tuple[int, str, Dict[int, int]]:
    """Load a histogram JSON. Returns (mode, kind, {value: count}).

    `kind` defaults to "hint_h" for backward compatibility with histograms
    produced before the kind tag was added.
    """
    with open(path) as f:
        data = json.load(f)
    mode = data["mode"]
    kind = data.get("kind", "hint_h")
    h = {int(k): int(v) for k, v in data["histogram"].items()}
    return mode, kind, h


# Map distribution kind -> (macro infix, var infix, human label)
KIND_META = {
    "hint_h":  ("HINT", "hint", "hint h"),
    "z1_high": ("Z1",   "z1",   "z[1..L] HighBits"),
}


def symmetrize(h: Dict[int, int]) -> Dict[int, int]:
    """Average h(v) and h(-v) to enforce a zero-mean symmetric distribution."""
    out = {}
    seen = set()
    for v, c in h.items():
        if v in seen:
            continue
        c_pos = h.get(v, 0)
        c_neg = h.get(-v, 0)
        avg = (c_pos + c_neg) // 2
        if v == 0:
            out[0] = h.get(0, 0)
        else:
            if avg > 0:
                out[v] = avg
                out[-v] = avg
        seen.add(v)
        seen.add(-v)
    # Rebalance: total may have lost 1-2 from integer division.
    return out


def quantize(h: Dict[int, int], prob_bits: int) -> Tuple[List[int], List[int]]:
    """Quantize a histogram to integer frequencies summing to 2^prob_bits.

    Returns (syms, freqs) lists aligned in ascending symbol order. Every
    observed symbol gets at least frequency 1 (so rANS can represent it).
    """
    total_count = sum(h.values())
    prob_total = 1 << prob_bits
    syms = sorted(h.keys())
    real_fs = [h[s] / total_count * prob_total for s in syms]
    # Floor, then top up to reach prob_total, preserving "at least 1" per seen symbol.
    freqs = [max(1, int(round(f))) for f in real_fs]
    diff = prob_total - sum(freqs)
    if diff != 0:
        # Adjust largest bucket to absorb the slack.
        max_idx = max(range(len(freqs)), key=lambda i: freqs[i])
        freqs[max_idx] += diff
        if freqs[max_idx] < 1:
            raise RuntimeError("quantization produced non-positive frequency after slack fix")
    assert sum(freqs) == prob_total, (sum(freqs), prob_total)
    return syms, freqs


def shannon_entropy(h: Dict[int, int]) -> float:
    total = sum(h.values())
    return sum(-c/total * math.log2(c/total) for c in h.values() if c > 0)


def rans_rate(syms: List[int], freqs: List[int], prob_bits: int) -> float:
    """Expected rANS code length assuming p(s) = freqs[s] / 2^prob_bits.
    Note: this may differ from the empirical entropy because of quantization."""
    prob_total = 1 << prob_bits
    return -sum(f / prob_total * math.log2(f / prob_total) for f in freqs if f > 0)


def emit_c_table(mode: int, kind: str, syms: List[int], freqs: List[int],
                 prob_bits: int, file) -> None:
    """Write a C snippet for one mode/kind rANS table.

    For kind="hint_h" we keep the legacy namespace
        shuttle{mode}_rans_{syms,freqs}     / SHUTTLE{MODE}_RANS_{...}
    so existing consumers (shuttle_rans.c, test_rans.c) stay untouched.
    For other kinds we use the kind infix:
        shuttle{mode}_rans_{infix}_{syms,freqs} / SHUTTLE{MODE}_RANS_{INFIX}_{...}
    """
    if kind not in KIND_META:
        raise ValueError(f"unknown kind {kind!r}")
    macro_infix, var_infix, label = KIND_META[kind]
    prob_total = 1 << prob_bits
    if kind == "hint_h":
        macro_pref = f"SHUTTLE{mode}_RANS"
        var_pref = f"shuttle{mode}_rans"
    else:
        macro_pref = f"SHUTTLE{mode}_RANS_{macro_infix}"
        var_pref = f"shuttle{mode}_rans_{var_infix}"
    print(f"/* ================================================================ */", file=file)
    print(f"/* rANS {label} frequency table for SHUTTLE-{mode}                   */", file=file)
    print(f"/* symbols: [{syms[0]}, {syms[-1]}]  total {len(syms)}, sum of freqs = {prob_total} */",
          file=file)
    print(f"/* ================================================================ */", file=file)
    print(f"#define {macro_pref}_PROB_BITS {prob_bits}", file=file)
    print(f"#define {macro_pref}_SYM_MIN  ({syms[0]})", file=file)
    print(f"#define {macro_pref}_SYM_MAX  ({syms[-1]})", file=file)
    print(f"#define {macro_pref}_NUM_SYMS {len(syms)}", file=file)
    print(f"", file=file)
    print(f"static const int16_t {var_pref}_syms[{len(syms)}] = {{", file=file)
    for i in range(0, len(syms), 12):
        line = ", ".join(f"{s:4d}" for s in syms[i:i+12])
        print(f"    {line},", file=file)
    print(f"}};", file=file)
    print(f"", file=file)
    print(f"static const uint16_t {var_pref}_freqs[{len(syms)}] = {{", file=file)
    for i in range(0, len(freqs), 12):
        line = ", ".join(f"{f:5d}" for f in freqs[i:i+12])
        print(f"    {line},", file=file)
    print(f"}};", file=file)
    print(f"", file=file)


def print_text_summary(path: str, mode: int, kind: str, syms: List[int], freqs: List[int],
                       prob_bits: int, raw_hist: Dict[int, int]):
    prob_total = 1 << prob_bits
    print(f"=== {path} (MODE-{mode}, KIND={kind}) ===")
    print(f"Shannon entropy:  {shannon_entropy(raw_hist):6.3f} bits/symbol")
    print(f"rANS rate (fit):  {rans_rate(syms, freqs, prob_bits):6.3f} bits/symbol")
    print(f"Naive fixed:      {math.ceil(math.log2(2 * max(abs(s) for s in syms) + 1)):6d} bits/symbol "
          f"(covers [{min(syms)}, {max(syms)}])")
    print(f"Quantization: {len(syms)} symbols, total freqs = {prob_total}")
    print(f"Top entries:")
    sorted_ = sorted(zip(syms, freqs), key=lambda t: -t[1])
    for v, f in sorted_[:12]:
        print(f"  sym={v:+4d}  freq={f:5d}  p={f/prob_total:.4f}")


def main():
    ap = argparse.ArgumentParser(description="rANS frequency table generator")
    ap.add_argument("histograms", nargs="+",
                    help="JSON files produced by sample_hints.py")
    ap.add_argument("--prob-bits", type=int, default=12,
                    help="log2(total frequency). default 12 (freqs sum to 4096).")
    ap.add_argument("--format", choices=["text", "c"], default="text")
    ap.add_argument("--out", default=None,
                    help="output path (default stdout; only meaningful for --format c)")
    args = ap.parse_args()

    outf = sys.stdout if args.out is None else open(args.out, "w")

    if args.format == "c":
        print("/* SHUTTLE rANS frequency tables. Auto-generated by "
              "SHUTTLE/tools/gen_rans_tables.py — do not edit. */", file=outf)
        print("#ifndef SHUTTLE_RANS_TABLES_H", file=outf)
        print("#define SHUTTLE_RANS_TABLES_H", file=outf)
        print("#include <stdint.h>", file=outf)
        print("", file=outf)

    for path in args.histograms:
        mode, kind, raw = load_histogram(path)
        sym = symmetrize(raw)
        syms, freqs = quantize(sym, args.prob_bits)
        if args.format == "text":
            print_text_summary(path, mode, kind, syms, freqs, args.prob_bits, raw)
        else:
            emit_c_table(mode, kind, syms, freqs, args.prob_bits, outf)

    if args.format == "c":
        print("#endif /* SHUTTLE_RANS_TABLES_H */", file=outf)

    if outf is not sys.stdout:
        outf.close()
        print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
