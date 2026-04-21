#!/usr/bin/env python3
"""Generate static rANS frequency tables from empirical hint histograms.

Reads a hint histogram JSON produced by sample_hints.py, quantizes to
integer frequencies summing to 2^prob_bits (typically 2^12), and emits
a C header snippet embedding the table.

Quantization procedure (safe against zero-frequency symbols):

  1. Enforce symmetry: replace p(v) with (p(v) + p(-v)) / 2.
  2. Optionally extend the vocabulary to [-M, M] for a target per-block
     OOV probability p_OOV (see SHUTTLE_rANS_analysis.md eq. 23).
  3. Compute exact real frequencies  f(v) = p(v) * PROB_TOTAL.
  4. Round each f(v) to nearest integer, minimum 1 for observed symbols.
  5. Adjust the single largest bucket's frequency up or down so the sum
     equals PROB_TOTAL exactly.

Output format is controlled by --format:
  --format text : human-readable table + stats
  --format c    : emits C header text on stdout

Usage:
  python3 gen_rans_tables.py hints_mode128.json --prob-bits 12 --format c

Combined-mode usage (emit a single C header for both modes):
  python3 gen_rans_tables.py hints_mode128.json hints_mode256.json \
      --prob-bits 12 --p-oov 4.77e-7 --format c --out rans_tables.h
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Dict, List, Optional, Tuple


# Per-kind block size (#samples per signature round) for OOV-probability math.
# n_block enters eq. (23) as: M >= sigma_eff * Phi^{-1}(1 - p/(2 n_b)).
#   hint_h : M * n  per round
#   z1_high: L * n  per round
#   z0     :     n  per round
def _n_block(kind: str, params: dict) -> int:
    n = params["n"]
    if kind == "hint_h":
        return params["M"] * n
    if kind == "z1_high":
        return params["L"] * n
    if kind == "z0":
        return n
    raise ValueError(f"unknown kind {kind!r}")


def load_histogram(path: str) -> Tuple[int, str, Dict[int, int], dict]:
    """Load a histogram JSON. Returns (mode, kind, {value: count}, params).

    `kind` defaults to "hint_h" for backward compatibility with histograms
    produced before the kind tag was added. `params` holds the recorded
    MODE_PARAMS dict (n/L/M/tau/...), used for OOV-probability math.
    """
    with open(path) as f:
        data = json.load(f)
    mode = data["mode"]
    kind = data.get("kind", "hint_h")
    h = {int(k): int(v) for k, v in data["histogram"].items()}
    params = data.get("params", {})
    return mode, kind, h, params


def _empirical_sigma(h: Dict[int, int]) -> float:
    """Estimate the effective standard deviation of the histogram (zero-mean)."""
    total = sum(h.values())
    if total == 0:
        return 0.0
    mean = sum(v * c for v, c in h.items()) / total
    var = sum((v - mean) ** 2 * c for v, c in h.items()) / total
    return math.sqrt(var)


def _normal_quantile(p: float) -> float:
    """Phi^{-1}(p) via Beasley-Springer-Moro (no scipy dependency)."""
    if not (0.0 < p < 1.0):
        raise ValueError(f"quantile out of range: {p}")
    a = [-3.969683028665376e+01,  2.209460984245205e+02,
         -2.759285104469687e+02,  1.383577518672690e+02,
         -3.066479806614716e+01,  2.506628277459239e+00]
    b = [-5.447609879822406e+01,  1.615858368580409e+02,
         -1.556989798598866e+02,  6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
          4.374664141464968e+00,  2.938163982698783e+00]
    d = [ 7.784695709041462e-03,  3.224671290700398e-01,
          2.445134137142996e+00,  3.754408661907416e+00]
    plow = 0.02425
    phigh = 1 - plow
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p <= phigh:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q / \
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)


def target_vocab_radius(h: Dict[int, int], kind: str, params: dict,
                        p_oov: float) -> int:
    """Required half-width M to meet a per-block OOV probability target.

    Uses the Gaussian-approximation of the distribution with empirical sigma
    and block size n_b = _n_block(kind, params). See analysis doc eq. (23):

        M >= sigma_eff * Phi^{-1}(1 - p_oov / (2 n_b))

    Falls back to the observed max |v| if params are missing (old JSONs).
    """
    observed = max(abs(v) for v in h) if h else 0
    if not params or p_oov is None or p_oov <= 0:
        return observed
    n_b = _n_block(kind, params)
    sigma = _empirical_sigma(h)
    if sigma <= 0:
        return observed
    q = _normal_quantile(1.0 - p_oov / (2.0 * n_b))
    return max(observed, int(math.ceil(sigma * q)))


# Map distribution kind -> (macro infix, var infix, human label)
KIND_META = {
    "hint_h":  ("HINT", "hint", "hint h"),
    "z1_high": ("Z1",   "z1",   "z[1..L] HighBits"),
    "z0":      ("Z0",   "z0",   "Z_0 = CompressY(z^(0))"),
}


def symmetrize(h: Dict[int, int], target_radius: Optional[int] = None) -> Dict[int, int]:
    """Build a zero-mean symmetric, contiguous-vocabulary distribution.

    Procedure:
      1. Average h[v] and h[-v] (ties on rare tail symbols default to 1).
      2. Ensure every integer in [-max_abs, max_abs] appears, filling gaps
         with count 1. Keeps the symbol list contiguous so the C decoder's
         `sym - sym_min` indexing stays valid and avoids false OOV rejections
         on values that happen to be unseen in the training window but land
         inside the observed range.
      3. If `target_radius` is given and larger than the observed max |v|,
         extend the vocabulary symmetrically with count=1 padding. This is
         how the generator meets a per-block p_OOV target: eq. (23) gives
         the minimum radius, and unseen-but-reserved tail symbols land on
         freq=1 after quantization (entropy cost < 0.01 bit/coef).
    """
    out = {}
    seen = set()
    for v in h:
        if v in seen:
            continue
        c_pos = h.get(v, 0)
        c_neg = h.get(-v, 0)
        avg = (c_pos + c_neg) // 2
        if v == 0:
            out[0] = h.get(0, 0)
        else:
            if c_pos > 0 or c_neg > 0:
                out[v] = max(1, avg)
                out[-v] = max(1, avg)
        seen.add(v)
        seen.add(-v)

    if not out:
        return out
    max_abs = max(abs(v) for v in out)
    if target_radius is not None and target_radius > max_abs:
        max_abs = target_radius
    for v in range(-max_abs, max_abs + 1):
        out.setdefault(v, 1)
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
    ap.add_argument("--p-oov", type=float, default=None,
                    help="target per-block OOV probability. When given, the "
                         "vocabulary is extended symmetrically with freq=1 "
                         "padding until M >= sigma_eff * Phi^{-1}(1-p/(2 n_b)). "
                         "Typical: 4.77e-7 (~2^{-21}) for p_block = 2^{-20}.")
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
        mode, kind, raw, params = load_histogram(path)
        target_m = target_vocab_radius(raw, kind, params, args.p_oov) \
            if args.p_oov is not None else None
        sym = symmetrize(raw, target_radius=target_m)
        syms, freqs = quantize(sym, args.prob_bits)
        if args.format == "text":
            print_text_summary(path, mode, kind, syms, freqs, args.prob_bits, raw)
            if target_m is not None:
                observed = max(abs(v) for v in raw) if raw else 0
                print(f"OOV target: p={args.p_oov:.2e}, observed M={observed}, "
                      f"required M={target_m}, final M={max(abs(s) for s in syms)}")
        else:
            emit_c_table(mode, kind, syms, freqs, args.prob_bits, outf)

    if args.format == "c":
        print("#endif /* SHUTTLE_RANS_TABLES_H */", file=outf)

    if outf is not sys.stdout:
        outf.close()
        print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
