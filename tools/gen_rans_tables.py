#!/usr/bin/env python3
"""Generate static rANS frequency tables from the *theoretical* discrete
Gaussian PMF — no empirical histogram required.

Background (see agent/rANS/SHUTTLE_rANS.tex §2.3 and §3):

  SHUTTLE has at most two rANS contexts per mode:

    z-hi  : sigma = r/alpha_r,   M_voc = ceil((11*r + tau*eta)/alpha_r)
    hint  : sigma = 2r/alpha_h,  M_voc = floor(2*(11*r + tau*eta)/alpha_h) + 1

  SampleY hard-truncates |y_k| <= 11*sigma_y, so the bounds above hold
  with *probability 1*. The vocabulary covers [-M_voc, M_voc] in full,
  hence OOV is mathematically impossible.

  For mode-128 the coincidence alpha_h = 2*alpha_r gives sigma_hint =
  sigma_zhi, so a single table doubles as both z-hi and hint context.
  For mode-256/512 the two scales differ; two tables are emitted.

Quantization procedure:

  1. PMF eval: p_tilde(s) = exp(-s^2 / (2 sigma^2)), s in [-M_voc, M_voc].
     Normalize to a probability mass function over the finite alphabet.
  2. Round: g(s) = max(1, round(p(s) * 2^t)) with t = 10. Every alphabet
     symbol receives at least one slot.
  3. Slack: adjust the largest bucket so sum(g) = 2^t exactly. The PMF
     is symmetric around 0 so the largest bucket is unique (= s=0).

Usage:
  python3 gen_rans_tables.py --out ../ref/rans_tables.h
  python3 gen_rans_tables.py --mode 128 --format text   # human inspection
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import Dict, List, Tuple


# SHUTTLE-Spec/main.tex Table 2 parameters that drive the table geometry.
MODE_PARAMS: Dict[int, dict] = {
    128: {"n": 256,  "lenS": 3, "lenE": 2, "eta": 1, "tau":  30,
          "r": 101,  "alpha_1":  8, "alpha_h":  128},
    256: {"n": 512,  "lenS": 3, "lenE": 2, "eta": 1, "tau":  58,
          "r": 149,  "alpha_1": 16, "alpha_h": 1024},
    512: {"n": 1024, "lenS": 3, "lenE": 2, "eta": 1, "tau": 115,
          "r": 202,  "alpha_1": 16, "alpha_h": 2048},
}

ALPHA_R = 64         # cross-mode constant (tex Tab 1)
T_SIGMA = 11         # SampleY truncation multiplier
PROB_BITS = 10       # quantization bits (tex §2.2)


def theoretical_pmf(sigma: float, M_voc: int) -> Dict[int, float]:
    """Discrete Gaussian PMF on the finite alphabet [-M_voc, M_voc]."""
    weights = {k: math.exp(-k * k / (2 * sigma * sigma))
               for k in range(-M_voc, M_voc + 1)}
    Z = sum(weights.values())
    return {k: w / Z for k, w in weights.items()}


def quantize(pmf: Dict[int, float], prob_bits: int) -> Tuple[List[int], List[int]]:
    """Quantize a PMF to integer frequencies summing to 2^prob_bits.

    Every symbol receives at least one slot (so it remains representable).
    The largest bucket absorbs the rounding slack.
    """
    prob_total = 1 << prob_bits
    syms = sorted(pmf.keys())
    real_fs = [pmf[s] * prob_total for s in syms]
    freqs = [max(1, int(round(f))) for f in real_fs]
    diff = prob_total - sum(freqs)
    if diff != 0:
        max_idx = max(range(len(freqs)), key=lambda i: freqs[i])
        freqs[max_idx] += diff
        if freqs[max_idx] < 1:
            raise RuntimeError("quantization produced non-positive freq")
    assert sum(freqs) == prob_total
    return syms, freqs


def zhi_table(params: dict) -> Tuple[List[int], List[int], float, int]:
    """z-hi context: sigma = r/alpha_r, M_voc = ceil((11r + tau*eta)/alpha_r).

    The tight bound covers both HighBits_{alpha_0'}(z^(0)) (Cor 3) and
    HighBits_{alpha_r}(z^(i)) for i >= 1 (Thm 5). Empirically Cor 3 may
    add an alpha_1-quantization +1 to z^(0); we take the max here so the
    shared vocabulary subsumes both cases.
    """
    r = params["r"]
    tau, eta = params["tau"], params["eta"]
    sigma = r / ALPHA_R
    M_voc = math.ceil((T_SIGMA * r + tau * eta) / ALPHA_R)
    # Add the +1 / alpha_0' robustness margin from Cor 3:
    M_voc = max(M_voc, math.ceil((math.ceil(T_SIGMA * r / params["alpha_1"]) + 1)
                                 / (ALPHA_R // params["alpha_1"])))
    syms, freqs = quantize(theoretical_pmf(sigma, M_voc), PROB_BITS)
    return syms, freqs, sigma, M_voc


def hint_table(params: dict) -> Tuple[List[int], List[int], float, int]:
    """hint context: sigma = 2r/alpha_h, M_voc per Thm 7."""
    r = params["r"]
    alpha_h = params["alpha_h"]
    tau, eta = params["tau"], params["eta"]
    sigma = 2 * r / alpha_h
    M_voc = (2 * (T_SIGMA * r + tau * eta)) // alpha_h + 1
    M_voc = max(M_voc, 1)
    syms, freqs = quantize(theoretical_pmf(sigma, M_voc), PROB_BITS)
    return syms, freqs, sigma, M_voc


def shares_table(params: dict) -> bool:
    """True iff sigma_hint == sigma_zhi (mode-128 coincidence)."""
    return params["alpha_h"] == 2 * ALPHA_R


def emit_table_block(mode: int, infix_upper: str, infix_lower: str,
                     syms: List[int], freqs: List[int],
                     prob_bits: int, sigma: float, file) -> None:
    """Write a labeled rANS frequency block to file.

    infix_upper = "ZHI" / "HINT" / "UNIFIED"     (macro infix)
    infix_lower = "zhi" / "hint" / "unified"     (variable infix)
    """
    prob_total = 1 << prob_bits
    macro_pref = f"SHUTTLE{mode}_RANS_{infix_upper}"
    var_pref = f"shuttle{mode}_rans_{infix_lower}"
    print(f"/* ================================================================ */",
          file=file)
    print(f"/* SHUTTLE-{mode} rANS {infix_lower} frequency table.                 */",
          file=file)
    print(f"/*   sigma = {sigma:.4f}, M_voc = {syms[-1]} (alphabet [-M_voc, M_voc]) */",
          file=file)
    print(f"/*   prob_bits = {prob_bits}, sum(freqs) = {prob_total}.                */",
          file=file)
    print(f"/* Generated by tools/gen_rans_tables.py from the theoretical PMF.    */",
          file=file)
    print(f"/* See SHUTTLE_rANS.tex §2.3 (PMF) and §3 (M_voc tight bound).        */",
          file=file)
    print(f"/* ================================================================ */",
          file=file)
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


def emit_header(modes: List[int], file) -> None:
    print("/* SHUTTLE rANS frequency tables.", file=file)
    print(" * Auto-generated by SHUTTLE/tools/gen_rans_tables.py.", file=file)
    print(" *", file=file)
    print(" * Tables come from the theoretical discrete Gaussian PMF.", file=file)
    print(" * The alphabet is the |y|_inf <= 11*sigma tight bound (Cor 3 / Thm 5", file=file)
    print(" * / Thm 7 of SHUTTLE_rANS.tex), so no signature coefficient can fall", file=file)
    print(" * outside the vocabulary; OOV is mathematically impossible.", file=file)
    print(" *", file=file)
    print(" * Two contexts per mode in general:", file=file)
    print(" *   z-hi  : carries HighBits_{alpha_0'}(z^(0)) and HighBits_{alpha_r}(z^(1..lenS))", file=file)
    print(" *   hint  : carries the MakeHint output", file=file)
    print(" * For mode-128 alpha_h == 2*alpha_r, so sigma_hint == sigma_zhi and a", file=file)
    print(" * single shared table is emitted (named *_unified) as well as redirect", file=file)
    print(" * macros that point the {zhi,hint} names at the same data.", file=file)
    print(" *", file=file)
    print(" * Do not edit by hand — regenerate with:", file=file)
    print(" *   python3 SHUTTLE/tools/gen_rans_tables.py --out SHUTTLE/ref/rans_tables.h", file=file)
    print(" */", file=file)
    print("#ifndef SHUTTLE_RANS_TABLES_H", file=file)
    print("#define SHUTTLE_RANS_TABLES_H", file=file)
    print("#include <stdint.h>", file=file)
    print("", file=file)


def emit_footer(file) -> None:
    print("#endif /* SHUTTLE_RANS_TABLES_H */", file=file)


def emit_unified_aliases(mode: int, file) -> None:
    """For mode-128: alias ZHI and HINT macros to UNIFIED."""
    print(f"/* mode-{mode} coincidence: sigma_hint == sigma_zhi, both contexts", file=file)
    print(f" * share the unified table above. */", file=file)
    for sub in ("PROB_BITS", "SYM_MIN", "SYM_MAX", "NUM_SYMS"):
        print(f"#define SHUTTLE{mode}_RANS_ZHI_{sub}  SHUTTLE{mode}_RANS_UNIFIED_{sub}",
              file=file)
    for sub in ("PROB_BITS", "SYM_MIN", "SYM_MAX", "NUM_SYMS"):
        print(f"#define SHUTTLE{mode}_RANS_HINT_{sub} SHUTTLE{mode}_RANS_UNIFIED_{sub}",
              file=file)
    print(f"#define shuttle{mode}_rans_zhi_syms   shuttle{mode}_rans_unified_syms",
          file=file)
    print(f"#define shuttle{mode}_rans_zhi_freqs  shuttle{mode}_rans_unified_freqs",
          file=file)
    print(f"#define shuttle{mode}_rans_hint_syms  shuttle{mode}_rans_unified_syms",
          file=file)
    print(f"#define shuttle{mode}_rans_hint_freqs shuttle{mode}_rans_unified_freqs",
          file=file)
    print("", file=file)


def text_summary(mode: int, params: dict) -> None:
    print(f"=== SHUTTLE-{mode} ===")
    z_syms, z_freqs, z_sigma, z_M = zhi_table(params)
    h_syms, h_freqs, h_sigma, h_M = hint_table(params)
    print(f"z-hi : sigma = {z_sigma:.4f}, M_voc = {z_M}, alphabet size = {len(z_syms)}")
    print(f"hint : sigma = {h_sigma:.4f}, M_voc = {h_M}, alphabet size = {len(h_syms)}")
    if shares_table(params):
        print("(mode-128: hint table aliased to z-hi table — sigmas coincide)")
    print(f"z-hi top buckets:")
    for s, f in sorted(zip(z_syms, z_freqs), key=lambda t: -t[1])[:8]:
        print(f"  s={s:+4d} g={f:5d} p={f/1024:.4f}")
    print(f"hint top buckets:")
    for s, f in sorted(zip(h_syms, h_freqs), key=lambda t: -t[1])[:8]:
        print(f"  s={s:+4d} g={f:5d} p={f/1024:.4f}")


def main() -> None:
    ap = argparse.ArgumentParser(description="rANS theoretical-PMF table generator")
    ap.add_argument("--mode", type=int, action="append", choices=[128, 256, 512],
                    help="restrict to one mode; default emits all three")
    ap.add_argument("--format", choices=["text", "c"], default="c")
    ap.add_argument("--out", default=None,
                    help="output path (default stdout; only meaningful for --format c)")
    args = ap.parse_args()

    modes = sorted(set(args.mode)) if args.mode else [128, 256, 512]

    if args.format == "text":
        for m in modes:
            text_summary(m, MODE_PARAMS[m])
        return

    outf = sys.stdout if args.out is None else open(args.out, "w")
    emit_header(modes, outf)
    for m in modes:
        params = MODE_PARAMS[m]
        if shares_table(params):
            # mode-128: emit unified table + aliases for ZHI/HINT.
            syms, freqs, sigma, _ = zhi_table(params)
            emit_table_block(m, "UNIFIED", "unified", syms, freqs,
                             PROB_BITS, sigma, outf)
            emit_unified_aliases(m, outf)
        else:
            z_syms, z_freqs, z_sigma, _ = zhi_table(params)
            h_syms, h_freqs, h_sigma, _ = hint_table(params)
            emit_table_block(m, "ZHI", "zhi", z_syms, z_freqs,
                             PROB_BITS, z_sigma, outf)
            emit_table_block(m, "HINT", "hint", h_syms, h_freqs,
                             PROB_BITS, h_sigma, outf)
    emit_footer(outf)
    if outf is not sys.stdout:
        outf.close()
        print(f"Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
