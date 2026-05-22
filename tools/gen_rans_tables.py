#!/usr/bin/env python3
"""Generate static rANS frequency tables from the *theoretical* discrete
Gaussian PMF — no empirical histogram required.

Background (see agent/rANS/SHUTTLE_rANS.tex §2.3, §3, and the hint
errata in §4.5):

  SHUTTLE has at most two rANS contexts per mode:

    z-hi  : "smooth" discrete Gaussian sigma_zhi = r/alpha_r,
            M_voc = ceil((11*r + tau*eta)/alpha_r)
    hint  : output of MakeHint = round(w/alpha_h) - round((w-2z)/alpha_h),
            *not* a discrete Gaussian. See `hint_pmf_correct` below for
            the actual distribution.
            M_voc = floor(2*(11*r + tau*eta)/alpha_h) + 1

  SampleY hard-truncates |y_k| <= 11*sigma_y, so the bounds above hold
  with *probability 1*. The vocabulary covers [-M_voc, M_voc] in full,
  hence OOV is mathematically impossible.

  For mode-128 we cannot share a single table any more (z-hi and hint
  have different shapes even when their scales coincide — hint is wider
  due to the bucket-crossing effect described in the errata). Each mode
  emits two distinct tables: ZHI and HINT.

Quantization procedure (same for both contexts):

  1. PMF eval:
       z-hi: p(s) ~ exp(-s^2 / (2 sigma^2))  on [-M_voc, M_voc].
       hint: bucket-crossing PMF over the same domain (see
             hint_pmf_correct).
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
    """Discrete Gaussian PMF on the finite alphabet [-M_voc, M_voc].

    Used for the z-hi context (which IS approximately discrete Gaussian).
    """
    weights = {k: math.exp(-k * k / (2 * sigma * sigma))
               for k in range(-M_voc, M_voc + 1)}
    Z = sum(weights.values())
    return {k: w / Z for k, w in weights.items()}


def hint_pmf_correct(r: float, alpha_h: int, M_voc: int,
                     T_sigma: float = 11.0) -> Dict[int, float]:
    """Exact MakeHint output PMF, derived from the bucket-crossing model.

    Derivation (SHUTTLE_rANS.tex §4.5, eq:hint-bernoulli):

      hint_k = round(comY_k / alpha_h) - round((comY_k - 2*z_2,k) / alpha_h)
             = round(u) - round(u - delta)        with  u = comY_k/alpha_h,
                                                       delta = 2*z_2,k/alpha_h.

      For comY_k uniform in [0, 2q) (NTT-mixed commitment), u mod 1 is
      uniform in [0,1) and independent of delta. Decompose
      delta = k + f with k = round_half_up(delta), f in [-1/2, 1/2).
      For u ~ Uniform[0,1):

          P(hint = k | delta)              = 1 - |f|
          P(hint = k + sign(f) | delta)    = |f|

      Averaging over z_2,k ~ D_{Z, r} truncated to |z| <= T_sigma * r
      (the SampleY truncation propagates to z_2 via the IRS).

    This is the CORRECT model. The previous documentation's "hint ~
    discrete Gaussian D_{Z, 2r/alpha_h}" model severely under-predicts
    the nonzero rate in narrow-sigma modes (mode-256/512). See the
    errata in SHUTTLE_rANS.tex §4.5.
    """
    M_z = int(T_sigma * r)
    # Truncated discrete Gaussian on z.
    z_weights = {z: math.exp(-z * z / (2 * r * r))
                 for z in range(-M_z, M_z + 1)}
    Z = sum(z_weights.values())
    z_pmf = {z: w / Z for z, w in z_weights.items()}

    pmf: Dict[int, float] = {h: 0.0 for h in range(-M_voc, M_voc + 1)}
    for z, pz in z_pmf.items():
        delta = 2.0 * z / alpha_h
        # round-half-up nearest integer (matches our highbits_mod_2q convention).
        k = int(math.floor(delta + 0.5))
        f = delta - k
        # Clamp to vocabulary; the bound from Thm 7 guarantees this is safe.
        def add(h: int, mass: float) -> None:
            if -M_voc <= h <= M_voc:
                pmf[h] += mass
            else:
                # Should be unreachable for |z| <= T_sigma*r and the
                # M_voc from Thm 7; flag if it happens.
                raise ValueError(f"hint value {h} out of vocab [-{M_voc}, {M_voc}]")
        add(k, pz * (1.0 - abs(f)))
        if abs(f) > 1e-15:
            sign = 1 if f > 0 else -1
            add(k + sign, pz * abs(f))

    # Numerical hygiene: re-normalize (rounding errors aside).
    total = sum(pmf.values())
    return {h: p / total for h, p in pmf.items()}


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
    """hint context: bucket-crossing PMF per SHUTTLE_rANS.tex §4.5.

    Returns (syms, freqs, sigma_nominal, M_voc). The "sigma_nominal" is
    the legacy 2r/alpha_h value reported for compatibility; it is NOT the
    sigma of the actual PMF (which isn't a discrete Gaussian).
    """
    r = params["r"]
    alpha_h = params["alpha_h"]
    tau, eta = params["tau"], params["eta"]
    sigma_nominal = 2 * r / alpha_h
    M_voc = (2 * (T_SIGMA * r + tau * eta)) // alpha_h + 1
    M_voc = max(M_voc, 1)
    pmf = hint_pmf_correct(r, alpha_h, M_voc)
    syms, freqs = quantize(pmf, PROB_BITS)
    return syms, freqs, sigma_nominal, M_voc


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
    print(" * Two contexts per mode (always two distinct tables):", file=file)
    print(" *   z-hi  : carries HighBits_{alpha_0'}(z^(0)) and HighBits_{alpha_r}(z^(1..lenS))", file=file)
    print(" *           Distribution: discrete Gaussian D_{Z, r/alpha_r}.", file=file)
    print(" *   hint  : carries MakeHint output.", file=file)
    print(" *           Distribution: bucket-crossing PMF (SHUTTLE_rANS.tex §4.5),", file=file)
    print(" *           NOT a discrete Gaussian. The earlier 'mode-128 hint shares", file=file)
    print(" *           the z-hi table' alias has been removed: even though the", file=file)
    print(" *           nominal sigma_h coincides with sigma_zhi for mode-128, the", file=file)
    print(" *           shapes differ (hint is wider).", file=file)
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


def text_summary(mode: int, params: dict) -> None:
    print(f"=== SHUTTLE-{mode} ===")
    z_syms, z_freqs, z_sigma, z_M = zhi_table(params)
    h_syms, h_freqs, h_sigma, h_M = hint_table(params)
    print(f"z-hi : sigma = {z_sigma:.4f}, M_voc = {z_M}, alphabet size = {len(z_syms)}")
    print(f"hint : sigma_nominal = {h_sigma:.4f}, M_voc = {h_M}, "
          f"alphabet size = {len(h_syms)}")
    print(f"       (hint PMF is bucket-crossing, NOT discrete Gaussian)")
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
