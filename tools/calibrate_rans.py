#!/usr/bin/env python3
"""Calibrate SHUTTLE rANS reservation budgets from empirical encode lengths.

Runs K simulated signing rounds per mode, encodes the three rANS blocks
(Z_0, hi(z[1..L]), hint h) with the current static frequency tables, and
reports the per-block length distribution (mean, std, max, empirical
(1 - p_ovf) quantile).

Recommended `SHUTTLE_*_RESERVED_BYTES` per block is:

    R = max( empirical_quantile(1 - p_ovf),
             ceil(mu + z * sigma) )

with z = Phi^{-1}(1 - p_ovf). A safety cushion of +4 B is added (rANS
flush adds up to ~4 B beyond the coef entropy sum). Results are rounded
up to a 2-byte boundary and emitted as a C snippet for direct paste
into `ref/params.h`.

Usage:

    python3 calibrate_rans.py --mode 128 --mode 256 \
        --trials 50000 --p-ovf 4.77e-7 \
        --hist128 hints_mode128.json --hist256 hints_mode256.json \
        --z1hist128 z1_mode128.json  --z1hist256 z1_mode256.json \
        --z0hist128 z0_mode128.json  --z0hist256 z0_mode256.json \
        --out calibration.json

Histograms define the active frequency tables (same inputs as
gen_rans_tables.py); the target p_OOV for vocab padding should match
the value passed to gen_rans_tables.py so the calibration uses the
same vocabulary that will ship.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

from rans import RansTable, encode                              # noqa: E402
from sample_hints import MODE_PARAMS, simulate_one_round        # noqa: E402
from sample_z0 import sample_z0_one_round                       # noqa: E402
from sample_z1 import sample_z1_one_round                       # noqa: E402
from gen_rans_tables import (load_histogram, symmetrize,        # noqa: E402
                             quantize, target_vocab_radius,
                             _normal_quantile)


def build_table(histogram_path: str, p_oov: float,
                prob_bits: int = 12) -> Tuple[RansTable, dict]:
    """Load a histogram and return the rANS table + its params dict."""
    _, kind, raw, params = load_histogram(histogram_path)
    target_m = target_vocab_radius(raw, kind, params, p_oov)
    sym = symmetrize(raw, target_radius=target_m)
    syms, freqs = quantize(sym, prob_bits)
    return RansTable(syms=syms, freqs=freqs, prob_bits=prob_bits), params


def encode_length(table: RansTable, syms: List[int]) -> int:
    """Return the rANS-encoded byte length for a symbol list."""
    return len(encode(syms, table))


def _align_up(x: int, step: int) -> int:
    return ((x + step - 1) // step) * step


def quantile(arr: List[int], q: float) -> float:
    """Empirical (1 - ...) quantile by np."""
    return float(np.quantile(np.asarray(arr), q, method="higher"))


def calibrate_mode(mode: int, trials: int, p_ovf: float, p_oov: float,
                   hist_paths: Dict[str, str], prob_bits: int = 12,
                   seed: int = 0x5CA1AB1E) -> dict:
    """Run `trials` signing rounds for one mode; return a stats dict."""
    params = MODE_PARAMS[mode]
    rng = np.random.default_rng(seed ^ mode)

    tbl_z0,   _ = build_table(hist_paths["z0"],   p_oov, prob_bits)
    tbl_z1hi, _ = build_table(hist_paths["z1"],   p_oov, prob_bits)
    tbl_h,    _ = build_table(hist_paths["hint"], p_oov, prob_bits)

    oov_z0 = oov_z1 = oov_h = 0
    lens_z0, lens_z1, lens_h = [], [], []

    progress_every = max(1, trials // 20)
    for t in range(trials):
        z0 = sample_z0_one_round(rng, params)
        z1_hi = sample_z1_one_round(rng, params)
        h = simulate_one_round(rng, params)

        if any(int(v) not in tbl_z0.sym_index for v in z0):
            oov_z0 += 1
        else:
            lens_z0.append(encode_length(tbl_z0, [int(v) for v in z0]))

        z1_flat = [int(v) for v in z1_hi.ravel()]
        if any(v not in tbl_z1hi.sym_index for v in z1_flat):
            oov_z1 += 1
        else:
            lens_z1.append(encode_length(tbl_z1hi, z1_flat))

        h_flat = [int(v) for v in h.ravel()]
        if any(v not in tbl_h.sym_index for v in h_flat):
            oov_h += 1
        else:
            lens_h.append(encode_length(tbl_h, h_flat))

        if (t + 1) % progress_every == 0:
            print(f"  mode-{mode} [{t + 1}/{trials}]")

    z = _normal_quantile(1.0 - p_ovf)
    result = {"mode": mode, "trials": trials, "p_ovf": p_ovf, "p_oov": p_oov,
              "z_ovf": z, "blocks": {}}

    for name, lens, oov in [("z0", lens_z0, oov_z0),
                            ("z1", lens_z1, oov_z1),
                            ("h",  lens_h,  oov_h)]:
        if not lens:
            raise RuntimeError(f"mode-{mode} block {name}: every sample hit OOV")
        arr = np.asarray(lens)
        mu = float(arr.mean())
        sigma = float(arr.std(ddof=0))
        mx = int(arr.max())
        q_emp = quantile(lens, 1.0 - p_ovf) if len(arr) > 0 else mx
        r_clt = int(math.ceil(mu + z * sigma)) + 4   # +4 safety margin
        r = max(mx, int(math.ceil(q_emp)), r_clt)
        r_aligned = _align_up(r, 2)
        result["blocks"][name] = {
            "mean": mu, "std": sigma, "max": mx,
            "quantile": q_emp, "r_clt": r_clt,
            "R": r_aligned,
            "oov": oov, "samples": len(lens),
            "vocab_M": max(abs(s) for s in (
                tbl_z0 if name == "z0" else
                tbl_z1hi if name == "z1" else tbl_h).syms),
        }
    return result


def emit_params_snippet(results: Dict[int, dict], out) -> None:
    print("/* ============================================================", file=out)
    print(" * rANS reservation budgets. Auto-calibrated by", file=out)
    print(" *   SHUTTLE/tools/calibrate_rans.py", file=out)
    print(" * Per-block failure target: p_block = p_ovf + p_OOV.", file=out)
    for mode in sorted(results):
        r = results[mode]
        print(f" *   mode-{mode}: p_ovf = {r['p_ovf']:.2e}, "
              f"p_OOV = {r['p_oov']:.2e}, trials = {r['trials']}",
              file=out)
        for name in ("z0", "z1", "h"):
            b = r["blocks"][name]
            print(f" *     {name}: mean {b['mean']:.1f} B, "
                  f"std {b['std']:.2f} B, max {b['max']} B, "
                  f"empirical q {b['quantile']:.1f} B -> R = {b['R']} B",
                  file=out)
    print(" * ============================================================ */",
          file=out)
    first = True
    for mode in sorted(results):
        blk = results[mode]["blocks"]
        guard = "#if" if first else "#elif"
        first = False
        print(f"{guard} SHUTTLE_MODE == {mode}", file=out)
        print(f"#  define SHUTTLE_HINT_RESERVED_BYTES    {blk['h']['R']}",
              file=out)
        print(f"#  define SHUTTLE_Z1_RANS_RESERVED_BYTES {blk['z1']['R']}",
              file=out)
        print(f"#  define SHUTTLE_Z0_RANS_RESERVED_BYTES {blk['z0']['R']}",
              file=out)
    print("#endif", file=out)


def main():
    ap = argparse.ArgumentParser(description="Calibrate rANS reservation budgets")
    ap.add_argument("--mode", type=int, action="append", choices=[128, 256],
                    required=True, help="repeatable")
    ap.add_argument("--trials", type=int, default=50000)
    ap.add_argument("--p-ovf", type=float, default=4.77e-7,
                    help="per-block overflow probability target; "
                         "4.77e-7 ~= 2^{-21} (for p_block = 2^{-20})")
    ap.add_argument("--p-oov", type=float, default=4.77e-7,
                    help="per-block OOV probability target "
                         "(must match --p-oov used for gen_rans_tables.py)")
    ap.add_argument("--prob-bits", type=int, default=12)
    ap.add_argument("--hist128",   default=None)
    ap.add_argument("--hist256",   default=None)
    ap.add_argument("--z1hist128", default=None)
    ap.add_argument("--z1hist256", default=None)
    ap.add_argument("--z0hist128", default=None)
    ap.add_argument("--z0hist256", default=None)
    ap.add_argument("--out", default=None,
                    help="write JSON stats + C snippet to this path; "
                         "also echoes the C snippet to stdout")
    ap.add_argument("--seed", type=int, default=0x5CA1AB1E)
    args = ap.parse_args()

    paths = {
        128: {"hint": args.hist128, "z1": args.z1hist128, "z0": args.z0hist128},
        256: {"hint": args.hist256, "z1": args.z1hist256, "z0": args.z0hist256},
    }

    results: Dict[int, dict] = {}
    for mode in sorted(set(args.mode)):
        missing = [k for k, v in paths[mode].items() if v is None]
        if missing:
            raise SystemExit(f"missing histogram args for mode-{mode}: {missing}")
        print(f"Calibrating mode-{mode} with {args.trials} trials ...")
        results[mode] = calibrate_mode(
            mode, args.trials, args.p_ovf, args.p_oov,
            paths[mode], args.prob_bits, args.seed,
        )

    emit_params_snippet(results, sys.stdout)

    if args.out is not None:
        with open(args.out, "w") as f:
            json.dump({str(m): r for m, r in results.items()}, f, indent=2)
        print(f"\nStats JSON written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
