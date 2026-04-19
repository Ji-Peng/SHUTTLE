#!/usr/bin/env python3
"""Collect empirical Z_0 = CompressY(z^{(0)}) histograms for Phase 6d rANS training.

Z_0 derivation (see docs/NGCC_Sign/SHUTTLE_draft.md, section "Z_0 coefficient
bound"):

   z^{(0)} = y^{(0)} + c_eff * sk_full[0]   (IRS output)
           = y^{(0)} + alpha_1 * c_eff      (because sk_full[0] is the constant
                                             polynomial alpha_1)
   Z_0 = round(z^{(0)} / alpha_1)
       = round(y^{(0)} / alpha_1) + c_eff   (exact because alpha_1 divides)

Each Z_0 coefficient is dominated by round(y / alpha_1) (discrete Gaussian with
effective width sigma/alpha_1), with a +/-1 shift on the tau positions where
c_eff is non-zero. The resulting distribution has low entropy (sigma_eff ~10-13),
so full-symbol rANS gives a large saving over the current 11-bit fixed-width
polyz0_pack.

Usage:
    python3 sample_z0.py --mode 128 --trials 300 --out z0_mode128.json
"""

from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np

from sample_hints import MODE_PARAMS, sample_gaussian


def sample_c_eff(rng: np.random.Generator, n: int, tau: int) -> np.ndarray:
    """Signed sparse challenge: tau positions with +/-1, rest 0.

    Mirrors the effect of SampleInBall (+1 at tau positions) combined with
    IRS signs (+/-1 per non-zero token).
    """
    positions = rng.choice(n, size=tau, replace=False)
    signs = rng.choice([-1, 1], size=tau)
    c = np.zeros(n, dtype=np.int64)
    c[positions] = signs
    return c


def sample_z0_one_round(rng: np.random.Generator, params: dict) -> np.ndarray:
    """Draw one Z_0 polynomial (length n) matching the signer's post-CompressY z^{(0)}."""
    n, alpha_1, tau = params['n'], params['alpha_1'], params['tau']
    sigma, gauss_bound = params['sigma'], params['gauss_bound']
    y0 = sample_gaussian(rng, n, sigma, gauss_bound)
    c_eff = sample_c_eff(rng, n, tau)
    # Round-half-up division for both signs, then add c_eff.
    y0_comp = (y0 + (alpha_1 // 2)) // alpha_1
    return y0_comp + c_eff


def collect_z0(mode: int, trials: int, seed: int):
    params = MODE_PARAMS[mode]
    rng = np.random.default_rng(seed)
    counter: Counter = Counter()
    progress_every = max(1, trials // 10)
    for t in range(trials):
        z0 = sample_z0_one_round(rng, params)
        vals, counts = np.unique(z0, return_counts=True)
        for v, c in zip(vals, counts):
            counter[int(v)] += int(c)
        if (t + 1) % progress_every == 0:
            print(f"  [{t + 1}/{trials}] samples so far: {sum(counter.values())}")
    return counter, params


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", type=int, choices=[128, 256], required=True)
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260421)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = MODE_PARAMS[args.mode]
    print(f"Sampling Z_0 for mode-{args.mode} "
          f"(n={params['n']}, sigma={params['sigma']}, alpha_1={params['alpha_1']}, "
          f"tau={params['tau']}) over {args.trials} trials ...")
    counter, params = collect_z0(args.mode, args.trials, args.seed)

    n_tot = sum(counter.values())
    min_v, max_v = min(counter), max(counter)
    mean = sum(v * c for v, c in counter.items()) / n_tot
    var = sum((v - mean) ** 2 * c for v, c in counter.items()) / n_tot
    print(f"Collected {n_tot} samples, range [{min_v}, {max_v}], "
          f"mean={mean:.3f}, std={var ** 0.5:.3f}")

    data = {
        "mode": args.mode,
        "kind": "z0",
        "params": params,
        "trials": args.trials,
        "total_samples": n_tot,
        "histogram": {str(k): v for k, v in sorted(counter.items())},
    }
    with open(args.out, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
