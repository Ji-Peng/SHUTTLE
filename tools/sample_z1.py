#!/usr/bin/env python3
"""Collect empirical HighBits(z[1..L]) histograms for Phase 6c rANS training.

z[i] for i in 1..L satisfies z[i] = y[i] + c_eff * s[i] (mod Z-operations), with
|z[i]|_inf <= 11*sigma + tau*eta. For signature compression, we split each
coefficient into:
    z_high = round(z / alpha_h)
    z_low  = z - alpha_h * z_high        # in (-alpha_h/2, alpha_h/2]
and encode `z_low` directly at log2(alpha_h) bits/coef while rANS-compressing
`z_high`. The resulting per-coef rate is close to the Gaussian entropy
log2(sigma * sqrt(2*pi*e)).

Usage:
    python3 sample_z1.py --mode 128 --trials 300 --out z1_mode128.json
"""

from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np

from sample_hints import (MODE_PARAMS, sample_gaussian, sample_cbd,
                           sample_challenge, sample_uniform_modq,
                           poly_mul_modq, _poly_mul_signed)


def sample_z1_one_round(rng: np.random.Generator, params: dict) -> np.ndarray:
    """Simulate one signing round to produce z[1..L] HighBits (L, N) int array."""
    n, q, L, M = params['n'], params['q'], params['L'], params['M']
    eta, tau = params['eta'], params['tau']
    sigma, alpha_h = params['sigma'], params['alpha_h']
    gauss_bound = params['gauss_bound']

    # Small secrets.
    s = np.stack([sample_cbd(rng, n, eta) for _ in range(L)])

    # y[1..L] only; we don't need y[0] or y[L+1..] or public matrix for this.
    y_mid = np.stack([sample_gaussian(rng, n, sigma, gauss_bound) for _ in range(L)])

    # Challenge + IRS-idealized z[i] = y[i] + c * s[i] (signed).
    c = sample_challenge(rng, n, tau)
    z_mid = np.empty((L, n), dtype=np.int64)
    for i in range(L):
        z_mid[i] = y_mid[i] + _poly_mul_signed(c, s[i])

    # Split via alpha_h.
    z_hi = (z_mid + (alpha_h // 2)) // alpha_h
    return z_hi


def collect_z1_highs(mode: int, trials: int, seed: int):
    params = MODE_PARAMS[mode]
    rng = np.random.default_rng(seed)
    counter = Counter()
    progress_every = max(1, trials // 10)
    for t in range(trials):
        z_hi = sample_z1_one_round(rng, params)
        vals, counts = np.unique(z_hi.ravel(), return_counts=True)
        for v, c in zip(vals, counts):
            counter[int(v)] += int(c)
        if (t + 1) % progress_every == 0:
            print(f"  [{t + 1}/{trials}] samples so far: {sum(counter.values())}")
    return counter, params


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", type=int, choices=[128, 256], required=True)
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260420)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = MODE_PARAMS[args.mode]
    print(f"Sampling z[1..L] HighBits for mode-{args.mode} "
          f"(L={params['L']}, n={params['n']}, sigma={params['sigma']}, "
          f"alpha_h={params['alpha_h']}) over {args.trials} trials ...")
    counter, params = collect_z1_highs(args.mode, args.trials, args.seed)

    n = sum(counter.values())
    min_v, max_v = min(counter), max(counter)
    mean = sum(v * c for v, c in counter.items()) / n
    var = sum((v - mean) ** 2 * c for v, c in counter.items()) / n
    print(f"Collected {n} samples, range [{min_v}, {max_v}], "
          f"mean={mean:.3f}, std={var ** 0.5:.3f}")

    data = {
        "mode": args.mode,
        "kind": "z1_high",
        "params": params,
        "trials": args.trials,
        "total_samples": n,
        "histogram": {str(k): v for k, v in sorted(counter.items())},
    }
    with open(args.out, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
