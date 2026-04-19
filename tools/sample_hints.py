#!/usr/bin/env python3
"""Simulate the SHUTTLE signing pipeline to collect empirical hint histograms.

The simulation matches the NGCC-Signature Alg 2 logic at a polynomial level,
using schoolbook arithmetic (no NTT):

    1. Sample a_gen, A_gen uniform mod q.
    2. Sample s, e from CBD(eta).
    3. Compute b = a_gen + A_gen * s + e (mod q).
    4. Sample y component-wise from discrete Gaussian (sigma), truncated at
       +- 11*sigma per the SHUTTLE sampler.
    5. CompressY: Y_0 = round(y[0] / alpha_1).
    6. Compute comY = hat_A * Y (mod 2q).
    7. Sample challenge c: tau random +1 positions.
    8. IRS output (idealized): z = y + c * sk_stretched; then CompressY on z[0].
    9. z_2 = z[L+1..L+M] = y[L+1..L+M] + c * e.
   10. Compute h = MakeHint(z_2, comY).

The simulation is statistical; run with --trials 1000 or more for stable
histograms. Output is a JSON file with per-coefficient hint histograms for
both M-many component positions (or aggregated if --aggregate is set).

Usage:
    python3 sample_hints.py --mode 128 --trials 1000 --out hints_mode128.json
    python3 sample_hints.py --mode 256 --trials 1000 --out hints_mode256.json

Dependencies: numpy (for fast polynomial multiplication).
"""

from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np


# ---------------------------------------------------------------------------
# Mode-dependent parameters (matching params.h exactly).
# ---------------------------------------------------------------------------

MODE_PARAMS = {
    128: dict(n=256, q=13313, L=3, M=2, eta=1, tau=30,
              sigma=101, alpha_h=128, alpha_1=8,
              gauss_bound=11*101),
    256: dict(n=512, q=13313, L=3, M=2, eta=1, tau=58,
              sigma=149, alpha_h=256, alpha_1=16,
              gauss_bound=11*149),
}


# ---------------------------------------------------------------------------
# Sampling helpers
# ---------------------------------------------------------------------------

def sample_gaussian(rng: np.random.Generator, n: int, sigma: float,
                    bound: int) -> np.ndarray:
    """Discrete Gaussian approximation: continuous Gaussian + round, truncated."""
    while True:
        xs = rng.normal(0.0, sigma, size=n).round().astype(np.int64)
        if np.all(np.abs(xs) <= bound):
            return xs
        # Resample out-of-bound coefficients only.
        out = np.abs(xs) > bound
        if not out.any():
            return xs
        xs[out] = 0
        # Slow path: re-sample individual coefficients until in bound.
        for i in np.where(out)[0]:
            while True:
                v = int(round(rng.normal(0.0, sigma)))
                if -bound <= v <= bound:
                    xs[i] = v
                    break
        return xs


def sample_cbd(rng: np.random.Generator, n: int, eta: int) -> np.ndarray:
    """Centered binomial distribution: sum of eta bits - sum of eta bits."""
    a = rng.integers(0, 2, size=(n, eta), dtype=np.int64).sum(axis=1)
    b = rng.integers(0, 2, size=(n, eta), dtype=np.int64).sum(axis=1)
    return a - b


def sample_challenge(rng: np.random.Generator, n: int, tau: int) -> np.ndarray:
    """Sparse +1 challenge: exactly tau coefficients set to +1, rest 0."""
    positions = rng.choice(n, size=tau, replace=False)
    c = np.zeros(n, dtype=np.int64)
    c[positions] = 1
    return c


def sample_uniform_modq(rng: np.random.Generator, n: int, q: int) -> np.ndarray:
    """Uniform [0, q) polynomial."""
    return rng.integers(0, q, size=n, dtype=np.int64)


# ---------------------------------------------------------------------------
# Polynomial arithmetic (schoolbook negacyclic convolution in R_q = Z_q[x]/(x^n+1)).
# ---------------------------------------------------------------------------

def poly_mul_modq(a: np.ndarray, b: np.ndarray, q: int) -> np.ndarray:
    """Negacyclic convolution mod q. a, b are length-n arrays.

    Uses numpy convolution + wrap. Specifically, c[k] = sum_{i+j=k} a[i]*b[j]
    (mod x^n+1). For i+j >= n, the contribution is -a[i]*b[j] to coefficient
    (i+j-n).
    """
    n = len(a)
    # Full linear convolution (length 2n-1).
    full = np.convolve(a, b, mode='full')
    # Wrap the upper half with a sign flip.
    wrapped = full[:n].copy()
    wrapped[:len(full) - n] -= full[n:]
    return wrapped % q


def polyveck_sub_modq(u: np.ndarray, v: np.ndarray, q: int) -> np.ndarray:
    return (u - v) % q


# ---------------------------------------------------------------------------
# mod 2q arithmetic
# ---------------------------------------------------------------------------

def lift_to_2q(u_mod_q: np.ndarray, parity_source: np.ndarray, q: int) -> np.ndarray:
    """Lift from mod q to mod 2q using parity of `parity_source`.

    out = (2 * u_mod_q + q * (parity_source & 1)) mod 2q."""
    parity = parity_source & 1
    r = 2 * u_mod_q + q * parity
    return r % (2 * q)


def highbits_mod_2q(w: np.ndarray, alpha_h: int, hint_max: int) -> np.ndarray:
    """round(w / alpha_h), wrapped to [0, hint_max)."""
    hb = (w + (alpha_h >> 1)) // alpha_h
    hb = np.where(hb >= hint_max, hb - hint_max, hb)
    return hb


def reduce_mod_2q(a: np.ndarray, q: int) -> np.ndarray:
    return a % (2 * q)


# ---------------------------------------------------------------------------
# One SHUTTLE signing simulation round -> returns hint polyveck.
# ---------------------------------------------------------------------------

def simulate_one_round(rng: np.random.Generator, params: dict) -> np.ndarray:
    """Execute one full signing round. Returns hint array of shape (M, N)."""
    n, q, L, M = params['n'], params['q'], params['L'], params['M']
    eta, tau = params['eta'], params['tau']
    sigma, alpha_h, alpha_1 = params['sigma'], params['alpha_h'], params['alpha_1']
    gauss_bound = params['gauss_bound']
    dq = 2 * q
    hint_max = (2 * (q - 1)) // alpha_h

    # ---- key material
    a_gen = np.stack([sample_uniform_modq(rng, n, q) for _ in range(M)])  # (M,N)
    A_gen = np.stack([
        np.stack([sample_uniform_modq(rng, n, q) for _ in range(M)])
        for _ in range(L)
    ])                                                                    # (L,M,N)
    s = np.stack([sample_cbd(rng, n, eta) for _ in range(L)])             # (L,N)
    e = np.stack([sample_cbd(rng, n, eta) for _ in range(M)])             # (M,N)

    # b = a_gen + A_gen * s + e  (mod q)
    b = np.empty_like(a_gen)                                              # (M,N)
    for i in range(M):
        acc = np.zeros(n, dtype=np.int64)
        for j in range(L):
            acc = (acc + poly_mul_modq(A_gen[j, i], s[j], q)) % q
        b[i] = (a_gen[i] + acc + e[i]) % q

    # ---- sample y and compute CompressY, comY
    y = np.stack([
        sample_gaussian(rng, n, sigma, gauss_bound) for _ in range(1 + L + M)
    ])                                                                    # (1+L+M, N)

    Y_0 = (y[0] + (alpha_1 // 2)) // alpha_1         # CompressY on slot 0
    Y = y.copy()
    Y[0] = Y_0

    # U = (a_gen - b) * Y_0 + A_gen * Y[1..L] + Y[1+L..1+L+M]  (mod q)
    U = np.zeros((M, n), dtype=np.int64)
    for i in range(M):
        ab = (a_gen[i] - b[i]) % q
        acc = poly_mul_modq(ab, Y[0], q)
        for j in range(L):
            acc = (acc + poly_mul_modq(A_gen[j, i], Y[1 + j], q)) % q
        U[i] = (acc + Y[1 + L + i]) % q

    # comY = lift(U, Y_0)  (mod 2q)
    comY = np.stack([lift_to_2q(U[i], Y[0], q) for i in range(M)])

    # ---- challenge c, IRS (idealized): z = y + c * sk_stretched
    # For z_2 (used by MakeHint) we want signed integer values (NOT mod q).
    # c and s, e are small (CBD +- eta), so c*s, c*e are small signed polys.
    c = sample_challenge(rng, n, tau)
    z = np.empty_like(y)
    z[0] = y[0] + alpha_1 * c    # sk_stretched[0] = alpha_1 (scalar)
    for i in range(L):
        z[1 + i] = y[1 + i] + _poly_mul_signed(c, s[i])
    for i in range(M):
        z[1 + L + i] = y[1 + L + i] + _poly_mul_signed(c, e[i])

    # ---- CompressY on IRS output: z[0] = round(z'[0] / alpha_1)
    # (Not strictly needed for hint; MakeHint uses z_2 which is z[1+L..].)

    # z_2 = z[1+L..1+L+M-1] (signed, int)
    z_2 = z[1 + L:1 + L + M]                                              # (M,N)

    # ---- MakeHint
    # tilde_w = (comY - 2 * z_2) mod 2q
    tilde_w = reduce_mod_2q(comY - 2 * z_2, q)
    w_h = highbits_mod_2q(comY, alpha_h, hint_max)
    tilde_w_h = highbits_mod_2q(tilde_w, alpha_h, hint_max)
    h = w_h - tilde_w_h                                                   # (M,N)

    # Normalize to centered representative in (-HINT_MAX/2, HINT_MAX/2].
    # This merges the "wrap" cluster (e.g., -207 == +1 mod HINT_MAX) with the
    # near-zero cluster, giving a single-peak distribution amenable to rANS.
    h = np.where(h >  hint_max // 2, h - hint_max, h)
    h = np.where(h <= -hint_max // 2, h + hint_max, h)

    return h


def _poly_mul_signed(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Negacyclic convolution over Z (no mod). Used for small coefficient polys."""
    n = len(a)
    full = np.convolve(a, b, mode='full')
    wrapped = full[:n].copy()
    wrapped[:len(full) - n] -= full[n:]
    return wrapped


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def collect_hints(mode: int, trials: int, seed: int):
    """Run `trials` signing rounds and return a Counter over all hint values."""
    params = MODE_PARAMS[mode]
    rng = np.random.default_rng(seed)

    counter = Counter()
    progress_every = max(1, trials // 10)
    for t in range(trials):
        h = simulate_one_round(rng, params)
        # Flatten the M*N hint coefficients and count.
        vals, counts = np.unique(h.ravel(), return_counts=True)
        for v, c in zip(vals, counts):
            counter[int(v)] += int(c)
        if (t + 1) % progress_every == 0:
            print(f"  [{t + 1}/{trials}] total samples so far: {sum(counter.values())}")
    return counter, params


def main():
    ap = argparse.ArgumentParser(description="SHUTTLE hint histogram sampler")
    ap.add_argument("--mode", type=int, choices=[128, 256], required=True)
    ap.add_argument("--trials", type=int, default=200,
                    help="number of full signing rounds (more = stabler "
                         "histogram; 200 gives ~100-250k samples).")
    ap.add_argument("--seed", type=int, default=20260419)
    ap.add_argument("--out", required=True, help="output JSON path")
    args = ap.parse_args()

    print(f"Sampling mode-{args.mode} ({MODE_PARAMS[args.mode]['n']} coef, "
          f"sigma={MODE_PARAMS[args.mode]['sigma']}, "
          f"alpha_h={MODE_PARAMS[args.mode]['alpha_h']}, "
          f"alpha_1={MODE_PARAMS[args.mode]['alpha_1']}) over {args.trials} trials ...")

    counter, params = collect_hints(args.mode, args.trials, args.seed)

    data = {
        "mode": args.mode,
        "kind": "hint_h",
        "params": params,
        "trials": args.trials,
        "total_samples": sum(counter.values()),
        "histogram": {str(k): v for k, v in sorted(counter.items())},
    }
    with open(args.out, "w") as f:
        json.dump(data, f, indent=2)

    # Summary.
    n_samples = sum(counter.values())
    min_v = min(counter)
    max_v = max(counter)
    mean = sum(v * c for v, c in counter.items()) / n_samples
    var = sum((v - mean) ** 2 * c for v, c in counter.items()) / n_samples
    print(f"Done. {n_samples} hint samples, range [{min_v}, {max_v}], "
          f"mean={mean:.3f}, std={var**0.5:.3f}")
    print(f"Histogram written to {args.out}")


if __name__ == "__main__":
    main()
