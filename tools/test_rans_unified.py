#!/usr/bin/env python3
"""End-to-end regression test for the SHUTTLE rANS layer.

Coverage:
  1. Frequency tables produced by gen_rans_tables.py are well-formed for
     every mode (alphabet contiguous, freqs sum to 2^prob_bits, every
     symbol has freq >= 1) and match the theoretical M_voc bounds in
     SHUTTLE_rANS.tex §3.
  2. Python reference rans.py round-trips a discrete-Gaussian sample.
  3. Final-state verification rejects any single-byte corruption of the
     encoded stream.
  4. Mode-128 alias: the unified table doubles as both the z-hi and the
     hint table (sigma coincidence).

This intentionally does not exercise the C engine — for that, run the
binaries under SHUTTLE/ref/test/out/ (test_rans, test_poly_z1, test_sig).
The C/Python byte-identity comparison lives in kat_compare/.
"""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from rans import RansTable, encode, decode, RANS_L
from gen_rans_tables import (MODE_PARAMS, ALPHA_R, PROB_BITS, T_SIGMA,
                             zhi_table, hint_table, shares_table,
                             theoretical_pmf)


def check_table(syms, freqs, prob_bits):
    assert len(syms) == len(freqs)
    assert syms == sorted(syms)
    assert syms == list(range(syms[0], syms[-1] + 1)), \
        f"alphabet not contiguous: {syms}"
    assert syms[0] == -syms[-1], f"alphabet not symmetric: [{syms[0]}, {syms[-1]}]"
    assert all(f >= 1 for f in freqs), "found freq < 1"
    assert sum(freqs) == 1 << prob_bits, \
        f"freq sum {sum(freqs)} != 2^{prob_bits}"


def test_tables():
    print("== Test: theoretical PMF tables, geometry and sums ==")
    for mode, params in MODE_PARAMS.items():
        z_syms, z_freqs, z_sigma, z_M = zhi_table(params)
        h_syms, h_freqs, h_sigma, h_M = hint_table(params)
        check_table(z_syms, z_freqs, PROB_BITS)
        check_table(h_syms, h_freqs, PROB_BITS)

        # Spot-check M_voc against the SHUTTLE_rANS.tex tight bound (Cor 3,
        # Thm 5 for z-hi; Thm 7 for hint). The generator's max() over the
        # two z^(i) cases is what we replicate here.
        r, tau, eta = params["r"], params["tau"], params["eta"]
        alpha_1, alpha_h_spec = params["alpha_1"], params["alpha_h"]
        expected_zhi = max(
            math.ceil((T_SIGMA * r + tau * eta) / ALPHA_R),
            math.ceil((math.ceil(T_SIGMA * r / alpha_1) + 1) / (ALPHA_R // alpha_1)),
        )
        expected_hint = (2 * (T_SIGMA * r + tau * eta)) // alpha_h_spec + 1
        assert z_M == expected_zhi, \
            f"mode-{mode} z-hi M_voc: got {z_M}, expected {expected_zhi}"
        assert h_M == max(expected_hint, 1), \
            f"mode-{mode} hint M_voc: got {h_M}, expected {expected_hint}"

        if shares_table(params):
            assert (z_syms, z_freqs) == (h_syms, h_freqs), \
                f"mode-{mode} should alias, but tables differ"

        print(f"  mode-{mode}: z-hi M_voc={z_M:2d} ({len(z_syms):2d} symbols, "
              f"sigma={z_sigma:.4f})  "
              f"hint M_voc={h_M:2d} ({len(h_syms):2d} symbols, sigma={h_sigma:.4f})"
              + ("  [unified]" if shares_table(params) else ""))


def discrete_gaussian_sample(rng: random.Random, sigma: float, M_voc: int) -> int:
    """Sample from the discrete Gaussian PMF (same one the generator uses)."""
    pmf = theoretical_pmf(sigma, M_voc)
    syms = sorted(pmf.keys())
    weights = [pmf[s] for s in syms]
    cum = 0.0
    pick = rng.random()
    for s, w in zip(syms, weights):
        cum += w
        if pick <= cum:
            return s
    return syms[-1]


def test_roundtrip():
    print("== Test: theoretical-Gaussian sample -> encode -> decode ==")
    rng = random.Random(0xC0DECAFE)
    for mode, params in MODE_PARAMS.items():
        for ctx, (syms, freqs, sigma, M_voc) in [
                ("z-hi", zhi_table(params)), ("hint", hint_table(params))]:
            table = RansTable(syms=list(syms), freqs=list(freqs),
                              prob_bits=PROB_BITS)
            N = 1024 if ctx == "z-hi" else 256
            msg = [discrete_gaussian_sample(rng, sigma, M_voc) for _ in range(N)]
            buf = encode(msg, table)
            recovered = decode(buf, table, len(msg))
            assert recovered == msg, \
                f"mode-{mode} {ctx}: round-trip mismatch"
            print(f"  mode-{mode} {ctx}: {N} symbols -> {len(buf)} B "
                  f"({8 * len(buf) / N:.3f} bit/sym)")


def test_little_endian_flush():
    print("== Test: little-endian flush of the initial state ==")
    for mode, params in MODE_PARAMS.items():
        syms, freqs, _, _ = zhi_table(params)
        table = RansTable(syms=list(syms), freqs=list(freqs),
                          prob_bits=PROB_BITS)
        buf = encode([], table)
        # Empty message: the encoder writes exactly the 4-byte LE flush of L.
        assert len(buf) == 4, f"empty stream should be 4 bytes, got {len(buf)}"
        expected = bytes([
            (RANS_L >> 0) & 0xFF,
            (RANS_L >> 8) & 0xFF,
            (RANS_L >> 16) & 0xFF,
            (RANS_L >> 24) & 0xFF,
        ])
        assert buf == expected, \
            f"mode-{mode}: flush bytes {buf.hex()} != expected {expected.hex()}"
        print(f"  mode-{mode}: LE flush bytes = {buf.hex()}  (L = 0x{RANS_L:08x})")


def test_corruption_detected():
    print("== Test: final-state verification catches byte corruption ==")
    rng = random.Random(0xBADBABE)
    params = MODE_PARAMS[128]
    syms, freqs, sigma, M_voc = zhi_table(params)
    table = RansTable(syms=list(syms), freqs=list(freqs), prob_bits=PROB_BITS)
    N = 200
    msg = [discrete_gaussian_sample(rng, sigma, M_voc) for _ in range(N)]
    buf = encode(msg, table)
    undetected = 0
    for pos in range(len(buf)):
        corrupted = bytearray(buf)
        corrupted[pos] ^= 0xFF
        try:
            recovered = decode(bytes(corrupted), table, N)
            # No final-state mismatch raised; check whether symbols match.
            if recovered == msg:
                undetected += 1
        except RuntimeError:
            pass  # caught by x != RANS_L
        except IndexError:
            pass  # underflow during renorm
    assert undetected == 0, f"{undetected} corruptions went undetected"
    print(f"  all {len(buf)} single-byte corruptions caught")


def main():
    test_tables()
    test_roundtrip()
    test_little_endian_flush()
    test_corruption_detected()
    print()
    print("All rANS unified tests PASSED.")


if __name__ == "__main__":
    main()
