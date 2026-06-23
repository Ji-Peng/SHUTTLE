#!/usr/bin/env python3
"""Generate the accept-branch exp table for the SHUTTLE BLISS sampler.

The improved proof only needs acceptance-probability relative error.  A direct
Q64 table is faster and easier to audit than an arithmetic approximation: the
runtime path is an indexed load, and the verifier enumerates every possible
(x,y) pair.
"""
import argparse
import math
from pathlib import Path

import mpmath as mp

mp.mp.dps = 100

R = 825
K = 256
X_MAX = 36
Y_MAX = 255
TARGET_BITS = 53.0
UINT64_MAX = (1 << 64) - 1
SCALE = 1 << 64


def exponent(x: int, y: int) -> mp.mpf:
    return -mp.mpf(y * (y + 2 * K * x)) / mp.mpf(2 * R * R)


def q64_threshold(x: int, y: int) -> int:
    if y == 0:
        return UINT64_MAX
    p = exponent(x, y)
    return int(mp.floor(mp.exp(p) * mp.mpf(SCALE) + mp.mpf("0.5")))


def write_header(path: Path) -> None:
    rows = [[q64_threshold(x, y) for y in range(Y_MAX + 1)] for x in range(X_MAX + 1)]
    lines = []
    lines.append("#ifndef SHUTTLE_APPROX_EXP_ACCEPT_TABLE_H")
    lines.append("#define SHUTTLE_APPROX_EXP_ACCEPT_TABLE_H")
    lines.append("")
    lines.append("#include <stdint.h>")
    lines.append("")
    lines.append("#define SHUTTLE_EXP_ACCEPT_X_MAX 36")
    lines.append("#define SHUTTLE_EXP_ACCEPT_Y_MAX 255")
    lines.append("#define SHUTTLE_EXP_ACCEPT_QBITS 64")
    lines.append("")
    lines.append("static const uint64_t kShuttleExpAcceptQ64[37][256] = {")
    for x, row in enumerate(rows):
        lines.append(f"    /* x = {x} */")
        lines.append("    {")
        for i in range(0, len(row), 4):
            chunk = row[i:i + 4]
            values = ", ".join(f"UINT64_C({v})" for v in chunk)
            comma = "," if i + 4 < len(row) else ""
            lines.append(f"        {values}{comma}")
        comma = "," if x < X_MAX else ""
        lines.append(f"    }}{comma}")
    lines.append("};")
    lines.append("")
    lines.append("static inline uint64_t shuttle_exp_accept_q64(int x, int y)")
    lines.append("{")
    lines.append("    return kShuttleExpAcceptQ64[x][y];")
    lines.append("}")
    lines.append("")
    lines.append("#endif")
    path.write_text("\n".join(lines) + "\n")


def write_log(path: Path) -> None:
    min_p = mp.mpf(0)
    min_arg = (0, 0)
    max_rel = mp.mpf(0)
    worst = (0, 0)
    for x in range(X_MAX + 1):
        for y in range(Y_MAX + 1):
            p = exponent(x, y)
            if p < min_p:
                min_p = p
                min_arg = (x, y)
            if y == 0:
                rel = 0.0
            else:
                got = mp.mpf(q64_threshold(x, y)) / mp.mpf(SCALE)
                ref = mp.exp(p)
                rel = abs(got - ref) / ref
            if rel > max_rel:
                max_rel = rel
                worst = (x, y)
    lines = []
    lines.append("SHUTTLE accept-branch exp approximation")
    lines.append(f"r = {R}, k = {K}")
    lines.append("x in {0,...,36}, y in {0,...,255}")
    lines.append("p = -y*(y+2*k*x)/(2*r^2)")
    lines.append(f"minimum p = {mp.nstr(min_p, 30)} at x={min_arg[0]} y={min_arg[1]}")
    lines.append("scheme = direct Q64 threshold table, y=0 accepted unconditionally")
    lines.append(f"table bytes = {(X_MAX + 1) * (Y_MAX + 1) * 8}")
    lines.append("runtime arithmetic = indexed load only; no 64x64 product in exp approximation")
    lines.append(f"mpmath precheck max relative error = {mp.nstr(max_rel, 30)}")
    lines.append(f"mpmath precheck precision = {mp.nstr(-mp.log(max_rel, 2), 20)} bits at x={worst[0]} y={worst[1]}")
    lines.append(f"target precision = {TARGET_BITS:.0f} bits")
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    out_dir = args.out_dir
    log_dir = out_dir / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    write_header(out_dir / "approx_exp_accept_table.h")
    write_log(log_dir / "approx_exp_accept_generation.txt")


if __name__ == "__main__":
    main()
