#!/usr/bin/env python3
"""Generate and audit the BLISS convolution exp approximation.

The audited C scheme uses a high-half integer Chebyshev/Hermite polynomial for
normal points and a small-p Taylor branch.  The Taylor branch is needed because
for y=1,x=0 the rejection denominator is about 2^-20; a 64-bit-only final
probability cannot by itself provide 2^-52 rejection-relative precision.
"""
import math
from pathlib import Path

R = 825
K = 256
TARGET = 52.0
TAYLOR_CUTOFF = 1.0 / 32.0
COEFF = [
    6025547929741064358,
    2635708275570594412,
    477776027135403486,
    73356906713080848,
    9723463780491741,
    1130674987138628,
    116898973677229,
    10866892872815,
    916890868687,
    70782307249,
    5033905542,
    331766219,
    20367973,
    1170091,
    63152,
    3214,
    155,
    7,
]

def exponent(x: int, y: int) -> float:
    return -float(y * (y + 2 * K * x)) / float(2 * R * R)


def main() -> None:
    out = Path(__file__).resolve().parent / "log" / "approx_exp_bliss_generation.txt"
    max_p = 0.0
    min_p = 0.0
    arg_min = (0, 0)
    small = 0
    for x in range(37):
        for y in range(256):
            p = exponent(x, y)
            if p < min_p:
                min_p = p
                arg_min = (x, y)
            if abs(p) < TAYLOR_CUTOFF:
                small += 1
    lines = []
    lines.append("BLISS convolution exp approximation generation record")
    lines.append(f"r = {R}, k = {K}, x in [0,36], y in [0,255]")
    lines.append("p = -y*(y+2*k*x)/(2*r^2)")
    lines.append(f"range = [{min_p:.18e}, {max_p:.18e}], minimum at {arg_min}")
    lines.append(f"target combined precision = {TARGET:.0f} bits")
    lines.append("normal branch: P(p)=1+p+p^2*R(z), z=(p+1.755)/1.755")
    lines.append("R(z)=2^-64*sum c_i*T_i(z), evaluated by integer Clenshaw")
    lines.append("64x64 products use signed high-half extraction (product >> 64); Q rescaling is a small left shift after high-half extraction")
    lines.append(f"small branch: |p| < {TAYLOR_CUTOFF:.18e}, Taylor exp(p) in __float128")
    lines.append(f"number of discrete points in small branch = {small}")
    lines.append("Chebyshev coefficients c_i:")
    for i, c in enumerate(COEFF):
        lines.append(f"  c[{i:2d}] = {c}")
    out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nlog written to {out}")


if __name__ == "__main__":
    main()
