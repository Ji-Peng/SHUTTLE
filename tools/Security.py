"""Security: Divergence exponent estimates for each SHUTTLE-NGCC-SUF security level.

Converted from Security.ipynb. All output is redirected to log/Security.txt.
Run: python Security.py
"""

import math
import os
import sys
from contextlib import redirect_stdout


def calculate_powers(lam, q_exp, q_bs):
    # --- Formula 1 ---
    # denominator = sqrt(2) * (2 * lambda + 1) * Q_exp
    denominator1 = math.sqrt(2 * (2 * lam - 1) * q_exp)
    # Get the power of 2.
    power1 = math.log2(denominator1)
    # --- Formula 2 ---
    # denominator = 4 * Q_bs
    denominator2 = 4 * q_bs
    # Get the power of 2.
    power2 = math.log2(denominator2)
    return -power1, -power2


def main():
    ## For SHUTTLE-NGCC-SUF 128
    # SIS security level
    LAMBDA = 166
    N = 256
    L = 3
    M = 3
    # Number of exponential calls
    Q_EXP = (2**80 * N * (L + M + 1)) / 0.8
    # Number of base-sampling calls
    Q_BS = Q_EXP
    res1, res2 = calculate_powers(LAMBDA, Q_EXP, Q_BS)
    print("\nSHUTTLE-NGCC-SUF 128:")
    print(f"First formula approximation: 2^{{{res1:.2f}}}")
    print(f"Second formula approximation: 2^{{{res2:.2f}}}")

    ## For SHUTTLE-NGCC-SUF 256
    # SIS security level
    LAMBDA = 267
    N = 512
    L = 3
    M = 2
    # Number of exponential calls
    Q_EXP = (2**80 * N * (L + M + 1)) / 0.8
    # Number of base-sampling calls
    Q_BS = Q_EXP
    res1, res2 = calculate_powers(LAMBDA, Q_EXP, Q_BS)
    print("\nSHUTTLE-NGCC-SUF 256:")
    print(f"First formula approximation: 2^{{{res1:.2f}}}")
    print(f"Second formula approximation: 2^{{{res2:.2f}}}")

    ## For SHUTTLE-NGCC-SUF 512
    # SIS security level
    LAMBDA = 523
    N = 1024
    L = 3
    M = 2
    # Number of exponential calls
    Q_EXP = (2**80 * N * (L + M + 1)) / 0.8
    # Number of base-sampling calls
    Q_BS = Q_EXP
    res1, res2 = calculate_powers(LAMBDA, Q_EXP, Q_BS)
    print("\nSHUTTLE-NGCC-SUF 512:")
    print(f"First formula approximation: 2^{{{res1:.2f}}}")
    print(f"Second formula approximation: 2^{{{res2:.2f}}}")


if __name__ == "__main__":
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "log")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "Security.txt")
    with open(log_path, "w", encoding="utf-8") as f:
        with redirect_stdout(f):
            main()
    print(f"Output written to {log_path}", file=sys.stderr)
