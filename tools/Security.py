"""Security: SHUTTLE-NGCC-SUF 各安全级别的散度幂次估计。

由 Security.ipynb 转换而来。所有输出重定向到 log/Security.txt。
运行: python Security.py
"""

import math
import os
import sys
from contextlib import redirect_stdout


def calculate_powers(lam, q_exp, q_bs):
    # --- 公式 1 ---
    # 分母 = sqrt(2) * (2 * lambda + 1) * Q_exp
    denominator1 = math.sqrt(2 * (2 * lam - 1) * q_exp)
    # 得到 2 的幂次
    power1 = math.log2(denominator1)
    # --- 公式 2 ---
    # 分母 = 4 * Q_bs
    denominator2 = 4 * q_bs
    # 得到 2 的幂次
    power2 = math.log2(denominator2)
    return -power1, -power2


def main():
    ## For SHUTTLE-NGCC-SUF 128
    # SIS security level
    LAMBDA = 166
    N = 256
    L = 3
    M = 3
    # 指数调用次数
    Q_EXP = (2**80 * N * (L + M + 1)) / 0.8
    # 基础采样调用次数
    Q_BS = Q_EXP
    res1, res2 = calculate_powers(LAMBDA, Q_EXP, Q_BS)
    print("\nSHUTTLE-NGCC-SUF 128:")
    print(f"第一个公式近似结果: 2^{{{res1:.2f}}}")
    print(f"第二个公式近似结果: 2^{{{res2:.2f}}}")

    ## For SHUTTLE-NGCC-SUF 256
    # SIS security level
    LAMBDA = 267
    N = 512
    L = 3
    M = 2
    # 指数调用次数
    Q_EXP = (2**80 * N * (L + M + 1)) / 0.8
    # 基础采样调用次数
    Q_BS = Q_EXP
    res1, res2 = calculate_powers(LAMBDA, Q_EXP, Q_BS)
    print("\nSHUTTLE-NGCC-SUF 256:")
    print(f"第一个公式近似结果: 2^{{{res1:.2f}}}")
    print(f"第二个公式近似结果: 2^{{{res2:.2f}}}")

    ## For SHUTTLE-NGCC-SUF 512
    # SIS security level
    LAMBDA = 523
    N = 1024
    L = 3
    M = 2
    # 指数调用次数
    Q_EXP = (2**80 * N * (L + M + 1)) / 0.8
    # 基础采样调用次数
    Q_BS = Q_EXP
    res1, res2 = calculate_powers(LAMBDA, Q_EXP, Q_BS)
    print("\nSHUTTLE-NGCC-SUF 512:")
    print(f"第一个公式近似结果: 2^{{{res1:.2f}}}")
    print(f"第二个公式近似结果: 2^{{{res2:.2f}}}")


if __name__ == "__main__":
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "log")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "Security.txt")
    with open(log_path, "w", encoding="utf-8") as f:
        with redirect_stdout(f):
            main()
    print(f"输出已写入 {log_path}", file=sys.stderr)
