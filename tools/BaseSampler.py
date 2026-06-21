"""BaseSampler: RCDT 基础采样器表生成与 Rényi 散度分析。

由 BaseSampler.ipynb 转换而来。所有输出重定向到 log/BaseSampler.txt。
运行: python BaseSampler.py
"""

import math
import os
import sys
from contextlib import redirect_stdout
from decimal import Decimal, getcontext

# 设置全局极高精度，避免 a=509 次方时发生下溢或精度丢失
getcontext().prec = 300


def compute_ideal_gaussian(sigma, max_val):
    """
    计算理想的离散半高斯分布 D_{Z+, sigma} 的真实概率。
    max_val 应足够大（如 100 * sigma），以逼近无穷级数的和。
    """
    sigma_d = Decimal(str(sigma))
    weights = [Decimal(0)] * max_val
    for x in range(max_val):
        x_d = Decimal(x)
        # exp(-x^2 / (2 * sigma^2))
        weights[x] = Decimal.exp(-(x_d**2) / (Decimal(2) * sigma_d**2))

    total_S = sum(weights)
    return [w / total_S for w in weights]


def compute_base_sampler_tables(sigma, w, theta):
    """
    根据论文算法生成 PDT, CDT, RCDT。
    :param sigma: 标准差
    :param w: 截断边界 (Tailcut)
    :param theta: 绝对精度位宽
    :return: pdt_int, cdt_int, rcdt_int (均为乘以 2^theta 后的整数数组)
    """
    sigma_d = Decimal(str(sigma))
    two_theta = Decimal(2) ** theta

    # 步骤 1: 计算截断在 {0, ..., w} 上的受限理想分布 D_[w],sigma
    weights_w = [Decimal.exp(-(Decimal(x)**2) / (Decimal(2) * sigma_d**2)) for x in range(w + 1)]
    S_w = sum(weights_w)
    D_w = [weight / S_w for weight in weights_w]

    # 步骤 2: 计算定点精度概率密度表 PDT
    pdt_int = [0] * (w + 1)
    sum_pdt_tail = 0

    # 对于 z >= 1，向下取整截断
    for z in range(1, w + 1):
        # PDT(z) = 2^-theta * floor(2^theta * D_w(z))
        val_int = int(D_w[z] * two_theta)
        pdt_int[z] = val_int
        sum_pdt_tail += val_int

    # 对于 z = 0，用 1 减去剩余部分以保证总和严格为 1
    pdt_int[0] = int(two_theta) - sum_pdt_tail

    # 步骤 3: 从 PDT 派生 CDT 和 RCDT
    cdt_int = [0] * (w + 1)
    rcdt_int = [0] * (w + 1)

    current_sum = 0
    for i in range(w + 1):
        current_sum += pdt_int[i]
        cdt_int[i] = current_sum

    total_space = cdt_int[-1]  # 必然等于 2^theta
    for i in range(w + 1):
        rcdt_int[i] = total_space - cdt_int[i]

    return pdt_int, cdt_int, rcdt_int


def calculate_renyi_divergence(pdt_int, theta, sigma, a):
    """
    计算截断舍入后的分布与理想分布之间的 a 阶瑞利散度 R_a
    R_a(P || Q) = ( sum( P(x)^a / Q(x)^(a-1) ) ) ^ (1/(a-1))
    """
    a_d = Decimal(a)
    two_theta = Decimal(2) ** theta

    # 获取理想的无限域分布（计算到 100 * sigma 保证精度）
    Q_ideal = compute_ideal_gaussian(sigma, max_val=100 * math.ceil(sigma))

    divergence_sum = Decimal(0)

    # 只需遍历 PDT 的支撑集 {0, ..., w}，因为超出该范围 P(x) = 0
    for x in range(len(pdt_int)):
        P_x = Decimal(pdt_int[x]) / two_theta
        if P_x == 0:
            continue

        Q_x = Q_ideal[x]

        # P(x)^a / Q(x)^(a-1)
        term = (P_x ** a_d) / (Q_x ** (a_d - Decimal(1)))
        divergence_sum += term

    # (sum)^(1/(a-1))
    r_a = divergence_sum ** (Decimal(1) / (a_d - Decimal(1)))
    return r_a


def calculate_renyi_divergence_inf(pdt_int, theta, sigma):
    """
    计算无穷阶 Rényi 散度 R_inf(P || Q) = max_x P(x) / Q(x)

    无穷阶 Rényi 散度是有限阶的极限：当 a -> inf 时，
    R_a 退化为 P 和 Q 逐点比值的最大值。

    直觉：有限阶 R_a 对 P(x)/Q(x) 做"加权幂平均"，阶数越高越偏向最大比值；
    无穷阶就是直接取最大值。
    """
    two_theta = Decimal(2) ** theta

    Q_ideal = compute_ideal_gaussian(sigma, max_val=100 * max(math.ceil(sigma), 2))

    max_ratio = Decimal(0)
    for x in range(len(pdt_int)):
        P_x = Decimal(pdt_int[x]) / two_theta
        if P_x == 0:
            continue
        Q_x = Q_ideal[x]
        ratio = P_x / Q_x
        if ratio > max_ratio:
            max_ratio = ratio

    return max_ratio


def format_divergence(r_a):
    """将散度格式化为易读的 1 + 2^(-x) 形式"""
    diff = r_a - Decimal(1)
    if diff <= 0:
        return "1.0 (完美拟合)"

    # 计算 x = -log2(R_a - 1)
    bits = -float(diff.ln() / Decimal(2).ln())
    return f"1 + 2^(-{bits:.2f})"


def format_signed_power_of_two(x):
    """将带符号差值格式化为 +/- 2^(-bits)。"""
    x = Decimal(x)
    if x == 0:
        return "0"
    sign = "-" if x < 0 else "+"
    bits = -float(abs(x).ln() / Decimal(2).ln())
    return f"{sign}2^(-{bits:.2f})"


def calculate_rcdt_sampler_stddev(rcdt_int, theta, ideal_sigma=1, zero_fold=True):
    """
    计算实际 RCDT 表访问采样器诱导出的标准差。

    RCDT 采样器逻辑为 z = sum_i [u < RCDT[i]], u 在 {0,...,2^theta-1} 上均匀。
    对 keygen 的 sigma=1 噪声采样器，还会随机赋号，并对 z=0 执行 1/2 的 zero-fold 拒绝；
    zero_fold=True 时返回的就是最终有符号输出分布的标准差。
    """
    if theta <= 0:
        raise ValueError("theta 必须为正整数")

    total = Decimal(2) ** theta
    thresholds = [int(x) for x in rcdt_int]
    while thresholds and thresholds[-1] == 0:
        thresholds.pop()

    if not thresholds:
        return {
            "support_max": 0,
            "magnitude_pmf": [Decimal(1)],
            "accept_mass": Decimal(1) if zero_fold else None,
            "variance": Decimal(0),
            "stddev": Decimal(0),
            "ideal_sigma": Decimal(str(ideal_sigma)),
            "abs_gap": -Decimal(str(ideal_sigma)),
            "rel_gap": -Decimal(1),
        }

    for i, t in enumerate(thresholds):
        if t < 0 or t > int(total):
            raise ValueError(f"RCDT[{i}] 超出 [0, 2^theta] 范围: {t}")
        if i > 0 and thresholds[i - 1] < t:
            raise ValueError("RCDT 表必须单调不增")

    # Pr[z >= k] = RCDT[k-1] / 2^theta，因此可由相邻 tail 差分恢复 PMF。
    tail = [Decimal(t) / total for t in thresholds] + [Decimal(0)]
    magnitude_pmf = [Decimal(1) - tail[0]]
    for k in range(1, len(tail)):
        magnitude_pmf.append(tail[k - 1] - tail[k])

    second_moment = sum(Decimal(k * k) * p for k, p in enumerate(magnitude_pmf))
    if zero_fold:
        accept_mass = Decimal(1) - magnitude_pmf[0] / Decimal(2)
        if accept_mass <= 0:
            raise ValueError("zero-fold 接受概率必须为正")
        variance = second_moment / accept_mass
    else:
        accept_mass = None
        variance = second_moment

    stddev = variance.sqrt()
    ideal_sigma_d = Decimal(str(ideal_sigma))
    abs_gap = stddev - ideal_sigma_d
    rel_gap = abs_gap / ideal_sigma_d if ideal_sigma_d != 0 else Decimal(0)

    return {
        "support_max": len(magnitude_pmf) - 1,
        "magnitude_pmf": magnitude_pmf,
        "accept_mass": accept_mass,
        "variance": variance,
        "stddev": stddev,
        "ideal_sigma": ideal_sigma_d,
        "abs_gap": abs_gap,
        "rel_gap": rel_gap,
    }


def print_rcdt_table(rcdt_values, total_bits, base_bits, table_name, drop_terminal_zero=True):
    """
    将 rcdt 表按 base_bits 分割成若干 limbs，并打印为 C 数组。
    规则：
    - limbs = ceil(total_bits / base_bits)
    - limb0 是最低位分组（little-endian limb 顺序）
    - base_bits <= 32 用 uint32_t / U
    - base_bits > 32 用 uint64_t / ULL
    - 默认去掉末尾所有值为 0 的 RCDT 表项，使输出更贴近常见实现表格式
    """
    if base_bits <= 0:
        raise ValueError("base_bits 必须为正整数")

    values = list(rcdt_values)
    if drop_terminal_zero:
        while len(values) > 0 and values[-1] == 0:
            values = values[:-1]

    limb_count = (total_bits + base_bits - 1) // base_bits
    if limb_count <= 0:
        raise ValueError("limb_count 计算错误")

    c_type = "uint32_t" if base_bits <= 32 else "uint64_t"
    suffix = "U" if base_bits <= 32 else "ULL"

    # 固定十六进制宽度：每 limb 显示满宽，方便对齐阅读
    hex_width = (base_bits + 3) // 4
    limb_mask = (1 << base_bits) - 1

    print(f"\n/* {total_bits}-bit RCDT -> {limb_count}x{base_bits}-bit limbs */")
    print(f"const {c_type} {table_name}[{len(values)}][{limb_count}] = {{")

    limit = 1 << total_bits
    for i, v in enumerate(values):
        if v < 0:
            raise ValueError(f"RCDT[{i}] 为负数: {v}")
        if v >= limit:
            raise ValueError(f"RCDT[{i}] 超出 {total_bits} bit: {v}")

        limbs = []
        for j in range(limb_count):
            limb = (v >> (j * base_bits)) & limb_mask
            limbs.append(f"0x{limb:0{hex_width}X}{suffix}")

        comma = "," if i != len(values) - 1 else ""
        print(f"    {{{', '.join(limbs)}}}{comma}")

    print("};")


def run_pk_gaussian_sampler():
    """SHUTTLE-NGCC-SUF 128,256,512 的 pk 高斯采样器表与散度分析。"""
    SIGMA = 825 / 256
    THETA = 96
    # For SHUTTLE-NGCC-SUF 512
    SIS_Security = 523
    A_ORDER = SIS_Security * 2 - 1
    target_w = math.ceil(11 * SIGMA)

    print(f"=== SHUTTLE-NGCC-SUF 128/256/512: sigma={SIGMA}, theta={THETA} bits, a={A_ORDER} ===")

    pdt, cdt, rcdt = compute_base_sampler_tables(SIGMA, target_w, THETA)
    r_a = calculate_renyi_divergence(pdt, THETA, SIGMA, A_ORDER)
    r_inf = calculate_renyi_divergence_inf(pdt, THETA, SIGMA)
    print_rcdt_table(
        rcdt_values=rcdt,
        total_bits=96,
        base_bits=32,
        table_name="GAUSS0_96_3x32"
    )
    print(f"\nRenyi divergence R_{A_ORDER}: {format_divergence(r_a)}")
    print(f"Renyi divergence R_inf:  {format_divergence(r_inf)}")
    for w_test in range(math.ceil(9 * SIGMA), math.ceil(12 * SIGMA)):
        test_pdt, _, _ = compute_base_sampler_tables(SIGMA, w_test, THETA)
        test_ra = calculate_renyi_divergence(test_pdt, THETA, SIGMA, A_ORDER)
        test_rinf = calculate_renyi_divergence_inf(test_pdt, THETA, SIGMA)
        print(f"w = {w_test:2d} | R_{A_ORDER}: {format_divergence(test_ra)} | R_inf: {format_divergence(test_rinf)}")


def run_sk_gaussian_sampler():
    """SHUTTLE-NGCC-SUF sk 高斯采样器表、散度分析与实际标准差。"""
    SIGMA = [0.85, 0.9, 1]
    THETA = 93
    # We did not care about the Renyi divergence for the sk Gaussian sampler.
    A_ORDER = 512 * 2 - 1

    for each_sigma in SIGMA:
        target_w = math.ceil(11 * each_sigma)
        print(f"=== SHUTTLE-NGCC-SUF sk Gaussian Sampler: sigma={each_sigma}, theta={THETA} bits, a={A_ORDER} ===")

        pdt, cdt, rcdt = compute_base_sampler_tables(each_sigma, target_w, THETA)
        r_a = calculate_renyi_divergence(pdt, THETA, each_sigma, A_ORDER)
        r_inf = calculate_renyi_divergence_inf(pdt, THETA, each_sigma)
        print_rcdt_table(
            rcdt_values=rcdt,
            total_bits=93,
            base_bits=31,
            table_name="GAUSS0_93_3x31"
        )
        print(f"\nRenyi divergence R_{A_ORDER}: {format_divergence(r_a)}")
        print(f"Renyi divergence R_inf:  {format_divergence(r_inf)}")
        for w_test in range(math.ceil(9 * each_sigma), math.ceil(12 * each_sigma)):
            test_pdt, _, _ = compute_base_sampler_tables(each_sigma, w_test, THETA)
            test_ra = calculate_renyi_divergence(test_pdt, THETA, each_sigma, A_ORDER)
            test_rinf = calculate_renyi_divergence_inf(test_pdt, THETA, each_sigma)
            print(f"w = {w_test:2d} | R_{A_ORDER}: {format_divergence(test_ra)} | R_inf: {format_divergence(test_rinf)}")

        std_stats = calculate_rcdt_sampler_stddev(
            rcdt_int=rcdt,
            theta=THETA,
            ideal_sigma=each_sigma,
            zero_fold=True,
        )
        print("\n--- RCDT 实际标准差（keygen sigma=1 噪声采样器）---")
        print(f"RCDT 支撑范围: [-{std_stats['support_max']}, {std_stats['support_max']}]")
        print(f"zero-fold 接受概率: {std_stats['accept_mass']:.30f}")
        print(f"实际方差: {std_stats['variance']:.30f}")
        print(f"实际标准差: {std_stats['stddev']:.30f}")
        print(f"相对理想标准差 {each_sigma} 的绝对差: {std_stats['abs_gap']:.6E} ({format_signed_power_of_two(std_stats['abs_gap'])})")
        print(f"相对理想标准差 {each_sigma} 的相对差: {std_stats['rel_gap']:.6E} ({format_signed_power_of_two(std_stats['rel_gap'])})")


def main():
    run_pk_gaussian_sampler()
    print()
    run_sk_gaussian_sampler()


if __name__ == "__main__":
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "log")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "BaseSampler.txt")
    with open(log_path, "w", encoding="utf-8") as f:
        with redirect_stdout(f):
            main()
    print(f"输出已写入 {log_path}", file=sys.stderr)
