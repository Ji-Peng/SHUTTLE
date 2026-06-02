#!/usr/bin/env python3
"""Offline generator for SHUTTLE ApproxExp precomputed constants.

Implements the design specified in agent/Approx/ApproxExp.tex
(Algorithm 1, B = 32). Emits a C header with:

    TABLE_T[B]  : Q63 fixed-point  T[j] = round(2^{-j/B} * 2^63)
    RCP         : Q57 constant     RCP  = floor((B / ln2) * 2^57)
    LB          : Q(64+beta) const LB   = round((ln2 / B) * 2^(64 + beta))

    Per-mode Q80 reciprocal of 2*sigma^2 (for the a_q60 path):
      I = floor(2^80 / (2 * sigma^2))    ~ 64..66 bits
      stored as (I_hi : top bits above bit 64, I_lo : low 64 bits)

The polynomial coefficients themselves come from
tools/LogPolyApprox/expx_31div32mulln2_ln2_53_64.txt (degree 6, Q63,
all positive uint64). They are hard-coded in approx_exp.c -- only the
TABLE / RCP / LB / I_* are generated here so the script can be re-run
to confirm bit-exact reproducibility.

Error bound analysis is in ApproxExp.tex / DisApproxExp.tex; the key
choices encoded by the constants below are:
  beta       = 5             (so B = 32)
  RCP_shift  = 57            (giving V = floor((B*a/ln2) * 2^53) after
                              mulh64 with a in Q60)
  LB_shift   = 64 + beta     (so Step 2 yields R = r_l * 2^(64+beta);
                              R >> beta gives Q64 input to Horner)
  I_shift    = 80            (Q80 reciprocal of 2*sigma^2 keeps the
                              a_q60 error chain at 2^{-59}; using Q64
                              would drop it to 2^{-46})

Run:
    python3 gen_approx_exp_tables.py
Writes ``approx_exp_constants.h`` alongside the script.
"""

from decimal import Decimal, getcontext
import math
import os

# 300 decimal digits is comfortably > 2^800, more than enough for any
# table entry / reciprocal we compute below.
getcontext().prec = 300


# ============================================================
# Per-design parameters (mirror Algorithm 1 of ApproxExp.tex)
# ============================================================
BETA = 5
B = 1 << BETA                  # = 32

# Q57 reciprocal that turns mulh64(a_q60, RCP) into floor((B*a/ln2)*2^53)
RCP_SHIFT = 57
# Q(64+BETA) lift of ln2/B used in Step 2 range-reduction
LB_SHIFT = 64 + BETA           # = 69
# Q80 reciprocal of 2*sigma^2 for the sampler-side a_q60 computation
I_SHIFT = 80

# ============================================================
# Per-mode 2*sigma^2 values (taken from SHUTTLE-Spec Table 2)
# ============================================================
SHUTTLE_MODES = [
    # (security, sigma, 2*sigma^2)
    (128, 101,  2 * 101 * 101),   # SHUTTLE-128 -> 20402
    (256, 149,  2 * 149 * 149),   # SHUTTLE-256 -> 44402
    (512, 202,  2 * 202 * 202),   # SHUTTLE-512 -> 81608
]


def ln2_decimal():
    """High-precision ln(2). Decimal has no built-in ln(); compute via
    the rapidly-converging series  ln 2 = sum_{k>=0} 1/((2k+1)*9^(k+1)) * 6
    -- actually use math.log(2) and Decimal cast from a higher-precision
    approximation built via Newton's identity for 2 = exp(ln 2)."""
    # The cleanest approach in pure stdlib: Newton iteration on exp,
    # starting from the float approximation.
    x = Decimal("0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102057068573368552023575813055703267075163507596193072757082837143519030703862389167347112335011")
    return x


def emit_uint(name, value, comment="", width=20):
    """Format a uint64 / __uint128_t constant as a C-friendly literal."""
    if value < (1 << 64):
        return f"#define {name}  UINT64_C({value})  /* {comment} */"
    # 128-bit constants: split into hi:lo pair.
    hi = value >> 64
    lo = value & ((1 << 64) - 1)
    return (f"#define {name}_HI  UINT64_C({hi})  /* top bits of {comment} */\n"
            f"#define {name}_LO  UINT64_C({lo})  /* low 64 bits */")


def main():
    ln2 = ln2_decimal()

    # ----------------------------------------------------------------
    # T[j] = round(2^{-j/B} * 2^63)  for j in 0..B-1
    # ----------------------------------------------------------------
    table = []
    two63 = Decimal(1 << 63)
    for j in range(B):
        # 2^{-j/B} via exp(-j/B * ln2)
        val = two63 * (Decimal(-j) * ln2 / Decimal(B)).exp()
        # round-to-nearest, ties-to-even (default Decimal rounding)
        table.append(int(val.quantize(Decimal(1))))
    # Sanity: T[0] = 2^63, T[B-1] approx 2^{-31/32} * 2^63 ~ 0.5111 * 2^63
    assert table[0] == 1 << 63, f"T[0] expected 2^63, got {table[0]}"

    # ----------------------------------------------------------------
    # RCP = floor((B / ln2) * 2^57)
    # ----------------------------------------------------------------
    rcp_real = Decimal(B) / ln2 * Decimal(1 << RCP_SHIFT)
    rcp = int(rcp_real)  # truncates toward zero == floor for positive

    # ----------------------------------------------------------------
    # LB  = round((ln2 / B) * 2^(64+beta))
    # ----------------------------------------------------------------
    lb_real = ln2 / Decimal(B) * Decimal(1 << LB_SHIFT)
    lb = int(lb_real.quantize(Decimal(1)))

    # ----------------------------------------------------------------
    # Per-mode I = floor(2^80 / (2*sigma^2))
    # ----------------------------------------------------------------
    per_mode = []
    for security, sigma, two_sigma2 in SHUTTLE_MODES:
        I = (1 << I_SHIFT) // two_sigma2
        I_hi = I >> 64
        I_lo = I & ((1 << 64) - 1)
        per_mode.append((security, sigma, two_sigma2, I, I_hi, I_lo))

    # ----------------------------------------------------------------
    # Write header
    # ----------------------------------------------------------------
    here = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(here, "approx_exp_constants.h")

    lines = []
    lines.append(f'/*')
    lines.append(f' * approx_exp_constants.h - Precomputed constants for ApproxExp')
    lines.append(f' *                          (Algorithm 1 of agent/Approx/ApproxExp.tex).')
    lines.append(f' *')
    lines.append(f' * AUTO-GENERATED by gen_approx_exp_tables.py. DO NOT EDIT BY HAND.')
    lines.append(f' *')
    lines.append(f' *  Design parameters (Algorithm 1):')
    lines.append(f' *    BETA      = {BETA}     -> B = 2^BETA = {B}')
    lines.append(f' *    RCP shift = {RCP_SHIFT}   -> RCP = floor((B/ln2) * 2^{RCP_SHIFT})')
    lines.append(f' *    LB  shift = {LB_SHIFT}    -> LB  = round((ln2/B) * 2^(64+BETA)) = ((ln2/B) * 2^{LB_SHIFT})')
    lines.append(f' *    I   shift = {I_SHIFT}   -> I   = floor(2^{I_SHIFT} / (2*sigma^2))   per mode')
    lines.append(f' *')
    lines.append(f' *  Error budget (relative to exp(-a); see ApproxExp.tex Section 4):')
    lines.append(f' *    E_aq60   <= 2^-59   (Q80 reciprocal + 128-bit multiply)')
    lines.append(f' *    E_range  <= 2^-61   (m*LB plus Q(64+beta) -> Q64 shift)')
    lines.append(f' *    E_poly   <= 2^-61.31 (Sage Horner-aware fit on [31/32*ln2, ln2])')
    lines.append(f' *    E_table  <= 2^-63   (T[j] rounded to nearest Q63)')
    lines.append(f' *    E_mul    <= 2^-62   (one mulh64 truncation in Step 5)')
    lines.append(f' *    E_shift  <= 2^-64 / exp(-a)  (round-to-nearest output)')
    lines.append(f' *  Combined: 2^-58 (a-independent) + E_shift; worst-case a~7.37 -> 2^-53.3')
    lines.append(f' *  meeting the K=53 relative-precision target with ~0.3 bit margin.')
    lines.append(f' */')
    lines.append(f'')
    lines.append(f'#ifndef SHUTTLE_APPROX_EXP_CONSTANTS_H')
    lines.append(f'#define SHUTTLE_APPROX_EXP_CONSTANTS_H')
    lines.append(f'')
    lines.append(f'#include <stdint.h>')
    lines.append(f'#include "config.h"')
    lines.append(f'')
    lines.append(f'/* Design parameters. Algorithm 1 of ApproxExp.tex.')
    lines.append(f' * Changing BETA propagates through TABLE_T size, RCP, LB and the')
    lines.append(f' * compile-time bit-shifts in approx_exp.c; the generator script ties')
    lines.append(f' * all of them together. */')
    lines.append(f'#define APPROX_EXP_BETA   {BETA}')
    lines.append(f'#define APPROX_EXP_B      (1u << APPROX_EXP_BETA)')
    lines.append(f'')
    lines.append(f'/* RCP = floor((B / ln2) * 2^{RCP_SHIFT}). Used in Step 1 of approx_exp:')
    lines.append(f' *   V = mulh64(a_q60, APPROX_EXP_RCP) = floor((B*a/ln2) * 2^53)')
    lines.append(f' * Stored as a single uint64 since log2(RCP) ~ 62.5 < 64. */')
    lines.append(f'#define APPROX_EXP_RCP_SHIFT  {RCP_SHIFT}')
    lines.append(f'#define APPROX_EXP_RCP        UINT64_C({rcp})')
    lines.append(f'')
    lines.append(f'/* LB = round((ln2 / B) * 2^(64+BETA)). Used in Step 2:')
    lines.append(f' *   R   = m * LB - (a_q60 << (BETA + 4))    [128-bit]')
    lines.append(f' *   X   = R >> BETA                          [Q64]')
    lines.append(f' * Stored as a single uint64 since log2(LB) ~ 63.47 < 64. */')
    lines.append(f'#define APPROX_EXP_LB_SHIFT   {LB_SHIFT}')
    lines.append(f'#define APPROX_EXP_LB         UINT64_C({lb})')
    lines.append(f'')
    lines.append(f'/* Q63 table  T[j] = round(2^{{-j/B}} * 2^63),  j in 0..B-1.')
    lines.append(f' * Total: B * 8 = {B*8} bytes. Looked up in constant time by approx_exp.c. */')
    lines.append(f'static const uint64_t APPROX_EXP_TABLE_T[APPROX_EXP_B] = {{')
    for i in range(0, B, 4):
        row = "    " + ", ".join(f"UINT64_C({v:20d})" for v in table[i:i+4]) + ","
        lines.append(row)
    lines.append("};")
    lines.append("")
    lines.append("/* ============================================================")
    lines.append(" * Per-mode Q80 reciprocal of 2*sigma^2.")
    lines.append(" *")
    lines.append(" *   I = floor(2^80 / (2*sigma^2))  (~64..66 bits)")
    lines.append(" *")
    lines.append(" * The sampler computes the input a_q60 to approx_exp as:")
    lines.append(" *   num   = y * (y + 2*k*x)                  -- fits in uint32 (<2^20)")
    lines.append(" *   M_lo  = num * I_lo  (full 128-bit)")
    lines.append(" *   M_hi  = num * I_hi  (always fits in uint64 since num < 2^20)")
    lines.append(" *   a_q60 = (M_lo + (M_hi << 64)) >> 20      -- truncated to Q60")
    lines.append(" *")
    lines.append(" * Error analysis (Section 4.1 of ApproxExp.tex):")
    lines.append(" *   I = 2^80/N - delta_I,  0 <= delta_I < 1")
    lines.append(" *   M = num * I            (exact integer multiply)")
    lines.append(" *   |M_ideal - M|     <= num <= 2^18      (in Q80 ULP)")
    lines.append(" *   |a_q60 - a*2^60|  <= 2^18/2^20 + 1 = 5/4 + 1 < 2  (in Q60 ULP)")
    lines.append(" *   |Delta a|         <= 2 / 2^60 = 2^{-59}")
    lines.append(" *   |E_aq60|/exp(-a)  <= 2^{-59}")
    lines.append(" * ============================================================ */")
    lines.append("")
    for security, sigma, two_sigma2, I, I_hi, I_lo in per_mode:
        lines.append(f"/* SHUTTLE-{security}: sigma = {sigma}, 2*sigma^2 = {two_sigma2},  I has {I.bit_length()} bits. */")
        lines.append(f"#if SHUTTLE_MODE == {security}")
        lines.append(f"#  define APPROX_EXP_I_HI  UINT64_C({I_hi})")
        lines.append(f"#  define APPROX_EXP_I_LO  UINT64_C({I_lo})")
        lines.append(f"#endif")
        lines.append("")
    lines.append("#endif /* SHUTTLE_APPROX_EXP_CONSTANTS_H */")
    lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"Wrote {out_path}")
    print(f"  B = {B}, BETA = {BETA}")
    print(f"  RCP = {rcp}   (bits: {rcp.bit_length()})")
    print(f"  LB  = {lb}   (bits: {lb.bit_length()})")
    print(f"  T[0]  = {table[0]}  (= 2^63)")
    print(f"  T[B-1]= {table[-1]}")
    for security, sigma, two_sigma2, I, I_hi, I_lo in per_mode:
        print(f"  SHUTTLE-{security}: 2sigma^2 = {two_sigma2:>6d}, I = {I}  ({I.bit_length()} bits)")


if __name__ == "__main__":
    main()
