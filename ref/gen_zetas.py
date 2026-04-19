#!/usr/bin/env python3
"""Offline generator for SHUTTLE NTT zeta tables and Montgomery constants.

Produces drop-in constants for reduce.h and ntt.c, parameterized by (q, N).
The layout follows the Dilithium-style forward-NTT (full-N zetas, int32_t
coefficients, R = 2^32 Montgomery base).

Run:  python3 gen_zetas.py

For the shipping configuration (q = 13313) this prints MONT, QINV, and two
zeta tables (N = 256 for SHUTTLE-128 and N = 512 for SHUTTLE-256), as well
as the final invntt scale factor f = MONT^2 / N mod q for each.
"""

def mod_order(g, q):
    """Multiplicative order of g mod prime q."""
    n = q - 1
    # q - 1 = 2^10 * 13 for q=13313
    for p in [2, 13]:
        while n % p == 0 and pow(g, n // p, q) == 1:
            n //= p
    return n

def primitive_root(q, order_want):
    for g in range(2, q):
        if mod_order(g, q) == order_want:
            return g
    raise ValueError(f"No element of order {order_want} mod {q}")

def bitrev(x, bits):
    r = 0
    for i in range(bits):
        r = (r << 1) | ((x >> i) & 1)
    return r

def signed(x, q):
    x %= q
    if x > q // 2:
        x -= q
    return x

def gen(q, N):
    # Montgomery constants
    MONT = pow(2, 32, q)
    QINV = pow(q, -1, 1 << 32)
    # Primitive 2N-th root of unity
    rho = primitive_root(q, 2 * N)
    # Full-N zetas, bit-reversed
    bits = N.bit_length() - 1
    zetas = [signed((MONT * pow(rho, bitrev(i, bits), q)) % q, q) for i in range(N)]
    # Invntt final scale
    f = signed((MONT * MONT * pow(N, -1, q)) % q, q)
    # Barrett V for reduce32
    V = round((1 << 32) / q)
    return MONT, QINV, V, rho, zetas, f

def fmt_table(name, zs):
    out = [f"static const int32_t {name}[{len(zs)}] = {{"]
    for i in range(0, len(zs), 8):
        out.append("    " + ", ".join(f"{v:8d}" for v in zs[i:i+8]) + ",")
    out.append("};")
    return "\n".join(out)

def main():
    q = 13313
    print(f"q = {q}")
    for N in (256, 512):
        MONT, QINV, V, rho, zetas, f = gen(q, N)
        print(f"\n=== N = {N} ===")
        print(f"  MONT = {MONT}U")
        print(f"  QINV = {QINV}U")
        print(f"  Barrett V = {V}")
        print(f"  primitive 2N-th root = {rho}")
        print(f"  invntt scale f = {f}")
        print()
        print(fmt_table("zetas", zetas))

if __name__ == "__main__":
    main()
