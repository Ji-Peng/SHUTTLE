# 16-way AVX512 SM3 / KDF-SM3 (pseudoXOF)

The AVX-512 sibling of `../avx2`. Same design, twice the width: each call evaluates **16 independent SM3 instances at once** by packing one 32-bit SM3 word per lane of a 512-bit register. Output is **byte-identical** to the scalar reference (`../ref/auxfunc.c`) run 16 times.

## API

```c
int pseudoXOF_avx512(unsigned long long output_len_bits,
                     const unsigned char *const msg[16],
                     unsigned long long msg_len_bits,
                     unsigned char *const output[16]);

int sm3hash_avx512(const unsigned char *const msg[16],
                   unsigned long long msg_len_bits,
                   unsigned char *const digest[16]);   /* 32-byte digests */
```

`output_len_bits` and `msg_len_bits` are shared across the 16 lanes; `msg[]` and `output[]` are per-lane. Semantics are exactly those of the scalar `pseudoXOF` (KDF-SM3, GB/T 32918.4-2016 §5.4.3).

## Differences from the AVX2 build

- **16×16 dword transpose.** A 64-byte block of each message is a single ZMM load; the 16 loads are transposed to the message words $W[0..15]$ by a four-stage `unpack32 → unpack64 → shuffle_i32x4 → shuffle_i32x4` sequence. The stage configuration and the resulting column order were verified against simulated intrinsic semantics before use.
- **`ternarylogic` boolean functions.** The XOR-of-three folds ($FF_1, GG_1, P_0, P_1$, and both message-expansion XORs) become a single `_mm512_ternarylogic_epi32` with immediate `0x96`; majority ($FF_2$) uses `0xE8`; $GG_2 = ((f\oplus g)\,\&\,e)\oplus g$ uses `0xCA`. The truth tables are derived in the source comments.
- **Native rotate.** Rotations use `_mm512_rol_epi32` directly instead of the `slli | srli` pair.
- **Digest extraction.** The 8 byte-swapped state words are split into their low (lanes 0–7) and high (lanes 8–15) 256-bit halves and transposed with two 8×8 routines — the digest size matches a 256-bit half exactly, so unlike a zero-padded 16×16 transpose nothing is computed and discarded.

The scalar trailing-block / counter / padding logic and the constant-prefix precompute are identical to the AVX2 build, and reproduce `../ref/auxfunc.c` byte-for-byte.

## DRNG (16-way SM3 Hash-DRBG)

`drng_avx512.{c,h}` mirror `../avx2/drng_avx2.{c,h}` at 16 lanes, built on `sm3hash_avx512`: `init_random_number_avx512` / `get_random_number_avx512` with a 16-lane `DRNG_ctx_avx512`. `seed_len_bytes` / `random_number_len_bits` are shared; `seed[]` / `random_number[]` are per-lane. Byte-identical (output and state) to the scalar `../ref/drng.c` run 16 times.

## Build & test

```sh
make test          # XOF (vs ../ref/auxfunc.c) + DRNG (vs ../ref/drng.c)
make check-const   # audit SM3_T[]/SM3_IV[] against gen_sm3_const.py
```

`sm3_const.h` and `gen_sm3_const.py` are relative symlinks into `../avx2` (single source of truth for the SM3 constants). On an i7-11700K the small-message XOF case runs ~17–19× faster than 16 scalar calls; large messages are far faster thanks to the prefix precompute.

## Requirements

Needs `AVX512F`, `AVX512BW`, `AVX512DQ`, `AVX512VL` (the Makefile passes the matching `-m` flags). The i7-11700K (Rocket Lake) used for development has all of these.
