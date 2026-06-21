# 8-way AVX2 SM3 / KDF-SM3 (pseudoXOF)

Vectorised drop-in companions for the scalar auxiliary functions in `../ref/auxfunc.c`. Each call evaluates **8 independent SM3 instances at once** by packing one 32-bit SM3 word per lane of a 256-bit register. The output is **byte-identical** to running the scalar reference 8 times — this is what keeps the future signature scheme KAT-compatible regardless of which lane width actually runs.

## Why 8-way (and not, say, 4-way)

SM3's state and message words are 32 bit. AVX2 holds eight 32-bit lanes per `__m256i`, so the natural width is 8. The design note in `agent/Prompts/SHUTTLE-Impl-Design.md` fixes the algorithm flow at an 8-way XOF on purpose: an 8-way batch is the common denominator that both AVX2 (one register) and a future AVX-512 8-of-16 path can serve without changing the pseudo-random byte sequence the signature consumes.

## API

```c
int pseudoXOF_avx2(unsigned long long output_len_bits,
                   const unsigned char *const msg[8],
                   unsigned long long msg_len_bits,
                   unsigned char *const output[8]);

int sm3hash_avx2(const unsigned char *const msg[8],
                 unsigned long long msg_len_bits,
                 unsigned char *const digest[8]);   /* 32-byte digests */
```

`output_len_bits` and `msg_len_bits` are **shared** across the 8 lanes; `msg[]` and `output[]` are **per-lane** pointers. This matches the brief: in `pseudoXOF_avx2` the two lengths are unified for all 8 lanes while the messages and outputs are 8 separate buffers. `pseudoXOF` follows GB/T 32918.4-2016 §5.4.3 (KDF-SM3): output block $i$ is $\mathrm{SM3}(\text{msg} \parallel ct)$ with the 32-bit counter $ct = i+1$ appended MSB-first.

## How it works

- **Load + transpose.** A 64-byte block of each of the 8 messages is loaded as two 256-bit vectors and turned into the 16 message words $W[0..15]$ by an 8×8 dword transpose (`unpacklo/hi` + `permute2x128`) followed by a per-dword byte swap, because SM3 reads words big-endian. The transpose sequence was verified against simulated intrinsic semantics before being committed.
- **Rounds.** Message expansion and all 64 rounds are a direct lane-wise translation of the reference; rotations are synthesised as `slli | srli`.
- **Constant-prefix precompute.** Because `msg` is identical across the $ct$ iterations, the complete 512-bit blocks that consist purely of message bits are compressed **once** into a cached state; only the trailing block(s) carrying $ct$ + padding are rebuilt per counter value. For a short seed (single block) this is a no-op; for multi-block messages it removes the reference's redundant re-hashing of the message prefix.
- **Trailing block.** The counter insertion, `normalize`, the `0x80` padding bit and the 64-bit length field are reproduced with the **exact** byte arithmetic of `../ref/auxfunc.c` (including its length-encoding split), so the result is bit-for-bit equal even for non-byte-aligned message lengths.

## Build & test

```sh
make test          # builds against ../ref/auxfunc.c and checks 392 vectors
make check-const   # re-derives SM3_T[]/SM3_IV[] from the spec and diffs the header
```

`make test` compares `pseudoXOF_avx2` / `sm3hash_avx2` against the scalar reference across a sweep of message lengths (empty, sub-byte, byte-aligned, block-boundary, multi-block, and the two-final-block window where $\text{msg\_tail\_bits} > 415$) and output lengths, then prints cycle counts. On an i7-11700K (Rocket Lake) the small-message case is ~8× faster than 8 scalar calls (pure SIMD width) and large messages are much faster still thanks to the prefix precompute.

## Auditable constants

Every constant lives in `sm3_const.h` and is reproduced by `gen_sm3_const.py`:

- `SM3_IV[8]` — the GB/T 32905-2016 initial value.
- `SM3_T[64]` — per-round constants $T_j = \mathrm{ROTL32}(T_{\text{base}}, i \bmod 32)$, with $T_{\text{base}} = \mathtt{0x79CC4519}$ for rounds 0–15 and $\mathtt{0x7A879D8A}$ for 16–63. These are precomputed only so the vector round can broadcast `SM3_T[i]`; the reference computes the identical value on the fly. Run `python3 gen_sm3_const.py --check` to verify.
