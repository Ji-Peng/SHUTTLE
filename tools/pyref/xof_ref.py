#!/usr/bin/env python3
"""xof_ref.py -- the SHUTTLE XOF compatibility layer + 16-stream fill engines
for the Python reference.  Byte-exact to ref/symmetric.c + ref/xof.h and the
ref/polyvec.c stream structs (uniform_stream / gauss_stream).

Two MODEs (ref/xof.h):
  SHA3_MODE : xof128 = SHAKE128 (rate 168), xof256 = SHAKE256 (rate 136),
              each init+absorb_once over the whole (tag||seed||nonce) buffer,
              then incremental rate-buffered squeeze.  hashlib provides both.
  NGCC_MODE : both xof128 and xof256 collapse onto the SM3 Hash-DRBG
              (drng_ref.DRNG): init = instantiate(seed), squeeze = generate.

A `XofCtx` wraps the chosen primitive so the producers (ExpandA / ExpandS /
SampleY / SampleC / IRS) are written once against init+squeeze.
"""
import hashlib

import drng_ref


# ===================================================================== #
#  SHAKE incremental squeeze (mirror Keccak rate-buffered squeeze)      #
# ===================================================================== #
# hashlib.shake_*.digest(n) re-runs from scratch each call, but a continuous
# Keccak squeeze is a single stream.  We emulate the continuous stream by
# squeezing a growing prefix and slicing: digest(total) is byte-stable as the
# squeeze count grows (SHAKE is a prefix-stable stream), so consecutive
# squeezes off the same ctx chain correctly.
class _ShakeCtx:
    def __init__(self, which, seed):
        self.which = which                       # 128 or 256
        self.seed = bytes(seed)
        self.pos = 0                             # bytes already consumed

    def squeeze(self, n):
        total = self.pos + n
        h = hashlib.shake_128() if self.which == 128 else hashlib.shake_256()
        h.update(self.seed)
        out = h.digest(total)[self.pos:total]
        self.pos = total
        return out


# ===================================================================== #
#  NGCC SM3 Hash-DRBG ctx (mirror drng.c via drng_ref)                  #
# ===================================================================== #
class _DrngCtx:
    def __init__(self, seed):
        self.d = drng_ref.DRNG.instantiate(bytes(seed))

    def squeeze(self, n):
        return drng_ref.get_random_number(self.d, n * 8)   # BYTES -> BITS shim


class XofCtx:
    """Unified XOF ctx.  `which` is 128 or 256 (selects the primitive in
    SHA3_MODE; ignored in NGCC_MODE where both collapse onto SM3)."""

    def __init__(self, mode_sha3, which, seed):
        self.mode_sha3 = mode_sha3
        if mode_sha3:
            self.impl = _ShakeCtx(which, seed)
        else:
            self.impl = _DrngCtx(seed)

    def squeeze(self, n):
        return self.impl.squeeze(n)


def xof256_init(mode_sha3, seed):
    return XofCtx(mode_sha3, 256, seed)


def xof128_init(mode_sha3, seed):
    return XofCtx(mode_sha3, 128, seed)


# ===================================================================== #
#  uniform_stream (ExpandA): one-squeeze-per-fill, UNIFORM_BLOCK bytes  #
# ===================================================================== #
UNIFORM_BLOCK = 4096


class UniformStream:
    """Mirror ref/polyvec.c uniform_stream (xof128).  nonce = tag||seed||
    LE16(lane)||LE16(refill); each fill draws UNIFORM_BLOCK in one squeeze."""

    def __init__(self, mode_sha3, tag, seed, lane):
        self.mode_sha3 = mode_sha3
        self.tag = tag
        self.seed = bytes(seed)
        self.lane = lane
        self.refill = 0
        self.buf = b""
        self.pos = 0
        self.avail = 0
        self._fill()

    def _nonce(self):
        return (bytes([self.tag]) + self.seed
                + self.lane.to_bytes(2, "little")
                + self.refill.to_bytes(2, "little"))

    def _fill(self):
        ctx = xof128_init(self.mode_sha3, self._nonce())
        self.buf = ctx.squeeze(UNIFORM_BLOCK)
        self.pos = 0
        self.avail = UNIFORM_BLOCK

    def next_candidate(self, nbytes, mask):
        """us_next_candidate: pull nbytes (BQ), refill across the block edge."""
        if self.pos + nbytes > self.avail:
            self.refill += 1
            self._fill()
        x = 0
        for i in range(nbytes):
            x |= self.buf[self.pos + i] << (8 * i)
        self.pos += nbytes
        return x & mask


# ===================================================================== #
#  gauss_stream (ExpandS / SampleY): GAUSS_STREAM_BLOCK one-squeeze fill #
# ===================================================================== #
class GaussStream:
    """Mirror ref/polyvec.c gauss_stream + gs_fill/gs_ensure.  The buffer is
    a persistent byte array; gs_ensure memmoves the leftover to the front and
    appends a fresh GAUSS_STREAM_BLOCK squeeze (refill counter bumped)."""

    def __init__(self, mode_sha3, tag, seed, lane, block):
        self.mode_sha3 = mode_sha3
        self.tag = tag
        self.seed = bytes(seed)
        self.lane = lane
        self.block = block               # GAUSS_STREAM_BLOCK
        self.refill = 0
        self.buf = bytearray()
        self.pos = 0
        self.avail = 0

    def _nonce(self):
        return (bytes([self.tag]) + self.seed
                + self.lane.to_bytes(2, "little")
                + self.refill.to_bytes(2, "little"))

    def _fill(self):
        # gs_fill: draw GAUSS_STREAM_BLOCK bytes appended after avail.
        ctx = xof256_init(self.mode_sha3, self._nonce())
        chunk = ctx.squeeze(self.block)
        # ensure buf has room: buf is logically [0, avail); append at avail.
        if len(self.buf) < self.avail + self.block:
            self.buf.extend(bytes(self.avail + self.block - len(self.buf)))
        self.buf[self.avail:self.avail + self.block] = chunk
        self.avail += self.block

    def ensure(self, need):
        if self.avail - self.pos >= need:
            return
        left = self.avail - self.pos
        if left:
            self.buf[0:left] = self.buf[self.pos:self.pos + left]
        self.pos = 0
        self.avail = left
        self.refill += 1
        self._fill()
