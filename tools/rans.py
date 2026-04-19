"""Pure-Python reference rANS (range Asymmetric Numeral Systems) encoder/decoder.

Used in Phase 6b-4 to:
  1. Validate frequency tables generated from empirical hint histograms.
  2. Cross-check the C-side rANS implementation (Phase 6b-5) against a
     well-tested reference.

The design mirrors Fabian Giesen's rANS_byte engine (rygorous/ryg_rans),
which is what HAETAE ships. Byte-wise renormalization, 32-bit state,
encode proceeds right-to-left, decode left-to-right over a reversed stream.

Frequency table convention
--------------------------

A frequency table is a pair (syms, freqs) with aligned lists:
  - syms  : the symbols in a canonical order (ints).
  - freqs : the integer frequency counts, summing to prob_total = 2**prob_bits.

Each symbol also has a cumulative start, cdf[i] = sum(freqs[0..i-1]), used
by the encoder. The decoder reverses cdf lookup via a table
  sym_lookup[c] = index of symbol whose bucket contains c,
for c in [0, prob_total).

Only positive frequencies are allowed (every symbol in syms has freq >= 1);
symbols that would map to freq=0 in quantization must be handled up front
by the caller (either merged into "escape" or dropped via rejection).

References:
    Duda, J., "Asymmetric numeral systems: entropy coding combining speed
        of Huffman coding with compression rate of arithmetic coding",
        arXiv:1311.2540, 2013.
    Giesen, F., "ryg_rans", https://github.com/rygorous/ryg_rans, 2014.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence


# ---------------------------------------------------------------------------
# Constants (chosen to match HAETAE's rans_byte.h defaults).
# ---------------------------------------------------------------------------

RANS_L = 1 << 23                # low threshold; below this, refill from stream
BYTE_BITS = 8
BYTE_SIZE = 1 << BYTE_BITS

STATE_MASK = (1 << 32) - 1      # keep state in uint32 range


# ---------------------------------------------------------------------------
# Frequency-table helpers
# ---------------------------------------------------------------------------

@dataclass
class RansTable:
    """Immutable container for an rANS frequency table."""
    syms: List[int]            # symbol values (ints), in canonical order
    freqs: List[int]           # frequency counts; sum = 2^prob_bits
    prob_bits: int             # log2(total frequency)

    # Derived tables, filled by _build_cdf / _build_sym_lookup.
    cdf: List[int] = None
    sym_index: dict = None     # symbol -> index in syms list
    sym_lookup: List[int] = None   # [prob_total] -> index in syms list

    def __post_init__(self):
        total = 1 << self.prob_bits
        if sum(self.freqs) != total:
            raise ValueError(
                f"freqs sum to {sum(self.freqs)}, expected 2^{self.prob_bits} = {total}"
            )
        if len(self.syms) != len(self.freqs):
            raise ValueError("syms and freqs length mismatch")
        if any(f <= 0 for f in self.freqs):
            raise ValueError("all frequencies must be >= 1")

        # Build cumulative distribution and symbol lookup table.
        self.sym_index = {s: i for i, s in enumerate(self.syms)}

        self.cdf = [0] * (len(self.syms) + 1)
        for i, f in enumerate(self.freqs):
            self.cdf[i + 1] = self.cdf[i] + f

        self.sym_lookup = [0] * total
        for i, f in enumerate(self.freqs):
            start = self.cdf[i]
            for k in range(start, start + f):
                self.sym_lookup[k] = i

    @property
    def prob_total(self) -> int:
        return 1 << self.prob_bits


# ---------------------------------------------------------------------------
# Encoder
# ---------------------------------------------------------------------------

class RansEncoder:
    """Byte-wise rANS encoder.

    Call .put(symbol) for each symbol, IN REVERSE ORDER (last symbol first);
    rANS is last-in-first-out. After all symbols are pushed, call .finalize()
    to get the final byte stream (flush state, then reverse).
    """

    def __init__(self, table: RansTable):
        self.table = table
        self.state = RANS_L
        self.out_bytes: List[int] = []   # built in reverse; reversed by finalize()

    def put(self, sym) -> None:
        t = self.table
        if sym not in t.sym_index:
            raise KeyError(f"symbol {sym!r} not in frequency table")
        idx = t.sym_index[sym]
        freq = t.freqs[idx]
        start = t.cdf[idx]

        # Renormalize: if x would exceed freq * (RANS_L >> prob_bits) << 8 after
        # encoding, emit a byte and shift.
        # Equivalent guard: while x >= freq * (RANS_L >> prob_bits) * BYTE_SIZE.
        x_max = (RANS_L >> t.prob_bits) * BYTE_SIZE * freq
        x = self.state
        while x >= x_max:
            self.out_bytes.append(x & 0xFF)
            x >>= BYTE_BITS

        # Encode: x' = (x // freq) << prob_bits + start + (x % freq).
        self.state = ((x // freq) << t.prob_bits) + start + (x % freq)
        self.state &= STATE_MASK

    def finalize(self) -> bytes:
        """Flush the state into the output stream and return a byte buffer."""
        x = self.state
        for _ in range(4):
            self.out_bytes.append(x & 0xFF)
            x >>= BYTE_BITS
        # out_bytes was built with "most-recently-emitted byte last"; rANS
        # decode reads bytes in the opposite order, so reverse now.
        self.out_bytes.reverse()
        data = bytes(self.out_bytes)
        # Reset encoder for reuse.
        self.state = RANS_L
        self.out_bytes = []
        return data


# ---------------------------------------------------------------------------
# Decoder
# ---------------------------------------------------------------------------

class RansDecoder:
    """Byte-wise rANS decoder that mirrors RansEncoder.

    Initialize with (table, byte_stream). Call .get() to pop one symbol.
    """

    def __init__(self, table: RansTable, stream: bytes):
        self.table = table
        self.stream = stream
        self.pos = 0

        # Prime the 32-bit state from the first 4 bytes (they were the
        # last 4 bytes written during finalize()).
        if len(stream) < 4:
            raise ValueError("stream too short to decode rANS state")
        x = 0
        for i in range(4):
            x = (x << 8) | stream[i]
        self.state = x
        self.pos = 4

    def get(self):
        t = self.table
        prob_mask = t.prob_total - 1

        # Extract bucket index c in [0, prob_total).
        c = self.state & prob_mask
        idx = t.sym_lookup[c]
        freq = t.freqs[idx]
        start = t.cdf[idx]

        # Update state: x = freq * (x >> prob_bits) + c - start.
        self.state = freq * (self.state >> t.prob_bits) + (c - start)

        # Renormalize: pull bytes while state is below L.
        while self.state < RANS_L:
            if self.pos >= len(self.stream):
                raise RuntimeError("rANS decode underflow (stream exhausted)")
            self.state = ((self.state << BYTE_BITS) | self.stream[self.pos]) & STATE_MASK
            self.pos += 1

        return t.syms[idx]


# ---------------------------------------------------------------------------
# Convenience API
# ---------------------------------------------------------------------------

def encode(symbols: Sequence[int], table: RansTable) -> bytes:
    """Encode a sequence of symbols with the given frequency table."""
    enc = RansEncoder(table)
    for s in reversed(symbols):
        enc.put(s)
    return enc.finalize()


def decode(stream: bytes, table: RansTable, n: int) -> List[int]:
    """Decode n symbols from the byte stream."""
    dec = RansDecoder(table, stream)
    return [dec.get() for _ in range(n)]


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _selftest():
    # Toy distribution: 4 symbols, skewed.
    syms = [-2, -1, 0, 1, 2]
    raw_freqs = [80, 300, 600, 300, 80]   # sum 1360; need to normalize to 2^bits.
    # Scale to 2^12 = 4096.
    import math
    scale = 4096 / sum(raw_freqs)
    freqs = [max(1, round(f * scale)) for f in raw_freqs]
    diff = 4096 - sum(freqs)
    freqs[len(freqs) // 2] += diff
    assert sum(freqs) == 4096

    table = RansTable(syms=syms, freqs=freqs, prob_bits=12)

    # Encode a random symbol sequence, decode, compare.
    import random
    random.seed(0xC0DE)
    msg = [random.choices(syms, weights=freqs, k=1)[0] for _ in range(1000)]

    buf = encode(msg, table)
    recovered = decode(buf, table, len(msg))

    assert recovered == msg, "rANS round-trip failed"
    print(f"rans.py self-test: OK. {len(msg)} symbols -> {len(buf)} bytes "
          f"= {8 * len(buf) / len(msg):.3f} bits/symbol "
          f"(Shannon entropy ~= {sum((-f/4096 * math.log2(f/4096)) for f in freqs):.3f} bits/symbol).")


if __name__ == "__main__":
    _selftest()
