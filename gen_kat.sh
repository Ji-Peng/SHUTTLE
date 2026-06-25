#!/usr/bin/env bash
#
# gen_kat.sh -- generate the NGCC Known-Answer-Test (KAT) signature vectors
# for every SHUTTLE parameter set.
#
# For each instance (SHUTTLE-128/256/512) this builds the ICCS-provided KAT
# driver (KAT_SIG.c) against the SHUTTLE NGCC adapter (SIG_AlgorithmInstance.c)
# with OUTPUT_BLANK_TEST_VECTORS=0 (real vectors, fixed in the committed
# header), runs it, and copies the produced KAT_SIG_<instance>.txt into the
# output directory.
#
# Usage:
#   gen_kat.sh [--out <dir>] [--sha3] [--keep-build]
#
#   --out <dir>    Destination directory for the KAT_SIG_*.txt files
#                  (default: <repo>/SHUTTLE/dist/Test_Vectors).
#   --sha3         Build with the SHAKE (SHA3) XOF backend instead of the
#                  default NGCC SM3 Hash-DRBG backend, and emit the vectors
#                  into <out>/SHA3 so they do not overwrite the competition
#                  set.  (The harness always seeds via the SM3 Hash-DRBG; the
#                  flag only changes the scheme-internal XOF.)
#   --keep-build   Leave the temporary build directory in place (for debugging).
#
# The default (no flags) emits the NGCC competition vectors used by the
# submission.  The driver is built straight from the reference (pure C99)
# sources, so this script needs only gcc + a POSIX shell.

set -euo pipefail

# --- locate the SHUTTLE tree (this script lives at its root) ----------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
REF="${SCRIPT_DIR}/ref"
TOOLS="${SCRIPT_DIR}/tools"

# --- defaults ---------------------------------------------------------------
OUT="${SCRIPT_DIR}/dist/Test_Vectors"
XOF_MODE="NGCC"     # NGCC (SM3 Hash-DRBG) | SHA3 (SHAKE)
KEEP_BUILD=0
CC="${CC:-gcc}"

while [ $# -gt 0 ]; do
    case "$1" in
        --out)        OUT="$2"; shift 2 ;;
        --out=*)      OUT="${1#*=}"; shift ;;
        --sha3)       XOF_MODE="SHA3"; shift ;;
        --keep-build) KEEP_BUILD=1; shift ;;
        -h|--help)
            sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "gen_kat.sh: unknown argument '$1'" >&2; exit 2 ;;
    esac
done

# A SHA3 set lands in a side directory so it never shadows the NGCC vectors.
if [ "$XOF_MODE" = "SHA3" ]; then
    OUT="${OUT%/}/SHA3"
fi

# --- instance -> (mode, q-set) mapping --------------------------------------
# 128 -> q15361n256, 256 -> q61441n512, 512 -> q59393n1024.
modes="128 256 512"
qset_for() {
    case "$1" in
        128) echo "q15361n256" ;;
        256) echo "q61441n512" ;;
        512) echo "q59393n1024" ;;
        *)   echo "gen_kat.sh: bad mode '$1'" >&2; exit 2 ;;
    esac
}

# --- reference compile flags (NGCC mandated; pure ISO C99) ------------------
CFLAGS="-std=c99 -Wpedantic -Wall -Wextra -O2"

# XOF backend selection and its extra link sources.  The KAT harness ALWAYS
# seeds via the SM3 Hash-DRBG (drng.c + auxfunc.c), independent of the XOF.
if [ "$XOF_MODE" = "SHA3" ]; then
    MODEFLAG="-DSHA3_MODE"
    XOF_SRCS="symmetric.c fips202.c drng.c auxfunc.c"
else
    MODEFLAG=""
    XOF_SRCS="symmetric.c drng.c auxfunc.c"
fi

# Scheme + NGCC-infra source list (mirrors ref/Makefile's test_kat list, but
# with the ICCS KAT_SIG.c driver -- which carries main() -- instead of the
# internal test harness).
SCHEME_SRCS="SIG_AlgorithmInstance.c sign.c polyvec.c sampler.c sampler_u.c \
irs.c rounding.c packing.c poly.c poly_ntt.c reduce.c rans.c \
approx_exp.c approx_log.c"

mkdir -p "$OUT"
BUILD_DIR="$(mktemp -d "${TMPDIR:-/tmp}/shuttle_kat.XXXXXX")"
cleanup() { [ "$KEEP_BUILD" -eq 1 ] || rm -rf "$BUILD_DIR"; }
trap cleanup EXIT

echo "gen_kat.sh: XOF backend = $XOF_MODE, out = $OUT"

for m in $modes; do
    inst="SHUTTLE-${m}"
    qs="$(qset_for "$m")"
    bdir="${BUILD_DIR}/${inst}"
    mkdir -p "$bdir"

    # Build straight from ref/ (the C99 source of truth).  -I<ntt> resolves the
    # per-set scalar oracle header that poly_ntt.h pulls in by relative path.
    # shellcheck disable=SC2086
    "$CC" $CFLAGS $MODEFLAG -DSHUTTLE_MODE="$m" \
        -I"$REF" -I"$REF/ntt/$qs" -I"$TOOLS" \
        "$REF/KAT_SIG.c" \
        $(for s in $SCHEME_SRCS $XOF_SRCS; do printf '%s ' "$REF/$s"; done) \
        "$REF/ntt/$qs/ntt_ref.c" \
        -o "$bdir/kat_sig" \
        || { echo "gen_kat.sh: build failed for $inst" >&2; exit 1; }

    # The driver writes output/KAT_SIG_<instance>.txt under its own cwd.
    ( cd "$bdir" && ./kat_sig >/dev/null ) \
        || { echo "gen_kat.sh: run failed for $inst" >&2; exit 1; }

    src="${bdir}/output/KAT_SIG_${inst}.txt"
    [ -s "$src" ] || { echo "gen_kat.sh: empty KAT for $inst" >&2; exit 1; }
    cp "$src" "${OUT}/KAT_SIG_${inst}.txt"

    bytes=$(wc -c < "${OUT}/KAT_SIG_${inst}.txt")
    counts=$(grep -c '^Count = ' "${OUT}/KAT_SIG_${inst}.txt" || true)
    echo "  ok  ${inst}: ${counts} records, ${bytes} bytes -> ${OUT}/KAT_SIG_${inst}.txt"
done

echo "gen_kat.sh: done ($XOF_MODE)"
