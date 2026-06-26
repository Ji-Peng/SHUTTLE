#!/usr/bin/env bash
#
# pack_ngcc.sh -- assemble the NGCC submission tree for the SHUTTLE signature
# scheme.
#
# Produces, under <out>/:
#
#   Implementations/
#     Reference_Implementation/SHUTTLE-{128,256,512}/   (pure ISO C99)
#     Optimized_Implementation/SHUTTLE-{128,256,512}/    (AVX2)
#     Additional_Implementation/SHUTTLE-{128,256,512}/   (AVX-512)
#     README
#   Test_Vectors/
#     KAT_SIG_SHUTTLE-{128,256,512}.txt
#
# Each implementation folder is SELF-CONTAINED: every symlink the working tree
# uses (the backends symlink shared scalar sources from ref/) is dereferenced
# into a real file, the four generated constant headers (normally reached via
# -I../tools) are copied in, the per-set scalar NTT oracle is vendored under
# ntt/<qset>/, and a folder-local Makefile builds the ICCS KAT driver with the
# mandated per-backend compile flags.  Nothing outside the folder is needed.
#
# Usage:
#   pack_ngcc.sh [--out <dir>] [--no-sanity]
#
#   --out <dir>    Destination root (default: <repo>/SHUTTLE/dist).
#   --no-sanity    Skip the standalone self-containment build/KAT check.
#
# The script does a clean rebuild of <out>, generates the Test_Vectors with
# gen_kat.sh, and (unless --no-sanity) proves self-containment by copying one
# Reference folder into an isolated temp dir with NO access to the original
# tree, building it there, and confirming its KAT matches Test_Vectors/.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
REF="${SCRIPT_DIR}/ref"
AVX2="${SCRIPT_DIR}/avx2"
AVX512="${SCRIPT_DIR}/avx512"
TOOLS="${SCRIPT_DIR}/tools"

OUT="${SCRIPT_DIR}/dist"
DO_SANITY=1
CC="${CC:-gcc}"

while [ $# -gt 0 ]; do
    case "$1" in
        --out)        OUT="$2"; shift 2 ;;
        --out=*)      OUT="${1#*=}"; shift ;;
        --no-sanity)  DO_SANITY=0; shift ;;
        -h|--help)    sed -n '2,34p' "$0"; exit 0 ;;
        *) echo "pack_ngcc.sh: unknown argument '$1'" >&2; exit 2 ;;
    esac
done

modes="128 256 512"
qset_for() {
    case "$1" in
        128) echo "q15361n256" ;;
        256) echo "q61441n512" ;;
        512) echo "q59393n1024" ;;
    esac
}

# ---------------------------------------------------------------------------
# File manifests.
#
# Scalar sources + headers are the SAME across all three backends; they resolve
# (possibly through a symlink) inside whichever backend directory we copy from,
# and `cp -L` dereferences the symlink into a real file.
# ---------------------------------------------------------------------------

# NGCC fixed-infrastructure files (ICCS-provided; ship verbatim).
INFRA_FILES="KAT_SIG.c SIG_AlgorithmInstance.c SIG_AlgorithmInstance.h \
auxfunc.c auxfunc.h drng.c drng.h README.txt"

# Scheme scalar C sources.  fips202.c is carried so a SHA3-mode rebuild also
# links; it is not on the default (NGCC SM3) link line.
SCHEME_C="sign.c polyvec.c sampler.c sampler_u.c irs.c rounding.c packing.c \
poly.c poly_ntt.c reduce.c rans.c approx_exp.c approx_log.c symmetric.c fips202.c"

# Every scheme/config header (copy the whole set -- small, and guarantees the
# #include closure is complete under either XOF mode).
SCHEME_H="api.h config.h namespace.h params.h poly.h poly_ntt.h polyvec.h \
reduce.h rounding.h packing.h rans.h sampler.h sampler_u.h sign.h irs.h \
approx_exp.h approx_log.h symmetric.h xof.h fips202.h"

# Generated constant headers normally found via -I../tools; vendored to root.
TOOLS_H="approx_exp_poly.h approx_log_poly.h rcdt_tables.h rounding_consts.h"

# Profiling probe header included (as a no-op) by the scheme sources; vendored
# under test/ so its relative include path "test/prof.h" still resolves.
PROF_H="test/prof.h"

# Per-backend SIMD extras (the N-way XOF + the vendored NTT live per <qset>).
AVX2_SIMD_C="symmetric_avx2.c drng_avx2.c auxfunc_avx2.c fips202x4.c f1600x4.S"
AVX2_SIMD_H="auxfunc_avx2.h drng_avx2.h fips202x4.h sm3_const.h simd_red.h"
AVX512_SIMD_C="symmetric_avx512.c drng_avx512.c auxfunc_avx512.c fips202x8.c keccakf1600x8.c"
AVX512_SIMD_H="auxfunc_avx512.h drng_avx512.h fips202x8.h sm3_const.h simd_red.h"

# ---------------------------------------------------------------------------
# Per-folder Makefile emitter.
#
#   $1 backend folder root (the assembled instance dir)
#   $2 mode (128|256|512)
#   $3 q-set (q15361n256 ...)
#   $4 backend label (reference|optimized|additional)
# ---------------------------------------------------------------------------
emit_makefile() {
    local dir="$1" mode="$2" qs="$3" backend="$4"
    local cflags simd_c ntt_simd inc_self

    inc_self='-I. -Intt/$(QSET)'
    case "$backend" in
        reference)
            cflags='-std=c99 -Wpedantic -Wall -Wextra -O2'
            simd_c=''
            ntt_simd=''
            ;;
        optimized)
            cflags='-O3 -march=x86-64 -mavx2 -mtune=native -flto -fomit-frame-pointer -std=c99 -Wpedantic -Wall -Wextra -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -DUSE_AVX2_XOF_NWAY'
            simd_c="symmetric_avx2.c drng_avx2.c auxfunc_avx2.c"
            ntt_simd='$(QSET)/ntt.S $(QSET)/ntt_consts.c'
            ;;
        additional)
            cflags='-O3 -march=x86-64 -mavx2 -mtune=native -flto -fomit-frame-pointer -std=c99 -Wpedantic -Wall -Wextra -mbmi2 -mpopcnt -mavx512f -mavx512bw -mavx512dq -mavx512vl -mavx512vbmi2 -DUSE_AVX512_NTT -DUSE_AVX512_SAMPLER -DUSE_AVX512_SHAKE8X -DUSE_AVX512_XOF_NWAY'
            simd_c="symmetric_avx512.c drng_avx512.c auxfunc_avx512.c"
            ntt_simd='$(QSET)/ntt_avx512.S $(QSET)/ntt_consts_avx512.c'
            ;;
    esac

    {
        cat <<EOF
# Self-contained NGCC build for ${backend^} Implementation of SHUTTLE-${mode}.
#
#   make            build the KAT generator (kat_sig) -> ./kat_sig
#   make vectors    build + run it; writes output/KAT_SIG_SHUTTLE-${mode}.txt
#   make clean
#
# This folder is self-contained: no symlinks, no include paths outside it.
# The compile flags are the NGCC-mandated set for this implementation tier.

CC      ?= ${CC}
MODE     = ${mode}
QSET     = ${qs}
INSTANCE = SHUTTLE-${mode}

CFLAGS  = ${cflags}
CPPFLAGS = ${inc_self} -DSHUTTLE_MODE=\$(MODE)

# Scheme + ICCS-infra scalar sources (the SM3 Hash-DRBG XOF backend).
SCHEME_SRCS = \\
  SIG_AlgorithmInstance.c sign.c polyvec.c sampler.c sampler_u.c irs.c \\
  rounding.c packing.c poly.c poly_ntt.c reduce.c rans.c approx_exp.c \\
  approx_log.c symmetric.c drng.c auxfunc.c \\
  ntt/\$(QSET)/ntt_ref.c
EOF
        if [ -n "$simd_c" ]; then
            cat <<EOF

# Backend SIMD kernels (vectorized NTT + N-way SM3 XOF).
SIMD_SRCS = ${ntt_simd} ${simd_c}
EOF
        else
            printf '\nSIMD_SRCS =\n'
        fi
        cat <<'EOF'

ALL_SRCS = $(SCHEME_SRCS) $(SIMD_SRCS)

.PHONY: all vectors clean

all: kat_sig

kat_sig: KAT_SIG.c $(ALL_SRCS)
	$(CC) $(CFLAGS) $(CPPFLAGS) KAT_SIG.c $(ALL_SRCS) -o $@

vectors: kat_sig
	./kat_sig
	@echo "wrote output/KAT_SIG_$(INSTANCE).txt"

clean:
	rm -f kat_sig
	rm -rf output
EOF
    } > "${dir}/Makefile"
}

# ---------------------------------------------------------------------------
# Assemble one instance folder for one backend.
# ---------------------------------------------------------------------------
assemble() {
    local backend="$1" srcdir="$2" mode="$3" destroot="$4"
    local qs dir f
    qs="$(qset_for "$mode")"
    dir="${destroot}/SHUTTLE-${mode}"
    mkdir -p "${dir}/ntt/${qs}" "${dir}/test" "${dir}/${qs}"

    # NGCC infra + scalar scheme sources/headers (dereferenced).
    for f in $INFRA_FILES $SCHEME_C $SCHEME_H; do
        cp -L "${srcdir}/${f}" "${dir}/${f}"
    done

    # config.h MUST come from the matching backend (its namespace tail differs:
    # _ref / _avx2 / _avx512).  It is a real fork already copied via SCHEME_H
    # from $srcdir, so nothing more to do.

    # Generated constant headers (vendored to root; included by basename).
    for f in $TOOLS_H; do
        cp "${TOOLS}/${f}" "${dir}/${f}"
    done

    # Profiling probe header (no-op; relative include "test/prof.h").
    cp -L "${REF}/${PROF_H}" "${dir}/test/prof.h"

    # Per-set scalar NTT oracle (relative include "ntt/<qset>/ntt_ref.h").
    cp -L "${REF}/ntt/${qs}/ntt_ref.c" "${dir}/ntt/${qs}/ntt_ref.c"
    cp -L "${REF}/ntt/${qs}/ntt_ref.h" "${dir}/ntt/${qs}/ntt_ref.h"

    # Backend SIMD kernels (only for optimized/additional).
    case "$backend" in
        optimized)
            for f in $AVX2_SIMD_C $AVX2_SIMD_H; do
                cp -L "${AVX2}/${f}" "${dir}/${f}"
            done
            cp -L "${AVX2}/${qs}/ntt.S"        "${dir}/${qs}/ntt.S"
            cp -L "${AVX2}/${qs}/ntt_consts.c" "${dir}/${qs}/ntt_consts.c"
            ;;
        additional)
            for f in $AVX512_SIMD_C $AVX512_SIMD_H; do
                cp -L "${AVX512}/${f}" "${dir}/${f}"
            done
            cp -L "${AVX512}/${qs}/ntt_avx512.S"        "${dir}/${qs}/ntt_avx512.S"
            cp -L "${AVX512}/${qs}/ntt_consts_avx512.c" "${dir}/${qs}/ntt_consts_avx512.c"
            ;;
        reference)
            rmdir "${dir}/${qs}" 2>/dev/null || true
            ;;
    esac

    emit_makefile "$dir" "$mode" "$qs" "$backend"
}

# ---------------------------------------------------------------------------
# Top-level Implementations/README.
# ---------------------------------------------------------------------------
write_impl_readme() {
    local impl="$1"
    cat > "${impl}/README" <<'EOF'
SHUTTLE signature -- NGCC Implementations
=========================================

Three implementation tiers, each covering all parameter sets
(SHUTTLE-128, SHUTTLE-256, SHUTTLE-512):

  Reference_Implementation/   Pure ISO C99 (no assembly, intrinsics, or
                              compiler extensions).  Builds warning-clean
                              under: gcc -std=c99 -Wpedantic -Wall -Wextra -O2.
                              This is the correctness reference and the KAT
                              oracle.

  Optimized_Implementation/   AVX2-optimized for mainstream 64-bit x86 PCs.
                              Vectorized NTT + N-way SM3 XOF.  Built with:
                              gcc -O3 -march=x86-64 -mavx2 -mtune=native -flto
                              -fomit-frame-pointer -std=c99 -Wpedantic -Wall
                              -Wextra (a pure-AVX2 path, no AVX-512).

  Additional_Implementation/  AVX-512-optimized variant (extra tier).  Same as
                              Optimized but with the AVX-512 subsets enabled
                              (-mavx512f/bw/dq/vl/vbmi2).  Requires an AVX-512
                              host to run.

Each SHUTTLE-<set> folder is self-contained and ships with a Makefile:

  make           build the ICCS KAT generator (kat_sig).
  make vectors   build + run it; writes output/KAT_SIG_SHUTTLE-<set>.txt.
  make clean

Per-folder file inventory:

  KAT_SIG.c                 ICCS-provided KAT driver (SHALL NOT be modified).
  SIG_AlgorithmInstance.{c,h}  NGCC programming-interface adapter for SHUTTLE
                            (sig_keygen / sig_sign / sig_verify) and the
                            ALGORITHM_INSTANCE / OUTPUT_BLANK_TEST_VECTORS
                            macros.
  drng.{c,h}, auxfunc.{c,h} ICCS-provided SM3 Hash-DRBG + auxiliary functions
                            (SHALL NOT be modified).
  README.txt                ICCS instructions for the infrastructure files.
  sign.c, polyvec.c, sampler.c, sampler_u.c, irs.c, rounding.c, packing.c,
  poly.c, poly_ntt.c, reduce.c, rans.c, approx_exp.c, approx_log.c,
  symmetric.c, fips202.c    SHUTTLE scheme sources.
  *.h                       SHUTTLE scheme + generated-constant headers
                            (approx_exp_poly.h, approx_log_poly.h,
                            rcdt_tables.h, rounding_consts.h are reproducible
                            constant tables).
  ntt/<qset>/ntt_ref.{c,h}  Per-set scalar NTT.
  <qset>/ntt*.{S,c}         Per-set vectorized NTT (Optimized/Additional only).
  test/prof.h               No-op profiling-probe header.

The matching Known-Answer-Test vectors are in ../Test_Vectors/.
EOF
}

# ===========================================================================
# Drive.
# ===========================================================================
echo "pack_ngcc.sh: out = $OUT"
echo "pack_ngcc.sh: clean rebuild"
rm -rf "$OUT"
IMPL="${OUT}/Implementations"
TV="${OUT}/Test_Vectors"
mkdir -p "${IMPL}/Reference_Implementation" \
         "${IMPL}/Optimized_Implementation" \
         "${IMPL}/Additional_Implementation" \
         "${TV}"

for m in $modes; do
    assemble reference  "$REF"    "$m" "${IMPL}/Reference_Implementation"
    assemble optimized  "$AVX2"   "$m" "${IMPL}/Optimized_Implementation"
    assemble additional "$AVX512" "$m" "${IMPL}/Additional_Implementation"
    echo "  assembled SHUTTLE-${m} (reference / optimized / additional)"
done

write_impl_readme "$IMPL"

# Test vectors (default NGCC SM3 backend).
echo "pack_ngcc.sh: generating Test_Vectors"
bash "${SCRIPT_DIR}/gen_kat.sh" --out "$TV"

# ---------------------------------------------------------------------------
# Sanity: build every assembled folder in place (proves the flags compile),
# then prove SELF-CONTAINMENT by building one Reference folder in an isolated
# copy with no access to the original tree, and matching its KAT.
# ---------------------------------------------------------------------------
build_in_place() {
    local dir="$1" label="$2"
    if ( cd "$dir" && make -s kat_sig >/dev/null 2>build.log ); then
        echo "  build OK   : ${label}"
        rm -f "${dir}/build.log" "${dir}/kat_sig"
    else
        echo "  build FAIL : ${label}"
        sed 's/^/      /' "${dir}/build.log" | head -20
        return 1
    fi
}

host_has_avx512() {
    grep -qm1 'avx512f' /proc/cpuinfo 2>/dev/null
}

if [ "$DO_SANITY" -eq 1 ]; then
    echo "pack_ngcc.sh: sanity -- in-place build of every assembled folder"
    for m in $modes; do
        build_in_place "${IMPL}/Reference_Implementation/SHUTTLE-${m}" "reference SHUTTLE-${m}"
        build_in_place "${IMPL}/Optimized_Implementation/SHUTTLE-${m}" "optimized SHUTTLE-${m}"
        if host_has_avx512; then
            build_in_place "${IMPL}/Additional_Implementation/SHUTTLE-${m}" "additional SHUTTLE-${m} (avx512)"
        else
            echo "  build SKIP : additional SHUTTLE-${m} (no AVX-512 on this host)"
        fi
    done

    echo "pack_ngcc.sh: sanity -- isolated self-containment build + KAT match"
    ISO="$(mktemp -d "${TMPDIR:-/tmp}/shuttle_iso.XXXXXX")"
    trap 'rm -rf "$ISO"' EXIT
    for m in $modes; do
        inst="SHUTTLE-${m}"
        # Copy ONLY the one folder into a fresh location -- nothing from the
        # original tree is reachable from there.
        cp -a "${IMPL}/Reference_Implementation/${inst}" "${ISO}/${inst}"
        if ! ( cd "${ISO}/${inst}" && make -s vectors >/dev/null 2>iso.log ); then
            echo "  FAIL: ${inst} did not build standalone in isolation"
            sed 's/^/      /' "${ISO}/${inst}/iso.log" | head -20
            exit 1
        fi
        iso_kat="${ISO}/${inst}/output/KAT_SIG_${inst}.txt"
        ref_kat="${TV}/KAT_SIG_${inst}.txt"
        if cmp -s "$iso_kat" "$ref_kat"; then
            echo "  PASS: ${inst} standalone build + KAT matches Test_Vectors/"
        else
            echo "  FAIL: ${inst} standalone KAT differs from Test_Vectors/"
            exit 1
        fi
    done
fi

# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------
echo
echo "==================== assembled tree ===================="
if command -v tree >/dev/null 2>&1; then
    tree -L 3 --noreport "$IMPL"
else
    ( cd "$OUT" && find Implementations -maxdepth 3 -type d | sort | sed 's/^/  /' )
fi
echo
echo "Test_Vectors:"
for m in $modes; do
    f="${TV}/KAT_SIG_SHUTTLE-${m}.txt"
    printf '  %-32s %8s bytes  (%s records)\n' \
        "KAT_SIG_SHUTTLE-${m}.txt" "$(wc -c < "$f")" "$(grep -c '^Count = ' "$f")"
done
echo
echo "Reference folder sizes:"
for m in $modes; do
    d="${IMPL}/Reference_Implementation/SHUTTLE-${m}"
    printf '  %-28s %s\n' "SHUTTLE-${m}" "$(du -sh "$d" | cut -f1)"
done
echo
echo "pack_ngcc.sh: done -> $OUT"
