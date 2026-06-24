#!/bin/sh
# mklinks.sh - materialize the avx2/ and avx512/ backend trees by relative
# symlinking the SHARED scalar source from ref/, skipping the per-backend FORK
# set (Overview 3 / P01-T10).
#
# ref/ is the single source of truth.  In avx2/ and avx512/ the REAL FORKS are:
#   config.h, polyvec.c, sampler.c, irs.c, rounding.c, sign.c
#   + the AVX-only files (fips202x4/x8.*, f1600x4.S/keccakf1600x8.c, the SM3
#     N-way auxfunc_avx2/512.* + drng_avx2/512.*, the vendored NTT asm, and
#     ntt_avx_decls.h).
# Everything else is a relative symlink `ln -s ../ref/<f>`.
#
# This script is IDEMPOTENT and re-runnable: it (re)creates only the symlinks
# for shared files that actually exist in ref/ right now and are not in the
# fork/skip set, and never clobbers a real (non-symlink) file in a backend dir.
# As later plans add shared sources to ref/, re-running this picks them up.
#
# Run from anywhere:  sh tools/mklinks.sh   (or: cd tools && sh mklinks.sh)

set -eu

# Resolve the SHUTTLE root from this script's location (tools/..).
here=$(CDPATH= cd "$(dirname "$0")" && pwd)
root=$(CDPATH= cd "$here/.." && pwd)
ref="$root/ref"

# --- FORK / SKIP set: real files per backend, never symlinked from ref/. ---
# (Files in ref/ whose vectorized twin lives in the backend dir, plus NGCC
#  fixed-infra that the backend builds via -I../ref rather than symlinking.)
is_fork() {
    case "$1" in
    config.h | polyvec.c | sampler.c | irs.c | rounding.c | sign.c) return 0 ;;
    # P03: the NTT shim has a per-backend FORK body (avx2/avx512 poly_ntt.c
    # hard-code the two asm signature families + signed canonicalization); the
    # scalar ref/poly_ntt.c is NOT symlinked into the backends.
    poly_ntt.c) return 0 ;;
    # AVX-only / N-way XOF files are real in the backend dirs (P02).
    fips202x4.* | fips202x8.* | f1600x4.* | keccakf1600x8.* | ntt_avx_decls.h) return 0 ;;
    auxfunc_avx2.* | auxfunc_avx512.* | drng_avx2.* | drng_avx512.*) return 0 ;;
    sm3_const.h | gen_sm3_const.py) return 0 ;;
    # NGCC fixed infra: consumed via -I../ref, not symlinked.
    SIG_AlgorithmInstance.* | KAT_SIG.c | README.* | drng.* | auxfunc.*) return 0 ;;
    # build/dir artifacts and the test dir are not symlinked file-by-file.
    Makefile | out | test | explore | .gitignore) return 0 ;;
    *) return 1 ;;
    esac
}

link_into() {
    backend="$1"
    dir="$root/$backend"
    mkdir -p "$dir"
    n=0
    for path in "$ref"/*; do
        f=$(basename "$path")
        # Only link regular files (skip directories like test/ out/).
        [ -f "$path" ] || continue
        if is_fork "$f"; then
            continue
        fi
        target="../ref/$f"
        dest="$dir/$f"
        if [ -L "$dest" ]; then
            # already a symlink: refresh it (idempotent).
            ln -sf "$target" "$dest"
            n=$((n + 1))
        elif [ -e "$dest" ]; then
            # a REAL file lives here (a fork we didn't list, or stale): leave it.
            echo "  skip (real file present): $backend/$f"
        else
            ln -s "$target" "$dest"
            n=$((n + 1))
        fi
    done
    # Mirror the test/ source dir: backends build `test/<t>.c` (resolved in the
    # backend dir), so symlink each ref/test/*.c|*.h into <backend>/test/.  The
    # backends fork no test sources, so all of them are symlinked.
    if [ -d "$ref/test" ]; then
        mkdir -p "$dir/test"
        for tpath in "$ref"/test/*; do
            [ -f "$tpath" ] || continue
            tf=$(basename "$tpath")
            tdest="$dir/test/$tf"
            if [ -L "$tdest" ]; then
                ln -sf "../../ref/test/$tf" "$tdest"
            elif [ -e "$tdest" ]; then
                echo "  skip (real file present): $backend/test/$tf"
            else
                ln -s "../../ref/test/$tf" "$tdest"
            fi
        done
    fi
    echo "$backend: linked/refreshed $n shared file(s) + test/ from ref/"
}

link_into avx2
link_into avx512
echo "mklinks.sh: done"
