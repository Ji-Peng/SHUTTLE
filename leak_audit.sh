#!/bin/sh
# Memory-leak / memory-safety gate.
#
# Builds + runs test/diff_fuzz.c (keygen -> sign -> verify -> tamper, the path
# that exercises the N-way gauss bulk-preset malloc/free in expand_s/sample_y)
# under AddressSanitizer + LeakSanitizer on each backend x mode x xof.  A leak
# (e.g. a missed free() on the bulk buffer or the gss[] scratch, on the success
# path OR the malloc-failure fallback) or any heap error aborts with a nonzero
# exit and a FAIL line.
#
# This is the leak counterpart to integration_audit_matrix.sh's differential
# fuzz: that one proves the backends agree byte-for-byte, this one proves none
# of them leak while doing it.  Run both before shipping a sampler change.
#
# Knobs (env): LEAK_AUDIT_BACKENDS, LEAK_AUDIT_MODES, LEAK_AUDIT_XOF (space-
# separated, default "ngcc sha3"), LEAK_FUZZ_N (cases per run, default 8).
set -eu
cd "$(dirname "$0")"
. ./tools/audit_lib.sh

mkdir -p test
REPORT=${LEAK_AUDIT_REPORT:-test/leak_audit.txt}
BACKENDS=${LEAK_AUDIT_BACKENDS:-ref avx2 avx512}
MODES=${LEAK_AUDIT_MODES:-128 256 512}
XOFS=${LEAK_AUDIT_XOF:-ngcc sha3}
LEAK_FUZZ_N=${LEAK_FUZZ_N:-8}

# -fsanitize=address pulls in LeakSanitizer at exit; ASAN_OPTIONS forces a
# nonzero exit on any leak/error so the gate actually fails.  -g keeps the
# leak report's stack readable.  We deliberately do NOT add UBSan: the crypto
# core relies on defined unsigned wraparound that UBSan would flag as noise.
ASAN_CFLAGS="-fsanitize=address -fno-omit-frame-pointer -g"
ASAN_RUNOPTS="detect_leaks=1:abort_on_error=1:exitcode=99"

overall=PASS
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT INT HUP TERM

{
    printf '=== SHUTTLE memory-leak audit (ASan/LSan) ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\nmodes: %s\nxof: %s\nLEAK_FUZZ_N: %s\n\n' \
        "$BACKENDS" "$MODES" "$XOFS" "$LEAK_FUZZ_N"
} > "$REPORT"

for xof in $XOFS; do
    for backend in $BACKENDS; do
        for mode in $MODES; do
            bin="out/leak_fuzz_${mode}"
            log="$tmpdir/$backend.$mode.$xof"
            if ! audit_build "$backend" "$mode" "$xof" test/diff_fuzz.c "$bin" \
                    -DDIFF_FUZZ_N="$LEAK_FUZZ_N" $ASAN_CFLAGS \
                    > "$log.build" 2>&1; then
                printf 'BUILDFAIL backend=%s mode=%s xof=%s\n' \
                    "$backend" "$mode" "$xof" | tee -a "$REPORT"
                cat "$log.build"
                overall=FAIL
                continue
            fi
            if ASAN_OPTIONS="$ASAN_RUNOPTS" "$backend/$bin" \
                    > "$log.run" 2>&1; then
                printf 'PASS backend=%s mode=%s xof=%s cases=%s\n' \
                    "$backend" "$mode" "$xof" "$LEAK_FUZZ_N" | tee -a "$REPORT"
            else
                printf 'LEAK/ASAN-FAIL backend=%s mode=%s xof=%s\n' \
                    "$backend" "$mode" "$xof" | tee -a "$REPORT"
                cat "$log.run" | tee -a "$REPORT"
                overall=FAIL
            fi
        done
    done
done

printf '\nOVERALL: %s\n' "$overall" | tee -a "$REPORT"
[ "$overall" = PASS ]
