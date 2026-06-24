#!/bin/sh
# End-to-end SHUTTLE gate-of-gates (P13-T11). Thin fail-fast wrapper around the
# Makefile targets (owned by P01/P10) + the P13 audit matrices. The Makefile
# target BODIES have one owner; this script wires their AUDIT SEMANTICS in the
# fail-fast order from 13-Security-Audit.md.
#
# Some Makefile targets named in the plan are owned by other plans and may not
# be wired yet in the in-flight tree (check-rans-failrate, the per-backend
# ct-scan/dudect-sign targets). This wrapper runs every gate it can and treats a
# *missing Makefile target* as a SKIP (loud), while a *failing* gate is fatal.
# NOTE: SHUTTLE ships NO opt-in float-exp variant, so check-kat-avx2exp /
# check-kat-avx512ifma and AVX512IFMAEXP=1 are intentionally absent (P13-T12).
set -eu
cd "$(dirname "$0")"

run() {
    printf '\n== %s ==\n' "$*"
    "$@"
}

# run_make BACKEND TARGET [extra make args] -- run only if the target exists;
# else print SKIP and continue (the target's owner has not wired it yet).
run_make() {
    backend=$1; target=$2; shift 2
    if make -C "$backend" -n "$target" >/dev/null 2>&1; then
        printf '\n== make -C %s %s %s ==\n' "$backend" "$target" "$*"
        make -C "$backend" "$target" "$@"
    else
        printf '\n== make -C %s %s: SKIP (target not wired by its owner yet) ==\n' "$backend" "$target"
    fi
}

run_make ref tables                  # materialize generated tables (P01..P10)
run_make ref check                   # correctness + KAT-accept, all TESTS x MODES
run_make ref check-rans-failrate     # rANS overflow-restart rate <= 2^-35 (P10)
run ./ct_scan_matrix.sh              # MANDATORY static CT scan (HARD gate)
run ./security_audit_matrix.sh       # fault injection + parser-negative (+ trace-KAT)
run ./integration_audit_matrix.sh    # cross-backend diff-fuzz + perf-checklist

run_make avx2 check-consts           # reproducible-constants drift gate
run_make avx2 check                  # correctness + ref==avx2 byte equality
run_make avx2 check-kat              # shared KAT (ref==avx2==avx512)
run_make avx2 check-rans-failrate
run_make avx2 ct-scan                # single-backend CT scan -> avx2/test/ct_scan_avx2.txt
run_make avx2 check-no-avx512        # prove avx2 build has NO zmm/EVEX/AVX-512 opcode
run_make avx2 check-refill           # Gaussian refill correctness (P07)

run_make avx512 check-consts
run_make avx512 check
run_make avx512 check-kat
run_make avx512 check-rans-failrate
run_make avx512 ct-scan
run_make avx512 check-refill         # NB: no check-no-avx512 for avx512

# Optional KyberSlash/TIMECOP dynamic cross-check (one-time patched-Valgrind
# build). Off by default; the variable-latency property is already gated
# statically by ct_scan and statistically by run_dudect_all.sh.
if [ "${TIMECOP:-0}" = 1 ]; then
    run sh tools/timecop/run_all.sh
fi

printf '\nrun_tests.sh: all wired gates passed\n'
