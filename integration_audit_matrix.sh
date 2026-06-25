#!/bin/sh
# Full integration audit gate: cross-backend differential fuzz plus
# the perf-review-checklist plumbing. Catches "functional-clean but backend-
# divergent on some input" (the ML-DSA AABBCC bug class).
#
# Builds + runs test/diff_fuzz.c on each backend x mode (DIFF_FUZZ_N cases),
# normalizes the CASE / PASS lines, and cmps every backend's normalized output
# against the first backend; any mismatch -> FAIL with a diff -u. Also checks
# the perf-review checklist files exist and the README links them.
set -eu
cd "$(dirname "$0")"
. ./tools/audit_lib.sh

mkdir -p test
DIFF_REPORT=${DIFF_FUZZ_REPORT:-test/diff_fuzz_matrix.txt}
CHECKLIST_REPORT=${PERF_REVIEW_REPORT:-test/perf_review_checklist_status.txt}
MATRIX_REPORT=${INTEGRATION_AUDIT_REPORT:-test/integration_audit_matrix.txt}
BACKENDS=${INTEGRATION_AUDIT_BACKENDS:-ref avx2 avx512}
MODES=${INTEGRATION_AUDIT_MODES:-128 256 512}
XOF=${INTEGRATION_AUDIT_XOF:-ngcc}
DIFF_FUZZ_N=${DIFF_FUZZ_N:-32}

overall=PASS
: > "$DIFF_REPORT"
{
    printf '=== SHUTTLE differential fuzz matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'xof: %s\n' "$XOF"
    printf 'DIFF_FUZZ_N: %s\n\n' "$DIFF_FUZZ_N"
} >> "$DIFF_REPORT"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT INT HUP TERM

for backend in $BACKENDS; do
    out="$tmpdir/$backend.out"
    {
        printf '## backend=%s target=diff_fuzz\n' "$backend"
        for mode in $MODES; do
            bin="out/diff_fuzz_${mode}"
            if ! audit_build "$backend" "$mode" "$XOF" test/diff_fuzz.c "$bin" \
                    -DDIFF_FUZZ_N="$DIFF_FUZZ_N" > "$tmpdir/build.$backend.$mode" 2>&1; then
                printf 'BACKEND %s mode %s BUILDFAIL\n' "$backend" "$mode"
                cat "$tmpdir/build.$backend.$mode"
                exit 1
            fi
            "$backend/$bin" || { printf 'BACKEND %s mode %s RUNFAIL\n' "$backend" "$mode"; exit 1; }
        done
        printf '\n'
    } > "$out" 2>&1 || overall=FAIL
    cat "$out" >> "$DIFF_REPORT"
    # Normalize: keep CASE / PASS diff_fuzz lines (backend-independent transcript).
    awk '/^(CASE|PASS diff_fuzz)/ { print }' "$out" > "$tmpdir/$backend.norm"
done

ref_backend=
for backend in $BACKENDS; do ref_backend=$backend; break; done
if [ -n "$ref_backend" ]; then
    for backend in $BACKENDS; do
        if ! cmp -s "$tmpdir/$ref_backend.norm" "$tmpdir/$backend.norm"; then
            printf 'DIFF-CONSISTENCY backend=%s vs=%s FAIL\n' "$backend" "$ref_backend" >> "$DIFF_REPORT"
            diff -u "$tmpdir/$ref_backend.norm" "$tmpdir/$backend.norm" >> "$DIFF_REPORT" || true
            overall=FAIL
        else
            lines=$(wc -l < "$tmpdir/$backend.norm" | tr -d ' ')
            printf 'DIFF-CONSISTENCY backend=%s vs=%s PASS normalized_lines=%s\n' "$backend" "$ref_backend" "$lines" >> "$DIFF_REPORT"
        fi
    done
fi
printf 'OVERALL: %s\n' "$overall" >> "$DIFF_REPORT"

check_status=PASS
{
    printf '=== SHUTTLE performance-review checklist status ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    for path in docs/design-notes/PERF_PATCH_REVIEW.md .github/pull_request_template.md; do
        if [ -f "$path" ]; then
            printf 'PATH %s PASS\n' "$path"
        else
            printf 'PATH %s FAIL missing\n' "$path"
            check_status=FAIL
        fi
    done
    if [ -f README.md ] && grep -q 'PERF_PATCH_REVIEW.md' README.md; then
        printf 'README_LINK PERF_PATCH_REVIEW.md PASS\n'
    else
        printf 'README_LINK PERF_PATCH_REVIEW.md FAIL\n'
        check_status=FAIL
    fi
    printf 'OVERALL: %s\n' "$check_status"
} > "$CHECKLIST_REPORT"
[ "$check_status" = PASS ] || overall=FAIL

{
    printf '=== SHUTTLE integration audit matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'DIFF_FUZZ_N: %s\n' "$DIFF_FUZZ_N"
    printf 'diff_report: %s\n' "$DIFF_REPORT"
    printf 'checklist_report: %s\n\n' "$CHECKLIST_REPORT"
    printf '## %s\n' "$DIFF_REPORT"
    tail -n 16 "$DIFF_REPORT"
    printf '\n## %s\n' "$CHECKLIST_REPORT"
    cat "$CHECKLIST_REPORT"
    printf '\nOVERALL: %s\n' "$overall"
} > "$MATRIX_REPORT"
cat "$MATRIX_REPORT"
[ "$overall" = PASS ]
