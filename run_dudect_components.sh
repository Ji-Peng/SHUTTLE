#!/bin/sh
# Component-level dudect smoke (P13-T6 / CT-5, SMOKE, non-gating). Builds
# tools/dudect/dudect_components.c per backend x mode and runs it with
# DUDECT_COMPONENT_N samples. The isochronous primitive cdt_scan96 /
# sampler_sigma2 (K11 branchless full-table scan) is the gating probe;
# approx_exp / approx_log are provably branchless integer kernels (ct_scan owns
# their CT property) so a high |t| on those few-cycle kernels is uarch noise,
# not a leak. WARN is allowed; this never hard-gates run_tests.sh.
set -eu
cd "$(dirname "$0")"
. ./tools/audit_lib.sh

mkdir -p test
REPORT=${DUDECT_COMPONENTS_REPORT:-test/dudect_components_matrix.txt}
BACKENDS=${DUDECT_COMPONENTS_BACKENDS:-ref avx2 avx512}
MODES=${DUDECT_COMPONENTS_MODES:-128 256 512}
DUDECT_COMPONENT_N=${DUDECT_COMPONENT_N:-512}
XOF=${DUDECT_COMPONENTS_XOF:-ngcc}

status=PASS
{
    printf '=== SHUTTLE optional dudect components matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'DUDECT_COMPONENT_N: %s\n\n' "$DUDECT_COMPONENT_N"
} > "$REPORT"

# The components harness only needs the sampler/approx sources, not the whole
# sign spine -- build it with a trimmed source list for speed.
build_components() {
    backend=$1; mode=$2
    cc=${CC:-gcc}
    cflags=$(audit_cflags "$backend" "$mode" "$XOF")
    srcs="sampler.c reduce.c approx_exp.c approx_log.c sampler_u.c symmetric.c"
    if [ "$XOF" = sha3 ]; then srcs="$srcs fips202.c drng.c auxfunc.c"; else srcs="$srcs drng.c auxfunc.c"; fi
    ( cd "$backend" && mkdir -p out && \
      # shellcheck disable=SC2086
      $cc $cflags -DDUDECT_COMPONENT_N="$DUDECT_COMPONENT_N" \
        "$AUDIT_ROOT/tools/dudect/dudect_components.c" $srcs -lm \
        -o "out/dudect_components_${mode}" )
}

for backend in $BACKENDS; do
    for mode in $MODES; do
        {
            printf '## backend=%s mode=%s target=dudect_components\n' "$backend" "$mode"
            if ! build_components "$backend" "$mode"; then
                printf 'BUILDFAIL backend=%s mode=%s\n' "$backend" "$mode"
                exit 1
            fi
            "$backend/out/dudect_components_${mode}"
            printf '\n'
        } >> "$REPORT" 2>&1 || status=FAIL
    done
done
# WARN if a probe crossed the |t| threshold (smoke; non-gating).
if grep -q 'WARN dudect_components' "$REPORT"; then
    [ "$status" = FAIL ] || status=WARN
fi
printf 'OVERALL: %s\n' "$status" >> "$REPORT"
cat "$REPORT"
[ "$status" != FAIL ]
