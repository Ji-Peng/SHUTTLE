#!/usr/bin/env bash
#
# run.sh - Build SHUTTLE/ref and SHUTTLE/avx2 kat_dump binaries for modes
# 128 and 256, run each against the same deterministic DRBG seed, and
# verify the resulting hex dumps are byte-identical.
#
# Exit 0 iff every (ref,avx2) pair diffs empty.
#
set -euo pipefail

cd "$(dirname "$0")"

NROUNDS=${NROUNDS:-10}
MODES=(128 256)
IMPLS=(ref avx2)

mkdir -p out

echo "Building kat_dump binaries..."
# Clean once up front, then build all targets into the same out/ dir.
make -s clean >/dev/null 2>&1 || true
for impl in "${IMPLS[@]}"; do
  for mode in "${MODES[@]}"; do
    make -s SHUTTLE_IMPL=$impl SHUTTLE_MODE=$mode all
  done
done

echo "Running kat_dump for every (impl, mode) against the same master seed..."
for impl in "${IMPLS[@]}"; do
  for mode in "${MODES[@]}"; do
    out_file="out/kat_${impl}${mode}.txt"
    ./out/kat_dump_${impl}${mode} "$NROUNDS" > "$out_file"
    echo "  wrote $out_file ($(wc -l < "$out_file") lines, $(wc -c < "$out_file") bytes)"
  done
done

echo "Comparing ref vs avx2 for each mode..."
fail=0
for mode in "${MODES[@]}"; do
  ref_file="out/kat_ref${mode}.txt"
  avx_file="out/kat_avx2${mode}.txt"
  # Strip the first 'algorithm=' line (contains algname, which differs
  # only in the internal naming but not in the cryptographic output we
  # care about).
  ref_body="out/kat_ref${mode}.body"
  avx_body="out/kat_avx2${mode}.body"
  tail -n +2 "$ref_file" > "$ref_body"
  tail -n +2 "$avx_file" > "$avx_body"

  if diff -q "$ref_body" "$avx_body" > /dev/null; then
    echo "  MODE=$mode: IDENTICAL"
  else
    echo "  MODE=$mode: DIFFER"
    diff "$ref_body" "$avx_body" | head -20
    fail=1
  fi
done

if [ $fail -eq 0 ]; then
  echo "ALL KATS IDENTICAL (modes: ${MODES[*]})"
  exit 0
else
  echo "KAT MISMATCH"
  exit 1
fi
