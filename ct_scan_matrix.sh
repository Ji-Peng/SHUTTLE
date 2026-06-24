#!/bin/sh
# Full machine-code constant-time scan matrix (P13 CT-1, MANDATORY HARD GATE).
# Thin wrapper over tools/ct_scan.py. Scans every SHUTTLE secret-handling object
# across {ref,avx2,avx512} x {128,256,512} x {gcc,clang} x {-O3,-Os} x {ngcc,sha3}
# for forbidden CT mnemonics (div/idiv/gather/scatter/sqrt/cvt...), allowlisting
# only documented public-length divisions / public SHAKE lane gathers.
set -eu
cd "$(dirname "$0")"
REPORT="${CT_SCAN_REPORT:-test/ct_scan_matrix.txt}"
exec python3 tools/ct_scan.py \
  --backends "${CT_BACKENDS:-all}" \
  --modes "${CT_MODES:-all}" \
  --ccs "${CT_CCS:-gcc clang}" \
  --opts="${CT_OPTS:--O3 -Os}" \
  --xofs "${CT_XOFS:-ngcc sha3}" \
  --include-optional \
  --report "$REPORT"
