# SHUTTLE-NGCC

SHUTTLE is a lattice-based digital signature scheme built on Iterated Rejection Sampling (IRS). This tree holds three byte-exact backends -- `ref/` (pure C99 scalar reference), `avx2/`, and `avx512/` -- across three parameter sets, SHUTTLE-128 / -256 / -512, each selected at build time via `-DSHUTTLE_MODE=128|256|512`. The default symmetric backend is the NGCC SM3 Hash-DRBG; `-DSHA3_MODE` selects the SHAKE backend.

## Build & test

```sh
make -C ref check          # build + run every test for every mode
make -C ref check-kat      # recorded KAT hashes (ref == avx2 == avx512)
make -C ref check-consts   # reproducible-constants drift gate
```

## Security & constant-time audit (P13)

The whole audit harness is the gate-of-gates `./run_tests.sh` (functional + KAT + reproducible-constants + the audit matrices, fail-fast). Two defense families:

- **KyberSlash (variable-latency / secret division)** -- owned by the mandatory static machine-code scan `./ct_scan_matrix.sh` (-> `tools/ct_scan.py`), plus statistical `./run_dudect_all.sh` and the opt-in dynamic taint `TIMECOP=1 ./run_tests.sh`.
- **ML-DSA-bug class (correctness/leakage bugs that pass functional tests)** -- owned by the independent Python reference + deterministic KAT + cross-backend byte equality + `./integration_audit_matrix.sh` (diff-fuzz) + `./security_audit_matrix.sh` (fault injection + negative parser corpus + trace-KAT).

Key audit documents:

- `SECRET_PUBLIC_AUDIT.md` -- the secret/public/public-length/declassified classification baseline; every `ct_scan.py` allowlist exemption cites a section here, and the K1-K15 -> gate traceability table lives at the end.
- `docs/design-notes/PERF_PATCH_REVIEW.md` -- the six-question checklist for any perf/sampler/packing/parser/constant patch (embedded in `.github/pull_request_template.md`).
- `docs/design-notes/TIMECOP.md` -- the patched-Valgrind variable-latency methodology (why avx512 is SKIP, the `--undef-value-errors` / `Memcheck:Cond`-vs-`Value` pitfalls).

Run before approving any change: `./run_tests.sh` (hard gate) and `./integration_audit_matrix.sh` (integration). See `docs/design-notes/PERF_PATCH_REVIEW.md`.
