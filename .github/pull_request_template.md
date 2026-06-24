## Change Summary

Describe the patch and the intended performance or maintainability effect.

## Performance Patch Review

If this PR touches signing, verification, sampling, packing, parsing, generated constants, vector code, or approximation code, fill out `docs/design-notes/PERF_PATCH_REVIEW.md` and summarize the answers here.

- Distribution changed: yes/no, with evidence.
- Randomness consumption changed (SM3-DRBG / XOF count/bytes/ordering/retry): yes/no, with trace evidence.
- Parser rejection conditions changed (rANS K15 / hint-range K14 / length / padding): yes/no, with negative-parser evidence.
- Machine-code CT envelope changed: yes/no, with `./ct_scan_matrix.sh` evidence.
- Default byte-for-byte KAT unchanged (ref==avx2==avx512, NGCC + SHA3): yes/no, with backend KAT evidence.
- New public-data-only allowlist needed: yes/no, citing `SECRET_PUBLIC_AUDIT.md`.

## Commands Run

- `./run_tests.sh`
- `./integration_audit_matrix.sh`
- Optional: `./run_dudect_components.sh`, `TIMECOP=1 ./run_tests.sh`, `./run_dudect_all.sh`
