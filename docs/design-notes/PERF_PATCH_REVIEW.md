# Performance Patch Review Checklist (SHUTTLE-NGCC)

Use this checklist for any optimization, vectorization, approximation, packing, parser, sampler, or generated-constant patch that can affect the signing or verification implementation. A reviewer should be able to answer each item from the pull request description, the linked reports, or a short code comment near the changed invariant.

## Required Statements

- **Distribution.** State whether the patch changes any sampled distribution, rejection distribution, acceptance bound, IRS transition, or retry distribution; if it does, link the design note and the new KAT evidence, and re-run the Python reference (`tools/pyref/`).
- **Randomness consumption.** State whether the patch changes `randombytes()` / SM3-DRBG / XOF call count, byte count, ordering, or retry-dependent consumption; if it does, update the trace-KAT report. (SHUTTLE's IRS consumes a fixed `18*tau` bytes per attempt — any change there is load-bearing.)
- **Parser rejection conditions.** State whether the patch changes signature, key, rANS canonical-decode, padding, length, or the non-power-of-2 hint range (`>= H_h` must REJECT) rejection conditions; if it does, add a negative parser case to `test/parser_negative.c`.
- **Machine-code CT envelope.** State whether the patch changes secret-dependent arithmetic, memory access, branches, division, modulo, gather/scatter, floating-point, or table lookups; if it does, link the `./ct_scan_matrix.sh` report and justify any public-data-only allowlist.
- **Byte-for-byte KAT status.** State whether the default build remains byte-for-byte identical to the recorded KAT hashes (`ref == avx2 == avx512`, both NGCC and SHA3 MODEs); if not, explain why the change is intentional and update ALL backend KAT records together.
- **Public-data-only allowlist.** State whether new divisions, modulo operations, gathers, variable shifts, early exits, or table lookups operate only on public data; cite `SECRET_PUBLIC_AUDIT.md` for the classification.

## Reviewer Gate

Before approval, run or inspect the relevant retained evidence. The default hard-gate command is `./run_tests.sh`, and the integration-audit command is `./integration_audit_matrix.sh`. Optional timing evidence for deterministic components is `./run_dudect_components.sh`; treat `WARN` as a request for a quieter rerun or a stronger static argument, not as an automatic rejection. The whole-sign statistical timing gate is `./run_dudect_all.sh` (`DUDECT_N=200000`), recorded in `test/dudect_all.txt`. The optional dynamic taint cross-check is `TIMECOP=1 ./run_tests.sh` (see `TIMECOP.md`).
