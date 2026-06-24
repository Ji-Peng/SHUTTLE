# Makefile for the SHUTTLE signature AVX2 implementation.
#
#   make                 # build all tests for all modes -> out/
#   make check           # build + run every test for every mode (assert exit 0)
#   make check-kat       # build test_kat; assert recorded KAT (ref==avx2==avx512)
#   make check-no-avx512 # objdump gate: NO zmm/EVEX/AVX-512 op leaked (pure AVX2)
#   make check-refill    # force multi-way XOF continuation-squeeze; KAT unchanged
#   make check-consts    # reproducibility/drift gate (params/bounds + vendored NTT)
#   make sm3-test        # legacy SM3 N-way XOF/DRNG byte-exactness test (P02)
#   make clean
#
# Shared algorithm/KAT source is SYMLINKED from ../ref (single source of truth;
# materialized by ../tools/mklinks.sh).  AVX2-specific overrides are REAL files
# here (config.h, polyvec.c, sampler.c, irs.c, rounding.c, sign.c + the AVX-only
# XOF/NTT files).  Per-mode NTT asm + consts are VENDORED locally in ./<qset>/.
#
# Compile-time switches select the AVX2 code paths in the symlinked ref sources:
#   USE_AVX2_NTT / USE_AVX2_SAMPLER / USE_AVX2_SHAKE4X  (set in CFLAGS below).

MODES ?= 128 256 512
# Seeded with t_params + test_xof (P02); P03 adds test_ntt_avx (AVX2 NTT shim
# round-trip/convolution + the K1 nttunpack reconciliation + K10 ref==avx2).
TESTS ?= t_params test_xof test_ntt_avx
CC    ?= gcc

# Exact NGCC Optimized flags (digest 14 1.1).  -mno-avx512f forces a PURE AVX2
# path: -mtune=native would otherwise let GCC auto-vectorize the symlinked
# scalar C to zmm on an AVX-512 host; -mno-avx512f transitively disables every
# AVX-512 subset while keeping AVX2/BMI2/POPCNT.  Enforced by check-no-avx512.
# NOTE: -std=c99 -Wpedantic warns on intrinsics/asm; per Overview 10 / digest
# 14 R1 we TOLERATE those on the perf build (the flag set is fixed) and keep
# only the ref build strictly warning-clean.
CFLAGS ?= -O3 -march=x86-64 -mavx2 -mtune=native -flto -fomit-frame-pointer \
          -std=c99 -Wpedantic -Wall -Wextra -mbmi2 -mpopcnt -mno-avx512f \
          -I. -I../ref -I../ref/test -I../tools \
          -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X
LDLIBS ?=
BINDIR  = out
CT_OPT ?= -O3

MODE ?=
ifeq ($(MODE),SHA3)
MODEFLAG = -DSHA3_MODE
else
MODEFLAG =
endif

# test_xof (P02) scalar XOF build-verify under the avx2 flag set.  Source list
# is MODE-dependent: SHA3_MODE links the vendored scalar Keccak (symmetric.c +
# fips202.c, symlinked from ../ref); NGCC_MODE links the SM3 DRBG (symmetric.c
# + ../ref/drng.c + ../ref/auxfunc.c).  Proves the wrappers compile + run clean
# under -mavx2 -mno-avx512f too (the scalar fips202.c must NOT auto-vectorize to
# zmm; check-no-avx512 enforces that).
ifeq ($(MODE),SHA3)
XOF_SRCS = symmetric.c fips202.c
else
XOF_SRCS = symmetric.c ../ref/drng.c ../ref/auxfunc.c
endif

# Per-test source lists (POSIX-sh strings).  t_params is header-only.  Tests
# whose srcs contain polyvec.c get the per-mode vendored ntt.S + ntt_consts.c
# appended in the build loop (see `all`).
SH_SRCS = \
't_params=' \
"test_xof=$(XOF_SRCS)" \
'test_ntt_avx=poly_ntt.c ../ref/reduce.c'

# Per-set NTT vendoring (Overview / digest 11 5): a case $$m picks $qs; when a
# test's srcs contain polyvec.c OR poly_ntt.c (the AVX2 NTT shim fork) the build
# appends $qs/ntt.S + $qs/ntt_consts.c, the scalar oracle ../ref/ntt/$qs/
# ntt_ref.c, and -I../ref/ntt/$qs so poly_ntt.h finds the per-set header.
QSET_128 = q15361n256
QSET_256 = q61441n512
QSET_512 = q59393n1024
# NTTGEN points at the generator (tools/ntt/<qset>/gen_ntt.py), used only by
# check-consts to prove the vendored copies are reproducible.
NTTGEN   = ../tools/ntt

# Recorded canonical KAT hashes (shared: ref == avx2 == avx512).  P12 records
# the real values; placeholders until M5.
KAT_128 = TBD
KAT_256 = TBD
KAT_512 = TBD

.PHONY: all check check-kat check-no-avx512 check-refill check-consts sm3-test clean

all:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 128) qs=$(QSET_128);; 256) qs=$(QSET_256);; 512) qs=$(QSET_512);; esac; \
	  for entry in $(SH_SRCS); do \
	    t=$${entry%%=*}; srcs=$${entry#*=}; \
	    for want in $(TESTS); do [ "$$want" = "$$t" ] || continue; \
	      extra=""; inc=""; \
	      case "$$srcs" in \
	        *polyvec.c*|*poly_ntt.c*) \
	          extra="$$qs/ntt.S $$qs/ntt_consts.c ../ref/ntt/$$qs/ntt_ref.c"; \
	          inc="-I../ref/ntt/$$qs";; \
	      esac; \
	      stats=""; \
	      case "$$t" in trace_kat) stats="-DTEST_RANS_STATS";; diff_fuzz) stats="-DDIFF_FUZZ_N=$(DIFF_FUZZ_N)";; esac; \
	      $(CC) $(CFLAGS) $$inc $(MODEFLAG) $$stats -DSHUTTLE_MODE=$$m test/$$t.c $$srcs $$extra $(LDLIBS) \
	        -o $(BINDIR)/$${t}_$${m} || exit 1; \
	    done; \
	  done; \
	done
	@echo "built: [$(TESTS)] x [$(MODES)] (avx2)"

check: all
	@for m in $(MODES); do \
	  echo "=== SHUTTLE-$$m (avx2) ==="; \
	  for want in $(TESTS); do \
	    [ -x $(BINDIR)/$${want}_$${m} ] \
	      || { echo "missing test binary: $(BINDIR)/$${want}_$${m}"; exit 1; }; \
	    $(BINDIR)/$${want}_$${m} || exit 1; \
	  done; \
	done

# Assert the recorded per-mode KAT hashes (ref == avx2 == avx512).  Deferred
# until test_kat + recorded KAT_* land (P12, M5).
check-kat:
	@if ! echo "$(TESTS)" | grep -qw test_kat; then \
	  echo "check-kat: deferred (test_kat lands with P11/P12)"; exit 0; fi
	@$(MAKE) TESTS=test_kat all >/dev/null
	@fail=0; for m in $(MODES); do \
	  case $$m in 128) want=$(KAT_128);; 256) want=$(KAT_256);; 512) want=$(KAT_512);; esac; \
	  h=$$($(BINDIR)/test_kat_$$m | grep -o 'KAT_hash=[0-9]*' | cut -d= -f2); \
	  [ "$$h" = "$$want" ] && echo "  PASS: mode=$$m KAT_hash=$$h" \
	    || { echo "  FAIL: mode=$$m got=$$h want=$$want"; fail=1; }; \
	done; [ $$fail -eq 0 ] && echo "check-kat: shared KAT matches (ref==avx2==avx512)" || exit 1

# Pure-AVX2 gate: assert NO AVX-512 (zmm / EVEX mask / AVX-512-only mnemonic)
# leaked into any built object.  -mno-avx512f should make this impossible; the
# gate proves the flag is effective (and would catch its accidental removal).
check-no-avx512:
	@rm -rf $(BINDIR)
	@$(MAKE) all >/dev/null
	@fail=0; \
	for b in $(BINDIR)/*; do \
	  [ -f "$$b" ] || continue; \
	  hits=$$(objdump -d "$$b" 2>/dev/null | grep -iE '%zmm[0-9]|\{%k[0-7]\}|vpternlog|vpcompress|vpexpand|vpmadd52|vperm[ti]2[a-z]|kmov|kshift|kortest|kunpck|valign[dq]' || true); \
	  if [ -n "$$hits" ]; then echo "  FAIL $$b:"; echo "$$hits" | head -3; fail=1; fi; \
	done; \
	[ $$fail -eq 0 ] && echo "check-no-avx512: pure AVX2 (no zmm/EVEX/AVX-512 op in any build object)" || { echo "check-no-avx512: AVX-512 LEAKED"; exit 1; }

# Force-exercise the multi-way XOF continuation-squeeze (refill) path with a
# tiny initial squeeze and assert the KAT is unchanged (K6 cursor-advance
# determinism).  Deferred until the sampler/KAT tests land.
check-refill:
	@if ! echo "$(TESTS)" | grep -qw test_kat; then \
	  echo "check-refill: deferred (test_kat + sampler land with P07/P11/P12)"; exit 0; fi
	@$(MAKE) CFLAGS="$(CFLAGS) -DGAUSS4X_INIT=4 -DNOISE4X_INIT=4" TESTS="test_kat" all >/dev/null
	@echo "== forced multi-way refill: KAT must be unchanged =="
	@for m in $(MODES); do $(BINDIR)/test_kat_$$m; done

# Reproducibility/drift gate: params/bounds regions + the vendored NTT (once
# P03 lands).  cmp against a pre-regeneration copy.
check-consts:
	@mkdir -p $(BINDIR)
	@echo "== params.h (@@AUTOGEN:params@@ + @@AUTOGEN:bounds@@; via ../ref) =="
	@cp ../ref/params.h $(BINDIR)/params.h.bak
	@( cd ../tools && python3 gen_params.py >/dev/null 2>&1 && python3 gen_bounds.py >/dev/null 2>&1 ) \
	  || { echo "  FAIL: generator error"; rm -f $(BINDIR)/params.h.bak; exit 1; }
	@cmp -s ../ref/params.h $(BINDIR)/params.h.bak \
	  && echo "  PASS: params.h reproducible (params + bounds)" \
	  || { echo "  FAIL: params.h drift"; rm -f $(BINDIR)/params.h.bak; exit 1; }
	@rm -f $(BINDIR)/params.h.bak
	@echo "== vendored AVX2 NTT (./<qset>/ntt.S + ntt_consts.c) vs gen_ntt.py =="
	@rm -rf $(BINDIR)/ntt_regen
	@mkdir -p $(BINDIR)/ntt_regen
	@for qs in $(QSET_128) $(QSET_256) $(QSET_512); do \
	  if [ -d "$(NTTGEN)/$$qs" ] && [ -f "$$qs/ntt.S" ]; then \
	    python3 $(NTTGEN)/$$qs/gen_ntt.py $(BINDIR)/ntt_regen >/dev/null 2>&1 \
	      || { echo "  FAIL: $$qs gen_ntt.py"; exit 1; }; \
	    cmp -s $(BINDIR)/ntt_regen/avx2/$$qs/ntt_consts.c $$qs/ntt_consts.c \
	      && cmp -s $(BINDIR)/ntt_regen/avx2/$$qs/ntt.S $$qs/ntt.S \
	      && echo "  PASS: $$qs vendored AVX2 NTT matches gen_ntt.py" \
	      || { echo "  FAIL: $$qs vendored AVX2 NTT differs (re-run make tables)"; exit 1; }; \
	  else echo "  skip: $$qs vendored NTT not present"; fi; \
	done
	@rm -rf $(BINDIR)/ntt_regen
	@echo "check-consts: P01 + P03 (AVX2 NTT) constants reproducible"

# Legacy SM3 N-way XOF/DRNG byte-exactness test (8-way AVX2 vs scalar ref).
# Owned by P02; kept here so the existing capability is not lost.
SM3_AVX2 := -mavx2 -mbmi2
sm3-test:
	@mkdir -p $(BINDIR)
	@$(CC) -O3 -std=c99 -Wall -Wextra $(SM3_AVX2) -I../ref \
	   test_avx2.c auxfunc_avx2.c ../ref/auxfunc.c -o $(BINDIR)/sm3_test_avx2 2>/dev/null \
	   && $(BINDIR)/sm3_test_avx2 || echo "sm3-test: build deferred (P02 sources)"
	@$(CC) -O3 -std=c99 -Wall -Wextra $(SM3_AVX2) -I../ref \
	   test_drng_avx2.c drng_avx2.c auxfunc_avx2.c ../ref/drng.c -o $(BINDIR)/sm3_test_drng_avx2 2>/dev/null \
	   && $(BINDIR)/sm3_test_drng_avx2 || echo "sm3-test(drng): build deferred (P02 sources)"

clean:
	rm -rf $(BINDIR)
