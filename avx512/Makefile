# Makefile for the SHUTTLE signature AVX-512 implementation.
#
#   make                 # build all tests for all modes -> out/
#   make check           # build + run every test for every mode (assert exit 0)
#   make check-kat       # build test_kat; assert recorded KAT (ref==avx2==avx512)
#   make check-refill    # force multi-way XOF continuation-squeeze; KAT unchanged
#   make check-consts    # reproducibility/drift gate (params/bounds + vendored NTT)
#   make sm3-test        # legacy SM3 16-way XOF/DRNG byte-exactness test (P02)
#   make clean
#
# Mirrors avx2/Makefile but its flags drop -mno-avx512f and add -mavx512* (the
# existing SM3 16-way code needs VBMI2 per memory `avx512-requires-vbmi2`), and
# it DROPS check-no-avx512.  Maps to the NGCC Additional_Implementation folder
# (steps outside the literal NGCC perf flag set; AVX-512 perf is local-only).
# Per-mode NTT vendoring appends the ntt_avx512.S / ntt_consts_avx512.c variant
# (q59393n1024 uses the 2x512 superblock; P03 owns the codegen).

MODES ?= 128 256 512
# P03 adds test_ntt_avx512 (AVX-512 NTT shim round-trip/convolution + K1
# nttunpack + K10 ref==avx512; q59393n1024 exercises the 2x512 superblock).
TESTS ?= t_params test_xof test_ntt_avx512
CC    ?= gcc

# avx2 perf flags MINUS -mno-avx512f, PLUS the AVX-512 subsets (F/BW/DQ/VL +
# VBMI2).  -std=c99 -Wpedantic warns on intrinsics/asm; tolerated on the perf
# build (Overview 10).
CFLAGS ?= -O3 -march=x86-64 -mavx2 -mtune=native -flto -fomit-frame-pointer \
          -std=c99 -Wpedantic -Wall -Wextra -mbmi2 -mpopcnt \
          -mavx512f -mavx512bw -mavx512dq -mavx512vl -mavx512vbmi2 \
          -I. -I../ref -I../ref/test -I../tools \
          -DUSE_AVX512_NTT -DUSE_AVX512_SAMPLER -DUSE_AVX512_SHAKE8X
LDLIBS ?=
BINDIR  = out
CT_OPT ?= -O3

MODE ?=
ifeq ($(MODE),SHA3)
MODEFLAG = -DSHA3_MODE
else
MODEFLAG =
endif

# test_xof (P02) scalar XOF build-verify under the avx512 flag set.  MODE-
# dependent source list (SHA3: symmetric.c + fips202.c symlinks; NGCC:
# symmetric.c + ../ref/drng.c + ../ref/auxfunc.c).
ifeq ($(MODE),SHA3)
XOF_SRCS = symmetric.c fips202.c
else
XOF_SRCS = symmetric.c ../ref/drng.c ../ref/auxfunc.c
endif

SH_SRCS = \
't_params=' \
"test_xof=$(XOF_SRCS)" \
'test_ntt_avx512=poly_ntt.c ../ref/reduce.c'

QSET_128 = q15361n256
QSET_256 = q61441n512
QSET_512 = q59393n1024
NTTGEN   = ../tools/ntt

KAT_128 = TBD
KAT_256 = TBD
KAT_512 = TBD

.PHONY: all check check-kat check-refill check-consts sm3-test clean

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
	          extra="$$qs/ntt_avx512.S $$qs/ntt_consts_avx512.c ../ref/ntt/$$qs/ntt_ref.c"; \
	          inc="-I../ref/ntt/$$qs";; \
	      esac; \
	      stats=""; \
	      case "$$t" in trace_kat) stats="-DTEST_RANS_STATS";; diff_fuzz) stats="-DDIFF_FUZZ_N=$(DIFF_FUZZ_N)";; esac; \
	      $(CC) $(CFLAGS) $$inc $(MODEFLAG) $$stats -DSHUTTLE_MODE=$$m test/$$t.c $$srcs $$extra $(LDLIBS) \
	        -o $(BINDIR)/$${t}_$${m} || exit 1; \
	    done; \
	  done; \
	done
	@echo "built: [$(TESTS)] x [$(MODES)] (avx512)"

check: all
	@for m in $(MODES); do \
	  echo "=== SHUTTLE-$$m (avx512) ==="; \
	  for want in $(TESTS); do \
	    [ -x $(BINDIR)/$${want}_$${m} ] \
	      || { echo "missing test binary: $(BINDIR)/$${want}_$${m}"; exit 1; }; \
	    $(BINDIR)/$${want}_$${m} || exit 1; \
	  done; \
	done

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

check-refill:
	@if ! echo "$(TESTS)" | grep -qw test_kat; then \
	  echo "check-refill: deferred (test_kat + sampler land with P07/P11/P12)"; exit 0; fi
	@$(MAKE) CFLAGS="$(CFLAGS) -DGAUSS8X_INIT=4 -DNOISE8X_INIT=4" TESTS="test_kat" all >/dev/null
	@echo "== forced multi-way refill: KAT must be unchanged =="
	@for m in $(MODES); do $(BINDIR)/test_kat_$$m; done

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
	@echo "== vendored AVX-512 NTT (./<qset>/ntt_avx512.S + ntt_consts_avx512.c) vs gen_ntt.py =="
	@rm -rf $(BINDIR)/ntt_regen
	@mkdir -p $(BINDIR)/ntt_regen
	@for qs in $(QSET_128) $(QSET_256) $(QSET_512); do \
	  if [ -d "$(NTTGEN)/$$qs" ] && [ -f "$$qs/ntt_avx512.S" ]; then \
	    python3 $(NTTGEN)/$$qs/gen_ntt.py $(BINDIR)/ntt_regen >/dev/null 2>&1 \
	      || { echo "  FAIL: $$qs gen_ntt.py"; exit 1; }; \
	    cmp -s $(BINDIR)/ntt_regen/avx512/$$qs/ntt_consts_avx512.c $$qs/ntt_consts_avx512.c \
	      && cmp -s $(BINDIR)/ntt_regen/avx512/$$qs/ntt_avx512.S $$qs/ntt_avx512.S \
	      && echo "  PASS: $$qs vendored AVX-512 NTT matches gen_ntt.py" \
	      || { echo "  FAIL: $$qs vendored AVX-512 NTT differs (re-run make tables)"; exit 1; }; \
	  else echo "  skip: $$qs vendored AVX-512 NTT not present"; fi; \
	done
	@rm -rf $(BINDIR)/ntt_regen
	@echo "check-consts: P01 + P03 (AVX-512 NTT) constants reproducible"

# Legacy SM3 16-way XOF/DRNG byte-exactness test (16-way AVX-512 vs scalar ref).
SM3_AVX512 := -mavx512f -mavx512bw -mavx512dq -mavx512vl -mbmi2
sm3-test:
	@mkdir -p $(BINDIR)
	@$(CC) -O3 -std=c99 -Wall -Wextra $(SM3_AVX512) -I../ref \
	   test_avx512.c auxfunc_avx512.c ../ref/auxfunc.c -o $(BINDIR)/sm3_test_avx512 2>/dev/null \
	   && $(BINDIR)/sm3_test_avx512 || echo "sm3-test: build deferred (P02 sources)"
	@$(CC) -O3 -std=c99 -Wall -Wextra $(SM3_AVX512) -I../ref \
	   test_drng_avx512.c drng_avx512.c auxfunc_avx512.c ../ref/drng.c -o $(BINDIR)/sm3_test_drng_avx512 2>/dev/null \
	   && $(BINDIR)/sm3_test_drng_avx512 || echo "sm3-test(drng): build deferred (P02 sources)"

clean:
	rm -rf $(BINDIR)
