CC ?= gcc
CFLAGS_BASE = -Wall -Wextra -Wpedantic -Wmissing-prototypes -O3 -std=c99 -fPIC
LDFLAGS =
LDLIBS = -lm

# ============================================================
# Source sets
# ============================================================
SAMPLER_SOURCES = fips202.c symmetric-shake.c approx_exp.c \
                  sampler.c sampler_4x.c randombytes.c

SIGN_SOURCES    = fips202.c symmetric-shake.c approx_exp.c approx_log.c \
                  sampler.c sampler_4x.c rejsample.c \
                  reduce.c ntt.c poly.c polyvec.c rounding.c \
                  packing.c sign.c randombytes.c

NTT_SOURCES     = fips202.c randombytes.c reduce.c ntt.c poly.c rounding.c

MOD2Q_SOURCES   = fips202.c randombytes.c reduce.c ntt.c poly.c polyvec.c rounding.c

RANS_SOURCES    = shuttle_rans.c

POLY_Z1_SOURCES = fips202.c randombytes.c reduce.c ntt.c poly.c rounding.c shuttle_rans.c

SIG_V2_SOURCES  = fips202.c randombytes.c reduce.c ntt.c poly.c polyvec.c rounding.c shuttle_rans.c packing_v2.c

SIG_V3_SOURCES  = fips202.c randombytes.c reduce.c ntt.c poly.c polyvec.c rounding.c shuttle_rans.c packing_v3.c

SIGN_V2_SOURCES = fips202.c symmetric-shake.c approx_exp.c approx_log.c \
                  sampler.c sampler_4x.c rejsample.c \
                  reduce.c ntt.c poly.c polyvec.c rounding.c \
                  packing.c packing_v2.c shuttle_rans.c sign_v2.c randombytes.c

SIGN_V3_SOURCES = fips202.c symmetric-shake.c approx_exp.c approx_log.c \
                  sampler.c sampler_4x.c rejsample.c \
                  reduce.c ntt.c poly.c polyvec.c rounding.c \
                  packing.c packing_v3.c shuttle_rans.c sign_v3.c randombytes.c

TEST_COMMON     = test/cpucycles.c test/speed_print.c

OUTDIR = test/out

# ============================================================
# Multi-mode matrix
#   SHUTTLE_MODE in {128, 256}:
#     default `all` target builds both
#   SHUTTLE_MODE = 512:
#     guarded behind SHUTTLE_BUILD_512 (expected to #error at compile
#     time until the SHUTTLE-512 parameters are filled in)
# ============================================================
SHUTTLE_BUILD_512 ?= 0
ifeq ($(SHUTTLE_BUILD_512),1)
MODES = 128 256 512
else
MODES = 128 256
endif

# Per-mode compile flag.
CFLAGS_MODE = -DSHUTTLE_MODE=$(1)

# ============================================================
# Phony targets
# ============================================================
.PHONY: all clean check \
        test-sampler test-ntt test-sign \
        speed-sign profile-sign \
        libs

# ============================================================
# Target lists (generated per mode)
# ============================================================
TEST_SAMPLER_BINS  = $(foreach m,$(MODES),$(OUTDIR)/test_sampler$(m))
SPEED_SAMPLER_BINS = $(foreach m,$(MODES),$(OUTDIR)/speed_sampler$(m))
TEST_NTT_BINS      = $(foreach m,$(MODES),$(OUTDIR)/test_ntt$(m))
TEST_MOD2Q_BINS    = $(foreach m,$(MODES),$(OUTDIR)/test_mod2q$(m))
TEST_RANS_BINS     = $(foreach m,$(MODES),$(OUTDIR)/test_rans$(m))
TEST_POLY_Z1_BINS  = $(foreach m,$(MODES),$(OUTDIR)/test_poly_z1$(m))
TEST_SIG_V2_BINS   = $(foreach m,$(MODES),$(OUTDIR)/test_sig_v2$(m))
TEST_SIG_V3_BINS   = $(foreach m,$(MODES),$(OUTDIR)/test_sig_v3$(m))
TEST_SIGN_V2_BINS  = $(foreach m,$(MODES),$(OUTDIR)/test_sign_v2$(m))
TEST_SIGN_V3_BINS  = $(foreach m,$(MODES),$(OUTDIR)/test_sign_v3$(m))
TEST_SIGN_BINS     = $(foreach m,$(MODES),$(OUTDIR)/test_sign$(m))
SPEED_SIGN_BINS    = $(foreach m,$(MODES),$(OUTDIR)/speed_sign$(m))
PROFILE_SIGN_BINS  = $(foreach m,$(MODES),$(OUTDIR)/profile_sign$(m))
LIB_BINS           = $(foreach m,$(MODES),$(OUTDIR)/libshuttle$(m)_ref.so)

all: $(TEST_SAMPLER_BINS) $(SPEED_SAMPLER_BINS) \
     $(TEST_NTT_BINS) $(TEST_MOD2Q_BINS) $(TEST_RANS_BINS) \
     $(TEST_POLY_Z1_BINS) \
     $(TEST_SIG_V2_BINS) $(TEST_SIG_V3_BINS) \
     $(TEST_SIGN_V2_BINS) $(TEST_SIGN_V3_BINS) $(TEST_SIGN_BINS) \
     $(SPEED_SIGN_BINS) $(PROFILE_SIGN_BINS) \
     $(LIB_BINS)

libs: $(LIB_BINS)

test-sampler: $(TEST_SAMPLER_BINS)
test-ntt:     $(TEST_NTT_BINS)
test-mod2q:   $(TEST_MOD2Q_BINS)
test-rans:    $(TEST_RANS_BINS)
test-poly-z1: $(TEST_POLY_Z1_BINS)
test-sig-v2:  $(TEST_SIG_V2_BINS)
test-sig-v3:  $(TEST_SIG_V3_BINS)
test-sign-v2: $(TEST_SIGN_V2_BINS)
test-sign-v3: $(TEST_SIGN_V3_BINS)
test-sign:    $(TEST_SIGN_BINS)
speed-sign:   $(SPEED_SIGN_BINS)
profile-sign: $(PROFILE_SIGN_BINS)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# ============================================================
# Per-mode build template
# ============================================================
define MODE_template
$(OUTDIR)/test_sampler$(1): test/test_sampler.c $(SAMPLER_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_sampler.c $(SAMPLER_SOURCES) $(LDLIBS)

$(OUTDIR)/speed_sampler$(1): test/speed_sampler.c $(TEST_COMMON) \
                              $(SAMPLER_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/speed_sampler.c $(TEST_COMMON) $(SAMPLER_SOURCES) $(LDLIBS)

$(OUTDIR)/test_ntt$(1): test/test_ntt.c $(NTT_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_ntt.c $(NTT_SOURCES) $(LDLIBS)

$(OUTDIR)/test_mod2q$(1): test/test_mod2q.c $(MOD2Q_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_mod2q.c $(MOD2Q_SOURCES) $(LDLIBS)

$(OUTDIR)/test_rans$(1): test/test_rans.c $(RANS_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_rans.c $(RANS_SOURCES) $(LDLIBS)

$(OUTDIR)/test_poly_z1$(1): test/test_poly_z1.c $(POLY_Z1_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_poly_z1.c $(POLY_Z1_SOURCES) $(LDLIBS)

$(OUTDIR)/test_sig_v2$(1): test/test_sig_v2.c $(SIG_V2_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_sig_v2.c $(SIG_V2_SOURCES) $(LDLIBS)

$(OUTDIR)/test_sig_v3$(1): test/test_sig_v3.c $(SIG_V3_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_sig_v3.c $(SIG_V3_SOURCES) $(LDLIBS)

$(OUTDIR)/test_sign_v2$(1): test/test_sign_v2.c $(SIGN_V2_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_sign_v2.c $(SIGN_V2_SOURCES) $(LDLIBS)

$(OUTDIR)/test_sign_v3$(1): test/test_sign_v3.c $(SIGN_V3_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_sign_v3.c $(SIGN_V3_SOURCES) $(LDLIBS)

$(OUTDIR)/test_sign$(1): test/test_sign.c $(SIGN_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/test_sign.c $(SIGN_SOURCES) $(LDLIBS)

$(OUTDIR)/speed_sign$(1): test/speed_sign.c $(TEST_COMMON) $(SIGN_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/speed_sign.c $(TEST_COMMON) $(SIGN_SOURCES) $(LDLIBS)

$(OUTDIR)/profile_sign$(1): test/profile_sign.c $(TEST_COMMON) $(SIGN_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) $(LDFLAGS) \
	    -o $$@ test/profile_sign.c $(TEST_COMMON) $(SIGN_SOURCES) $(LDLIBS)

$(OUTDIR)/libshuttle$(1)_ref.so: $(SIGN_SOURCES) | $(OUTDIR)
	$$(CC) $(CFLAGS_BASE) $(call CFLAGS_MODE,$(1)) -shared \
	    -o $$@ $(SIGN_SOURCES) $(LDLIBS)
endef

$(foreach m,$(MODES),$(eval $(call MODE_template,$(m))))

# ============================================================
# make check: run correctness tests for all active modes
# ============================================================
check: $(TEST_SAMPLER_BINS) $(TEST_NTT_BINS) $(TEST_MOD2Q_BINS) \
       $(TEST_RANS_BINS) $(TEST_POLY_Z1_BINS) \
       $(TEST_SIG_V2_BINS) $(TEST_SIG_V3_BINS) \
       $(TEST_SIGN_V2_BINS) $(TEST_SIGN_V3_BINS) $(TEST_SIGN_BINS)
	@set -e ; \
	for m in $(MODES) ; do \
	  echo "=== MODE=$$m: test_sampler ===" ; \
	  $(OUTDIR)/test_sampler$$m | tail -n 3 ; \
	  echo "=== MODE=$$m: test_ntt ===" ; \
	  $(OUTDIR)/test_ntt$$m | tail -n 5 ; \
	  echo "=== MODE=$$m: test_mod2q ===" ; \
	  $(OUTDIR)/test_mod2q$$m | tail -n 12 ; \
	  echo "=== MODE=$$m: test_rans ===" ; \
	  $(OUTDIR)/test_rans$$m | tail -n 6 ; \
	  echo "=== MODE=$$m: test_poly_z1 ===" ; \
	  $(OUTDIR)/test_poly_z1$$m | tail -n 10 ; \
	  echo "=== MODE=$$m: test_sig_v2 ===" ; \
	  $(OUTDIR)/test_sig_v2$$m | tail -n 6 ; \
	  echo "=== MODE=$$m: test_sig_v3 ===" ; \
	  $(OUTDIR)/test_sig_v3$$m | tail -n 10 ; \
	  echo "=== MODE=$$m: test_sign_v2 ===" ; \
	  $(OUTDIR)/test_sign_v2$$m | tail -n 10 ; \
	  echo "=== MODE=$$m: test_sign_v3 ===" ; \
	  $(OUTDIR)/test_sign_v3$$m | tail -n 10 ; \
	  echo "=== MODE=$$m: test_sign ===" ; \
	  $(OUTDIR)/test_sign$$m | tail -n 10 ; \
	done
	@echo "=== All mode check runs completed ==="

clean:
	rm -rf $(OUTDIR)
