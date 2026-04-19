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
TEST_SIGN_BINS     = $(foreach m,$(MODES),$(OUTDIR)/test_sign$(m))
SPEED_SIGN_BINS    = $(foreach m,$(MODES),$(OUTDIR)/speed_sign$(m))
PROFILE_SIGN_BINS  = $(foreach m,$(MODES),$(OUTDIR)/profile_sign$(m))
LIB_BINS           = $(foreach m,$(MODES),$(OUTDIR)/libshuttle$(m)_ref.so)

all: $(TEST_SAMPLER_BINS) $(SPEED_SAMPLER_BINS) \
     $(TEST_NTT_BINS) $(TEST_SIGN_BINS) \
     $(SPEED_SIGN_BINS) $(PROFILE_SIGN_BINS) \
     $(LIB_BINS)

libs: $(LIB_BINS)

test-sampler: $(TEST_SAMPLER_BINS)
test-ntt:     $(TEST_NTT_BINS)
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
check: $(TEST_SAMPLER_BINS) $(TEST_NTT_BINS) $(TEST_SIGN_BINS)
	@set -e ; \
	for m in $(MODES) ; do \
	  echo "=== MODE=$$m: test_sampler ===" ; \
	  $(OUTDIR)/test_sampler$$m | tail -n 3 ; \
	  echo "=== MODE=$$m: test_ntt ===" ; \
	  $(OUTDIR)/test_ntt$$m | tail -n 5 ; \
	  echo "=== MODE=$$m: test_sign ===" ; \
	  $(OUTDIR)/test_sign$$m | tail -n 10 ; \
	done
	@echo "=== All mode check runs completed ==="

clean:
	rm -rf $(OUTDIR)
