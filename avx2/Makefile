CC = gcc
AS = gcc
# Discrete Gaussian sampler mode (SHUTTLE_SIGMA ∈ {128, 101, 149}).
# Default 101 = SHUTTLE-128 parameter set. Override: `make SHUTTLE_SIGMA=149`.
SHUTTLE_SIGMA ?= 101
CFLAGS = -Wall -Wextra -Wpedantic -Wmissing-prototypes -O3 -std=c99 \
         -mavx2 -mbmi2 -mpopcnt -march=native -mtune=native \
         -DSHUTTLE_SIGMA=$(SHUTTLE_SIGMA)
LDFLAGS =
LDLIBS = -lm

# Source files: existing AVX2 sampler/shake + copied ref NTT/poly/sign
SOURCES = fips202.c fips202x4.c symmetric-shake.c approx_exp.c \
          sampler.c sampler_4x.c randombytes.c \
          reduce.c ntt.c poly.c polyvec.c rounding.c packing.c \
          sign.c rejsample.c approx_log.c
ASM_SOURCES = f1600x4.S

HEADERS = fips202.h fips202x4.h symmetric.h approx_exp.h sampler.h \
          params.h config.h randombytes.h \
          reduce.h ntt.h poly.h polyvec.h rounding.h packing.h \
          sign.h rejsample.h approx_log.h api.h

TEST_COMMON = test/cpucycles.c test/speed_print.c

OUTDIR = test/out

.PHONY: all clean

all: $(OUTDIR)/test_sampler $(OUTDIR)/speed_sampler \
     $(OUTDIR)/test_sign $(OUTDIR)/speed_sign \
     $(OUTDIR)/profile_sign

$(OUTDIR):
	mkdir -p $(OUTDIR)

# Existing sampler tests
$(OUTDIR)/test_sampler: test/test_sampler.c $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/test_sampler.c $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

$(OUTDIR)/speed_sampler: test/speed_sampler.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/speed_sampler.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

# New sign tests
$(OUTDIR)/test_sign: test/test_sign.c $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/test_sign.c $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

$(OUTDIR)/speed_sign: test/speed_sign.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/speed_sign.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

$(OUTDIR)/profile_sign: test/profile_sign.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/profile_sign.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

clean:
	rm -rf $(OUTDIR)
