CC = gcc
# Discrete Gaussian sampler mode (SHUTTLE_SIGMA ∈ {128, 101, 149}).
# Default 101 = SHUTTLE-128 parameter set. Override: `make SHUTTLE_SIGMA=149`.
SHUTTLE_SIGMA ?= 101
CFLAGS = -Wall -Wextra -Wpedantic -Wmissing-prototypes -O3 -std=c99 \
         -DSHUTTLE_SIGMA=$(SHUTTLE_SIGMA)
LDFLAGS =
LDLIBS = -lm

SOURCES = fips202.c symmetric-shake.c approx_exp.c sampler.c sampler_4x.c randombytes.c
HEADERS = fips202.h symmetric.h approx_exp.h sampler.h params.h randombytes.h
TEST_COMMON = test/cpucycles.c test/speed_print.c

OUTDIR = test/out

.PHONY: all clean

all: $(OUTDIR)/test_sampler $(OUTDIR)/speed_sampler

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(OUTDIR)/test_sampler: test/test_sampler.c $(SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/test_sampler.c $(SOURCES) $(LDLIBS)

$(OUTDIR)/speed_sampler: test/speed_sampler.c $(TEST_COMMON) $(SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/speed_sampler.c $(TEST_COMMON) $(SOURCES) $(LDLIBS)

clean:
	rm -rf $(OUTDIR)
