CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -Wmissing-prototypes -O3 -std=c99
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
