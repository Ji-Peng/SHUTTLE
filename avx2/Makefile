CC = gcc
AS = gcc
CFLAGS = -Wall -Wextra -Wpedantic -Wmissing-prototypes -O3 -std=c99 \
         -mavx2 -mbmi2 -mpopcnt -march=native -mtune=native
LDFLAGS =
LDLIBS = -lm

SOURCES = fips202.c fips202x4.c symmetric-shake.c approx_exp.c sampler.c sampler_4x.c randombytes.c
ASM_SOURCES = f1600x4.S
HEADERS = fips202.h fips202x4.h symmetric.h approx_exp.h sampler.h params.h randombytes.h
TEST_COMMON = test/cpucycles.c test/speed_print.c

OUTDIR = test/out

.PHONY: all clean

all: $(OUTDIR)/test_sampler $(OUTDIR)/speed_sampler

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(OUTDIR)/test_sampler: test/test_sampler.c $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/test_sampler.c $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

$(OUTDIR)/speed_sampler: test/speed_sampler.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(HEADERS) | $(OUTDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test/speed_sampler.c $(TEST_COMMON) $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

clean:
	rm -rf $(OUTDIR)
