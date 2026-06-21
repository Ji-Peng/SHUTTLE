# 8-way AVX2 SM3 XOF build.
#
# The reference object (../ref/auxfunc.c) is built without AVX2 flags so the
# benchmark compares a portable scalar baseline against the 8-way code.

CC      ?= gcc
CFLAGS  ?= -O3 -std=c99 -Wall -Wextra
AVX2    := -mavx2 -mbmi2
REFDIR  := ../ref

.PHONY: all test clean check-const

all: test_avx2

# vectorised translation unit
auxfunc_avx2.o: auxfunc_avx2.c auxfunc_avx2.h sm3_const.h
	$(CC) $(CFLAGS) $(AVX2) -c $< -o $@

# scalar reference (ground truth + baseline); no AVX2
auxfunc_ref.o: $(REFDIR)/auxfunc.c $(REFDIR)/auxfunc.h
	$(CC) $(CFLAGS) -I$(REFDIR) -c $< -o $@

test_avx2.o: test_avx2.c auxfunc_avx2.h $(REFDIR)/auxfunc.h
	$(CC) $(CFLAGS) $(AVX2) -I$(REFDIR) -c $< -o $@

test_avx2: test_avx2.o auxfunc_avx2.o auxfunc_ref.o
	$(CC) $(CFLAGS) $^ -o $@

test: test_avx2
	./test_avx2

# audit the embedded SM3 constants against the spec-derived generator
check-const:
	python3 gen_sm3_const.py --check

clean:
	rm -f *.o test_avx2
