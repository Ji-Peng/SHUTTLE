# 16-way AVX512 SM3 XOF build.

CC      ?= gcc
CFLAGS  ?= -O3 -std=c99 -Wall -Wextra
AVX512  := -mavx512f -mavx512bw -mavx512dq -mavx512vl -mbmi2
REFDIR  := ../ref

.PHONY: all test clean check-const

all: test_avx512

auxfunc_avx512.o: auxfunc_avx512.c auxfunc_avx512.h sm3_const.h
	$(CC) $(CFLAGS) $(AVX512) -c $< -o $@

auxfunc_ref.o: $(REFDIR)/auxfunc.c $(REFDIR)/auxfunc.h
	$(CC) $(CFLAGS) -I$(REFDIR) -c $< -o $@

test_avx512.o: test_avx512.c auxfunc_avx512.h $(REFDIR)/auxfunc.h
	$(CC) $(CFLAGS) $(AVX512) -I$(REFDIR) -c $< -o $@

test_avx512: test_avx512.o auxfunc_avx512.o auxfunc_ref.o
	$(CC) $(CFLAGS) $^ -o $@

test: test_avx512
	./test_avx512

check-const:
	python3 gen_sm3_const.py --check

clean:
	rm -f *.o test_avx512
