CC ?= gcc

SRC_DIR = ..

CFLAGS_BASE = -Wall -Wextra -Wpedantic -O2 -std=c99 \
              -I. -I$(SRC_DIR)
LDLIBS = -lm

# Same source set as the ref/ test-sign target, minus randombytes.c
# (we link against our deterministic rng.c instead).
SOURCES = $(SRC_DIR)/fips202.c $(SRC_DIR)/symmetric-shake.c \
          $(SRC_DIR)/approx_exp.c $(SRC_DIR)/approx_log.c \
          $(SRC_DIR)/sampler.c $(SRC_DIR)/sampler_4x.c \
          $(SRC_DIR)/rejsample.c \
          $(SRC_DIR)/reduce.c $(SRC_DIR)/ntt.c \
          $(SRC_DIR)/poly.c $(SRC_DIR)/polyvec.c \
          $(SRC_DIR)/rounding.c $(SRC_DIR)/shuttle_rans.c \
          $(SRC_DIR)/packing.c $(SRC_DIR)/sign.c

NROUNDS ?= 10

.PHONY: all build vectors vectors-128 vectors-256 clean clean-bin clean-vectors

# Default: build both modes and produce the vector files.
all: vectors

build: bin/gen_vectors_mode128 bin/gen_vectors_mode256

bin/gen_vectors_mode128: gen_vectors.c rng.c rng.h $(SOURCES) | bin
	$(CC) $(CFLAGS_BASE) -DSHUTTLE_MODE=128 \
	    -o $@ gen_vectors.c rng.c $(SOURCES) $(LDLIBS)

bin/gen_vectors_mode256: gen_vectors.c rng.c rng.h $(SOURCES) | bin
	$(CC) $(CFLAGS_BASE) -DSHUTTLE_MODE=256 \
	    -o $@ gen_vectors.c rng.c $(SOURCES) $(LDLIBS)

bin:
	mkdir -p bin

# Run each mode's generator. Output files land directly in this
# directory (ref/vector/) so they're easy to ship alongside the source.
vectors: vectors-128 vectors-256

vectors-128: bin/gen_vectors_mode128
	./bin/gen_vectors_mode128 . $(NROUNDS)

vectors-256: bin/gen_vectors_mode256
	./bin/gen_vectors_mode256 . $(NROUNDS)

clean: clean-bin

clean-bin:
	rm -rf bin

clean-vectors:
	rm -f keygen_mode128.kat sign_mode128.kat verify_mode128.kat \
	      keygen_mode256.kat sign_mode256.kat verify_mode256.kat
