CC ?= gcc
# SHUTTLE_IMPL: ref | avx2
SHUTTLE_IMPL ?= ref
# SHUTTLE_MODE: 128 | 256
SHUTTLE_MODE ?= 128

SRC_DIR = ../../$(SHUTTLE_IMPL)

CFLAGS_BASE = -Wall -Wextra -Wpedantic -O2 -std=c99 \
              -I. -I$(SRC_DIR) -DSHUTTLE_MODE=$(SHUTTLE_MODE)
LDLIBS = -lm

ifeq ($(SHUTTLE_IMPL),avx2)
CFLAGS = $(CFLAGS_BASE) -mavx2 -mbmi2 -mpopcnt -march=native -mtune=native
SOURCES = $(SRC_DIR)/fips202.c $(SRC_DIR)/symmetric-shake.c \
          $(SRC_DIR)/fips202x4.c $(SRC_DIR)/approx_exp.c \
          $(SRC_DIR)/sampler.c $(SRC_DIR)/sampler_4x.c \
          $(SRC_DIR)/approx_log.c $(SRC_DIR)/rejsample.c \
          $(SRC_DIR)/reduce.c $(SRC_DIR)/ntt.c \
          $(SRC_DIR)/poly.c $(SRC_DIR)/polyvec.c \
          $(SRC_DIR)/rounding.c $(SRC_DIR)/shuttle_rans.c \
          $(SRC_DIR)/packing.c $(SRC_DIR)/sign.c
ASM_SOURCES = $(SRC_DIR)/f1600x4.S
else
CFLAGS = $(CFLAGS_BASE)
SOURCES = $(SRC_DIR)/fips202.c $(SRC_DIR)/symmetric-shake.c \
          $(SRC_DIR)/approx_exp.c \
          $(SRC_DIR)/sampler.c $(SRC_DIR)/sampler_4x.c \
          $(SRC_DIR)/approx_log.c $(SRC_DIR)/rejsample.c \
          $(SRC_DIR)/reduce.c $(SRC_DIR)/ntt.c \
          $(SRC_DIR)/poly.c $(SRC_DIR)/polyvec.c \
          $(SRC_DIR)/rounding.c $(SRC_DIR)/shuttle_rans.c \
          $(SRC_DIR)/packing.c $(SRC_DIR)/sign.c
ASM_SOURCES =
endif

OUT = out/kat_dump_$(SHUTTLE_IMPL)$(SHUTTLE_MODE)

.PHONY: all clean
all: $(OUT)

$(OUT): kat_dump.c rng.c rng.h $(SOURCES) $(ASM_SOURCES) | out
	$(CC) $(CFLAGS) -o $@ kat_dump.c rng.c $(SOURCES) $(ASM_SOURCES) $(LDLIBS)

out:
	mkdir -p out

clean:
	rm -rf out
