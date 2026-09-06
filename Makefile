# Build the five programs.  C99, no libraries beyond libc; the Python side
# needs only numpy.
CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99 -D_POSIX_C_SOURCE=200809L
BIN      = bin
PROGS    = $(BIN)/enum $(BIN)/shards $(BIN)/orbits $(BIN)/pools $(BIN)/pack

all: $(PROGS)

$(BIN):
	mkdir -p $(BIN)

$(BIN)/enum:   src/enum.c   src/lines.h | $(BIN); $(CC) $(CFLAGS) -o $@ src/enum.c
$(BIN)/shards: src/shards.c             | $(BIN); $(CC) $(CFLAGS) -o $@ src/shards.c
$(BIN)/orbits: src/orbits.c src/lines.h | $(BIN); $(CC) $(CFLAGS) -o $@ src/orbits.c
$(BIN)/pools:  src/pools.c              | $(BIN); $(CC) $(CFLAGS) -o $@ src/pools.c
$(BIN)/pack:   src/pack.c               | $(BIN); $(CC) $(CFLAGS) -o $@ src/pack.c

# The order-8 controls.  Everything the order-9 argument runs through is
# exercised here at an order where the answer is independently known.
test: all
	python3 tests/test_n8.py

clean:
	rm -rf $(BIN)

.PHONY: all test clean
