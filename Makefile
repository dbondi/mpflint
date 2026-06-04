CC = cc
CFLAGS = -O2 -march=native -shared -fPIC
INCLUDES = -I/opt/homebrew/include/flint -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lflint -lgmp -lmpfr -lm

all: pslqm3.dylib

pslqm3.dylib: pslqm3.c pslq_arb.h pslq_sort.c pslq_matmul.c pslq_4level.c
	$(CC) $(CFLAGS) -o $@ pslqm3.c $(INCLUDES) $(LDFLAGS)

clean:
	rm -f pslqm3.dylib

.PHONY: all clean
