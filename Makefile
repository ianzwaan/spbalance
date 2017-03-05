CFLAGS=-std=c99 -pedantic-errors -Wall -Wextra -O2

MEX=mex
MEXFLAGS=CFLAGS='$$CFLAGS $(CFLAGS)' -largeArrayDims

TRGT=mexspbalance
DEPS=Makefile
SRCS=mexspbalance.c
LIBS=-lm

all: $(TRGT)

$(TRGT): $(SRCS) $(DEPS)
	$(MEX) $(MEXFLAGS) -output $@ $(SRCS) $(LIBS)
