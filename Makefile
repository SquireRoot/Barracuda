CC=nvcc

INCDIR=./include
OBJDIR=./build
SRCDIR=./src

CFLAGS=-I$(INCDIR)

_OBJS=main.o Kernels.o HelperFunctions.o
_DEPS=UserConfig.h Config.h Kernels.cuh HelperFunctions.h

OBJS=$(patsubst %,$(OBJDIR)/%,$(_OBJS))
DEPS=$(patsubst %,$(INCDIR)/%,$(_DEPS))

barracuda.out: DATA $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cu $(DEPS) build 
	$(CC) -c -o $@ $< $(CFLAGS)

build:
	mkdir $(OBJDIR)

DATA: cleanDATA
	mkdir DATA
	mkdir DATA/v DATA/curl DATA/dp

.PHONY: clean cleanDATA

clean: 
	rm -rf $(OBJDIR)/ barracuda.out DATA/

cleanDATA:
	rm -rf DATA/
