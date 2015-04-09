include ../Flags.mk
DEST = ../bin

PROGS = exons2underlays

all : $(PROGS)

exons2underlays :
	$(CC) $(CFLAGS) -O3 -I. exons2underlays.c exons.c util.c $(LFLAGS) -o $(DEST)/exons2underlays

