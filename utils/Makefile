include ../Flags.mk
DEST = ../bin

PROGS = conv_head remove_headers remove_lower_tri_and_diag

all : $(PROGS)

conv_head :  
	$(CC) $(CFLAGS) conv_head.c util.c $(LFLAGS) -o $(DEST)/conv_head

remove_headers :  
	$(CC) $(CFLAGS) remove_headers.c util.c $(LFLAGS) -o $(DEST)/remove_headers

remove_lower_tri_and_diag :  
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c mz_scores.c remove_lower_tri_and_diag.c -o $(DEST)/remove_lower_tri_and_diag
