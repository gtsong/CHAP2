include ../Flags.mk
DEST = ../bin

PROGS = expand_branches reformat_nonredundant gc-info

all : $(PROGS)

expand_branches :
	$(CC) $(CFLAGS) branch.c tree_op.c regions.c util_i.c util.c $(LFLAGS) -o $(DEST)/expand_branches

reformat_nonredundant :
	$(CC) $(CFLAGS) reformat.c regions.c util_i.c util.c $(LFLAGS) -o $(DEST)/reformat_nonredundant

gc-info :
	$(CC) $(CFLAGS) gc-info.c tree_op.c regions.c util_i.c util.c read_maf.c kd_tree.c util_gen.c $(LFLAGS) -o $(DEST)/gc-info

