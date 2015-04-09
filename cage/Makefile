include ../Flags.mk
DEST = ../bin

PROGS = install_mask install_chain install_ortho install_util install_fragment install_rechain install_best_ortho

all : $(PROGS)

install_mask :
	$(CC) $(CFLAGS) mask.d/*.c $(LFLAGS) -o $(DEST)/cage_mask

install_chain :
	$(CC) $(CFLAGS) chain.d/*.c $(LFLAGS) -o $(DEST)/cage_chain

install_ortho :
	$(CC) $(CFLAGS) ortho.d/*.c $(LFLAGS) -o $(DEST)/cage_ortho

install_util :
	$(CC) $(CFLAGS) util.d/*.c $(LFLAGS) -o $(DEST)/check_mask

install_fragment :
	$(CC) $(CFLAGS) fragment.d/*.c $(LFLAGS) -o $(DEST)/fragment_blocks

install_rechain :
	$(CC) $(CFLAGS) rechain.d/*.c $(LFLAGS) -o $(DEST)/rechain_blocks

install_best_ortho :
	$(CC) $(CFLAGS) convert_mto1.d/*.c $(LFLAGS) -o $(DEST)/best_orthologs

