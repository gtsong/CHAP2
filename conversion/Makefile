include ../Flags.mk
DEST = ../bin

PROGS = remove_lower_tri remove_repeats chain md_quadruplet_entire_region GC_history 3seq_2D_cluster separate_pvalue FDR_cluster FDR_test combine_all_pvalue_cluster show_GC_cluster check_species_tree 3seq_2D_cluster_fast

all : $(PROGS)

remove_lower_tri : util.h util.c seq.h seq.c nib.h nib.c charvec.h charvec.c multi_util.c multi_util.h maf.h maf.c mz_scores.h mz_scores.c remove_lower_tri.c
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c mz_scores.c remove_lower_tri.c $(LFLAGS) -o $(DEST)/remove_lower_tri

remove_repeats : util.h util.c seq.h seq.c nib.h nib.c charvec.h charvec.c multi_util.c multi_util.h maf.h maf.c mz_scores.h mz_scores.c  remove_repeats.h remove_repeats.c 
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c mz_scores.c remove_repeats.c $(LFLAGS) -o $(DEST)/remove_repeats

chain : util.h util.c seq.h seq.c nib.h nib.c charvec.h charvec.c multi_util.c multi_util.h maf.h maf.c mz_scores.h mz_scores.c chain.h chain.c
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c mz_scores.c chain.c $(LFLAGS) -o $(DEST)/chain

md_quadruplet_entire_region : util.h util.c seq.h seq.c nib.h nib.c charvec.h charvec.c contigs_op.h contigs_op.c multi_util.c multi_util.h maf.h maf.c mz_scores.h mz_scores.c md_quadruplet_entire_region.c
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c contigs_op.c multi_util.c maf.c mz_scores.c md_quadruplet_entire_region.c $(LFLAGS) -o $(DEST)/md_quadruplet_entire_region

GC_history : util.h util.c seq.h seq.c nib.h nib.c charvec.h charvec.c multi_util.c multi_util.h maf.h maf.c mz_scores.h mz_scores.c contigs_op.c GC_history.c
	$(CC) $(CFLAGS) util.c seq.c nib.c charvec.c multi_util.c maf.c mz_scores.c contigs_op.c GC_history.c $(LFLAGS) -o $(DEST)/GC_history

3seq_2D_cluster : 
	$(CC) $(CFLAGS) 3seq_2D_cluster.c $(LFLAGS) -o $(DEST)/3seq_2D_cluster

separate_pvalue :
	$(CC) $(CFLAGS) separate_pvalue.c $(LFLAGS) -o $(DEST)/separate_pvalue

FDR_cluster :
	$(CC) $(CFLAGS) FDR_cluster.c $(LFLAGS) -o $(DEST)/FDR_cluster

FDR_test :
	$(CC) $(CFLAGS) FDR_test.c $(LFLAGS) -o $(DEST)/FDR_test

combine_all_pvalue_cluster :
	$(CC) $(CFLAGS) combine_all_pvalue_cluster.c $(LFLAGS) -o $(DEST)/combine_all_pvalue_cluster

show_GC_cluster :
	$(CC) $(CFLAGS) show_GC_cluster.c $(LFLAGS) -o $(DEST)/show_GC_cluster

check_species_tree :
	$(CC) $(CFLAGS) check_species_tree.c $(LFLAGS) -o $(DEST)/check_species_tree

3seq_2D_cluster_fast :
	$(CC) $(CFLAGS) 3seq_2D_cluster_fast.c $(LFLAGS) -o $(DEST)/3seq_2D_cluster_fast

#gc-find : gc-find.c
#	$(CC) $(CFLAGS) gc-find.c lib.c $(LFLAGS) -o $(DEST)/gc-find

