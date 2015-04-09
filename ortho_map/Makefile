include ../Flags.mk
DEST = ../bin

PROGS = ortho_map adjust_symmetry get_stat count_ch find_ps ortho_view_info print_genes renumber adjust_events refine_maf

all : $(PROGS)

ortho_map :  
	$(CC) $(CFLAGS) ortho_map.c contigs_op.c regions.c read_maf.c util_gen.c util.c util_i.c util_ops.c util_algns.c id_ortho_conv.c cal_pid_conv.c sec_round.c redo_ops.c apply_ops.c tree_op.c const_graph.c find_gene_loss.c map_algns.c util_input.c read_algn.c id_ortho.c find_merging.c pred_regions.c kd_tree.c find_dup_copy.c ins_dup_copy.c extend_slist.c handle_tandem_dup.c update_init_algns.c pred_ops.c rollback.c pred_dels.c deal_gaps.c chain_pair_alg.c kd_op.c check_copy.c adjust_algn.c $(LFLAGS) -o $(DEST)/ortho_map

adjust_symmetry :  
	$(CC) $(CFLAGS) adjust_symmetry.c contigs_op.c read_maf.c util.c util_input.c util_i.c util_gen.c regions.c id_ortho_conv.c read_algn.c util_algns.c id_ortho.c kd_tree.c find_merging.c update_init_algns.c deal_gaps.c apply_ops.c util_ops.c $(LFLAGS) -o $(DEST)/adjust_symmetry

get_stat :   
	$(CC) $(CFLAGS) get_stat.c contigs_op.c read_maf.c util_input.c util.c util_i.c util_gen.c regions.c read_algn.c analysis.c util_ops.c $(LFLAGS) -o $(DEST)/get_stat

count_ch :   
	$(CC) $(CFLAGS) cal_init_algns.c contigs_op.c read_maf.c util_input.c util.c util_i.c util_gen.c regions.c read_algn.c util_ops.c $(LFLAGS) -o $(DEST)/count_ch

find_ps :  
	$(CC) $(CFLAGS) pseudogenes.c find_ps.c contigs_op.c id_ortho_conv.c id_ortho.c regions.c util_gen.c util.c util_i.c read_algn.c util_algns.c read_maf.c util_input.c kd_tree.c util_ops.c $(LFLAGS) -o $(DEST)/find_ps

ortho_view_info :
	$(CC) $(CFLAGS) ortho_view_info.c contigs_op.c map_genes.c id_ortho_conv.c id_ortho.c regions.c util_gen.c util.c util_i.c util_ops.c read_algn.c util_algns.c read_maf.c util_input.c kd_tree.c $(LFLAGS) -o $(DEST)/ortho_view_info

print_genes :
	$(CC) $(CFLAGS) print_genes.c contigs_op.c util_input.c util_i.c util_gen.c regions.c util.c util_ops.c $(LFLAGS) -o $(DEST)/print_genes

renumber :
	$(CC) $(CFLAGS) renumber.c contigs_op.c renum_util.c tree_op.c regions.c util_gen.c util.c util_i.c util_input.c kd_tree.c util_ops.c  $(LFLAGS) -o $(DEST)/renumber

adjust_events :
	$(CC) $(CFLAGS) adjust_ops.c contigs_op.c util_algns.c read_maf.c read_algn.c tree_op.c regions.c util_gen.c util.c util_i.c util_input.c util_ops.c $(LFLAGS) -o $(DEST)/adjust_events

refine_maf :
	$(CC) $(CFLAGS) refine_maf.c refine_algns.c contigs_op.c update_init_algns.c id_ortho_conv.c id_ortho.c regions.c util_gen.c util_ops.c util.c util_i.c read_algn.c util_algns.c read_maf.c util_input.c kd_tree.c find_merging.c deal_gaps.c apply_ops.c const_graph.c $(LFLAGS) -o $(DEST)/refine_maf

