/*
	A Feature of GAGE: Indentifying the Orthologous Mappings of a Gene Cluster between Two Species
	Input - two self-alignments and a pairwise alignment for all two pairs of sequences and the three alignments should be merged in a maf file.
	In the merged maf file, each alignment is required to be seperated by '#' indicator.
	The default is one-to-many
*/

#include "main.h"
#include "regions.h"
#include "pred_sp.h"
#include "pred_ops.h"
#include "id_inpar.h"
#include "map_algns.h"
#include "read_maf.h"
#include "util_gen.h"
#include "rollback.h"
#include "ins_dup_copy.h"
#include "util.h"
#include "util_i.h"
#include "update_init_algns.h"
#include "chain_pair_alg.h"
#include "id_ortho.h"
#include "apply_ops.h"

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

int main(int argc, char **argv)
{
	struct DotList *algns; // the sequence starts at 1 and each range of a segment is a form of [a, b), i.e. the nucleotide at 'a' is included and one at 'b' is not
	struct DotList *init_algns; 
	int num_a[3];
	int count = 0; // the number of the alignments in the initial dot-plot
	int algn_type = -1; // first self-alignment if SELF1(0), second self-alignment if SELF2(1), and pairwise alignment if PAIR(2)
	int *num_algns; // the number of local alignments in the dot plot
	int *num_init_algns; // the number of local alignments in the initial dot plot
	int num_cur_init = 0;

	int *size1, *size2;
	int size;
	bool *is_x;
	int prev_num = 0;
	int opt_id;

	int *num_suspend_pairs;
	struct ID_List *suspend_list; // the suspend list
	int num_ins_regs;
	float scaling_value;
	struct slist rm_list;
	
	int *num_id; // the number of pairs of repeats to be eliminated
	int *threshold;
	int i;
	int num_sp;
	int *sp_order;
	int *rm_sp;
	int *left_sp;
	int cur_num_sp = 0;
	FILE *f;
	int maf_mode, mode;
	char species[100], species2[100];
	int run_mode; // the default is one-to-one
	struct ops_list *ops, *ops_cur_pos;
	int *num_ops;
	int num_dup_ops = 0;
	struct n_pair *contigs1, *contigs2;
	int *num_contigs1, *num_contigs2;
	int *len_sum1, *len_sum2;
//	int ctg_id1 = -1, ctg_id2 = -1;
//	int index = -1;

	debug_mode = FALSE;
	if( argc == 3 ) { 
		if( strcmp(argv[2], "debug-mode") == 0 ) {
			debug_mode = TRUE;
//			run_mode = INF_DUP;
			run_mode = ONE_TO_MANY;
		}	
		else if( strcmp(argv[2], "inf-dup") == 0 ) run_mode = INF_DUP;
		else if( strcmp(argv[2], "many-to-one") == 0 ) run_mode = ONE_TO_MANY;
		else if( strcmp(argv[2], "one-to-many") == 0 ) run_mode = ONE_TO_MANY;
		else if( strcmp(argv[2], "many-to-many") == 0 ) run_mode = MANY_TO_MANY;
		else if( strcmp(argv[2], "one-to-one") == 0 ) run_mode = ONE_TO_ONE; // the ancestral alignment is reconstructed
		else {
			fatal("args: dots-file and mode (many-to-many, one-to-many, or one-to-one) ) \n");
//			run_mode = INC_CONV;
		}
	}
	else if( argc == 4 ) {
		if( strcmp(argv[3], "debug-mode") == 0 ) {
			debug_mode = TRUE;
			if( strcmp(argv[2], "inf-dup") == 0 ) run_mode = INF_DUP;
			else if( strcmp(argv[2], "many-to-many") == 0 ) run_mode = MANY_TO_MANY;
			else if( strcmp(argv[2], "many-to-one") == 0 ) run_mode = MANY_TO_MANY;
			else if( strcmp(argv[2], "one-to-many") == 0 ) run_mode = ONE_TO_MANY;
			else if( strcmp(argv[2], "one-to-one") == 0 ) run_mode = ONE_TO_ONE; // the ancestral alignment is reconstructed
			else {
				fatal("args: dots-file and mode (many-to-many, one-to-many, or one-to-one) ) \n");
//				run_mode = INC_CONV;
			}
		}	
		else {
			fatal("args: dots-file ( and conv-file or mode ) \n");
		}
	}
	else if(argc != 2 )
	{
		fatal("args: dots-file ( and conv-file or mode ) \n");
	}
//	else run_mode = INF_DUP; // default run-mode
	else run_mode = ONE_TO_MANY; // default run-mode

	is_x = (bool *) ckalloc(sizeof(bool));
  num_init_algns = (int *) ckalloc(sizeof(int));
	num_suspend_pairs = (int *) ckalloc(sizeof(int));
	num_id = (int *) ckalloc(sizeof(int));
	threshold = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_ops = (int *) ckalloc(sizeof(int));

	if( run_mode == INC_CONV ) {
	}

	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("file %s does not exist\n", argv[1]);
	}

	while(fgets(S, BIG, f)) {
		if( S[0] == '#' ) {
			while( S[0] == '#' ) {
				fgets(S, BIG, f);
				if( strncmp(S, "##maf", 5) == 0 ) {
					num_a[algn_type] = count;
					algn_type++;
					count = 0;
				}
			}

			if( algn_type != -1 ) {
				if( (algn_type == SELF1) || (algn_type == SELF2) ) num_a[algn_type] = count;
					
				else num_a[algn_type] = count;
			}
			algn_type++;
			count = 0;	
		}

  	if( S[0] == 'a' ) {
			if( (algn_type == -1) || (algn_type > PAIR) ) {
				fatal("The input is not a ##maf file\n");
			}
			count++;

			if( (algn_type == PAIR) && (count == 1) ) {
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", species) != 1) || (sscanf(T, "%*s %s %*s", species2) != 1)) {}
  		}
		}
  }
	fclose(f);

	if( algn_type != -1 ) {
		if( (algn_type == SELF1) || (algn_type == SELF2) ) num_a[algn_type] = count;
		else num_a[algn_type] = count;
	}
	else {
		fatal("The input is not a ##maf file\n");
	}

	if( algn_type == PAIR ) num_sp = 2;
	else if (algn_type == SELF1) num_sp = 1;
	else {
		fatalf("more than two species not supported yet algn_type: %d\n", algn_type);
	}

	count = 0;
	for( i = 0; i <= algn_type; i++ )
	{
		count = count + num_a[i]; // the number of alignments in the initial dot-plot
	}
	
	if( algn_type == PAIR ) maf_mode = C_MODE;
	else if( (algn_type == SELF1) || (algn_type == SELF2) ) maf_mode = D_MODE;
	else maf_mode = G_MODE;

	sp_order = (int *) ckalloc(num_sp * sizeof(int));
	rm_sp = (int *) ckalloc(num_sp * sizeof(int));
	left_sp = (int *) ckalloc(num_sp * sizeof(int));

	for( i = 0; i < num_sp; i++ ) {
		sp_order[i] = i;
		rm_sp[i] = num_sp - i - 1;
		left_sp[i] = i;
	}

	if( count > 0 )
	{
		algns = (struct DotList *) ckalloc(sizeof(struct DotList) * count);
		init_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * count);
		ops = (struct ops_list *) ckalloc(sizeof(struct ops_list) * count);
		ops_cur_pos = (struct ops_list *) ckalloc(sizeof(struct ops_list) * count);
		suspend_list = (struct ID_List *) ckalloc(sizeof(struct ID_List) * count);

		initialize_algns(algns, count);
		initialize_algns(init_algns, count);
	}

	for( i = 0; i < count; i++ ) { // initialization for ops
		ops[i].sign = 'n'; // nothing
		ops[i].id = -1;	
		ops[i].srcStart = 0;	
		ops[i].srcEnd = 0;	
		ops[i].dstStart = 0;	
		ops[i].dstEnd = 0;	
		ops[i].src_b = 0;	
		ops[i].src_e = 0;	
		ops[i].dst_b = 0;	
		ops[i].dst_e = 0;	
		ops[i].pid = 0;	
		ops[i].sp_id = 0;	
		ops[i].dir = 0;	
		
		suspend_list[i].is_x = true;
		suspend_list[i].m_id = -1;
		suspend_list[i].left_id = -1;
		suspend_list[i].right_id = -1;
		suspend_list[i].f_is_x = true;
		suspend_list[i].is_t_ins = true;
	}

	*num_ops = 0;
	*num_algns = 0;
	*num_init_algns = 0;
	rm_list.id = NO_EXIST; // initialization
	rm_list.val = 0;
	rm_list.val_red = 0;
	rm_list.sp_state = 0;
	rm_list.add_sp_state = 0;
	rm_list.is_x = true;
	*size1 = 0;
	*size2 = 0;

	if( count > 0 ) {
		read_maf(argv[1], maf_mode, algns, num_algns, size1, size2);
		read_maf(argv[1], maf_mode, init_algns, num_init_algns, size1, size2);
	}

	if( debug_mode == TRUE ) {
		for(i = 0; i < (*num_algns); i++) {
			printf("%d: %d-%d, %d-%d\n", algns[i].sp_id, algns[i].x.lower, algns[i].x.upper, algns[i].y.lower, algns[i].y.upper);
		}
	}

	if( algn_type == PAIR ) { // for the combined alignment of two self-alignments and a pairwise alignment, the position of the second species is added to the length of the first sequence

		num_contigs1 = (int *) ckalloc(sizeof(int));
		num_contigs2 = (int *) ckalloc(sizeof(int));

		if( (*num_algns) > 0 ) {
			contigs1 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * 2 * (*num_algns));
			contigs2 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * 2 * (*num_algns));

			init_n_pair(contigs1, 0, (2*(*num_algns))-1);
			init_n_pair(contigs2, 0, (2*(*num_algns))-1);
		}
		*num_contigs1 = 0;
		*num_contigs2 = 0;

		adjust_multi_contig_pos(algns, *num_algns, size1, size2, contigs1, num_contigs1, contigs2, num_contigs2);
		adjust_algn_pos(init_algns, *num_init_algns, contigs1, *num_contigs1, size1, contigs2, *num_contigs2, size2, CTG_NOT_ASSIGNED);

		if( (*num_contigs1) > 0 ) len_sum1 = (int *) ckalloc(sizeof(int) * (*num_contigs1));
    else len_sum1 = (int *) ckalloc(sizeof(int));

    if( (*num_contigs2) > 0 ) len_sum2 = (int *) ckalloc(sizeof(int) * (*num_contigs2));
    else len_sum2 = (int *) ckalloc(sizeof(int));

    cal_length_sum(len_sum1, contigs1, *num_contigs1);
    cal_length_sum(len_sum2, contigs2, *num_contigs2);

		for(i = 0; i < (*num_algns); i++) {
			if( algns[i].sp_id == SELF2 ) {
				algns[i].x = assign_I(algns[i].x.lower+(*size1), algns[i].x.upper+(*size1));
				algns[i].y = assign_I(algns[i].y.lower+(*size1), algns[i].y.upper+(*size1));
			}
			else if( algns[i].sp_id == PAIR ) {
				algns[i].y = assign_I(algns[i].y.lower+(*size1), algns[i].y.upper+(*size1));
			}
		}
	}

	if( debug_mode == TRUE ) {
		for(i = 0; i < (*num_algns); i++) {
			printf("%d: %d-%d, %d-%d\n", algns[i].sp_id, algns[i].x.lower, algns[i].x.upper, algns[i].y.lower, algns[i].y.upper);
		}
	}

	opt_id = 0;
	if( algn_type == PAIR ) {
		size = (*size1) + (*size2);
		mode = BEFORE_SP;
	}
	else {
		size = (*size1);
		mode = AFTER_SP;
	}

	f = fopen(argv[1], "r");
	while((*num_algns > 0) && (opt_id != NO_EXIST) && (mode == BEFORE_SP))
	{
		prev_num = *num_algns;

		rm_list = iden_inpar(num_algns, algns, size, STRICT, num_suspend_pairs, suspend_list, mode, init_algns, f);

		if( rm_list.id == NO_EXIST ) rm_list = iden_inpar(num_algns, algns, size, LOOSE, num_suspend_pairs, suspend_list, mode, init_algns, f);

		if( (rm_list.id != SP_EVENT) && ((rm_list.id != NO_EXIST) && (rm_list.id != DEL_EXIST)) && ((rm_list.sp_state == INSERTED_IN_X) || (rm_list.sp_state == INSERTED_IN_Y) || (rm_list.sp_state == INS_IN_X_TAN) || (rm_list.sp_state == INS_IN_Y_TAN)) )
		{
			num_ins_regs = *num_suspend_pairs;
			rollback_ins_dup(rm_list.sp_state, rm_list.id, algns, suspend_list, num_ins_regs, init_algns, f); // when the inferred duplication occurs in an inserted region of the suspend list
		}
		
		if( rm_list.id == NO_EXIST ) opt_id = NO_EXIST;
		else if( rm_list.id == DEL_EXIST ) 
		{
			size = size + rm_list.val;
		}
		else if( rm_list.id == SP_EVENT )
		{
			predict_sp_op(sp_order[cur_num_sp], rm_sp[cur_num_sp], left_sp[cur_num_sp], num_algns, algns, num_ops, ops);
			*num_ops = (*num_ops) + 1;
			cur_num_sp++;
		}
		else
		{
/*
			printf("%d alignments left\n", *num_algns);
			if( rm_list.id != NO_EXIST ) {
				index = algns[rm_list.id].index;
				ctg_id1 = init_algns[index].ctg_id1;
				ctg_id2 = init_algns[index].ctg_id2;
				if( (ctg_id1 != -1) && (ctg_id2 != -1) ) {
					if( init_algns[index].sp_id == SELF1 ) {
						printf("%d: %d-%d %d-%d in %d.%s-%s\n", init_algns[index].identity, init_algns[index].x.lower-len_sum1[ctg_id1], init_algns[index].x.upper-len_sum1[ctg_id1], init_algns[index].y.lower-len_sum1[ctg_id2], init_algns[index].y.upper-len_sum1[ctg_id2], init_algns[index].sp_id, contigs1[ctg_id1].name2, contigs1[ctg_id2].name2);
					}
					else if( init_algns[index].sp_id == SELF2 ) {
						printf("%d: %d-%d %d-%d in %d:%s-%s\n", init_algns[index].identity, init_algns[index].x.lower-len_sum2[ctg_id1], init_algns[index].x.upper-len_sum2[ctg_id1], init_algns[index].y.lower-len_sum2[ctg_id2], init_algns[index].y.upper-len_sum2[ctg_id2], init_algns[index].sp_id, contigs2[ctg_id1].name2, contigs2[ctg_id2].name2);
					}
					else {
						fatalf("Not self-alignment: %d-%d %d-%d in %d.%s-%s\n", init_algns[index].x.lower-len_sum1[ctg_id1], init_algns[index].x.upper-len_sum1[ctg_id1], init_algns[index].y.lower-len_sum2[ctg_id2], init_algns[index].y.upper-len_sum2[ctg_id2], init_algns[index].sp_id, contigs1[ctg_id1].name2, contigs2[ctg_id2].name2);
					}
				}
			}
*/
			rollback_init_dots(algns, rm_list.id, rm_list.is_x, *num_algns, init_algns, *num_init_algns, f, *num_ops, ops, SECOND_RUN, *size1);
			predict_op(rm_list.is_x, rm_list.id, num_algns, algns, rm_list.sp_state, *num_ops, ops);
			*num_ops = (*num_ops) + 1;
			if( prev_num <= (*num_algns) ) {
				fatalf("possiblity of infinite iterations %d-%d\n", prev_num, *num_algns);
			}
		}


		if( (mode == BEFORE_SP) && (rm_list.id == SP_EVENT) ) {
			mode = AFTER_SP;
		}
		else if( (mode == BEFORE_SP) && (rm_list.id == NO_EXIST) ) {
			if( run_mode == INF_DUP ) {
				rollback_sp(algns, num_algns, init_algns, *num_init_algns, *num_ops, ops, SELF2);
				*num_ops = (*num_ops) + 1;
				mode = AFTER_SP;
				opt_id = 0;
			}
			else {
				obtain_ortho_algns(algns, *num_algns, init_algns, *num_init_algns); // one-to-one putative orthologous alignments of ancestral common genomic regions
				if( run_mode == ONE_TO_ONE ) {
					num_cur_init = *num_init_algns;
					map_one_to_one(num_cur_init, init_algns, f);
				}
				else if( run_mode == ONE_TO_MANY ) {
					if( *num_ops > 0 ) {
						num_dup_ops = cal_cur_pos_ops(*num_ops, ops, ops_cur_pos, SELF1, 0);
						redo_dups_for_mtom(num_dup_ops, ops_cur_pos, *num_init_algns, init_algns, f, REF_SEQ);
					}
					num_cur_init = *num_init_algns;
					map_one_to_one(num_cur_init, init_algns, f);
				}
				else if( run_mode == MANY_TO_MANY ) {
					if( *num_ops >  0 ) {
						redo_dups_for_mtom(*num_ops, ops, *num_init_algns, init_algns, f, REF_SEQ);
					}
					num_cur_init = *num_init_algns;
					map_one_to_one(num_cur_init, init_algns, f);
					if( *num_ops > 0 ) {
						redo_dups_for_mtom(*num_ops, ops, *num_init_algns, init_algns, f, SELF2);
					}
				}
				if( (run_mode == ONE_TO_ONE) || (run_mode == ONE_TO_MANY) || (run_mode == MANY_TO_MANY) ) write_init_maf_stdout(init_algns, *num_init_algns, contigs1, contigs2, len_sum1, len_sum2, *num_contigs1, *num_contigs2, *size1, *size2, f, PAIR);
				rollback_sp(algns, num_algns, init_algns, *num_init_algns, *num_ops, ops, SELF2);
				mode = AFTER_SP;
				opt_id = 0;
			}
		}
	}

	scaling_value = (float)1;
 	if( run_mode == INF_DUP) output_ops(*num_ops, ops, *size2, scaling_value);

	if( (debug_mode == TRUE) || (run_mode == INF_DUP) ) {
		printf("-------------\n");
		printf("For debugging\n");
		for( i = ((*num_ops)-2); i >= 0; i-- ) {
//			printf("%d-%d is copied from %d-%d: %f\n", ops[i].dstStart, ops[i].dstEnd, ops[i].srcStart, ops[i].srcEnd, ops[i].pid);
			printf("%d-%d is copied from %d-%d: %d\n", ops[i].dstStart, ops[i].dstEnd, ops[i].srcStart, ops[i].srcEnd, ops[i].sp_id);
		}
	}

	fclose(f);

	if( count > 0 ) {
		free(algns);
		free(init_algns);
		free(ops_cur_pos);
		free(ops);
		free(suspend_list);
	}
	else if( count == 0 ) {
		if( (run_mode == ONE_TO_ONE) || (run_mode == ONE_TO_MANY) || (run_mode == MANY_TO_MANY) ) 
		{
			printf("##maf version=1 scoring=lastz-pid\n");
		}
	}
	
	if( algn_type == PAIR ) {
		if( (*num_algns) > 0 ) {
			free(contigs1);
			free(contigs2);
		}
		free(num_contigs1);
		free(num_contigs2);
	}
	free(num_ops);
	free(size2);
	free(size1);
	free(num_algns);
	free(num_init_algns);
	free(threshold);
	free(num_id);
	free(num_suspend_pairs);
	free(is_x);
	free(sp_order);
	free(rm_sp);
	free(left_sp);
	return EXIT_SUCCESS;
}

void output_ops(int num_ops, struct ops_list * ops, int len, float sc)
{
	int i = 0;
	char op;
	int a, b, c, d;
	int t_a, t_b, t_c, t_d;

	for( i = num_ops - 1; i >= 0; i-- )
	{
		op = ops[i].sign;
		a = ops[i].src_b;
		b = ops[i].src_e;
		c = ops[i].dst_b;
		d = ops[i].dst_e;
		t_a = (int)((float)a/sc);
		t_b = (int)((float)b/sc);
		t_c = (int)((float)c/sc);
		t_d = (int)((float)d/sc);
		if( (op == '+') || (op == '-') ) printf("dup %c %d %d %d %d\n", op, t_a, t_b, t_c, t_d - t_c);
		else if( (op == 'd') || (op == 'l') ) printf("del . %d %d\n", t_a, t_b);
		else if( op == 'i') printf("inv . %d %d\n", t_a, t_b);
		else if( op == 'c') printf("con + %d %d %d %d\n", t_a, t_b, t_c, t_d);
		else if( op == 'v') printf("con + %d %d %d %d\n", t_a, t_b, t_c, t_d);
		else if( op == 's') printf("spe\n");
	}
	printf("len %d\n", len);
}

bool is_current_sp( int code, int *list1, int *list2, int num_sp)
{
	int i;
	bool res = false;

	for( i = 0; i < num_sp; i++ )
	{
		if( list1[i] == code ) res = true;
		else if( list2[i] == code ) res = true;
	}

	return res;
}
