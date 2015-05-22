/*
	Indentifying the Orthologous Mappings of a Gene Cluster between Two Species, 
	following the concept of orthology-by-content
	Input - two self-alignments and a pairwise alignment for all two pairs of sequences and the three alignments should be merged in a maf file. 
	all.gc - gene conversion information
	In the merged maf file, each alignment is required to be seperated by '#' indicator.
	the coordinates of all alignments, annotations, and events are adjusted temporally during main steps and they should be converted to the original ones based on contigs.
*/

#include "main.h"
#include "regions.h"
#include "read_maf.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "util_ops.h"
#include "util_algns.h"
#include "id_ortho_conv.h"
#include "id_ortho.h"
#include "cal_pid_conv.h"
#include "sec_round.h"
#include "redo_ops.h"
#include "apply_ops.h"
#include "tree_op.h"
#include "const_graph.h"
#include "find_gene_loss.h"
#include "map_algns.h"
#include "util_input.h"
#include "update_init_algns.h"
#include "contigs_op.h"

int count_node;
int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

void output_contigs(char *fname, struct n_pair *contigs, int *old_len_sum, int num_contigs);
void output_ops(int num_ops, struct ops_list *ops, int len, float sc, int sp_id);
void assign_dot_list(struct DotList *old_algn, int num_pair, struct DotList *new_algn);
bool is_current_sp(int code, int *list1, int *list2, int num_sp);

int main(int argc, char **argv)
{
	struct DotList *algns; // the sequence starts at 1 and each range of a segment is a form of [a, b), i.e. the nucleotide at 'a' is included and one at 'b' is not
	struct DotList *init_algns; 
	struct DotList *temp_algns;
	int num_a[3];
	int count = 0; // the number of the alignments in the initial dot-plot
	int *num_alloc_blocks;
	int *num_alloc_del_ops;
	int algn_type = -1; // first self-alignment if SELF1(0), second self-alignment if SELF2(1), and pairwise alignment if PAIR(2)
	int *num_algns; // the number of local alignments in the dot plot
	int *num_init_algns; // the number of local alignments in the initial dot plot

	int *size1, *size2;
	int size = 0;
	bool *is_x;
	int opt_id = -1;
	int *num_suspend_pairs;
	struct ID_List *suspend_list; // the suspend list
	float scaling_value;
	
	int *num_id; // the number of pairs of repeats to be eliminated
	int *threshold;
	int i = 0, j = 0;
	int num_sp = 0;
	int *sp_order;
	int *rm_sp;
	int *left_sp;
	FILE *f, *g;
	int maf_mode = C_MODE;
	int mode = AFTER_SP;
	char species[LEN_NAME], species2[LEN_NAME];
	char out_species[LEN_NAME];
	char temp_name1[LEN_NAME], temp_name2[LEN_NAME], temp_name3[LEN_NAME];
	char sp_name1[LEN_NAME], sp_name2[LEN_NAME];
	char ctg_name1[LEN_NAME], ctg_name2[LEN_NAME];
	bool is_spname_determined = false;
	int run_mode = ONE_TO_ONE; // the default is one-to-one
	int is_gc_given = FALSE;
	int is_exons_given = FALSE;
	int is_ancestral = FALSE;
	struct ops_list *ops;
	struct ops_list *ops_cur_pos;
	int num_dup_ops = 0;
	int *num_ops;
	int ops_id = 0;
	struct ops_list *prev_ops;
	int num_prev_ops = 0;
	struct n_pair *contigs1, *contigs2;
	int *num_contigs1, *num_contigs2;
	int ctg_id = -1;
	
 	struct cv_list *cv, *init_cv; // the conversion events detected by Chih-Hao's program
	struct cv_list *init_dup_cv;
	int *num_cv, num_init_cv = 0, num_dup_cv = 0, num_alloc_cv = 0;
	int num_cur_cv = 0;
	char buf[1000], tree_line[1000];
 	struct exons_list *exons, *init_exons; // the conversion events detected by Chih-Hao's program
	struct exons_list *genes;
	int *num_exons, num_init_exons = 0;
	int *num_genes, num_init_genes = 0;
	int algn_len_x = 0, algn_len_y = 0;
	struct exons_list *skip_reg1, *skip_reg2;
	int num_skip1 = 0, num_skip2 = 0;
	float avg_pid = (float)0, cur_pid1 = (float)0, cur_pid2 = (float) 0;
	float *temp_val;
	int total_len = 0;
	struct slist *sorted1, *sorted2;
	char op_name[10], ori = '+';
	int c1 = 0, c2 = 0, c3 = 0, c4 = 0;
	struct p_tree *sp_tree = NULL;
	struct sp_list *sp_code = NULL;
	int num_sp_code = 0, temp_num_sp = 0;
	int sp1_code = -1, sp2_code = -1, out_code = -1, temp_code = -1;
	struct I src, dst;
	int *old_dups;
	int num_old_dups = 0;
	struct ops_list *del_ops;
	int *num_del_ops;
	int num_ortho_algns = 0;
	bool is_event_output_required = false;
	int *len_sum1, *len_sum2;
	int len1 = 0, len2 = 0;
	bool is_contig_printed = false;
	char *status;
//	struct bipar_node *graph;

	src = assign_I(0,1);
	dst = assign_I(0,1);
	debug_mode = FALSE;
	count_node = 0;
	is_gc_given = FALSE;
	is_exons_given = FALSE;
	is_ancestral = FALSE;
	strcpy(species, "sp_none");
	strcpy(species2, "sp2_none");
	strcpy(out_species, "out_none");
	strcpy(temp_name1, "");
	strcpy(temp_name2, "");
	strcpy(temp_name3, "");
	strcpy(op_name, "");
	strcpy(tree_line, "");
	strcpy(buf, "");

	if( argc >= 5 ) {
		if( strcmp(argv[4], "debug-mode") == 0 ) {
			if( argc == 5 ) {
				debug_mode = TRUE;
				run_mode = MANY_TO_MANY;
			}
			else {
				fatal("args: dots-file and mode (inf-dup, many-to-many, one-to-many, one-to-one, position-ortho, or content-ortho) )\n");
			}
		}
		else if( strcmp(argv[2], "inf-dup") == 0 ) run_mode = INF_DUP;
		else if( strcmp(argv[2], "many-to-one") == 0 ) run_mode = ONE_TO_MANY;
		else if( strcmp(argv[2], "one-to-many") == 0 ) run_mode = ONE_TO_MANY;
		else if( strcmp(argv[2], "many-to-many") == 0 ) run_mode = MANY_TO_MANY;
		else if( strcmp(argv[2], "one-to-one") == 0 ) run_mode = ONE_TO_ONE; // the ancestral alignment is reconstructed
		else if( strcmp(argv[2], "content-ortho") == 0 ) run_mode = CONTENT_ORTHO; // the ancestral alignment is reconstructed
		else if( strcmp(argv[2], "position-ortho") == 0 ) run_mode = POSITION_ORTHO; // the ancestral alignment is reconstructed
		else {
			fatal("args: dots-file and mode (inf-dup, many-to-many, one-to-many, one-to-one, position-ortho, or content-ortho) )\n");
		}

		if( argc == 5 ) {}
		else if( argc == 6 ) {
			if( strcmp(argv[5], "debug-mode") == 0 ) {
				debug_mode = TRUE;
			}
			else if( (f = ckopen(argv[5], "r")) != NULL ) { // conversion results
				fclose(f);
				is_gc_given = TRUE;
//				is_event_output_required = true;	
			}
			else {
				fatal("args: dots-file, mode (many-to-many, one-to-many, one-to-one, position-ortho, or content-ortho) ), conv-file, and debug-mode\n");
			}
		}
		else if( argc == 7 ) {
			if( strcmp(argv[6], "debug-mode") == 0 ) {
				is_gc_given = TRUE;
				debug_mode = TRUE;
			}
			else if( ((f = ckopen(argv[5], "r")) != NULL) && ((g = ckopen(argv[6], "r")) != NULL) ) { // conversion results and gene annotations
				is_gc_given = TRUE;
				is_exons_given = TRUE;
//				is_event_output_required = true;	
				fclose(g);
				fclose(f);
			}
			else {
				fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, (and debug-mode)\n");
			}
		}
		else if( argc == 8 ) {
			is_gc_given = TRUE;
			is_exons_given = TRUE;
			if( ((f = ckopen(argv[5], "r")) != NULL) && ((g = ckopen(argv[6], "r")) != NULL) ) 
			{ 
				if (strcmp(argv[7], "debug-mode") == 0)  {
					debug_mode = TRUE;
				}
				else {
					is_event_output_required = true;
				}
			}
			else {
				fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, (and debug-mode)\n");
			}
		}
		else if( argc == 9 ) {
			is_gc_given = TRUE;
			is_exons_given = TRUE;
			if( ((f = ckopen(argv[5], "r")) != NULL) && ((g = ckopen(argv[6], "r")) != NULL) ) {
				if ( strcmp(argv[8], "debug-mode") == 0 ) {
					debug_mode = TRUE;
				}
				else {
					fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, events-file (or debug-mode)\n");
				}

				is_event_output_required = true;	
				fclose(g);
				fclose(f);
			}
			else {
				fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, events-file (or debug-mode)\n");
			}
		}
		else if( argc == 10 ) {
			if( ((f = ckopen(argv[5], "r")) != NULL) && ((g = ckopen(argv[6], "r")) != NULL) ) {
				is_gc_given = TRUE;
				is_exons_given = TRUE;
				is_contig_printed = true;
//				debug_mode = TRUE;
				if( strcmp(argv[7], "no-events") == 0 ) {
					is_event_output_required = false;	
				}
				else {
					is_event_output_required = true;	
				}
			}
			else {
				fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, events-file, contigs_file contigs_file (or debug-mode/ancestral)\n");
			}
			fclose(f);
			fclose(g);
		}
		else if( argc == 11 ) {
			if( ((f = ckopen(argv[5], "r")) != NULL) && ((g = ckopen(argv[6], "r")) != NULL) ) {

				if( strcmp(argv[10], "debug-mode") == 0 ) {
					debug_mode = TRUE;
				}
				else {
					fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, events-file, contigs_file contigs_file (and debug-mode)\n");
				}

				if( strcmp(argv[7], "no-events") == 0 ) {
					is_event_output_required = false;	
				}
				else {
					is_event_output_required = true;	
				}

				is_gc_given = TRUE;
				is_exons_given = TRUE;
				is_contig_printed = true;
			}
			else {
				fatal("args: dots-file, mode (many-to-many, one-to-many, or one-to-one) ), conv-file, exons, events-file gene/exon/event-mode (and debug-mode)\n");
			}
		}
		else {
			fatal("args: dots-file, mode, conv-file, exons\n");
		}
	}
	else if((argc < 5) || (argc > 11))
	{
		fatal("args: dots-file and mode ( and conv-file and exons )\n");
	}

	strcpy(species, argv[3]);
	strcpy(species2, argv[4]);
	is_x = (bool *) ckalloc(sizeof(bool));
  num_init_algns = (int *) ckalloc(sizeof(int));
	num_suspend_pairs = (int *) ckalloc(sizeof(int));
	num_id = (int *) ckalloc(sizeof(int));
	threshold = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	num_exons = (int *) ckalloc(sizeof(int));
	num_genes = (int *) ckalloc(sizeof(int));
	num_cv = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_ops = (int *) ckalloc(sizeof(int));
//	sp_tree = (struct p_tree *) ckalloc(sizeof(struct p_tree));
//	init_tree(sp_tree);

	*num_exons = 0;
	*num_genes = 0;
	num_init_exons = 0;
	num_init_genes = 0;
	f = ckopen(argv[1], "r");
	while((status = fgets(S, BIG, f)) != NULL) {
		if( S[0] == '#' ) {
			while( (status != NULL) && (S[0] == '#') ) {
				if( strncmp(S, "##maf", 5) == 0 ) {
					if( algn_type >= 0 ) num_a[algn_type] = count;
					algn_type++;
					count = 0;
				}
				status = fgets(S, BIG, f);
			}

			if( algn_type != -1 ) {
				if( (algn_type == SELF1) || (algn_type == SELF2) ) num_a[algn_type] = count;
					
				else num_a[algn_type] = count;
			}
			count = 0;	
		}

  	if( (status != NULL) && (S[0] == 'a') ) {
			if( (algn_type == -1) || (algn_type > PAIR) ) {
				fatal("The input is not a ##maf file\n");
			}
			count++;

			if( (algn_type == PAIR) && (count == 1) ) {
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", temp_name1) != 1) || (sscanf(T, "%*s %s %*s", temp_name2) != 1)) {}
  		}
			else if( (algn_type == PAIR) && (is_spname_determined == false) ) {
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", temp_name1) != 1) || (sscanf(T, "%*s %s %*s", temp_name2) != 1)) {}
				is_spname_determined = true;
			}
		}
  }

	num_a[algn_type] = count;
	fclose(f);

	if( is_spname_determined == true ) {
		concat_ctg_name(temp_name1, species, temp_name3);
		concat_ctg_name(temp_name2, species2, temp_name3);
	}
	else {
		if( count == 1 ) {
			concat_ctg_name(temp_name1, species, temp_name3);
			concat_ctg_name(temp_name2, species2, temp_name3);
		}
		else if( count != 0 ) {
			fatal("species names not in the MAF file\n");
		}
	}

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

	if( count == 0 ) {
		algns = (struct DotList *) ckalloc(sizeof(struct DotList));
		init_algns = (struct DotList *) ckalloc(sizeof(struct DotList));
		temp_algns = (struct DotList *) ckalloc(sizeof(struct DotList));
		ops = (struct ops_list *) ckalloc(sizeof(struct ops_list));
		suspend_list = (struct ID_List *) ckalloc(sizeof(struct ID_List));
	}
	else {
		algns = (struct DotList *) ckalloc(sizeof(struct DotList) * count);
		init_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * count);
		temp_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * count);
		ops = (struct ops_list *) ckalloc(sizeof(struct ops_list) * count);
		suspend_list = (struct ID_List *) ckalloc(sizeof(struct ID_List) * count);
	}

	num_alloc_blocks = (int *) ckalloc(sizeof(int));
	num_alloc_del_ops = (int *) ckalloc(sizeof(int));
	*num_alloc_blocks = count;

	initialize_algns(algns, 0, count);
  initialize_algns(init_algns, 0, count);

	init_ops(ops, 0, count);	
  for( i = 0; i < count; i++ ) { // initialization for ops
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
  *size1 = 0;
  *size2 = 0;

	read_maf(argv[1], maf_mode, algns, num_algns, size1, size2);
	read_maf(argv[1], maf_mode, init_algns, num_init_algns, size1, size2);

	if( debug_mode == TRUE ) {
		for(i = 0; i < (*num_algns); i++) {
			printf("%d: %d-%d, %d-%d\n", algns[i].sp_id, algns[i].x.lower, algns[i].x.upper, algns[i].y.lower, algns[i].y.upper);
		}
	}

	num_contigs1 = (int *) ckalloc(sizeof(int));
	num_contigs2 = (int *) ckalloc(sizeof(int));

	if( algn_type == PAIR ) { // for the combined alignment of two self-alignments and a pairwise alignment, the position of the second species is added to the length of the first sequence

		if( (*num_init_algns) > 0 ) {
			contigs1 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * 2 * (*num_init_algns));
			contigs2 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * 2 * (*num_init_algns));
      init_n_pair(contigs1, 0, (2*(*num_algns))-1);
      init_n_pair(contigs2, 0, (2*(*num_algns))-1);
		}
		else {
			contigs1 = (struct n_pair *) ckalloc(sizeof(struct n_pair));
			contigs2 = (struct n_pair *) ckalloc(sizeof(struct n_pair));
      init_n_pair(contigs1, 0, 1);
      init_n_pair(contigs2, 0, 1);
		}

		*num_contigs1 = 0;
		*num_contigs2 = 0;
		adjust_multi_contig_pos(algns, *num_algns, size1, size2, contigs1, num_contigs1, contigs2, num_contigs2);
		adjust_algn_pos(init_algns, *num_init_algns, contigs1, *num_contigs1, size1, contigs2, *num_contigs2, size2, CTG_NOT_ASSIGNED);
		switch_xy_coordinates(algns, *num_algns);
		switch_xy_coordinates(init_algns, *num_init_algns);

		if( (*num_contigs1) > 0 ) len_sum1 = (int *) ckalloc(sizeof(int) * (*num_contigs1));
		else {
			len_sum1 = (int *) ckalloc(sizeof(int));
			strcpy(contigs1[0].name1, species);
			strcpy(contigs1[0].name2, species);
			contigs1[0].len = *size1;
			contigs1[0].id = 0;
			*num_contigs1 = 1;
//			fatalf("error: number of contigs - %d - must be at least 1\n", *num_contigs1); 
		}		

		if( (*num_contigs2) > 0 ) len_sum2 = (int *) ckalloc(sizeof(int) * (*num_contigs2));
		else {
			len_sum2 = (int *) ckalloc(sizeof(int));
			strcpy(contigs2[0].name1, species2);
			strcpy(contigs2[0].name2, species2);
			contigs2[0].len = *size2;
			contigs2[0].id = 0;
			*num_contigs2 = 1;
//			fatalf("error: number of contigs - %d - must be at least 1\n", *num_contigs1); 
		}

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

	if( is_gc_given == TRUE ) {
		f = ckopen(argv[5], "r");
    num_alloc_cv = 0;
    while( fgets(buf, 1000, f) ) if( (buf[0] != '#') && (buf[0] != '(') ) num_alloc_cv++;

    cv = (struct cv_list *) ckalloc(num_alloc_cv * sizeof(struct cv_list));
    init_cv = (struct cv_list *) ckalloc(num_alloc_cv * sizeof(struct cv_list));
    init_dup_cv = (struct cv_list *) ckalloc(num_alloc_cv * sizeof(struct cv_list));

		init_conv(cv, 0, num_alloc_cv-1);
		init_conv(init_cv, 0, num_alloc_cv-1);
		init_conv(init_dup_cv, 0, num_alloc_cv-1);

    fseek(f, 0, SEEK_SET);
    i = 0;
    while( fgets(buf, 1000, f) ) {
			if( buf[0] == '(' ) {
				strcpy(tree_line, buf);
				leave_only_taxons(tree_line);
				num_sp_code = 0;
				j = 0;
				while(tree_line[j] != '\0') {
					if( tree_line[j] == ',') num_sp_code++;
					j++;
				}
				num_sp_code++;
				sp_code = (struct sp_list *) ckalloc(sizeof(struct sp_list) * num_sp_code);
				for( j = 0; j < num_sp_code; j++ ) {
					strcpy(sp_code[j].name, "");
					sp_code[j].id = -1;
				}

				temp_num_sp = assign_sp_code(tree_line, sp_code, num_sp_code);
				if( temp_num_sp != num_sp_code ) {
					fatalf("Number of species mismatched in %d and %d", num_sp_code, temp_num_sp);
				}

				for( j = 0; j < num_sp_code; j++ ) {
					if( strcmp(sp_code[j].name, species) == 0 ) sp1_code = sp_code[j].id;
					if( strcmp(sp_code[j].name, species2) == 0 ) sp2_code = sp_code[j].id;
				}
				sp_tree = read_one_line(tree_line); //// assgin a number for each species in all.gc
				assign_sp_id(sp_tree, sp_code, 0, num_sp_code);
			}
			else if( buf[0] == '#' ) {}
			else {
      	sscanf(buf, "%*s %s %d %d %s %d %d %c %*s %*s %*s %*s %d %d %d %d %d %s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d", temp_name1, &init_cv[i].s1, &init_cv[i].s2, temp_name2, &init_cv[i].t1, &init_cv[i].t2, &init_cv[i].ori, &init_cv[i].a1, &init_cv[i].a2, &init_cv[i].b1, &init_cv[i].b2, &init_cv[i].dir, temp_name3, &init_cv[i].stat);
				concat_ctg_name(temp_name1, sp_name1, ctg_name1);
				concat_ctg_name(temp_name2, sp_name2, ctg_name2);
				strcpy(out_species, temp_name3);

				if( strcmp(sp_name1, sp_name2) != 0 ) {
					fatalf("Species different in a conversion event: %s - %s in ortho_map.c \n", sp_name1, sp_name2);
				}

				for( j = 0; j < num_sp_code; j++ ) {
					if( strcmp(sp_code[j].name, out_species) == 0 ) out_code = sp_code[j].id;
					if( strcmp(sp_code[j].name, sp_name1) == 0 ) temp_code = sp_code[j].id;
				}
				
				if( temp_code == sp1_code ) {
					init_cv[i].ctg_id1 = get_ctg_id(sp_name1, ctg_name1, contigs1, *num_contigs1);
					init_cv[i].ctg_id2 = get_ctg_id(sp_name2, ctg_name2, contigs1, *num_contigs1);
					adjust_pos_conv(init_cv, i, len_sum1, *num_contigs1);
				}
				else if( temp_code == sp2_code ) {
					init_cv[i].ctg_id1 = get_ctg_id(sp_name1, ctg_name1, contigs2, *num_contigs2);
					init_cv[i].ctg_id2 = get_ctg_id(sp_name2, ctg_name2, contigs2, *num_contigs2);
					adjust_pos_conv(init_cv, i, len_sum2, *num_contigs2);
				}

				init_cv[i].sp_code = temp_code;
				init_cv[i].out_code = out_code;
				src = assign_I(init_cv[i].a1, init_cv[i].a2);
				dst = assign_I(init_cv[i].b1, init_cv[i].b2);
				if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((init_cv[i].s2-init_cv[i].s1+1) >= 50) && ((init_cv[i].a2-init_cv[i].a1+1) >= 30) && (init_cv[i].stat != 4) && (!((proper_overlap(src, dst) == true) && (width(intersect(src, dst)) > 100)))) 
				{
					init_cv[i].s2 = init_cv[i].s2 + 1; // to make [a,b) interval
					init_cv[i].t2 = init_cv[i].t2 + 1; // to make [a,b) interval
					init_cv[i].a2 = init_cv[i].a2 + 1; // to make [a,b) interval
					init_cv[i].b2 = init_cv[i].b2 + 1; // to make [a,b) interval
					init_cv[i].fid = i;
					if( temp_code == sp1_code ) {
						init_cv[i].sp_id = SELF1;
					}
					else if( temp_code == sp2_code ) {
						init_cv[i].sp_id = SELF2;
					}
					else {
						fatalf("%d: species code not supported\n", temp_code);
					}
   		   	i++;
					if( i > num_alloc_cv ) {
						fatalf("cv allocation is short %d\n", num_alloc_cv);
					}
				}
				else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((init_cv[i].s2-init_cv[i].s1+1) >= 50) && ((init_cv[i].a2-init_cv[i].a1+1) >= 30) && (init_cv[i].stat == 4) && (!((proper_overlap(src, dst) == true) && (width(intersect(src, dst)) > 100)))) 
				{
					init_dup_cv[num_dup_cv] = assign_conv(init_cv[i]);
					init_dup_cv[num_dup_cv].s2 = init_cv[i].s2 + 1; // to make [a,b) interval
					init_dup_cv[num_dup_cv].t2 = init_cv[i].t2 + 1; // to make [a,b) interval
					init_dup_cv[num_dup_cv].a2 = init_cv[i].a2 + 1; // to make [a,b) interval
					init_dup_cv[num_dup_cv].b2 = init_cv[i].b2 + 1; // to make [a,b) interval
					init_dup_cv[num_dup_cv].fid = num_dup_cv;
					if( temp_code == sp1_code ) {
						init_dup_cv[num_dup_cv].sp_id = SELF1;
					}
					else if( temp_code == sp2_code ) {
						init_dup_cv[num_dup_cv].sp_id = SELF2;
					}
					else {
						fatalf("%d: species code not supported\n", temp_code);
					}
   		   	num_dup_cv++;

					init_cv[i].s2 = init_cv[i].s2 + 1; // to make [a,b) interval
					init_cv[i].t2 = init_cv[i].t2 + 1; // to make [a,b) interval
					init_cv[i].a2 = init_cv[i].a2 + 1; // to make [a,b) interval
					init_cv[i].b2 = init_cv[i].b2 + 1; // to make [a,b) interval
					init_cv[i].fid = i;
					if( temp_code == sp1_code ) {
						init_cv[i].sp_id = SELF1;
					}
					else if( temp_code == sp2_code ) {
						init_cv[i].sp_id = SELF2;
					}
					else {
						fatalf("%d: species code not supported\n", temp_code);
					}
					i++;
					if( i > num_alloc_cv ) {
						fatalf("cv allocation is short %d\n", num_alloc_cv);
					}
				}
			}
    }
		fclose(f);
		num_init_cv = i;
		*num_cv = num_init_cv;

		old_dups = (int *) ckalloc(sizeof(int) * num_dup_cv);
		for( i = 0; i < num_dup_cv; i++ ) {
			old_dups[i] = -1;
		}

		if( is_exons_given == TRUE ) {
    	num_init_exons = 0;
			num_init_genes = 0;
			f = ckopen(argv[6], "r");
 	 	  while( fgets(buf, 1000, f) ) {
				if( (buf[0] != '#') && (buf[0] != '<') && (buf[0] != '>') ) num_init_exons++;
				if( (buf[0] == '<') || (buf[0] == '>') ) num_init_genes++;
			}
			fclose(f);
			*num_exons = num_init_exons;
			*num_genes = num_init_genes;

    	exons = (struct exons_list *) ckalloc(num_init_exons * sizeof(struct exons_list));
    	genes = (struct exons_list *) ckalloc(num_init_genes * sizeof(struct exons_list));
    	init_exons = (struct exons_list *) ckalloc(num_init_exons * sizeof(struct exons_list));
    	skip_reg1 = (struct exons_list *) ckalloc((num_init_exons + (2*num_init_cv))  * sizeof(struct exons_list));
    	skip_reg2 = (struct exons_list *) ckalloc((num_init_exons + (2*num_init_cv))  * sizeof(struct exons_list));
			initialize_exons_list(exons, 0, num_init_exons);
			initialize_exons_list(genes, 0, num_init_genes);
			initialize_exons_list(init_exons, 0, num_init_exons);
			initialize_exons_list(skip_reg1, 0, num_init_exons + (2*num_init_cv));
			initialize_exons_list(skip_reg2, 0, num_init_exons + (2*num_init_cv));

			read_exons(init_exons, exons, num_exons, genes, num_genes, skip_reg1, &num_skip1, skip_reg2, &num_skip2, argv[6], sp_code, num_sp_code, sp1_code, sp2_code, contigs1, *num_contigs1, len_sum1, contigs2, *num_contigs2, len_sum2);
			num_init_genes = *num_genes;
			num_init_exons = *num_exons;
		}
		else {
    	exons = (struct exons_list *) ckalloc(sizeof(struct exons_list));
    	genes = (struct exons_list *) ckalloc(sizeof(struct exons_list));
    	init_exons = (struct exons_list *) ckalloc(sizeof(struct exons_list));
    	skip_reg1 = (struct exons_list *) ckalloc((2*num_init_cv)  * sizeof(struct exons_list));
    	skip_reg2 = (struct exons_list *) ckalloc((2*num_init_cv)  * sizeof(struct exons_list));
			initialize_exons_list(exons, 0, 1);
			initialize_exons_list(genes, 0, 1);
			initialize_exons_list(init_exons, 0, 1);
			initialize_exons_list(skip_reg1, 0, 2*num_init_cv);
			initialize_exons_list(skip_reg2, 0, 2*num_init_cv);
			num_init_genes = 0;
			num_init_exons = 0;
			*num_genes = 0;
			*num_exons = 0;
		}

		if( is_ancestral == TRUE ) {
			f = ckopen(argv[7], "r");
			while( fgets(buf, 1000, f) ) {
				sscanf(buf, "%s %*s", op_name);
				if( (strcmp(op_name, "dup") == 0 ) || (strcmp(op_name, "del") == 0) ) num_prev_ops++;
			}
			
			prev_ops = (struct ops_list *) ckalloc(sizeof(struct ops_list) * num_prev_ops);
			fseek(f, 0, SEEK_SET);
			i = num_prev_ops-1;
			while( fgets(buf, 1000, f) ) {
				sscanf(buf, "%s %*s", op_name);
				if( strcmp(op_name, "dup") == 0 ) {
					sscanf(buf, "%s %c %d %d %d %d", op_name, &ori, &c1, &c2, &c3, &c4);
					prev_ops[i].sign = ori;
					prev_ops[i].id = i;
					prev_ops[i].src_b = c1;
					prev_ops[i].src_e = c2;
					prev_ops[i].dst_b = c3;
					prev_ops[i].dst_e = c3+c4;
					prev_ops[i].sp_id = SELF1;
					i--;
				}
				else if( strcmp(op_name, "del") == 0 ) {
					sscanf(buf, "%s %*s %d %d", op_name, &c1, &c2);
					prev_ops[i].sign = 'd';
					prev_ops[i].id = i;
					prev_ops[i].src_b = c1;
					prev_ops[i].src_e = c2;
					i--;
				}
			}
			fclose(f);
		}
	
		for( i = 0; i < *num_exons; i++ ) {
/*
			ctg_id = exons[i].ctg_id;
			if( (*num_contigs1) > 1 ) {
				if( ctg_id != -1 ) {
					if( exons[i].sp_id == SELF1 ) {
						if( ctg_id > (*num_contigs1) ) {
							exons[i].reg = assign_I(exons[i].reg.lower + len_sum1[ctg_id], exons[i].reg.upper + len_sum1[ctg_id]);
							init_exons[i].reg = assign_I(init_exons[i].reg.lower + len_sum1[ctg_id], init_exons[i].reg.upper + len_sum1[ctg_id]);
							genes[i].reg = assign_I(genes[i].reg.lower + len_sum1[ctg_id], genes[i].reg.upper + len_sum1[ctg_id]);
						}
					}
				}
			}

			if( (*num_contigs2) > 1 ) {
				if( ctg_id != -1 ) {
					if( exons[i].sp_id == SELF2 ) {
						if( ctg_id > (*num_contigs1) ) {
							exons[i].reg = assign_I(exons[i].reg.lower + len_sum2[ctg_id], exons[i].reg.upper + len_sum2[ctg_id]);
							init_exons[i].reg = assign_I(init_exons[i].reg.lower + len_sum1[ctg_id], init_exons[i].reg.upper + len_sum1[ctg_id]);
							genes[i].reg = assign_I(genes[i].reg.lower + len_sum1[ctg_id], genes[i].reg.upper + len_sum1[ctg_id]);
						}
					}
				}
			}
*/
			if( exons[i].sp_id == SELF2 ) exons[i].reg = assign_I(exons[i].reg.lower + (*size1), exons[i].reg.upper + (*size1));
		}	

		for( i = 0; i < *num_genes; i++ ) {
			if( genes[i].sp_id == SELF2 ) genes[i].reg = assign_I(genes[i].reg.lower + (*size1), genes[i].reg.upper + (*size1));
		}
	}

	f = ckopen(argv[1], "r");

	free(num_suspend_pairs);
	free(suspend_list);
	free(temp_algns);

	if( is_gc_given == TRUE ) {
		j = 0;
		for( i = 0; i < (*num_init_algns); i++ ) {
			if( init_algns[i].sp_id != PAIR ) {
				init_algns[i].xl_diff = 0;
				init_algns[i].yl_diff = 0;
				init_algns[i].xr_diff = 0;
				init_algns[i].yr_diff = 0;
				init_algns[i].xl_offset = 0;
				init_algns[i].yl_offset = 0;
				init_algns[i].xr_offset = 0;
				init_algns[i].yr_offset = 0;
				init_algns[i].c_id = -1;
				init_algns[i].m_id = -1;
				init_algns[i].sign = init_algns[i].init_sign;
				assign_algn(algns, j, init_algns[i]);
				j++;
			} 
			else if( init_algns[i].sign != DELETED ) {
				algn_len_x = init_algns[i].x.upper - init_algns[i].xr_diff - init_algns[i].x.lower - init_algns[i].xl_diff;
				algn_len_y = init_algns[i].y.upper - init_algns[i].yr_diff - init_algns[i].y.lower - init_algns[i].yl_diff;
				if( (algn_len_x > MIN_LEN_FOR_RECAL_PID) && (algn_len_y > MIN_LEN_FOR_RECAL_PID) ) 
				{ 
					assign_algn(algns, j, init_algns[i]);
					j++;
				}
			}
		}
		*num_algns = j;

		cal_pid_conv(algns, *num_algns, init_cv, num_init_cv, f); // init_cv is sorted in this function call

		j = 0;
		for( i = 0; i < num_init_cv; i++ ) {
			if( ((int)init_cv[i].conv_pid) != -1 ) {
				init_cv[j] = assign_conv(init_cv[i]);
				j++;	
			}
		}
		num_init_cv = j;
		*num_cv = num_init_cv;

		for( i = 0; i < num_init_cv; i++ ) {
//			if( init_cv[i].stat == 4 ) {} // conversion covering the entire paralogous alignment
			if( init_cv[i].sp_id == SELF1 ) {
				skip_reg1[num_skip1].reg = assign_I(init_cv[i].a1, init_cv[i].a2);
				skip_reg1[num_skip1].val = (int) (init_cv[i].conv_pid + 0.5);
				skip_reg1[num_skip1].ctg_id = init_cv[i].ctg_id1;
				num_skip1++;
				skip_reg1[num_skip1].reg = assign_I(init_cv[i].b1, init_cv[i].b2);
				skip_reg1[num_skip1].val = (int) (init_cv[i].conv_pid + 0.5);
				skip_reg1[num_skip1].ctg_id = init_cv[i].ctg_id2;
				num_skip1++;
			}
			else if( init_cv[i].sp_id == SELF2 ) {
				skip_reg2[num_skip2].reg = assign_I(init_cv[i].a1, init_cv[i].a2);
				skip_reg2[num_skip2].val = (int) (init_cv[i].conv_pid + 0.5);
				skip_reg1[num_skip2].ctg_id = init_cv[i].ctg_id1;
				num_skip2++;
				skip_reg2[num_skip2].reg = assign_I(init_cv[i].b1, init_cv[i].b2);
				skip_reg2[num_skip2].val = (int) (init_cv[i].conv_pid + 0.5);
				skip_reg1[num_skip2].ctg_id = init_cv[i].ctg_id2;
				num_skip2++;
			}

			cv[i] = assign_conv(init_cv[i]); // sorted by similarity levels
//			if( cv[i].sp_id == SELF1 ) adjust_pos_conv(cv, i, len_sum1, *num_contigs1);

			if( cv[i].sp_id == SELF2 ) {
//				adjust_pos_conv(cv, i, len_sum2, *num_contigs2);
				cv[i].s1 = cv[i].s1 + (*size1);
				cv[i].s2 = cv[i].s2 + (*size1);
				cv[i].t1 = cv[i].t1 + (*size1);
				cv[i].t2 = cv[i].t2 + (*size1);
				cv[i].a1 = cv[i].a1 + (*size1);
				cv[i].a2 = cv[i].a2 + (*size1);
				cv[i].b1 = cv[i].b1 + (*size1);
				cv[i].b2 = cv[i].b2 + (*size1);
			}
		}

		sort_exons(skip_reg1, num_skip1);
		sort_exons(skip_reg2, num_skip2);

		// cv is a list of printing results, so any elements must not be deleted

		if( algn_type == PAIR ) { // for the combined alignment of two self-alignments and a pairwise alignment, the position of the second species is added to the length of the first sequence
			total_len = 0;
			avg_pid = (float) 0;
			for(i = 0; i < (*num_algns); i++) {
				if( algns[i].sp_id == SELF2 ) {
					algns[i].x = assign_I(algns[i].x.lower+(*size1), algns[i].x.upper+(*size1));
					algns[i].y = assign_I(algns[i].y.lower+(*size1), algns[i].y.upper+(*size1));
				}
				else if( algns[i].sp_id == PAIR ) {
					algns[i].y = assign_I(algns[i].y.lower+(*size1), algns[i].y.upper+(*size1));
				}

				cur_pid1 = -1;
				if( algns[i].sp_id != PAIR ) {
					cur_pid1 = pick_sim_level_algn(init_algns, algns[i].index, f, skip_reg1, num_skip1, skip_reg2, num_skip2, algns[i].sp_id); 
				}

				if( cur_pid1 != -1 ) {  
					algns[i].identity = cur_pid1;
				}
				else { // converted regions cover the entire alignment
				}

				if((algns[i].sp_id == PAIR) && (algns[i].sign != DELETED)) {
					total_len = total_len + width(algns[i].x);
					avg_pid = avg_pid + (float)(((float)algns[i].identity/(float)100)	 * width(algns[i].x));
				}
			}
			avg_pid = (avg_pid/(float)total_len) * (float)100;
		}

		for( i = 0; i < (*num_algns); i++ ) {
			if( algns[i].identity <= (((int)(avg_pid+0.5) - PID_DIFF)) ) {
				algns[i].sign = DELETED;
				if( init_algns[algns[i].index].sign != DELETED ) init_algns[algns[i].index].sign = DELETED;
			}
		}
		overwrite_dots(num_algns, algns);

		num_old_dups = make_old_dup_list(init_algns, *num_init_algns, init_dup_cv, num_dup_cv, old_dups);

		for( i = 0; i < num_old_dups; i++ ) {
			j = old_dups[i];			
			init_algns[j].lock = OLD_DUP;
		}

		if( is_ancestral == TRUE ) {
			redo_ops(prev_ops, num_prev_ops, num_algns, algns, cv, num_cv, size, init_algns, *num_init_algns, f, *size1);
			free(prev_ops);
		}
		
		if( *num_algns > 0 ) {
			sorted1 = (struct slist *) ckalloc(sizeof(struct slist) * (*num_algns));
			initialize_slist(sorted1, 0, (*num_algns));
			for( i = 0; i < (*num_algns); i++ ) sorted1[i].id = i;
			sort_by_pid(sorted1, algns, *num_algns);
			if( debug_mode == TRUE ) {
				for( i = 0; i < (*num_algns); i++ ) {
					printf("%d-%d, %d-%d: %d, %d\n", algns[sorted1[i].id].x.lower, algns[sorted1[i].id].x.upper, algns[sorted1[i].id].y.lower, algns[sorted1[i].id].y.upper, algns[sorted1[i].id].sp_id, algns[sorted1[i].id].identity);
				}
	
				for( i = 0; i < num_init_cv; i++ ) {
					printf("%d-%d, %d-%d: %f\n", init_cv[i].a1, init_cv[i].a2, init_cv[i].b1, init_cv[i].b2, init_cv[i].conv_pid);
				}
			}
	
			if( (run_mode == CONTENT_ORTHO) || (run_mode == POSITION_ORTHO) ) {
				num_cur_cv = num_init_cv;
				if( num_dup_cv > 0 ) {
					if( (num_cur_cv+num_dup_cv) > num_alloc_cv ) {
						cv = ckrealloc(cv, sizeof(struct cv_list) * (num_cur_cv + num_dup_cv));
					}

					for( i = 0; i < num_dup_cv; i++ ) {
						cv[num_cur_cv+i] = assign_conv(init_dup_cv[i]);
						if( init_dup_cv[i].sp_id == SELF2 ) {
							cv[num_cur_cv+i].s1 = cv[num_cur_cv+i].s1 + (*size1);        
							cv[num_cur_cv+i].s2 = cv[num_cur_cv+i].s2 + (*size1);        
							cv[num_cur_cv+i].t1 = cv[num_cur_cv+i].t1 + (*size1);        
							cv[num_cur_cv+i].t2 = cv[num_cur_cv+i].t2 + (*size1);        
							cv[num_cur_cv+i].a1 = cv[num_cur_cv+i].a1 + (*size1);        
							cv[num_cur_cv+i].a2 = cv[num_cur_cv+i].a2 + (*size1);        
							cv[num_cur_cv+i].b1 = cv[num_cur_cv+i].b1 + (*size1);        
							cv[num_cur_cv+i].b2 = cv[num_cur_cv+i].b2 + (*size1);
						}
					}
					num_cur_cv = num_cur_cv + num_dup_cv;
				}
				*num_cv = num_cur_cv;
			}

			reorder_dups(sorted1, num_algns, algns, size, STRICT, init_algns, *num_init_algns, f, avg_pid, cv, num_cv, ops, num_ops, exons, num_exons, genes, num_genes, old_dups, num_old_dups, run_mode, *size1); 
			free(sorted1);
		}
		ops[*num_ops].sign = 's';
		*num_ops = *num_ops + 1;
	} 

	scaling_value = (float)1;
 	if( run_mode == INF_DUP) output_ops(*num_ops, ops, *size2, scaling_value, -1);

	if( (debug_mode == TRUE) || (run_mode == INF_DUP) ) {
		printf("-------------\n");
		printf("For debugging\n");
		for( i = ((*num_ops)-2); i >= 0; i-- ) {
			printf("%d-%d is copied from %d-%d: %d\n", ops[i].dstStart, ops[i].dstEnd, ops[i].srcStart, ops[i].srcEnd, ops[i].sp_id);
		}
	}

	free_p_tree(sp_tree);
	free(num_exons);
	free(num_cv);
	free(num_genes);
	free(threshold);
	free(num_id);
	free(is_x);
	free(sp_order);
	free(rm_sp);
	free(left_sp);

	if( is_gc_given == TRUE ) {
		free(sp_code);
		free(skip_reg1);
		free(skip_reg2);
		free(old_dups);
		free(init_dup_cv);
		free(init_cv);
		free(exons);
		free(init_exons);
		free(genes);

		if( (*num_init_algns) > 0 ) {
			update_pid_init_algns(*num_init_algns, init_algns, f, avg_pid);
			const_graph(*num_init_algns, init_algns, f);
//			const_graph(*num_init_algns, init_algns, f, ops, *num_ops, *size1);
//			map_one_to_one(*num_init_algns, init_algns, f);
		}

		num_ortho_algns = 0;

		total_len = 0;
		avg_pid = (float)0;
		for( i = 0; i < (*num_init_algns); i++ ) {
			if( (init_algns[i].sign != DELETED) && (init_algns[i].sp_id == PAIR ) ) {
				num_ortho_algns++;
				total_len = total_len + width(init_algns[i].x);
				avg_pid = avg_pid + (float)(((float)init_algns[i].identity/(float)100) * width(init_algns[i].x));
			}
		}

		if( total_len > 0 ) {
			avg_pid = (avg_pid/(float)total_len) * (float)100;
			temp_val = (float *) ckalloc(sizeof(float));
			update_pid_init_algns(*num_init_algns, init_algns, f, avg_pid);

		  initialize_algns(algns, 0, count);
			*num_algns = assign_ortho_algns(algns, init_algns, (*num_init_algns), temp_val);
			free(temp_val);
			sorted1 = (struct slist *) ckalloc((*num_algns) * sizeof(struct slist));
			sorted2 = (struct slist *) ckalloc((*num_algns) * sizeof(struct slist));
			initialize_slist(sorted1, 0, (*num_algns));
			initialize_slist(sorted2, 0, (*num_algns));
    	sort_init_algns(sorted1, algns, *num_algns, INIT_SELF1);
    	sort_init_algns(sorted2, algns, *num_algns, INIT_SELF2);

			ops_id = *num_ops;
			if( run_mode == CONTENT_ORTHO ) {
				for( i = 0; i < num_cur_cv; i++ ) {
					if( (cv[i].ori != 'd') && (cv[i].conv_pid > (avg_pid+2)) ) {
			      if( cv[i].ori == '+' ) ops[ops_id].sign = 'c';
			      else if( cv[i].ori == '-' ) ops[ops_id].sign = 'v';
						else {
							fatalf("unsupported ops type: %c\n", cv[i].ori);
						}

						src = assign_I(cv[i].a1, cv[i].a2);
						dst = assign_I(cv[i].b1, cv[i].b2);
						cur_pid1 = (float) 0;
						cur_pid2 = (float) 0;

						cur_pid1 = cal_ortho_algn_pid(src, cv[i].sp_id, sorted1, algns, *num_algns, f);	
						cur_pid2 = cal_ortho_algn_pid(dst, cv[i].sp_id, sorted2, algns, *num_algns, f);	

						if( (cv[i].conv_pid >= cur_pid1) && (cv[i].conv_pid >= cur_pid2) ) {
  				    ops[ops_id].dir = cv[i].dir;
  			 	  	ops[ops_id].sp_id = cv[i].sp_id;
							ops[ops_id].pid = cv[i].conv_pid;
							ops[ops_id].id = cv[i].fid;
   				   if( cv[i].dir == 1 )  {
    				    ops[ops_id].src_b = cv[i].b1;
    				    ops[ops_id].srcStart = cv[i].b1;
      				  ops[ops_id].src_e = cv[i].b2;
      				  ops[ops_id].srcEnd = cv[i].b2;
								ops[ops_id].ctg_id1 = cv[i].ctg_id2;
      				  ops[ops_id].dst_b = cv[i].a1;
      				  ops[ops_id].dstStart = cv[i].a1;
       				  ops[ops_id].dst_e = cv[i].a2;
       				  ops[ops_id].dstEnd = cv[i].a2;
								ops[ops_id].ctg_id2 = cv[i].ctg_id1;
								ops[ops_id].dir = 1;
    				  }
      				else {
      				  ops[ops_id].src_b = cv[i].a1;
      				  ops[ops_id].srcStart = cv[i].a1;
       				  ops[ops_id].src_e = cv[i].a2;
       				  ops[ops_id].srcEnd = cv[i].a2;
								ops[ops_id].ctg_id1 = cv[i].ctg_id1;
       				  ops[ops_id].dst_b = cv[i].b1;
       				  ops[ops_id].dstStart = cv[i].b1;
       		 			ops[ops_id].dst_e = cv[i].b2;
       		 			ops[ops_id].dstEnd = cv[i].b2;
								ops[ops_id].ctg_id2 = cv[i].ctg_id2;
								ops[ops_id].dir = 2;
     				 }
							cv[i].ori = 'd';
							ops_id++;
						}
					}
				}
			}
			free(sorted1);
			free(sorted2);
			*num_ops = ops_id;
		}

		free(cv);
		free(algns);
		free(num_algns);

		if( num_ortho_algns > ALLOC_UNIT ) {
			*num_alloc_del_ops = num_ortho_algns;
		}
		else {
			*num_alloc_del_ops = ALLOC_UNIT;
		}
		del_ops = (struct ops_list *) ckalloc(sizeof(struct ops_list) * (*num_alloc_del_ops));
		init_ops(del_ops, 0, (*num_alloc_del_ops));
		num_del_ops = (int *) ckalloc(sizeof(int));
		*num_del_ops = 0;
		if( (*num_ops) > 0 ) {
			ops_cur_pos = (struct ops_list *) ckalloc(sizeof(struct ops_list) * (*num_ops));
			num_dup_ops = cal_cur_pos_ops(*num_ops, ops, ops_cur_pos, SELF1, 0);
//			iden_gene_loss(*num_init_algns, init_algns, num_dup_ops, ops_cur_pos, f, &del_ops, num_del_ops, num_alloc_del_ops);

			if( (run_mode == MANY_TO_MANY) || (run_mode == CONTENT_ORTHO) || (run_mode == POSITION_ORTHO)) {
				if( is_event_output_required == true ) {
					g = ckopen(argv[7], "w");
					fprintf(g, "# %s %s\n", species, species2);
					for( i = 0; i < num_dup_ops; i++ ) {
						strcpy(temp_name1, species);
						strcpy(temp_name2, species);
						len1 = 0;
						len2 = 0;

						ctg_id = ops_cur_pos[i].ctg_id1;
						if( (ctg_id < 0) || (ctg_id >= (*num_contigs1)) ) {
							fatalf("contig id not assigned %d in ortho_map.c\n", ctg_id);
						}
						else if( ctg_id != -1 ) {
							if( strcmp(species, contigs1[ctg_id].name2) != 0 ) {
								sprintf(temp_name1, "%s.%s", species, contigs1[ctg_id].name2);
							}
							len1 = len_sum1[ctg_id];
						}

						if( ops_cur_pos[i].sign != 'd' ) {
							ctg_id = ops_cur_pos[i].ctg_id2;
							if( (ctg_id < 0) || (ctg_id >= (*num_contigs1)) ) {
								fatalf("contig id not assigned %d in ortho_map.c\n", ctg_id);
							}
							else if( ctg_id != -1 ) {
								if( strcmp(species, contigs1[ctg_id].name2) != 0 ) {
									sprintf(temp_name2, "%s.%s", species, contigs1[ctg_id].name2);
								}
								len2 = len_sum1[ctg_id];
							}
						}

						if( ((ops_cur_pos[i].srcStart - len1) < 0) || ((ops_cur_pos[i].dstStart-len2) < 0) ) {
							fatalf("unexpected negative values: %s %d %d %s %d %d\n", temp_name1, ops_cur_pos[i].srcStart-len1, ops_cur_pos[i].srcEnd-len1, temp_name2, ops_cur_pos[i].dstStart-len2, ops_cur_pos[i].dstEnd-len2);
						}
						else {
							fprintf(g, "%c %s %d %d %s %d %d 0 %d %f\n", ops_cur_pos[i].sign, temp_name1, ops_cur_pos[i].srcStart-len1, ops_cur_pos[i].srcEnd-len1, temp_name2, ops_cur_pos[i].dstStart-len2, ops_cur_pos[i].dstEnd-len2, ops_cur_pos[i].dir, ops_cur_pos[i].pid);
						}
					}
				}

				*num_init_algns = redo_dups_for_mtom_inc_conv(num_dup_ops, ops_cur_pos, *num_init_algns, &init_algns, f, REF_SEQ, num_alloc_blocks);
			}
		}
		else {
//			iden_gene_loss(*num_init_algns, init_algns, num_dup_ops, ops_cur_pos, f, &del_ops, num_del_ops, num_alloc_del_ops);
		}

		if( ((*num_ops) > 0) && ((run_mode == MANY_TO_MANY) || (run_mode == ONE_TO_MANY) || (run_mode == CONTENT_ORTHO) || (run_mode == POSITION_ORTHO))) {
			num_dup_ops = cal_cur_pos_ops(*num_ops, ops, ops_cur_pos, SELF2, *size1);

			if( is_event_output_required == true ) {
//			fprintf(g, "## %s\n", species2);
				for( i = 0; i < num_dup_ops; i++ ) {
					strcpy(temp_name1, species2);
					strcpy(temp_name2, species2);
					len1 = 0;
					len2 = 0;

					ctg_id = ops_cur_pos[i].ctg_id1;
					if( (ctg_id < 0) || (ctg_id >= (*num_contigs2)) ) {
						fatalf("contig id not assigned %d in ortho_map.c\n", ctg_id);
					}
					else if( (ctg_id < 0) || (ctg_id != -1) ) {
						if( strcmp(species2, contigs2[ctg_id].name2) != 0 ) {
							sprintf(temp_name1, "%s.%s", species2, contigs2[ctg_id].name2);
						}
						len1 = len_sum2[ctg_id];
					}

					if( ops_cur_pos[i].sign != 'd' ) {
						ctg_id = ops_cur_pos[i].ctg_id2;
						if( (ctg_id < 0) || (ctg_id >= (*num_contigs2)) ) {
							fatalf("contig id not assigned %d in ortho_map.c\n", ctg_id);
						}
						else if( ctg_id != -1 ) {
							if( strcmp(species2, contigs2[ctg_id].name2) != 0 ) {
								sprintf(temp_name2, "%s.%s", species2, contigs2[ctg_id].name2);
							}
							len2 = len_sum2[ctg_id];
						}
					}

					if( ((ops_cur_pos[i].srcStart - len1) < 0) || ((ops_cur_pos[i].dstStart-len2) < 0) ) {
						fatalf("unexpected negative values: %s %d %d %s %d %d\n", temp_name1, ops_cur_pos[i].srcStart-len1, ops_cur_pos[i].srcEnd-len1, temp_name2, ops_cur_pos[i].dstStart-len2, ops_cur_pos[i].dstEnd-len2);
					}
					else {
						fprintf(g, "%c %s %d %d %s %d %d 1 %d %f\n", ops_cur_pos[i].sign, temp_name1, ops_cur_pos[i].srcStart-len1, ops_cur_pos[i].srcEnd-len1, temp_name2, ops_cur_pos[i].dstStart-len2, ops_cur_pos[i].dstEnd-len2, ops_cur_pos[i].dir, ops_cur_pos[i].pid);
					}
				}
				fclose(g);
			}	
			
			if( num_dup_ops > 0 ) {
				*num_init_algns = redo_dups_for_mtom_inc_conv(num_dup_ops, ops_cur_pos, *num_init_algns, &init_algns, f, SELF2, num_alloc_blocks);
			}
		}
		if( (*num_ops) > 0 ) {
			free(ops_cur_pos);
		}

		free(del_ops);
		free(num_del_ops);

		if( run_mode != INF_DUP ) {
//			write_init_maf(stdout, init_algns, *num_init_algns, species, species2, *size1, *size2, f, PAIR);
			write_init_maf_stdout(init_algns, *num_init_algns, contigs1, contigs2, *size1, *size2, f, PAIR, species, species2);
		}
	}
	else {
		free(algns);
		free(num_algns);
	}

	fclose(f);

	if( is_contig_printed == true ) { // argc == 9
		output_contigs(argv[8], contigs1, len_sum1, *num_contigs1);
		output_contigs(argv[9], contigs2, len_sum2, *num_contigs2);
	}
	
	if( (*num_init_algns) > 0 ) {
		free(contigs1);
		free(contigs2);
	}

	free(len_sum1);
	free(len_sum2);
	free(num_contigs1);
	free(num_contigs2);
	free(num_alloc_blocks);
	free(num_alloc_del_ops);
	free(num_ops);
	free(ops);
	free(size2);
	free(size1);
	free(init_algns);
	free(num_init_algns);
	return EXIT_SUCCESS;
}

void output_contigs(char *fname, struct n_pair *contigs, int *old_len_sum, int num_contigs)
{
	FILE *f = NULL;
	int count = 0, i = 0, j = 0;
	struct n_pair *all_contigs = NULL;
	int num_all_contigs = 0;
	int *len_sum = NULL;

	f = ckopen(fname, "a+");
	count = count_lines(f);
	if( (count + num_contigs) > 0 ) {
	/*
		for( i = 0; i < num_contigs; i++ ) {
			contigs[i].id = contigs[i].len;
			contigs[i].len = old_len_sum[i];
		}
	*/
		all_contigs = (struct n_pair *) ckalloc((count + num_contigs) * sizeof(struct n_pair));
		len_sum = (int *) ckalloc( (count + num_contigs) * sizeof(int));
		num_all_contigs = read_one_contig_file(fname, f, all_contigs, count);

		for( i = 0; i < num_all_contigs; i++ ) {
			len_sum[i] = all_contigs[i].len;
			all_contigs[i].len = all_contigs[i].id;
		}

		j = num_all_contigs;
		for( i = 0; i < num_contigs; i++ ) {
			if( is_ctg_in(contigs[i].name1, contigs[i].name2, all_contigs, j) == -1 ) 
			{
				strcpy(all_contigs[j].name1, contigs[i].name1);		
				strcpy(all_contigs[j].name2, contigs[i].name2);		
				all_contigs[j].len = contigs[i].len;
				all_contigs[j].id = contigs[i].id;
				len_sum[j] = old_len_sum[i];
				j++;
			}
		}
		fclose(f);
	}
	print_contig_list(f, all_contigs, len_sum, num_all_contigs, j);

	if( (count + num_contigs) > 0 ) {
		free(all_contigs);
		free(len_sum);
	}
}

void output_ops(int num_ops, struct ops_list * ops, int len, float sc, int sp_id)
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
		if( ( sp_id == ops[i].sp_id ) || (sp_id == -1) ) {
			if( (op == '+') || (op == '-') ) printf("dup %c %d %d %d %d\n", op, t_a, t_b, t_c, t_d - t_c);
			else if( (op == 'd') || (op == 'l') ) printf("del . %d %d\n", t_a, t_b);
			else if( op == 'i') printf("inv . %d %d\n", t_a, t_b);
			else if( op == 'c') {
				if( ops[i].dir == 0 ) printf("con + %d %d %d %d ?\n", t_a, t_b, t_c, t_d);
				else printf("con + %d %d %d %d\n", t_a, t_b, t_c, t_d);
			}
			else if( op == 'v') {
				if( ops[i].dir == 0 ) printf("con - %d %d %d %d ?\n", t_a, t_b, t_c, t_d);
				else printf("con - %d %d %d %d\n", t_a, t_b, t_c, t_d);
			}
			else if( op == 's') printf("spe\n");
		}
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
