#include "main.h"
#include "adjust_ops.h"
#include "regions.h"
#include "util.h"
#include "util_i.h"
#include "util_input.h"
#include "util_algns.h"
#include "read_maf.h"
#include "tree_op.h"
#include "util_ops.h"
#include "util_gen.h"
#include "contigs_op.h"

int debug_mode;
int count_node; 
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

int main(int argc, char **argv)
{
	char buf[1000];
	char tree_line[1000];
	char name_line[1000];
	struct p_tree *sp_tree = NULL;
	struct p_tree *t = NULL, *out_t = NULL;
	int root_id = 0;
	int i = 0, j = 0;
	int *count = 0;
	int num_sp = 0, temp_num = 0;
	struct sp_list *species;
	FILE *f, *g;
	int **temp_list;
	int num_list = 0, ref_id = 0;
	struct ops_list **ops;
	int *num_ops;
	struct ops_list *cur_ops;
	int num_cur_ops = 0, prev_num = 0;
	struct DotList **ortho_algns;
	struct slist **sorted;
	int *num_algns;
	char filename[LEN_NAME];
	char name1[LEN_NAME], name2[LEN_NAME];
	int *size1, *size2;
	int num_total_ops = 0;
	int num_sp_event = 0;
	int index = 0;
	struct n_pair **contigs;
	int *num_contigs; 
	int ctg_id1 = -1, ctg_id2 = -1;
	int len1 = 0, len2 = 0;
	
	count_node = 0;
	debug_mode = FALSE;
	strcpy(name_line, "");
	strcpy(tree_line, "");
	if( argc == 7 ) {
		if( strcmp(argv[6], "debug_mode") == 0 ) debug_mode = TRUE;
		else {
			fatal("args: non-redundant.gc ref_sp events_combined_file ortho_maf_directory contigs_list (debug_mode) \n");
		}
	}
	else if( argc != 6 ) {
		fatal("args: non-redundant.gc ref_sp events_combined_file ortho_maf_directory contigs_list\n");
	}

	if((f = ckopen(argv[1], "r")) != NULL ) {
		if( fgets(buf, 1000, f) == NULL ) {
			fatalf("nothing in file %s\n", argv[1]);
		}
		else if( buf[0] == '(' ) {
			strcpy(tree_line, buf);
			leave_only_taxons(tree_line);
			num_sp = 0;
			j = 0;
			while( tree_line[j] != '\0') {
				if( tree_line[j] == ',' ) num_sp++;
				j++;
			}
			num_sp++;
			if( num_sp > 0 ) {
				species = (struct sp_list *) ckalloc(num_sp * sizeof(struct sp_list));
				temp_list = (int **) ckalloc(num_sp * sizeof(int *));
				for( j = 0; j < num_sp; j++ ) {
					temp_list[j] = (int *) ckalloc((num_sp+1) * sizeof(int));
				}
        for( j = 0; j < num_sp; j++ ) {
          strcpy(species[j].name, "");
          species[j].id = -1;
					temp_list[j][0] = 0;
					for( i = 1; i <= num_sp; i++ ) {
						temp_list[j][i] = -1;
					}
        }
        temp_num = assign_sp_code(tree_line, species, num_sp);
        for( j = 0; j < num_sp; j++ ) {
          species[j].id = j;
        }
        if( temp_num != num_sp ) {
          fatalf("Number of species mismatched in %d and %d", num_sp, temp_num);
        }

				sp_tree = read_one_line(tree_line);
				root_id = sp_tree->nid;
				assign_sp_id(sp_tree, species, 0, num_sp);
				if( debug_mode == TRUE ) {
					print_tree(sp_tree, TREE_PRINT);
					printf("\n");
				}
			}
		}
		else {
			fatalf("not newick format: %s", buf);
		}
		fclose(f);
	}
	else {
		fatalf("file open error %s\n", argv[1]);
	}

	ref_id = -1;
	for(i = 0; i < num_sp; i++) {
		if( strcmp(species[i].name, argv[2]) == 0 ) {
			ref_id = i;
		}
	}

	if( ref_id == -1 ) {
		fatalf("reference species %s not found\n", argv[2]);
	}
	
	if( (f = ckopen(argv[3], "r")) != NULL ) {}
	else {
		fatalf("ops file open error: %s\n", argv[3]);
	}

	if( (g = ckopen(argv[5], "r")) != NULL ) {}
	else {
		fatalf("contigs file open error: %s\n", argv[5]);
	}

	if( num_sp > 0 ) {
		num_ops = (int * ) ckalloc(num_sp * sizeof(int));	

		num_algns = (int * ) ckalloc(num_sp * sizeof(int));	
		ops = (struct ops_list **) ckalloc(num_sp * sizeof(struct ops_list *));	
		ortho_algns = (struct DotList **) ckalloc(num_sp * sizeof(struct DotList *));	
		num_contigs = (int *) ckalloc(num_sp * sizeof(int));
		contigs = (struct n_pair **) ckalloc(num_sp * sizeof(struct n_pair *));
		sorted = (struct slist **) ckalloc(num_sp * sizeof(struct slist *));	
		for(i = 0; i < num_sp; i++ ) {
			num_ops[i] = 0;
			if( i != ref_id ) {
				num_ops[i] = count_ops(species[i].name, f, SELF1); // events in the reference sequence (SELF1 is 0)
				if( num_ops[i] > 0 ) {
					ops[i] = (struct ops_list *) ckalloc(num_ops[i] * sizeof(struct ops_list));
				}
			}
			num_contigs[i] = 0;
			num_contigs[i] = count_contigs(species[i].name, g);
			if( num_contigs[i] > 0 ) {
					contigs[i] = (struct n_pair *) ckalloc(num_contigs[i] * sizeof(struct n_pair));
			}
			num_total_ops = num_total_ops + num_ops[i] + 1;
		}

		if( num_total_ops > 0 ) {
			cur_ops = (struct ops_list *) ckalloc(num_total_ops * sizeof(struct ops_list));
		}

		for(i = 0; i < num_sp; i++ ) {
			if( num_contigs[i] > 0 ) {
				read_contigs_file(species[i].name, g, contigs[i], num_contigs[i]);
			}
		}

		for(i = 0; i < num_sp; i++ ) {
			if( num_ops[i] > 0 ) {
				read_ops_file(species[i].name, f, ops[i], num_ops[i], SELF1, contigs[ref_id], num_contigs[ref_id], species[ref_id].name);
			}
		}
	}

	fclose(g);
	fclose(f);

	if( num_sp > 0 ) {
		size1 = (int *) ckalloc(sizeof(int));
		size2 = (int *) ckalloc(sizeof(int));
		count = (int *) ckalloc(sizeof(int));
		for(i = 0; i < num_sp; i++ ) {
			num_algns[i] = 0;
			if( i != ref_id ) {
				sprintf(filename, "%s/%s.%s.maf", argv[4], species[ref_id].name, species[i].name);
				if( (g = ckopen(filename, "r")) != NULL ) {
					fclose(g);
				}
				else {
					fatalf("maf file open error: %s\n", argv[4]);
				}
				num_algns[i] = count_local_algns(filename, name1, name2);

				if( (strcmp(species[ref_id].name, name1) != 0) || (strcmp(species[i].name, name2) != 0) ) {
					fatalf("species names %s and %s mismatch %s and %s\n", species[ref_id].name, species[i].name, name1, name2);
				}
 
				if( num_algns[i] > 0 ) {
					ortho_algns[i] = (struct DotList *) ckalloc(num_algns[i] * sizeof(struct DotList));
					sorted[i] = (struct slist *) ckalloc(num_algns[i] * sizeof(struct slist));
					for( j = 0; j < num_algns[i]; j++ ) {
						initialize_slist(sorted[i], 0, num_algns[i]-1);
					}
					initialize_algns(ortho_algns[i], 0, num_algns[i]);
					read_maf(filename, O_MODE, ortho_algns[i], count, size1, size2);
					adjust_algn_pos(ortho_algns[i], *count, contigs[ref_id], num_contigs[ref_id], size1, contigs[i], num_contigs[i], size2, CTG_NOT_ASSIGNED_BUT_LEN);
					sort_init_algns(sorted[i], ortho_algns[i], num_algns[i], INIT_SELF1);
					species[i].id = *size2;  
					if( species[ref_id].id == ref_id ) {
						species[ref_id].id = *size1; // species[i].id field stores the sequence length of species species[i]
					}
				}
			}
		}
		free(count);
		free(size1);
		free(size2);
	}

	t = sp_tree;
	while( is_leaf_node(t) == false ) {
		if( is_in_list(ref_id, t->left->ch_sp, t->left->num_csp) == true ) {
			t = t->left;
		}
		else t = t->right;	
	}

	while( t != sp_tree ) { // sp_tree is a root
		t = t->parent;
		if( is_in_list(ref_id, t->left->ch_sp, t->left->num_csp) == true ) {
			out_t = t->right;
		}
		else {
			out_t = t->left;
		}
		num_list = out_t->num_csp;
		temp_list[num_sp_event][0] = num_list;
		for( i = 0; i < num_list; i++ ) {
			temp_list[num_sp_event][i+1] = out_t->ch_sp[i];
		}

		prev_num = num_cur_ops;
		num_cur_ops = adjust_ops_list(cur_ops, num_cur_ops, ops, num_ops, temp_list[num_sp_event], num_list+1);
		if( (num_cur_ops - prev_num) > 0 ) {
//			for( i = 0; i < (num_cur_ops - prev_num); i++ ) {
//				printf("%c %d %d %d %d 0 %d %d\n", cur_ops[prev_num+i].sign, cur_ops[prev_num+i].srcStart, cur_ops[prev_num+i].srcEnd, cur_ops[prev_num+i].dstStart, cur_ops[prev_num+i].dstEnd, cur_ops[prev_num+i].dir, (int)cur_ops[prev_num+i].pid);
//			}
		}

		if( num_list > 0 ) {
			cur_ops[num_cur_ops].sign = '#'; // speciation event
			cur_ops[num_cur_ops].srcStart = num_sp_event; // list index for the species split in this speciation event 
			num_cur_ops++;
			num_sp_event++;
//			printf("# sp");
//			for( i = 0; i < num_list; i++ ) {
//				printf(" %s", species[temp_list[i]].name);
//			}
//			printf("\n");
		}
	}

	adjust_direction(cur_ops, num_cur_ops);
	for( i = 0; i < num_cur_ops; i++ ) {
		if( cur_ops[i].sign == 'n' ) {}
		else if( cur_ops[i].sign == '#' ) {
			printf("# sp");
			index = cur_ops[i].srcStart;
			num_list = temp_list[index][0];
			for( j = 1; j <= num_list; j++ ) {
				printf(" %s", species[temp_list[index][j]].name);
			}
			printf("\n");
		}
		else {
			ctg_id1 = cur_ops[i].ctg_id1;
			ctg_id2 = cur_ops[i].ctg_id2;
			len1 = 0;
			len2 = 0;
			if( ctg_id1 >= num_contigs[ref_id] ) {
				fatalf("%d exceeds %d in adjust_ops.c\n", ctg_id1, num_contigs[ref_id]); 
			}
			else if( ctg_id1 == -1 ) {
				strcpy(name1, argv[2]);
			}
			else {
//				sprintf(name1, "%s.%s", contigs[ref_id][ctg_id1].name1, contigs[ref_id][ctg_id1].name2);
				strcpy(name1, contigs[ref_id][ctg_id1].name1);
				len1 = contigs[ref_id][ctg_id1].len;
			}

			if( ctg_id2 >= num_contigs[ref_id] ) {
				fatalf("%d exceeds %d in adjust_ops.c\n", ctg_id2, num_contigs[ref_id]); 
			}
			else if( ctg_id2 == -1 ) {
				strcpy(name2, argv[2]);
			}
			else {
//				sprintf(name2, "%s.%s", contigs[ref_id][ctg_id2].name1, contigs[ref_id][ctg_id2].name2);
				strcpy(name2, contigs[ref_id][ctg_id2].name1);
				len2 = contigs[ref_id][ctg_id2].len;
			}
			printf("%c %s %d %d %s %d %d 0 %d %d\n", cur_ops[i].sign, name1, cur_ops[i].srcStart-len1, cur_ops[i].srcEnd-len1, name2, cur_ops[i].dstStart-len2, cur_ops[i].dstEnd-len2, cur_ops[i].dir, (int)cur_ops[i].pid);
		}
	}

	if( num_sp > 0 ) {
		for( i = 0; i < num_sp; i++ ) {
			if( num_ops[i] > 0 ) {
				free(ops[i]);
			}
			if( num_algns[i] > 0 ) {
				free(ortho_algns[i]);
				free(sorted[i]);
			}
			if( num_contigs[i] > 0 ) {
				free(contigs[i]);
			}
		}

		if( num_total_ops > 0 ) {
			free(cur_ops);
		}

		if( num_sp > 0 ) {
			for( i = 0; i < num_sp; i++ ) {
				free(temp_list[i]);
			}
			free(temp_list);
		}
		free(num_contigs);
		free(contigs);
		free(num_ops);
		free(num_algns);
		free(ops);
		free(ortho_algns);
		free(sorted);
		free(species);
		free_p_tree(sp_tree);
	}
	return(EXIT_SUCCESS);
}	

//int adjust_ops_list(struct ops_list *cur_ops, int num_cur_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list, struct DotList **ortho_algns, int *num_algns, struct slist **sorted)
int adjust_ops_list(struct ops_list *cur_ops, int num_cur_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list)
{
	int i = 0, j = 0, cur_sp = 0;
	int count = 0;
	bool *is_dir_same; // is the dup direction same? 

	is_dir_same = (bool *) ckalloc(sizeof(bool));

	count = num_cur_ops;
	for( i = 1; i < num_list; i++ ) {
		cur_sp = temp_list[i];
		for( j = 0; j < num_ops[cur_sp]; j++ ) {
			if( is_already_in(ops[cur_sp][j], cur_ops, count) == true ) {} 
			else if( do_all_species_agree(ops[cur_sp][j], ops, num_ops, temp_list, num_list, is_dir_same) == true ) {
				cur_ops[count] = assign_ops(ops[cur_sp][j]);
				if( (*is_dir_same) == false ) {
					cur_ops[count].srcStart = ops[cur_sp][j].dstStart;
					cur_ops[count].srcEnd = ops[cur_sp][j].dstEnd;
					cur_ops[count].dstStart = ops[cur_sp][j].srcStart;
					cur_ops[count].dstEnd = ops[cur_sp][j].srcEnd;
					cur_ops[count].src_b = ops[cur_sp][j].dst_b;
					cur_ops[count].src_e = ops[cur_sp][j].dst_e;
					cur_ops[count].dst_b = ops[cur_sp][j].src_b;
					cur_ops[count].dst_e = ops[cur_sp][j].src_e;
				}
				cur_ops[count].sp_id = j;
				count++;
			}
		}	
	} 

	free(is_dir_same);
	return(count);
}

bool is_already_in(struct ops_list tmp_ops, struct ops_list *ops, int num_ops)
{
	int i = 0;
	bool is_in = false;
	struct I ops1_src, ops1_dst, ops2_src, ops2_dst;

	ops1_src = assign_I(tmp_ops.srcStart,tmp_ops.srcEnd);
	ops1_dst = assign_I(tmp_ops.dstStart,tmp_ops.dstEnd);
	ops2_src = assign_I(0,1);
	ops2_dst = assign_I(0,1);
	while( (i < num_ops) && (is_in == false) ) {
		ops2_src = assign_I(ops[i].srcStart,ops[i].srcEnd);
		ops2_dst = assign_I(ops[i].dstStart,ops[i].dstEnd);
		if(tmp_ops.sign == ops[i].sign) {
			if( ((strict_almost_equal(ops1_src, ops2_src) == true) && (strict_almost_equal(ops1_dst, ops2_dst) == true)) || ((strict_almost_equal(ops1_src, ops2_dst) == true) && (strict_almost_equal(ops1_dst, ops2_src) == true))) {
				is_in = true;
			}
		}
		i++;
	}
	return(is_in);
}

//bool do_all_species_agree(struct ops_list tmp_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list, struct DotList **ortho_algns, int *num_algns, struct slist **sorted, bool *is_dir_same)
bool do_all_species_agree(struct ops_list tmp_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list, bool *is_dir_same)
{
	bool res = false, is_in = false;
	struct I src, dst;
	struct I ops_src, ops_dst;
	int i = 0, j = 0, cur_sp = 0;
	int count_same = 0, count_opp = 0;
	int count = 0;

	ops_src = assign_I(0, 1);
	ops_dst = assign_I(0, 1);
	src = assign_I(tmp_ops.srcStart,tmp_ops.srcEnd);
	dst = assign_I(tmp_ops.dstStart,tmp_ops.dstEnd);
	
	i = 1;
//	while( ( i < num_list) && ( res == true ) ) {
	while( i < num_list ) {
		cur_sp = temp_list[i];
		is_in = false;
		j = 0;
		while((j < num_ops[cur_sp]) && (is_in == false)) {
			ops_src = assign_I(ops[cur_sp][j].srcStart,ops[cur_sp][j].srcEnd);
			ops_dst = assign_I(ops[cur_sp][j].dstStart,ops[cur_sp][j].dstEnd);
			if((tmp_ops.sign == ops[cur_sp][j].sign) && ( ((equal(src, ops_src) == true) && (equal(dst, ops_dst) == true)) || ((equal(src, ops_dst) == true) && (equal(dst, ops_src) == true)))) {
				is_in = true;
				if((equal(src, ops_src) == true) && (equal(dst, ops_dst) == true)){
					count_same++;
				}
				else {
					count_opp++;	
				}
			}
			j++;
		}	
		
		if( is_in == false ) {
/*
			is_one_to_one = false;
			is_one_to_one = check_one_to_one_mapping(tmp_ops.pid, src, dst, ortho_algns[cur_sp], num_algns[cur_sp], sorted[cur_sp]);
			if( is_one_to_one == true ) {
				res = false;
			}
			res = false;
*/
		}
		else {
			count++;
		}
		i++;
	} 

	if( count >= (int)(0.7*((float)num_list)+0.5) ) {
		res = true;
	}

	if( res == true ) {
		if( count_same > count_opp ) {
			*is_dir_same = true; 
		}
		else if( count_opp > count_same ) {
			*is_dir_same = false; 
		}
		else {
			if( src.lower < dst.lower ) {
				*is_dir_same = false;
			}
			else {
				*is_dir_same = true;
			}
		}
	}
	return(res);
}

bool check_one_to_one_mapping(int pid, struct I src, struct I dst, struct DotList *algns, int num_algns, struct slist *sorted)
{
	int s_loc = 0, e_loc = 0;
	int i = 0, j = 0, k = 0;
	bool res = false;
	struct I temp_reg;

	temp_reg = assign_I(0, 1);
	k = 0;
	while( (k < 2) && (res == false) ) {
		if( k == 0 ) {
			temp_reg = assign_I(src.lower, src.upper);
		}
		else {
			temp_reg = assign_I(dst.lower, dst.upper);
		}
		s_loc = search_range_b(sorted, algns, num_algns, temp_reg.lower, SELF1);
		e_loc = search_range_b(sorted, algns, num_algns, temp_reg.upper, SELF1);

		res = false;
		i = s_loc;
		while( (i <= e_loc) && (res == false)) {
			j = sorted[i].id;
			if( (proper_overlap(temp_reg, algns[j].x) == true) && (width(intersect(temp_reg, algns[j].x)) >= (int)(0.8 * (float)width(temp_reg))) ) {
				if( pid >= algns[j].identity ) {}
				else {
					res = true;
				}
			}
			i++;
		}
		k++;
	}
	return(res);
}

void adjust_direction(struct ops_list *ops, int num_ops)
{
	int i = num_ops-1, j = num_ops-1;
	struct I ops_src, ops_dst;
	struct I prev_ops_src, prev_ops_dst;
	bool is_adjusted = false;
	bool is_switched = false;
	int b = 0, e = 0;

	ops_src = assign_I(0, 1);
	ops_dst = assign_I(0, 1);
	prev_ops_src = assign_I(0, 1);
	prev_ops_dst = assign_I(0, 1);
	for( i = (num_ops-1); i >= 0; i-- ) { 
		if( (ops[i].sign == '+') || (ops[i].sign == '-') ) {
			ops_src = assign_I(ops[i].srcStart,ops[i].srcEnd);
			ops_dst = assign_I(ops[i].dstStart,ops[i].dstEnd);
			j = num_ops-1;
			is_adjusted = false;
			is_switched = false;
			while( (j > i) && (is_adjusted == false) ) {
				if( (ops[j].sign == '+') || (ops[j].sign == '-') ) {
					prev_ops_src = assign_I(ops[j].srcStart,ops[j].srcEnd);
					prev_ops_dst = assign_I(ops[j].dstStart,ops[j].dstEnd);
					if( subset(ops_dst, prev_ops_dst) == true ) {
						ops[i].sign = 'n'; 
						is_adjusted = true;
					}

					if( (subset(ops_dst, prev_ops_src) == true) || ((width(prev_ops_src) >= 3*ERR_TH) && (subset(prev_ops_src, ops_dst) == true)) ) {
						if( is_switched == false ) {
							b = ops[i].dstStart;
							e = ops[i].dstEnd;
							ops[i].dstStart = ops[i].srcStart;
							ops[i].dstEnd = ops[i].srcEnd;
							ops[i].srcStart = b;
							ops[i].srcEnd = e;
							ops_src = assign_I(ops[i].srcStart,ops[i].srcEnd);
							ops_dst = assign_I(ops[i].dstStart,ops[i].dstEnd);
							is_switched = true;
							j = num_ops-1;
						}
						else {
							ops[i].sign = 'n'; 
							is_adjusted = true;
						}
					}
				}	
				j--;
			}
		}
	}
}
