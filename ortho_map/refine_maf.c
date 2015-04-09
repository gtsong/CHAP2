#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_input.h"
#include "util_ops.h"
#include "util_gen.h"
#include "util_algns.h"
#include "util_i.h"
#include "regions.h"
#include "id_ortho_conv.h"
#include "read_algn.h"
#include "refine_algns.h"
#include "update_init_algns.h"
#include "const_graph.h"
#include "contigs_op.h"

#define ORTHO_CONTENT 1
#define ORTHO_POSITION 2
#define SIMPLE_CUT 3

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

void extend_algn_ends(struct DotList *position_ortho, int num_position, FILE *f);

int main(int argc, char **argv)
{
	struct DotList *position_ortho;
	struct DotList *content_ortho;
	struct DotList *algns; // alignment from chaining 
	int *num_position, *num_content, *num_algns;
	int *size1, *size2;
	int count = 0;
	char species[100], species2[100];
	char ctg_name[100], name[100];
	int i = 0, j = 0;
	FILE *f, *g;
	int ortho_mode = ORTHO_CONTENT;
	struct ops_list *ops1, *ops2;
	int num_ops1 = 0, num_ops2 = 0;
	struct ops_list *cur_ops1, *cur_ops2;
	int num_cur_ops1 = 0, num_cur_ops2 = 0;
	int *num_alloc_blocks1, *num_alloc_blocks2;
	struct I reg;
	int total_len = 0;
	float avg_pid = (float)0;
	struct n_pair *contigs1, *contigs2;
	int num_contigs1 = 0, num_contigs2 = 0;
	int *len_sum1, *len_sum2;

	reg = assign_I(0, 1);
	
	num_alloc_blocks1 = (int *) ckalloc(sizeof(int));
	num_alloc_blocks2 = (int *) ckalloc(sizeof(int));

	*num_alloc_blocks1 = 0;
	*num_alloc_blocks2 = 0;
	debug_mode = FALSE;
	if( argc == 9 ) {
		debug_mode = TRUE;
	}
	else if( argc != 8 ) {
		fatal("args: ortho_maf_file initial_maf_file ops_file1 ops_file2 output_position_maf output_content_maf contigs_list\n");
	}

	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_position = (int *) ckalloc(sizeof(int));
	num_content = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	*size1 = 0;
	*size2 = 0;
	*num_position = 0;
	*num_content = 0;
	*num_algns = 0;
	count = count_local_algns(argv[1], species, species2);
	strcpy(name, species);
	concat_ctg_name(name, species, ctg_name);
	strcpy(name, species2);
	concat_ctg_name(name, species2, ctg_name);

	f = ckopen(argv[7], "r");
  num_contigs1 = count_contigs(species, f);
  num_contigs2 = count_contigs(species2, f);

	if( num_contigs1 > 0 ) {
		contigs1 = (struct n_pair *) ckalloc(num_contigs1 * sizeof(struct n_pair));
		len_sum1 = (int *) ckalloc(num_contigs1 * sizeof(int));
		read_contigs_file(species, f, contigs1, num_contigs1);
		cal_length_sum_from_ctglist(len_sum1, contigs1, num_contigs1);
	}

  if( num_contigs2 > 0 ) {
		contigs2 = (struct n_pair *) ckalloc(num_contigs2 * sizeof(struct n_pair));
		len_sum2 = (int *) ckalloc(num_contigs2 * sizeof(int));
		read_contigs_file(species2, f, contigs2, num_contigs2);
		cal_length_sum_from_ctglist(len_sum2, contigs2, num_contigs2);
	}
	fclose(f);

	if( count > 0 ) {
		position_ortho = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(position_ortho, 0, count);

		if( ortho_mode == SIMPLE_CUT ) {
		  read_maf(argv[1], PAIR_MODE, position_ortho, num_position, size1, size2);
		}
		else if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
		  read_maf(argv[1], O_MODE, position_ortho, num_position, size1, size2);
		}
		adjust_algn_pos(position_ortho, *num_position, contigs1, num_contigs1, size1, contigs2, num_contigs2, size2, CTG_NOT_ASSIGNED_BUT_LEN);
	}

	if( count != (*num_position) ) {
		fatalf("counting error 1: %d vs. %d\n", count, *num_position);
	}

	count = count_local_algns(argv[2], species, species2);
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		position_ortho = (struct DotList *) ckrealloc(position_ortho, (count+(*num_position)) * sizeof(struct DotList));
		content_ortho = (struct DotList *) ckalloc((count+(*num_position)) * sizeof(struct DotList));
		*num_alloc_blocks1 = count + *num_position;
		*num_alloc_blocks2 = count + *num_position;
		initialize_algns(position_ortho, *num_position, count+(*num_position));
		initialize_algns(algns, 0, count);
		initialize_algns(content_ortho, 0, *num_alloc_blocks2);
	  read_maf(argv[2], PAIR_MODE, algns, num_algns, size1, size2);
		adjust_algn_pos(algns, *num_algns, contigs1, num_contigs1, size1, contigs2, num_contigs2, size2, CTG_NOT_ASSIGNED_BUT_LEN);
		
		if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
			*num_position = adjust_all_algns_diff(position_ortho, *num_position, algns, *num_algns);
		}
	}

	if( count != (*num_algns) ) {
		fatalf("counting error 3: %d vs. %d\n", count, *num_algns);
	}

	if( (f = ckopen(argv[3], "r")) != NULL ) {
		num_ops1 = count_final_ops(species2, f, SELF1);
		if( num_ops1 > 0 ) {
			ops1 = (struct ops_list *) ckalloc(num_ops1 * sizeof(struct ops_list));
			cur_ops1 = (struct ops_list *) ckalloc(num_ops1 * sizeof(struct ops_list));
			read_final_ops_file(species2, f, ops1, num_ops1, SELF1, contigs1, num_contigs1);
		}
		fclose(f);
	}
	else {
		fatalf("file open error: %s\n", argv[3]);
	}

	if( (f = ckopen(argv[4], "r")) != NULL ) {
		num_ops2 = count_ops(species2, f, SELF2);
		if( num_ops2 > 0 ) {
			ops2 = (struct ops_list *) ckalloc(num_ops2 * sizeof(struct ops_list));
			cur_ops2 = (struct ops_list *) ckalloc(num_ops2 * sizeof(struct ops_list));
			read_ops_file(species2, f, ops2, num_ops2, SELF2, contigs2, num_contigs2, species);
		}
		fclose(f);
	}
	else {
		fatalf("file open error: %s\n", argv[4]);
	}

	f = ckopen(argv[2], "r");

  for( i = 0; i < (*num_position); i++ ) {      
		if( position_ortho[i].sign != DELETED) {
      total_len = total_len + width(position_ortho[i].x);
      avg_pid = avg_pid + (float)(((float)position_ortho[i].identity/(float)100) * width(position_ortho[i].x));      
		}
  }

  if( total_len > 0 ) {
    avg_pid = (avg_pid/(float)total_len) * (float)100;
	}

	j = 0;
	for( i = 0; i < num_ops2; i++ ) {
		if( (ops2[i].sign == '+') || (ops2[i].sign == '-') ) {
			cur_ops2[j] = assign_ops(ops2[i]);
			j++;
			reg = assign_I(ops2[i].dstStart, ops2[i].dstEnd);
			undo_events(reg, position_ortho, *num_position, SELF2, f);
		}
	}
	num_cur_ops2 = j;

	j = 0;
	for( i = 0; i < num_ops1; i++ ) {
		if( (ops1[i].sign == '+') || (ops1[i].sign == '-') ) {
			cur_ops1[j] = assign_ops(ops1[i]);
			j++;
			reg = assign_I(ops1[i].dstStart, ops1[i].dstEnd);
			undo_events(reg, position_ortho, *num_position, SELF1, f);
		}
	}
	num_cur_ops1 = j;

  update_pid_init_algns(*num_position, position_ortho, f, avg_pid);
	const_graph(*num_position, position_ortho, f);

  for( i = 0; i < (*num_position); i++ ) {      
		if( position_ortho[i].sign != DELETED) {
      total_len = total_len + width(position_ortho[i].x);
      avg_pid = avg_pid + (float)(((float)position_ortho[i].identity/(float)100) * width(position_ortho[i].x));      
		}
  }

  if( total_len > 0 ) {
    avg_pid = (avg_pid/(float)total_len) * (float)100;
	}
  update_pid_init_algns(*num_position, position_ortho, f, avg_pid);
	remove_false_mtom(position_ortho, *num_position, ops1, num_ops1, f);

	j = 0;
	for( i = 0; i < *num_position; i++ ) {
		assign_algn(content_ortho, j, position_ortho[i]);
		content_ortho[j].index = j;
		content_ortho[j].sp_id = PAIR;
		j++;
	}
	*num_content = j;

  *num_position = redo_dups_for_mtom_inc_conv(num_cur_ops1, cur_ops1, *num_position, &position_ortho, f, REF_SEQ, num_alloc_blocks1);
  *num_position = redo_dups_for_mtom_inc_conv(num_cur_ops2, cur_ops2, *num_position, &position_ortho, f, SELF2, num_alloc_blocks1);
	g = ckopen(argv[5], "w");

	extend_algn_ends(position_ortho, *num_position, f);
	replace_contigs_len(contigs1, num_contigs1);
	replace_contigs_len(contigs2, num_contigs2);
	write_init_maf(g, position_ortho, *num_position, contigs1, contigs2, num_contigs1, num_contigs2, f, PAIR);
	fclose(g);

  *num_content = redo_dups_for_mtom_inc_conv(num_ops1, ops1, *num_content, &content_ortho, f, REF_SEQ, num_alloc_blocks2);
  *num_content = redo_dups_for_mtom_inc_conv(num_ops2, ops2, *num_content, &content_ortho, f, SELF2, num_alloc_blocks2);
	g = ckopen(argv[6], "w");
	extend_algn_ends(content_ortho, *num_content, f);
	write_init_maf(g, content_ortho, *num_content, contigs1, contigs2, num_contigs1, num_contigs2, f, PAIR);
	fclose(g);

	fclose(f);

	if( num_ops1 > 0 ) {
		free(ops1);
		free(cur_ops1);
	}

	if( num_contigs1 > 0 ) {
		free(contigs1);
		free(len_sum1);
	}

	if( num_ops2 > 0 ) {
		free(ops2);
		free(cur_ops2);
	}

	if( num_contigs2 > 0 ) {
		free(contigs2);
		free(len_sum2);
	}

	free(num_alloc_blocks1);
	free(num_alloc_blocks2);
	free(size1);
	free(size2);

	if( (*num_algns) > 0 ) {
		free(algns);
		free(position_ortho);
		free(content_ortho);
	}

	free(num_position);
	free(num_content);
	free(num_algns);
	return(EXIT_SUCCESS);
}

void extend_algn_ends(struct DotList *position_ortho, int num_position, FILE *f)
{
	int i = 0, j = 0;
	bool is_overlap = false;
	struct I reg1, reg2, cmp1, cmp2;
	cmp1 = assign_I(0, 1);
	cmp2 = assign_I(0, 1);
	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);
	int max_id = -1;
	int max_width = 0, cur_width_sp1 = 0, cur_width_sp2 = 0;
	int cur_width = 0, cur_sp_id = SELF1;
	int sp_id = SELF1;
	struct I ov;

	ov = assign_I(0,1);

	for( i = 0; i < num_position; i++ ) {
		if( (position_ortho[i].sign != DELETED) && ((position_ortho[i].xl_diff != 0) || (position_ortho[i].xr_diff != 0) || (position_ortho[i].yl_diff != 0) || (position_ortho[i].yr_diff != 0) ) ) {
			reg1 = assign_I(position_ortho[i].x.lower, position_ortho[i].x.upper);
			reg2 = assign_I(position_ortho[i].y.lower, position_ortho[i].y.upper);
			is_overlap = false;
			j = 0;
			max_id = -1;
			max_width = 0;
			while((j < num_position) && (is_overlap == false)) {
				if( i == j ) {}
				else if( position_ortho[j].sign != DELETED ) {
					cmp1 = assign_I(position_ortho[j].x.lower+position_ortho[j].xl_diff, position_ortho[j].x.upper-position_ortho[j].xr_diff);
					cmp2 = assign_I(position_ortho[j].y.lower+position_ortho[j].yl_diff, position_ortho[j].y.upper-position_ortho[j].yr_diff);
					if( (proper_overlap(reg1, cmp1) == true) || (proper_overlap(reg2, cmp2) == true ) ) {
						is_overlap = true;
						cur_width_sp1 = 0;
						cur_width_sp2 = 0;
						cur_width = 0;
						if( (f_loose_subset(reg1, cmp1, STRICT) == false) && (f_loose_subset(cmp1, reg1, STRICT) == false) && (proper_overlap(reg1, cmp1) == true) ) {
							cur_width_sp1 = width(intersect(reg1, cmp1));
						}

						if( (f_loose_subset(reg2, cmp2, STRICT) == false) && (f_loose_subset(cmp2, reg2, STRICT) == false) && proper_overlap(reg2, cmp2) == true ) {
							cur_width_sp2 = width(intersect(reg2, cmp2));
						}

						if( (cur_width_sp1 == 0) && (cur_width_sp2 == 0) ) {}
						else if( cur_width_sp1 >= cur_width_sp2 ) {
							cur_width = cur_width_sp1;
							cur_sp_id = SELF1;
						}
						else {
							cur_width = cur_width_sp2;
							cur_sp_id = SELF2;
						}

						if( cur_width > max_width ) {
							max_width = cur_width;
							sp_id = cur_sp_id;
							max_id = j;
						}
					}
				}
				j++;
			}
			if( is_overlap == false ) {
				position_ortho[i].xl_diff = 0;
				position_ortho[i].xr_diff = 0;
				position_ortho[i].yl_diff = 0;
				position_ortho[i].yr_diff = 0;
			}
			else if( (max_width > 0) && (sp_id == SELF1) ) {
				if( (reg1.lower < cmp1.lower) && (position_ortho[i].xr_diff > max_width) )	
				{
					position_ortho[i].xr_diff = 0;
					ov = intersect(reg1, cmp1);
					adjust_init_algn(position_ortho, i, ov, true, f, DUP, true);
				}
				else if( (reg1.lower > cmp1.lower) && (position_ortho[i].xl_diff > max_width) )	
				{
					position_ortho[i].xl_diff = 0;
					ov = intersect(reg1, cmp1);
					adjust_init_algn(position_ortho, i, ov, true, f, DUP, true);
				}
			}
			else if( (max_width > 0) && (sp_id == SELF2) ) {
				ov = intersect(reg2, cmp2);
				if( (reg2.lower < cmp2.lower) && (position_ortho[i].yr_diff > max_width) )	
				{
					position_ortho[i].yr_diff = 0;
					adjust_init_algn(position_ortho, i, ov, false, f, DUP, false);
				}
				else if( (reg2.lower > cmp2.lower) && (position_ortho[i].yl_diff > max_width) )	
				{
					position_ortho[i].yl_diff = 0;
					adjust_init_algn(position_ortho, i, ov, false, f, DUP, false);
				}
			}
		}
	}
}
