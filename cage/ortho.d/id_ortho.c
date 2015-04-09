#include "main.h"
#include "id_ortho.h"
#include "regions.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "read_algn.h"
#include "util_algns.h"

extern int debug_mode;

void obtain_ortho_algns(struct DotList *algns, int num_algns, struct DotList *init_algns, int num_init_algns)
{
	struct DotList *pair_algns;
	int *num_pair;
	int i, j;
	int max_x = 0, max_y = 0;
	struct kdnode *tree;
	struct perm_pt *p_pts;
	float avg_pid;
	int sum_int;
	int cur_len;

	j = 0;
	for( i = 0; i < num_algns; i++ ) {
		if( (algns[i].pair_self == PAIR) && (algns[i].sign != DELETED) ) {
			if( algns[i].x.upper > max_x ) max_x = algns[i].x.upper;
			if( algns[i].y.upper > max_y ) max_y = algns[i].y.upper;
			j++;
		}
	}

	num_pair = (int *) ckalloc(sizeof(int));
	*num_pair = j;
	pair_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * j);
	p_pts = (struct perm_pt *) ckalloc(sizeof(struct perm_pt) * j);

	j = 0;
	for( i = 0; i < num_algns; i++ ) {
		if( (algns[i].pair_self == PAIR) && (algns[i].sign != DELETED) ) {
			assign_algn(pair_algns, j, algns[i]);
			j++;
		}
	}

	assign_perm(p_pts, j, pair_algns, LEFT);
	tree = build_kd(p_pts, 0, j-1);

	for( i = 0; i < num_init_algns; i++ ) {
		if( init_algns[i].pair_self == PAIR ) {
			init_algns[i].c_id = -1;
			init_algns[i].m_id = -1;
		}

		if( (init_algns[i].x.upper <= init_algns[i].x.lower) || (init_algns[i].y.upper <= init_algns[i].y.lower) ) {
			init_algns[i].x = assign_I(init_algns[i].x.lower-init_algns[i].xl_offset, init_algns[i].x.upper-init_algns[i].xr_offset);
			init_algns[i].y = assign_I(init_algns[i].y.lower-init_algns[i].yl_offset, init_algns[i].y.upper-init_algns[i].yr_offset);
			init_algns[i].xl_offset = 0;
			init_algns[i].xr_offset = 0;
			init_algns[i].yl_offset = 0;
			init_algns[i].yr_offset = 0;
		}
	}

	iden_ortho_update_init(num_pair, pair_algns, num_init_algns, init_algns);

	avg_pid = (float)0;
	sum_int = 0;
	for( i = 0; i < num_init_algns; i++ ) {
		if( init_algns[i].sp_id == PAIR ) {
			if( init_algns[i].sign != DELETED ) {
				cur_len = (init_algns[i].x.upper - init_algns[i].xr_diff) - (init_algns[i].x.lower + init_algns[i].xl_diff);
				if( cur_len > 0 ) {
					avg_pid = avg_pid + (float)(((float)init_algns[i].identity)/((float)100)) * ((float) cur_len);
					sum_int = sum_int + cur_len;
				}
			}
		}
	}
	avg_pid = (avg_pid / ((float)(sum_int))) * ((float)100);

	free(num_pair);
	free(p_pts);
	free_kd(tree);
	free(pair_algns);
}

int assign_ortho_algns(struct DotList *ortho_algns, struct DotList *init_algns, int num_init_algns, float *temp_val)
{
	int i = 0, j = 0;
	int sum_int = 0;
	float avg_pid;
	int cur_len = 0;

  j = 0;
  avg_pid = (float)0;
  sum_int = 0;
  for( i = 0; i < num_init_algns; i++ ) {
    if( init_algns[i].sp_id == PAIR ) {
      cur_len = (init_algns[i].x.upper - init_algns[i].xr_diff) - (init_algns[i].x.lower + init_algns[i].xl_diff);
      if( init_algns[i].sign != DELETED ) {
        if( cur_len > 0 ) {
          avg_pid = avg_pid + (float)(((float)init_algns[i].identity)/((float)100)) * ((float) cur_len);
          sum_int = sum_int + cur_len;
        }
				else {
					init_algns[i].sign = DELETED;
				}
      }
      cur_len = (init_algns[i].x.upper - init_algns[i].x.lower);
			if( cur_len > DEL_TH ) {
      	assign_algn(ortho_algns, j, init_algns[i]);
      	j++;
			}
    }
  }

	if( sum_int > 0 ) {
		avg_pid = (avg_pid / ((float)(sum_int))) * ((float)100);	
		*temp_val = avg_pid;
	}
	else {
		*temp_val = (float)(-1);
	}

	return(j);
}

void add_ortho_intervals_dups(struct ops_list *ops, int id, struct slist *sorted, struct DotList *ortho_algns, int num_ortho_algns, int cur_sp_id, int *exc_list, FILE *f, float avg_pid)
{
	struct I query, tmp;
	int s_loc = 0, e_loc = 0, s_d = 0, e_d = 0, num_exc = 0;
	int j = 0;
	bool src_exist = true;
	bool *skip_untagging;

	skip_untagging = (bool *) ckalloc(sizeof(bool));

	tmp = assign_I(0, 1);
	query = assign_I(ops[id].srcStart, ops[id].srcEnd);
	s_loc = search_range_b(sorted, ortho_algns, num_ortho_algns, query.lower, cur_sp_id);
	if( debug_mode == TRUE ) printf("start loc: %d-%d,%d-%d for %d-%d\n", ortho_algns[sorted[s_loc].id].x.lower, ortho_algns[sorted[s_loc].id].x.upper, ortho_algns[sorted[s_loc].id].y.lower, ortho_algns[sorted[s_loc].id].y.upper, ops[id].srcStart, ops[id].srcEnd);
	e_loc = search_range_e(sorted, ortho_algns, num_ortho_algns, query.upper, cur_sp_id);
	if( debug_mode == TRUE ) printf("end loc: %d-%d,%d-%d for %d-%d\n", ortho_algns[sorted[e_loc].id].x.lower, ortho_algns[sorted[e_loc].id].x.upper, ortho_algns[sorted[e_loc].id].y.lower, ortho_algns[sorted[e_loc].id].y.upper, ops[id].srcStart, ops[id].srcEnd);
	query = assign_I(ops[id].dstStart, ops[id].dstEnd);
	num_exc = make_exc_list(query, ops, id, exc_list);

	if( debug_mode == TRUE ) {
		if( num_exc > 0 ) printf("For %d-%d\n", ops[id].dstStart, ops[id].dstEnd);
		for( j = 0; j < num_exc; j++ ) {
			printf("%d-%d is excluded\n", ops[exc_list[j]].dstStart, ops[exc_list[j]].dstEnd);
		}
	}
	s_d = search_range_b(sorted, ortho_algns, num_ortho_algns, query.lower, cur_sp_id);
	e_d = search_range_e(sorted, ortho_algns, num_ortho_algns, query.upper, cur_sp_id);

	if( (ops[id].sign == 'c') || (ops[id].sign == 'v') ) {
		src_exist = true;
	}
	else {
		*skip_untagging = false;
		src_exist = check_source_existence(ops, id, sorted, s_loc, e_loc, s_d, e_d, ortho_algns, cur_sp_id, f, avg_pid, skip_untagging);
	}

	if( src_exist == true ) {
		add_ortho_intervals(ops, id, sorted, s_loc, e_loc, s_d, e_d, ortho_algns, cur_sp_id, f, avg_pid);
	}
	else {
		tmp = assign_I(ops[id].srcStart, ops[id].srcEnd);
		ops[id].srcStart = ops[id].dstStart;
		ops[id].srcEnd = ops[id].dstEnd;
		ops[id].dstStart = tmp.lower;
		ops[id].dstEnd = tmp.upper;
		add_ortho_intervals(ops, id, sorted, s_d, e_d, s_loc, e_loc, ortho_algns, cur_sp_id, f, avg_pid);
	}

	free(skip_untagging);
}

void adjust_boundary_init(struct DotList *init_algns, struct DotList *ortho_algns, int id, int num_init_algns, FILE *f, int sp_id)	
{
	if( init_algns[ortho_algns[id].index].sign != ortho_algns[id].sign ) {
		if( (sp_id == REF_SEQ) && (ortho_algns[id].sign != DELETED) ) {
			adjust_boundary(id, ortho_algns, init_algns, num_init_algns, f);
		}
		init_algns[ortho_algns[id].index].sign = ortho_algns[id].sign;
	} 

	if( init_algns[ortho_algns[id].index].sign != DELETED ) {
		if( (init_algns[ortho_algns[id].index].xl_diff != ortho_algns[id].xl_diff )  || ( init_algns[ortho_algns[id].index].yl_diff != ortho_algns[id].yl_diff ) || ( init_algns[ortho_algns[id].index].xr_diff != ortho_algns[id].xr_diff ) || ( init_algns[ortho_algns[id].index].yr_diff != ortho_algns[id].yr_diff ) )
		{
			if( sp_id == REF_SEQ) adjust_boundary(id, ortho_algns, init_algns, num_init_algns, f);
			init_algns[ortho_algns[id].index].xl_diff = ortho_algns[id].xl_diff;
			init_algns[ortho_algns[id].index].yl_diff = ortho_algns[id].yl_diff;
			init_algns[ortho_algns[id].index].xr_diff = ortho_algns[id].xr_diff;
			init_algns[ortho_algns[id].index].yr_diff = ortho_algns[id].yr_diff;
		}
	}
}

void redo_dups_for_mtom(int num_ops, struct ops_list *ops, int num_init_algns, struct DotList *init_algns, FILE *f, int sp_id)
{
	struct DotList *ortho_algns;
	int num_ortho_algns = 0;
	float *temp_val; 
  struct slist *sorted;
  float avg_pid;
  int *exc_list; // the list of intervals which should be excluded in adding putative orthologous alignments
  int mode;
  int cur_sp_id;
	int i, j;

	temp_val = (float *) ckalloc(sizeof(float));

	ortho_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_init_algns);
	initialize_algns(ortho_algns, num_init_algns);

	if( sp_id == REF_SEQ ) cur_sp_id = SELF1;
	else cur_sp_id = sp_id;

	j = 0;
	avg_pid = (float)0;
	num_ortho_algns = assign_ortho_algns(ortho_algns, init_algns, num_init_algns, temp_val);
	avg_pid = *temp_val;
	free(temp_val);

	if( (num_ortho_algns <= 0) || (avg_pid < 0) ) {}
	else {
	  if( debug_mode == TRUE ) printf("Average pid: %f\n", avg_pid);
 		sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_ortho_algns);
  	if( cur_sp_id == SELF1 ) mode = INIT_SELF1;
  	else if( cur_sp_id == SELF2 ) mode = INIT_SELF2;
  	else {
    	fatalf("invalid species id: it must be either SELF1(0) or SELF2(1), but %d here\n", cur_sp_id);
  	}

  	sort_init_algns(sorted, ortho_algns, num_ortho_algns, mode); // sort by original coordinates in the initial dot plot
  	exc_list = (int *) ckalloc(sizeof(int) * num_ops);

	  for( i = 0; i < num_ops; i++ ) exc_list[i] = 0;

		for( i = num_ops-1; i >= 0; i-- ) {
	    if(((ops[i].sign == '+') || (ops[i].sign == '-')) && (init_algns[ops[i].id].sp_id == cur_sp_id)) {
				add_ortho_intervals_dups(ops, i, sorted, ortho_algns, num_ortho_algns, cur_sp_id, exc_list, f, avg_pid);
			}	
		}

		for( i = 0; i < num_ortho_algns; i++ ) {
			adjust_boundary_init(init_algns, ortho_algns, i, num_init_algns, f, sp_id);	
		}

		free(exc_list);
		free(sorted);
	}

	free(ortho_algns);
}

void iden_ortho_update_init(int *num_pair, struct DotList *pair_alg, int num_init_algns, struct DotList *init_algns)
{
	int i;
	int t_val = STRICT;
	struct slist *sorted;
	struct slist *ortho_sorted;
	int cur_id;
	int cur_loc;
	int temp_loc;
	int *num_init;

	sorted = (struct slist *) ckalloc(sizeof(struct slist) * (*num_pair));
	ortho_sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_init_algns);
	num_init = (int *) ckalloc(sizeof(int));

	for( i = 0; i < (*num_pair); i++ )
	{
		sorted[i].id = i;
	}
	sort_by_width( sorted, pair_alg, *num_pair );
	check_loc_for_ortho_one_side( sorted, num_pair, pair_alg, t_val );

	for( i = 0; i < (*num_pair); i++ )
	{
		sorted[i].id = i;
	}
	sort_init_algns( sorted, pair_alg, *num_pair, SELF1 );

	check_loc_one_side( sorted, num_pair, pair_alg, t_val );
	for( i = 0; i < (*num_pair); i++ )
	{
		cur_id = pair_alg[i].index;
		if( (pair_alg[i].sign == ORTHO) || (pair_alg[i].sign == ORTHO_CANDI) || (pair_alg[i].sign == 0)) 
		{		
			pair_alg[i].sign = 0;
			init_algns[cur_id].sign = ORTHO;
			if( init_algns[cur_id].c_id != -1 )
			{
				cur_loc = init_algns[cur_id].c_id;
				init_algns[cur_loc].sign = ORTHO;			
				while( init_algns[cur_loc].c_id != -1 )
				{
					temp_loc = cur_loc;
					cur_loc = init_algns[cur_loc].c_id;
					init_algns[cur_loc].sign = ORTHO;			
				}
			}
		}
		else if( (pair_alg[i].sign == ORTHO_COMP) || (pair_alg[i].sign == ORTHO_COMP_CANDI) || (pair_alg[i].sign == 1)) 
		{
			pair_alg[i].sign = 1;
			init_algns[cur_id].sign = ORTHO_COMP;			
			if( init_algns[cur_id].c_id != -1 )
			{
				cur_loc = init_algns[cur_id].c_id;
				init_algns[cur_loc].sign = ORTHO_COMP;			
				while( init_algns[cur_loc].c_id != -1 )
				{
					cur_loc = init_algns[cur_loc].c_id;
					init_algns[cur_loc].sign = ORTHO_COMP;			
				}
			}
		}
		else if( pair_alg[i].sign == OUT_PAR ) 
		{
			pair_alg[i].sign = DELETED;
			init_algns[cur_id].sign = DELETED;			
			if( init_algns[cur_id].c_id != -1 )
			{
				cur_loc = init_algns[cur_id].c_id;
				init_algns[cur_loc].sign = DELETED;			
				while( init_algns[cur_loc].c_id != -1 )
				{
					cur_loc = init_algns[cur_loc].c_id;
					init_algns[cur_loc].sign = DELETED;			
				}
			}
		}
		else if( pair_alg[i].sign == OUT_PAR_REV )
		{
			pair_alg[i].sign = DELETED;
			init_algns[cur_id].sign = DELETED;			
			if( init_algns[cur_id].c_id != -1 )
			{
				cur_loc = init_algns[cur_id].c_id;
				init_algns[cur_loc].sign = DELETED;			
				while( init_algns[cur_loc].c_id != -1 )
				{
					cur_loc = init_algns[cur_loc].c_id;
					init_algns[cur_loc].sign = DELETED;			
				}
			}
		}
	}

	for( i = 0; i < num_init_algns; i++ )
	{
		if( (init_algns[i].sign == ORTHO) || (init_algns[i].sign == 0) )
		{
			init_algns[i].sign = 0;
		}
		else if( (init_algns[i].sign == ORTHO_COMP) || (init_algns[i].sign == 1) )
		{
			init_algns[i].sign = 1;
		}
		else init_algns[i].sign = DELETED;	
	}
	for( i = 0; i < num_init_algns; i++ )
	{
		ortho_sorted[i].id = i;
	}

	sort_by_width( ortho_sorted, init_algns, num_init_algns );
	*num_init = num_init_algns;
	check_loc_for_ortho_one_side( ortho_sorted, num_init, init_algns, t_val );

	for( i = 0; i < num_init_algns; i++ )
	{
		if( (init_algns[i].sign == 0) || (init_algns[i].sign == ORTHO) || (init_algns[i].sign == ORTHO_CANDI)) 
		{
			init_algns[i].sign = 0;
		}
		else if( (init_algns[i].sign == 1) || (init_algns[i].sign == ORTHO_COMP) || (init_algns[i].sign == ORTHO_COMP_CANDI)) 
		{
			init_algns[i].sign = 1;
		}
		else init_algns[i].sign = DELETED;

		ortho_sorted[i].id = i;
	}

	sort_init_algns( ortho_sorted, init_algns, num_init_algns, SELF1 );
	*num_init = num_init_algns;
	check_loc_one_side( ortho_sorted, num_init, init_algns, t_val );
	for( i = 0; i < num_init_algns; i++ )
	{
		if( (init_algns[i].sign == 0) || (init_algns[i].sign == ORTHO) || (init_algns[i].sign == ORTHO_CANDI)) 
		{
			init_algns[i].sign = 0;
		}
		else if( (init_algns[i].sign == 1) || (init_algns[i].sign == ORTHO_COMP) || (init_algns[i].sign == ORTHO_COMP_CANDI)) 
		{
			init_algns[i].sign = 1;
		}
		else init_algns[i].sign = DELETED;
	}
	free(num_init);
	free(sorted);
	free(ortho_sorted);
}

void remove_non_ortho(int *num_pair, struct DotList *pair_alg)
{
	int i;
	int t_val = LOOSE;
	struct slist *sorted;

	sorted = (struct slist *) ckalloc((*num_pair)*sizeof(struct slist));

	for( i = 0; i < (*num_pair); i++ )
	{
		sorted[i].id = i;
	}

	sort_by_width( sorted, pair_alg, *num_pair );
	check_loc_for_ortho( sorted, num_pair, pair_alg, t_val );

	for( i = 0; i < (*num_pair); i++ )
	{
		if( pair_alg[i].sign == OUT_PAR ) 
		{
			pair_alg[i].sign = DELETED;
		}
		else if( pair_alg[i].sign == OUT_PAR_REV )
		{
			pair_alg[i].sign = DELETED;
		}
		else if( pair_alg[i].sign == ORTHO_CANDI ) pair_alg[i].sign = 0;
		else if( pair_alg[i].sign == ORTHO_COMP_CANDI ) pair_alg[i].sign = 1;
	}

	overwrite_dots(num_pair, pair_alg);

	free(num_pair);
}

void check_loc_for_ortho(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val)
{
	int i, j;
	bool is_x;

	for( i = ((*num_pair)-1); i >= 0; i-- )
	{
		if( (pair_alg[sorted[i].id].sign == OUT_PAR) || (pair_alg[sorted[i].id].sign == OUT_PAR_REV) ) 
		{
		}
		else if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == 1) )
		{
			if( pair_alg[sorted[i].id].sign == 0 )
			{
				pair_alg[sorted[i].id].sign = ORTHO_CANDI;
			}
			else if( pair_alg[sorted[i].id].sign == 1 )
			{
				pair_alg[sorted[i].id].sign = ORTHO_COMP_CANDI;
			}

			j = 0;
			while( (j < i) && ((pair_alg[sorted[i].id].sign == ORTHO_CANDI) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI)) )
			{
				if( sorted[j].id == sorted[i].id ) {}
				else
				{
					if((strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == true) || (strict_almost_equal(pair_alg[sorted[j].id].y, pair_alg[sorted[i].id].y) == true ))
					{
						if( pair_alg[sorted[j].id].identity > pair_alg[sorted[i].id].identity )
						{
							if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) )
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							if( strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == true )
							{
								is_x = true;
							}
							else is_x = false;
						}
					}
					else if((f_loose_subset(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true) || (f_loose_subset(pair_alg[sorted[i].id].y, pair_alg[sorted[j].id].y, t_val) == true ))
					{
						if( pair_alg[sorted[i].id].identity < pair_alg[sorted[j].id].identity)
						{
							if( (pair_alg[sorted[i].id].sign == 0 ) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							if( f_loose_subset(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true ) is_x = true;
							else is_x = false;
						}
					}
					else
					{
					}
				}
				j++;
			}
		}
	}
}

void check_loc_one_side(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val) // alignments are sorted by x.lower and width
{
	int i, j, pre_j, k;
	int prev_end;
	int d, min_d, id;
	int start;
	int cur_sign;
	int count;
	int num_algns;

	num_algns = *num_pair;

	i = 0;
	while( (i < num_algns) && ((pair_alg[sorted[i].id].sign == DELETED) || (pair_alg[sorted[i].id].sign == OUT_PAR) || (pair_alg[sorted[i].id].sign == OUT_PAR_REV) || (pair_alg[sorted[i].id].pair_self == SELF)) ) i++;

	start = i;
	if((pair_alg[sorted[start].id].sign == ORTHO) || (pair_alg[sorted[start].id].sign == ORTHO_CANDI) || (pair_alg[sorted[start].id].sign == 0)) prev_end = pair_alg[sorted[start].id].y.upper;
	else prev_end = pair_alg[sorted[start].id].y.lower;

	i = start+1;
	while( i < num_algns ) // (i-1)th is prev_end
	{
		if( (pair_alg[sorted[i].id].sign == DELETED) || (pair_alg[sorted[i].id].pair_self == SELF) || (pair_alg[sorted[i].id].sign == OUT_PAR) || (pair_alg[sorted[i].id].sign == OUT_PAR_REV) ) 
		{
			i++;
		}
		else if( (pair_alg[sorted[i].id].sign == ORTHO) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) || (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_COMP) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI) || (pair_alg[sorted[i].id].sign == 1)) 
		{
			j = i;

			if( (pair_alg[sorted[j].id].sign == ORTHO) || (pair_alg[sorted[j].id].sign == ORTHO_CANDI) || (pair_alg[sorted[j].id].sign == 0) ) min_d = abs(prev_end - pair_alg[sorted[j].id].y.lower);
			else if((pair_alg[sorted[j].id].sign == ORTHO_COMP) || (pair_alg[sorted[j].id].sign == ORTHO_COMP_CANDI) || (pair_alg[sorted[j].id].sign == 1)) min_d = abs(prev_end - pair_alg[sorted[j].id].y.upper);
			id = j;

			count = 1;
			pre_j = j;
			j++;

			while( (j < num_algns) && ( (pair_alg[sorted[j].id].pair_self == SELF) || (pair_alg[sorted[j].id].sign == DELETED) || (pair_alg[sorted[j].id].sign == OUT_PAR) || (pair_alg[sorted[j].id].sign == OUT_PAR_REV)) ) j++;

			while( (j < num_algns) && (((pair_alg[sorted[j].id].ctg_id1 == pair_alg[sorted[pre_j].id].ctg_id1) && (strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[pre_j].id].x) == true)) || ((pair_alg[sorted[j].id].ctg_id1 == pair_alg[sorted[i].id].ctg_id1) && (f_loose_subset(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true))) ) {
				cur_sign = pair_alg[sorted[j].id].sign;
				if( (cur_sign == ORTHO) || (cur_sign == ORTHO_CANDI) || (cur_sign == 0) ) cur_sign = 0;
				else if( (cur_sign == ORTHO_COMP) || (cur_sign == ORTHO_COMP_CANDI) || (cur_sign == 1) ) cur_sign = 1;

				if( cur_sign == 0 ) {
					d = abs(prev_end - pair_alg[sorted[j].id].y.lower);
					count++;
				}
				else if( cur_sign == 1 ) {
					d = abs(prev_end - pair_alg[sorted[j].id].y.upper);
					count++;
				}
				
				if( min_d > d ) {
					min_d = d;
					id = j;
				}
				pre_j = j;
				j++;

				while( (j < num_algns) && ((pair_alg[sorted[j].id].sign == DELETED) || (pair_alg[sorted[j].id].pair_self == SELF) || (pair_alg[sorted[j].id].sign == OUT_PAR) || (pair_alg[sorted[j].id].sign == OUT_PAR_REV)) ) j++;
			}

			if( count >= 2 ) {
				if( j > num_algns ) j = num_algns;
				for( k = i; k < j; k++ ) {
					cur_sign = pair_alg[sorted[k].id].sign;
					if( (cur_sign == ORTHO) || (cur_sign == ORTHO_CANDI) || (cur_sign == 0) ) cur_sign = 0;
					else if( (cur_sign == ORTHO_COMP) || (cur_sign == ORTHO_COMP_CANDI) || (cur_sign == 1) ) cur_sign = 1;

					if( k != id ) {
						if( cur_sign == 0 ) pair_alg[sorted[k].id].sign = OUT_PAR;
						else if( cur_sign == 1 ) pair_alg[sorted[k].id].sign = OUT_PAR_REV;
					}
				}					

				cur_sign = pair_alg[sorted[id].id].sign;
				if( (cur_sign == ORTHO) || (cur_sign == ORTHO_CANDI) || (cur_sign == 0) ) cur_sign = 0;
				else if( (cur_sign == ORTHO_COMP) || (cur_sign == ORTHO_COMP_CANDI) || (cur_sign == 1) ) cur_sign = 1;
				if( cur_sign == 0 ) prev_end = pair_alg[sorted[id].id].y.upper;
				else if( cur_sign == 1 ) prev_end = pair_alg[sorted[id].id].y.lower;
				i = j;
			}
			else {
				cur_sign = pair_alg[sorted[i].id].sign;
				if( (cur_sign == ORTHO) || (cur_sign == ORTHO_CANDI) || (cur_sign == 0) ) cur_sign = 0;
				else if( (cur_sign == ORTHO_COMP) || (cur_sign == ORTHO_COMP_CANDI) || (cur_sign == 1) ) cur_sign = 1;
				if( cur_sign == 0 ) prev_end = pair_alg[sorted[i].id].y.upper;
				else if( cur_sign == 1 ) prev_end = pair_alg[sorted[i].id].y.lower;
				i++;
			}
		}
		else i++;
	}
}

void check_loc_for_ortho_one_side(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val)
{
	int i, j;
	bool is_x;

	for( i = ((*num_pair)-1); i >= 0; i-- )
	{
		if( (pair_alg[sorted[i].id].pair_self == SELF) || (pair_alg[sorted[i].id].sign == DELETED) || (pair_alg[sorted[i].id].sign == OUT_PAR) || (pair_alg[sorted[i].id].sign == OUT_PAR_REV) ) 
		{
		}
		else if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == 1) )
		{
			if( pair_alg[sorted[i].id].sign == 0 )
			{
				pair_alg[sorted[i].id].sign = ORTHO_CANDI;
			}
			else if( pair_alg[sorted[i].id].sign == 1 )
			{
				pair_alg[sorted[i].id].sign = ORTHO_COMP_CANDI;
			}

			j = 0;
			while( (j < i) && ((pair_alg[sorted[i].id].sign == ORTHO_CANDI) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI)) )
			{
				if( sorted[j].id == sorted[i].id ) {}
				else if( pair_alg[sorted[j].id].sign == DELETED ) {}
				else if( pair_alg[sorted[j].id].pair_self == SELF ) {}
				else
				{
					if((pair_alg[sorted[j].id].ctg_id1 == pair_alg[sorted[i].id].ctg_id1) && (strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == true)) 
					{
						if( pair_alg[sorted[j].id].identity > (pair_alg[sorted[i].id].identity + 1) )
						{
							if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) )
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							if( strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == true )
							{
								is_x = true;
							}
							else is_x = false;
						}
					}
					else if((pair_alg[sorted[j].id].ctg_id1 == pair_alg[sorted[i].id].ctg_id1) && (f_loose_subset(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true))
					{
//					if( pair_alg[sorted[i].id].identity <= pair_alg[sorted[j].id].identity)
//						{
							if( (pair_alg[sorted[i].id].sign == 0 ) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							if( (pair_alg[sorted[j].id].ctg_id1 == pair_alg[sorted[i].id].ctg_id1) && (f_loose_subset(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true) ) is_x = true;
							else is_x = false;
//						}
					}
					else
					{
					}
				}
				j++;
			}
		}
	}
}

void check_pid_for_ortho(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val)
{
	int i, j;
	bool is_x;

	for( i = 0; i < (*num_pair); i++ )
	{
		if( (pair_alg[sorted[i].id].sign == OUT_PAR) || (pair_alg[sorted[i].id].sign == OUT_PAR_REV) ) 
		{
		}
		else if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == 1) )
		{
			if( pair_alg[sorted[i].id].sign == 0 )
			{
				pair_alg[sorted[i].id].sign = ORTHO_CANDI;
			}
			else if( pair_alg[sorted[i].id].sign == 1 )
			{
				pair_alg[sorted[i].id].sign = ORTHO_COMP_CANDI;
			}

			j = i+1;
			while( j < (*num_pair) )
			{
				if( sorted[j].id == sorted[i].id ) {}
				else
				{
					if((strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == true) || (strict_almost_equal(pair_alg[sorted[j].id].y, pair_alg[sorted[i].id].y) == true ))
					{
						if( pair_alg[sorted[j].id].identity > pair_alg[sorted[i].id].identity )
						{
							if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) )
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							if( strict_almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == true )
							{
								is_x = true;
							}
							else is_x = false;
						}
					}
					else if((check_contain(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true) || (check_contain(pair_alg[sorted[i].id].y, pair_alg[sorted[j].id].y, t_val) == true ))
					{
						if( pair_alg[sorted[i].id].identity <= ( pair_alg[sorted[j].id].identity - 1 ))
						{
							if( (pair_alg[sorted[i].id].sign == 0 ) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							if( check_contain(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true ) is_x = true;
							else is_x = false;
						}
					}
					else if( (check_contain(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true) || (check_contain(pair_alg[sorted[j].id].y, pair_alg[sorted[i].id].y, t_val) == true) )
					{
						if( width(pair_alg[sorted[j].id].x) <= ((int)(0.3 * ((float)(width(pair_alg[sorted[i].id].x))))) ) {}
						else if( width(pair_alg[sorted[j].id].x) <= ((int)(0.7 * ((float)(width(pair_alg[sorted[i].id].x))))) )
						{
							if( pair_alg[sorted[j].id].identity > (pair_alg[sorted[i].id].identity + 2) )
							{

								if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) )
								{
									pair_alg[sorted[i].id].sign = OUT_PAR;
								}
								else if( (pair_alg[sorted[i].id].sign == 1) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI) )
								{
									pair_alg[sorted[i].id].sign = OUT_PAR_REV;
								}

								if( check_contain(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true ) is_x = true;
								else is_x = false;
							}
						}
						else if( width(pair_alg[sorted[j].id].x) <= ((int)(0.85 *((float)(width(pair_alg[sorted[i].id].x))))) )
						{
							if( pair_alg[sorted[j].id].identity > (pair_alg[sorted[i].id].identity + 1) )
							{
								if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) )
								{
									pair_alg[sorted[i].id].sign = OUT_PAR;
								}
								else if( (pair_alg[sorted[i].id].sign == 1) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI) )
								{
									pair_alg[sorted[i].id].sign = OUT_PAR_REV;
								}

								if( check_contain(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true ) is_x = true;
								else is_x = false;
							}
						}
						else 
						{
							if( pair_alg[sorted[j].id].identity > pair_alg[sorted[i].id].identity )
							{
								if( (pair_alg[sorted[i].id].sign == 0) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI) )
								{
									pair_alg[sorted[i].id].sign = OUT_PAR;
								}
								else if( (pair_alg[sorted[i].id].sign == 1) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI) )
								{
									pair_alg[sorted[i].id].sign = OUT_PAR_REV;
								}

								if( check_contain(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true ) is_x = true;
								else is_x = false;
							}
						}
					}
					else if(f_loose_overlap(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x, t_val) == true)
					{
						if( (pair_alg[sorted[i].id].identity <= ( pair_alg[sorted[j].id].identity - 1 )) && (width(intersect(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x)) >= ((int)(0.3 * ((float)(width(pair_alg[sorted[i].id].x)))))))
						{
							if( (pair_alg[sorted[i].id].sign == 0 ) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}
							
							is_x = true;
						}
						else if( (pair_alg[sorted[i].id].identity > pair_alg[sorted[j].id].identity) && (width(intersect(pair_alg[sorted[i].id].x, pair_alg[sorted[j].id].x)) >= ((int)(0.2 * ((float)(width(pair_alg[sorted[i].id].x)))))))
						{
							if( (pair_alg[sorted[j].id].sign == 0 ) || (pair_alg[sorted[j].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[j].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[j].id].sign == 1 ) || (pair_alg[sorted[j].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[j].id].sign = OUT_PAR_REV;
							}
						}
					}
					else if(f_loose_overlap(pair_alg[sorted[i].id].y, pair_alg[sorted[j].id].y, t_val) == true) 
					{
						if( (pair_alg[sorted[i].id].identity <= ( pair_alg[sorted[j].id].identity - 1 )) && (width(intersect(pair_alg[sorted[i].id].y, pair_alg[sorted[j].id].y)) >= ((int)(0.3 * ((float)(width(pair_alg[sorted[i].id].y))))))) 
						{
							if( (pair_alg[sorted[i].id].sign == 0 ) || (pair_alg[sorted[i].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[i].id].sign == 1 ) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[i].id].sign = OUT_PAR_REV;
							}

							is_x = false;
						}
						else if( (pair_alg[sorted[i].id].identity > pair_alg[sorted[j].id].identity) && (width(intersect(pair_alg[sorted[i].id].y, pair_alg[sorted[j].id].y)) >= ((int)(0.2 * ((float)(width(pair_alg[sorted[i].id].y))))))) 
						{
							if( (pair_alg[sorted[j].id].sign == 0 ) || (pair_alg[sorted[j].id].sign == ORTHO_CANDI))
							{
								pair_alg[sorted[j].id].sign = OUT_PAR;
							}
							else if( (pair_alg[sorted[j].id].sign == 1 ) || (pair_alg[sorted[j].id].sign == ORTHO_COMP_CANDI))
							{
								pair_alg[sorted[j].id].sign = OUT_PAR_REV;
							}
						}
					}
					else
					{
					}
				}
				j++;
			}

			j = i+1;
			while( j < (*num_pair) )
			{
				if( sorted[j].id == sorted[i].id ) {}
				else
				{
					if((check_contain(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true) || (check_contain(pair_alg[sorted[j].id].y, pair_alg[sorted[i].id].y, t_val) == true ))
					{
						if( (pair_alg[sorted[i].id].sign == ORTHO) || (pair_alg[sorted[i].id].sign == ORTHO_COMP_CANDI) )
						{
							if( pair_alg[sorted[j].id].sign == 0 )
							{
								pair_alg[sorted[j].id].sign = OUT_PAR;
							}
							else if( pair_alg[sorted[j].id].sign == 1 )
							{
								pair_alg[sorted[j].id].sign = OUT_PAR_REV;
							}
						}
						else
						{
							if( (almost_equal(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x) == false) && (is_x == true))
							{
								if(check_contain(pair_alg[sorted[j].id].x, pair_alg[sorted[i].id].x, t_val) == true) 
								{
									if( pair_alg[sorted[j].id].identity <= pair_alg[sorted[i].id].identity )
									{
										if( pair_alg[sorted[j].id].sign == 0 )
										{
											pair_alg[sorted[j].id].sign = OUT_PAR;
										}
										else if( pair_alg[sorted[j].id].sign == 1 )
										{
											pair_alg[sorted[j].id].sign = OUT_PAR_REV;
										}
									}
								}
							}

							if((almost_equal(pair_alg[sorted[j].id].y, pair_alg[sorted[i].id].y) == false) && (is_x == false))
							{
								if(check_contain(pair_alg[sorted[j].id].y, pair_alg[sorted[i].id].y, t_val) == true) 
								{
									if( pair_alg[sorted[j].id].identity <= pair_alg[sorted[i].id].identity )
									{
										if( pair_alg[sorted[j].id].sign == 0 )
										{
											pair_alg[sorted[j].id].sign = OUT_PAR;
										}
										else if( pair_alg[sorted[j].id].sign == 1 )
										{
											pair_alg[sorted[j].id].sign = OUT_PAR_REV;
										}
									}
								}
							}
						}
					}
				}
				j++;
			}
		}
	}
}

bool check_contain(struct I reg1, struct I reg2, int t_val)
{
	bool res = false;
	int overlap_width;
	int min_width;

	min_width = width(reg1);

	if(f_loose_subset(reg1, reg2, t_val) == true)
	{
		res = true;	
	}
	else if(f_loose_overlap(reg1, reg2, t_val) == true)
	{
		overlap_width = width(intersect(reg1, reg2));	
		if( overlap_width >= (int)(0.90*((float)min_width)) )
		{
			res = true;
		}
	}

	return(res);
}

void add_ortho_intervals(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, int s_new, int e_new, struct DotList *algns, int sp_id, FILE *f, float avg_pid)
{
	int i = 0, j = 0, k = 0;
	int tmp1_id = 0, tmp2_id = 0;
	struct I cmp, org_cmp;
	int b = -1, e = -1, cur_b = -1, cur_e = -1;
	int cur_src = -1, cur_dst = -1;
	int left_diff = 0, right_diff = 0;
	int loc = 0, cur = 0, old = 0;
	struct short_alist *src_list;
	int num_list = 0;
	int count = 0;
	int low = 0, hi = 0;
	struct I temp, org_temp;
	struct I *aj;
	bool *is_in;
	float cur_pid = (float)0;

	if( (sp_id != SELF1) && (sp_id != SELF2) ) fatalf("species ID should be either SELF1 or SELF2, but %d here\n", sp_id); 

	src_list = (struct short_alist *) ckalloc(sizeof(struct short_alist) * (e_loc-s_loc+1));
	is_in = (bool *) ckalloc(sizeof(bool) * (e_new-s_new+1));
	aj = (struct I *) ckalloc(sizeof(struct I));

	for( i = 0; i < (e_new-s_new+1); i++ ) {
		is_in[i] = true;
	}

	for( i = 0; i < (e_loc-s_loc+1); i++ ) {
		src_list[i].id = 0;
		src_list[i].x = assign_I(0, 1);
		src_list[i].y = assign_I(0, 1);
		src_list[i].val = 0;
	}

	(*aj) = assign_I(0, 1);
	cmp = assign_I(0, 1);
	org_cmp = assign_I(0, 1);
	temp = assign_I(0, 1);
	org_temp = assign_I(0, 1);

	for( j = s_loc; j <= e_loc; j++ ) 
	{
		b = -1;
		e = -1;
		i = sorted[j].id;
		if( sp_id == SELF1 ) {
			low = algns[i].x.lower + algns[i].xl_diff;
			hi = algns[i].x.upper - algns[i].xr_diff;
		}
		else {
			low = algns[i].y.lower + algns[i].yl_diff;
			hi = algns[i].y.upper - algns[i].yr_diff;
		}

		if( hi > low ) cmp = assign_I(low, hi);

		if( hi <= low ) {}
		else if( algns[i].sign == DELETED ) {}
		else if( cmp.upper < ops[cur_id].srcStart ) {}
		else if( cmp.lower > ops[cur_id].srcEnd ) {}
		else {
			if( ((ops[cur_id].srcStart >= cmp.lower) && (ops[cur_id].srcStart <= cmp.upper)) ) 
			{
				if( sp_id == SELF1 ) {
					if( algns[i].sign == 0 ) b = find_yloc_one(algns[i], f, abs(ops[cur_id].srcStart - algns[i].x.lower), NO_GAP_INC);
					else if( algns[i].sign == 1 ) e = find_yloc_one(algns[i], f, abs(ops[cur_id].srcStart - algns[i].x.lower), NO_GAP_INC);
				}
				else {
					if( algns[i].sign == 0 ) {
						loc = find_yloc_one(algns[i], f, abs(ops[cur_id].srcStart - algns[i].y.lower), GAP_INC_IN_Y);
						b = find_xloc_one(algns[i], f, loc, GAP_INC);
					}
					else if( algns[i].sign == 1 ) {
						loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].srcStart), GAP_INC_IN_Y);
						e = find_xloc_one(algns[i], f, loc, GAP_INC);
					}
				}
			}

			if((ops[cur_id].srcEnd >= cmp.lower) && (ops[cur_id].srcEnd <= cmp.upper) ) 
			{
				if( sp_id == SELF1 ) {
					if( algns[i].sign == 0 ) e = find_yloc_one(algns[i], f, abs(ops[cur_id].srcEnd - algns[i].x.lower), NO_GAP_INC);
					else if( algns[i].sign == 1 ) b = find_yloc_one(algns[i], f, abs(ops[cur_id].srcEnd - algns[i].x.lower), NO_GAP_INC);
				}
				else {
					if( algns[i].sign == 0 ) {
						loc = find_yloc_one(algns[i], f, abs(ops[cur_id].srcEnd - algns[i].y.lower), GAP_INC_IN_Y);
						e = find_xloc_one(algns[i], f, loc, GAP_INC);
					}
					else if( algns[i].sign == 1 ) {
						loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].srcEnd), GAP_INC_IN_Y);
						b = find_xloc_one(algns[i], f, loc, GAP_INC);
					
					}
				}
			}

			if( (b == -1) && (e != -1) ) {
				if( sp_id == SELF1 ) b = algns[i].y.lower + algns[i].yl_diff;
				else b = algns[i].x.lower + algns[i].xl_diff;
			}
			else if( (b != -1) && (e == -1) ) {
				if( sp_id == SELF1 ) e = algns[i].y.upper - algns[i].yr_diff;
				else e = algns[i].x.upper - algns[i].xr_diff;
			}
			else {
				if( subset(cmp, assign_I(ops[cur_id].srcStart, ops[cur_id].srcEnd)) == true ) {
					if( sp_id == SELF1 ) {
						b = algns[i].y.lower + algns[i].yl_diff;
						e = algns[i].y.upper - algns[i].yr_diff;
					}
					else {
						b = algns[i].x.lower + algns[i].xl_diff;
						e = algns[i].x.upper - algns[i].xr_diff;
					}
				}
			}

			if( (b != -1) && (e != -1) && (b < e)) {
				src_list[num_list].id = i;
				src_list[num_list].x = assign_I(cmp.lower, cmp.upper);
				src_list[num_list].y = assign_I(b, e);	
				src_list[num_list].val = b;
				num_list++;
			}
		}
		i++;
	}

	count = 0;
	if( num_list > 1 ) {
		sort_by_loc_short_alist(src_list, num_list);
		for( j = 0; j < num_list; j++ ) src_list[j].val = algns[src_list[j].id].identity;
		for( j = 1; j < num_list; j++ ) {
			if( subset(src_list[j].y, src_list[count].y) == true ) {}
			else if( proper_overlap(src_list[j].y, src_list[count].y) == true ) {
				src_list[count].val = ((src_list[count].val * width(src_list[count].y)) + (src_list[j].val * width(src_list[j].y))) / (width(src_list[count].y) + width(src_list[j].y));

				if( src_list[j].y.upper > src_list[count].y.lower ) src_list[count].y = assign_I(src_list[count].y.lower, src_list[j].y.upper);
				else {
					fatalf("sorting error in %d-%d and %d-%d\n", src_list[count].y.lower, src_list[count].y.upper, src_list[j].y.lower, src_list[j].y.upper); 
				}

				if( width(src_list[count].x) < width(src_list[j].x) ) {
					src_list[count].x = assign_I(src_list[j].x.lower, src_list[j].x.upper);
					src_list[count].id = src_list[j].id;
				}
			}
			else {
				count++;
				src_list[count] = assign_alist(src_list[j]);
			}
		}
		count++;
	}
	else if( num_list == 1 ) {
		src_list[count].val = algns[src_list[0].id].identity;
		count++;
	}

	num_list = count;

	cmp = assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd);
	for( j = 0; j < (e_new-s_new+1); j++ ) {
		if( is_in[j] == true ) {
			tmp1_id = sorted[j+s_new].id;
			if( sp_id == SELF1 ) {
				org_temp = assign_I(algns[tmp1_id].x.lower + algns[tmp1_id].xl_diff, algns[tmp1_id].x.upper - algns[tmp1_id].xr_diff);
			}
			else if( sp_id == SELF2 ) {
				org_temp = assign_I(algns[tmp1_id].y.lower + algns[tmp1_id].yl_diff, algns[tmp1_id].y.upper - algns[tmp1_id].yr_diff);
			}

			if( org_temp.upper <= org_temp.lower ) {
				if( sp_id == SELF1 ) b = (-1) * width(algns[tmp1_id].x);
				else if( sp_id == SELF2 ) b = (-1) * width(algns[tmp1_id].y);
			}
			else if( proper_overlap(org_temp, cmp) == true ) {
				b = width(intersect(org_temp, cmp));
			}
			else {
				if( abs(org_temp.lower - cmp.upper) > abs(org_temp.upper - cmp.lower) ) 
				{
					b = (-1) * abs(org_temp.upper - cmp.lower);
				}
				else {
					b = (-1) * abs(org_temp.lower - cmp.upper);
				}
			}

			for( k = (j+1); k < (e_new-s_new+1); k++ ) {
				if( is_in[k] == false ) {}
				else {
					tmp2_id = sorted[k+s_new].id;
					if( sp_id == SELF1 ) {
						temp = assign_I(algns[tmp2_id].x.lower + algns[tmp2_id].xl_diff, algns[tmp2_id].x.upper - algns[tmp2_id].xr_diff);	
					}
					else if( sp_id == SELF2 ) {
						temp = assign_I(algns[tmp2_id].y.lower + algns[tmp2_id].yl_diff, algns[tmp2_id].y.upper - algns[tmp2_id].yr_diff);	
					}

					if( temp.upper <= temp.lower ) {
						if( sp_id == SELF1 ) e = (-1) * width(algns[tmp2_id].x);
						else if( sp_id == SELF2 ) e = (-1) * width(algns[tmp2_id].y);
					}
					else if( proper_overlap(temp, cmp) == true ) {
						e = width(intersect(temp, cmp));
					}
					else {
						if( abs(temp.lower - cmp.upper) > abs(temp.upper - cmp.lower) ) {
							e = (-1) * abs(temp.upper - cmp.lower);
						}
						else {
							e = (-1) * abs(temp.lower - cmp.upper);
						}
					}

					if( ( ( algns[tmp1_id].l_id == -1 ) && ( algns[tmp2_id].l_id == algns[tmp1_id].index ) ) || ( (algns[tmp1_id].l_id != -1) && ( (algns[tmp2_id].l_id == algns[tmp1_id].l_id) || (algns[tmp1_id].l_id == algns[tmp2_id].index) ) ) ){
						if( b >= e ) {
							is_in[k] = false;
						}
						else {
							is_in[j] = false;
						}
					}
				}
			}
		}
	}

	for( j = s_new; j <= e_new; j++ ) {
		b = -1;
		e = -1;
		i = sorted[j].id;
		if( sp_id == SELF1 ) {
			cmp = assign_I(algns[i].x.lower + algns[i].xl_diff, algns[i].x.upper - algns[i].xr_diff);
			org_cmp = assign_I(algns[i].x.lower, algns[i].x.upper);
		}
		else {
			cmp = assign_I(algns[i].y.lower + algns[i].yl_diff, algns[i].y.upper - algns[i].yr_diff);
			org_cmp = assign_I(algns[i].y.lower, algns[i].y.upper);
		}

		if( is_in[j-s_new] == false ) {}
		else if( org_cmp.upper < ops[cur_id].dstStart ) {}
		else if( org_cmp.lower > ops[cur_id].dstEnd ) {}
		else {
			if( ((ops[cur_id].dstStart >= org_cmp.lower) && (ops[cur_id].dstStart <= org_cmp.upper)) ) 
			{
				if( sp_id == SELF1 ) {
					if( algns[i].init_sign == 0 ) b = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].x.lower), NO_GAP_INC);
					else if( algns[i].init_sign == 1 ) e = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].x.lower), NO_GAP_INC);
				}
				else {
					if( algns[i].init_sign == 0 ) {
						loc = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].y.lower), GAP_INC_IN_Y);
						b = find_xloc_one(algns[i], f, loc, GAP_INC);
					}
					else if( algns[i].init_sign == 1 ) {
						loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].dstStart), GAP_INC_IN_Y);
						e	= find_xloc_one(algns[i], f, loc, GAP_INC);
					}
				}
			}

			if((ops[cur_id].dstEnd >= org_cmp.lower) && (ops[cur_id].dstEnd <= org_cmp.upper) ) 
			{
				if( sp_id == SELF1 ) {
					if( algns[i].init_sign == 0 ) e = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].x.lower), NO_GAP_INC);
					else if( algns[i].init_sign == 1 ) b = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].x.lower), NO_GAP_INC);
				}
				else {
					if( algns[i].init_sign == 0 ) {
						loc = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].y.lower), GAP_INC_IN_Y);
						e = find_xloc_one(algns[i], f, loc, GAP_INC);
					}
					else if( algns[i].init_sign == 1 ) {
						loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].dstEnd), GAP_INC_IN_Y);
						b = find_xloc_one(algns[i], f, loc, GAP_INC);
					}
				}
			}

			if( (b == -1) && (e != -1) ) {
				if( sp_id == SELF1 ) b = algns[i].y.lower;
				else b = algns[i].x.lower;
			}
			else if( (b != -1) && (e == -1) ) {
				if( sp_id == SELF1 ) e = algns[i].y.upper;
				else e = algns[i].x.upper;
			}
			else if( subset(org_cmp, assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd)) == true ) {
				if( sp_id == SELF1 ) {
					b = algns[i].y.lower; 
					e = algns[i].y.upper;
				}
				else {
					b = algns[i].x.lower;
					e = algns[i].x.upper;
				}
			}

			if( sp_id == SELF1 ) {
				cur_b = algns[i].x.lower + algns[i].xl_diff;
				cur_e = algns[i].x.upper - algns[i].xr_diff;
			}
			else {
				cur_b = algns[i].y.lower + algns[i].yl_diff;
				cur_e = algns[i].y.upper - algns[i].yr_diff;
			}

			if( proper_overlap(org_cmp, assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd)) == false ) {}
			else if( (b != -1) && (e != -1) && ((e - b) >= ERR_LG_TH) ) {
				org_temp = assign_I(b, e);
				if( sp_id == SELF1 ) {
					cur_pid = cal_pid_part_algn(algns, i, abs(b-algns[i].y.lower), abs(e-algns[i].y.upper), f, SELF2);
				}
				else {
					cur_pid = cal_pid_part_algn(algns, i, abs(b-algns[i].x.lower), abs(e-algns[i].x.upper), f, SELF1);
				}

				if( (is_within_bound(ops[cur_id], org_temp, src_list, num_list, aj, algns, i, sp_id, f) == true) ) {
//					if( ((algns[i].sign == DELETED) && (algns[i].identity >= (((int)avg_pid)-PID_DIFF))) || (algns[i].sign != DELETED) ) 

					if( ((algns[i].sign == DELETED) && ((int)(cur_pid+0.5) >= ((int)(avg_pid+0.5)-PID_DIFF))) || (algns[i].sign != DELETED) ) 
//					if( ((algns[i].sign == DELETED) && (((ops[cur_id].dir < 0) && ((int)(cur_pid+0.5) >= ((int)avg_pid+0.5)-PID_DIFF)) || ((int)(cur_pid+0.5) >= (int)(avg_pid+0.5)))) || (algns[i].sign != DELETED) ) 
					{
						if( algns[i].sign == DELETED ) {
							algns[i].sign = algns[i].init_sign;
							if( sp_id == SELF1 ) {
								left_diff = width(algns[i].x);
								right_diff = width(algns[i].x);
								cur_e = algns[i].x.lower;
								cur_b = algns[i].x.upper;
							}
							else {
								left_diff = width(algns[i].y);
								right_diff = width(algns[i].y);
								cur_e = algns[i].y.lower;
								cur_b = algns[i].y.upper;
							}
						}
						else {
							if( sp_id == SELF1 ) {
								left_diff = algns[i].xl_diff;
								right_diff = algns[i].xr_diff;
							}
							else {
								left_diff = algns[i].yl_diff;
								right_diff = algns[i].yr_diff;
							}
						}

						if( sp_id == SELF1 ) {
							cur_src = algns[i].x.lower;
							cur_dst = algns[i].x.upper;
						}
						else {
							cur_src = algns[i].y.lower;
							cur_dst = algns[i].y.upper;
						}
						
						if( cur_e <= cur_b ) {
							if( ops[cur_id].dstStart > cur_src )
							{
								if( left_diff > (ops[cur_id].dstStart - cur_src) ) {
									if( sp_id == SELF1 ) {
										algns[i].xl_diff = ops[cur_id].dstStart - algns[i].x.lower;
										cur = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].x.lower), NO_GAP_INC);
										if( algns[i].init_sign == 0 ) {
											old = algns[i].y.lower;
											if( cur > old ) algns[i].yl_diff = cur - old;
										}
										else if( algns[i].init_sign == 1 ) {
											old = algns[i].y.upper;
											if( cur < old ) algns[i].yr_diff = old - cur;
										}
									}
									else if( sp_id == SELF2 ) {
										algns[i].yl_diff = ops[cur_id].dstStart - algns[i].y.lower;
										if( algns[i].init_sign == 0 ) {
											loc = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].y.lower), GAP_INC_IN_Y);
											cur = find_xloc_one(algns[i], f, loc, GAP_INC);
											old = algns[i].x.lower;
											if( cur > old ) algns[i].xl_diff = cur - old;
										}
										else if( algns[i].init_sign == 1 ) {
											loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].dstStart), GAP_INC_IN_Y);
											cur = find_xloc_one(algns[i], f, loc, GAP_INC);
											old = algns[i].x.upper;
											if( cur < old ) algns[i].xr_diff = old - cur;
										}
									}
								}
							}
							else {
								if( sp_id == SELF1 ) {
									algns[i].xl_diff = 0;
									if( algns[i].init_sign == 0 ) {
										algns[i].yl_diff = 0;
									}
									else if( algns[i].init_sign == 1 ) {
										algns[i].yr_diff = 0;
									}
								}
								else if( sp_id == SELF2 ) {
									algns[i].yl_diff = 0;
									if( algns[i].init_sign == 0 ) {
										algns[i].xl_diff = 0;
									}
									else if( algns[i].init_sign == 1 ) {
										algns[i].xr_diff = 0;
									}
								}
							}
	
							if(ops[cur_id].dstEnd < cur_dst) 
							{
								if( right_diff > (cur_dst - ops[cur_id].dstEnd) ) {
									if( sp_id == SELF1 ) {
										algns[i].xr_diff = algns[i].x.upper - ops[cur_id].dstEnd;
										cur = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].x.lower), NO_GAP_INC);
										if( algns[i].init_sign == 0 ) {
											old = algns[i].y.upper;
											if( cur < old ) algns[i].yr_diff = old - cur;
										}
										else if( algns[i].init_sign == 1 ) {
											old = algns[i].y.lower;
											if( cur > old ) algns[i].yl_diff = cur - old;
										}
									}
									else if( sp_id == SELF2 ) {
										algns[i].yr_diff = algns[i].y.upper - ops[cur_id].dstEnd;
										if( algns[i].init_sign == 0 ) {
											loc = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].y.lower), GAP_INC_IN_Y);
											cur = find_xloc_one(algns[i], f, loc, GAP_INC);
											old = algns[i].x.upper;
											if( cur < old ) algns[i].xr_diff = old - cur;
										}
										else if( algns[i].init_sign == 1 ) {
											loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].dstEnd), GAP_INC_IN_Y);
											cur = find_xloc_one(algns[i], f, loc, GAP_INC);
											old = algns[i].x.lower;
											if( cur > old ) algns[i].xl_diff = cur - old;
										}
									}
								}	
							}
							else {
								if( sp_id == SELF1 ) {
									algns[i].xr_diff = 0;
									if( algns[i].init_sign == 0 ) {
										algns[i].yr_diff = 0;
									}
									else if( algns[i].init_sign == 1 ) {
										algns[i].yl_diff = 0;
									}
								}
								else if( sp_id == SELF2 ) {
									algns[i].yr_diff = 0;
									if( algns[i].init_sign == 0 ) {
										algns[i].xr_diff = 0;
									}
									else if( algns[i].init_sign == 1 ) {
										algns[i].xl_diff = 0;
									}
								}
							}
						}
						else {
							temp = assign_I(cur_b, cur_e);
							if( subset(assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd), temp) == true ) 
							{
							}
							else if ( cur_dst < cur_src ) {
								if( debug_mode == TRUE ) {
									printf("%d-%d; %d\n", cur_src, cur_dst, i);
								}
							}
							else if ( proper_overlap(assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd), assign_I(cur_src, cur_dst)) == true){
								if( ops[cur_id].dstStart <= cur_src )  
								{
									if( sp_id == SELF1 ) {
										algns[i].xl_diff = 0;
										if( algns[i].init_sign == 0 ) algns[i].yl_diff = 0;
										else if( algns[i].init_sign == 1 ) algns[i].yr_diff = 0;
									}
									else {
										algns[i].yl_diff = 0;
										if( algns[i].init_sign == 0 ) algns[i].xl_diff = 0;
										else if( algns[i].init_sign == 1 ) algns[i].xr_diff = 0;
									}
								}
								else 
								{
									if( left_diff > (ops[cur_id].dstStart - cur_src) ) 
									{
										if( sp_id == SELF1 ) {
											algns[i].xl_diff = ops[cur_id].dstStart - algns[i].x.lower;
											cur = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].x.lower), NO_GAP_INC);
											if( algns[i].init_sign == 0 ) {
												old = algns[i].y.lower;
												if( cur > old ) algns[i].yl_diff = cur - old;
											}
											else if( algns[i].init_sign == 1 ) {
												old = algns[i].y.upper;
												if( cur < old ) algns[i].yr_diff = old - cur;
											}
										}
										else if( sp_id == SELF2 ) {
											algns[i].yl_diff = ops[cur_id].dstStart - algns[i].y.lower;
											if( algns[i].init_sign == 0 ) {
												loc = find_yloc_one(algns[i], f, abs(ops[cur_id].dstStart - algns[i].y.lower), GAP_INC_IN_Y);
												cur = find_xloc_one(algns[i], f, loc, GAP_INC);
												old = algns[i].x.lower;
												if( cur > old ) algns[i].xl_diff = cur - old;
											}
											else if( algns[i].init_sign == 1 ) {
												loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].dstStart), GAP_INC_IN_Y);
												cur = find_xloc_one(algns[i], f, loc, GAP_INC);
												old = algns[i].x.upper;
												if( cur < old ) algns[i].xr_diff = old - cur;
											}
										}
									}
								}

								if( ops[cur_id].dstEnd >= cur_dst ) 
								{
									if( sp_id == SELF1 ) {
										algns[i].xr_diff = 0;
										if( algns[i].init_sign == 0 ) algns[i].yr_diff = 0;
										else if( algns[i].init_sign == 1 ) algns[i].yl_diff = 0;
									}
									else {
										algns[i].yr_diff = 0;
										if( algns[i].init_sign == 0 ) algns[i].xr_diff = 0;
										else if( algns[i].init_sign == 1 ) algns[i].xl_diff = 0;
									}
								}
								else 
								{
									if( right_diff > (cur_dst - ops[cur_id].dstEnd) ) {
										if( sp_id == SELF1 ) {
											algns[i].xr_diff = algns[i].x.upper - ops[cur_id].dstEnd;
											cur = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].x.lower), NO_GAP_INC);
											if( algns[i].init_sign == 0 ) {
												old = algns[i].y.upper;
												if( cur < old ) algns[i].yr_diff = old - cur;
											}
											else if( algns[i].init_sign == 1 ) {
												old = algns[i].y.lower;
												if( cur > old ) algns[i].yl_diff = cur - old;
											}
										}
										else if( sp_id == SELF2 ) {
											algns[i].yr_diff = algns[i].y.upper - ops[cur_id].dstEnd;
											if( algns[i].init_sign == 0 ) {
												loc = find_yloc_one(algns[i], f, abs(ops[cur_id].dstEnd - algns[i].y.lower), GAP_INC_IN_Y);
												cur = find_xloc_one(algns[i], f, loc, GAP_INC);
												old = algns[i].x.upper;
												if( cur < old ) algns[i].xr_diff = old - cur;
											}
											else if( algns[i].init_sign == 1 ) {
												loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops[cur_id].dstEnd), GAP_INC_IN_Y);
												cur = find_xloc_one(algns[i], f, loc, GAP_INC);
												old = algns[i].x.lower;
												if( cur > old ) algns[i].xl_diff = cur - old;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}	
	}

	free(is_in);
	free(aj);
	free(src_list);
}

bool is_within_bound(struct ops_list ops, struct I temp, struct short_alist *src_list, int num_list, struct I *aj, struct DotList *init_algns, int id, int sp_id, FILE *f)
{
	int b, e;
	int i;
	bool res = false;
	bool is_end = false;
	int hi, low;
	bool is_tandem = false;
	struct I reg1, reg2, temp_reg;
	int len1 = 0, len2 = 0, len = 0;
	int *res_b, *res_e;

	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));
	*res_b = -1;
	*res_e = -1;
	temp_reg = assign_I(0,1);
	reg1 = assign_I(ops.srcStart, ops.srcEnd);
	reg2 = assign_I(ops.dstStart, ops.dstEnd);

	if( ( abs(reg1.lower-reg2.upper) <= DEL_TH ) || (abs(reg1.upper - reg2.lower) <= DEL_TH) ) {
		is_tandem = true;
	}

	if( num_list >= 1 ) {
		low = src_list[0].y.lower;
		hi = src_list[0].y.upper;
	}
	else {
		low = 0;
		hi = 1;
		free(res_b);
		free(res_e);
		return res;
	}

	for( i = 1; i < num_list; i++ ) {
		if( hi < src_list[i].y.upper ) hi = src_list[i].y.upper;
		if( low > src_list[i].y.lower ) low = src_list[i].y.lower;
	}

	b = temp.upper+1;
	e = temp.lower-1;

	i = 0;
	while( (i < num_list)  && (is_end == false)) {
		if( src_list[i].y.lower >= src_list[i].y.upper ) {
			fatalf("id_ortho.c: empty interval (%d, %d)\n", src_list[i].y.lower, src_list[i].y.upper);
		}

		if( (proper_overlap(temp, src_list[i].y) == true) && (width(intersect(temp, src_list[i].y)) > ERR_LG_TH) ) {
			if( is_tandem == true ) {
				if( sp_id == SELF1 ) {
					find_overlapping_ends(src_list[i].y, temp, SELF2, init_algns, id, f, res_b, res_e);
				}
				else {
					find_overlapping_ends(src_list[i].y, temp, SELF1, init_algns, id, f, res_b, res_e);
				}

				len1 = width(intersect(temp, src_list[i].y));

				if( ((*res_b) == -1) || ((*res_e) == -1) || (((*res_e) - (*res_b)) <= DEL_TH)) {}
				else {
					temp_reg = assign_I(*res_b, *res_e);
//					if(proper_overlap(temp_reg, reg2) == true ) {
//						len2 = width(intersect(temp_reg, reg2));
						len2 = width(temp_reg);

						if( len1 > len2 ) len = len1;
						else len = len2;

						if( abs(len1 - len2) < (int)(0.3*(float)len) ) {
							res = true;
						}
					}
//				}
			}
			else {
				res = true;
			}
			
			if( res == false ) {}
			else if( subset(src_list[i].y, temp) == true ) {
				if( (b > src_list[i].y.lower) && (b != temp.lower) ) b = src_list[i].y.lower;
				if( (e < src_list[i].y.upper) && (e != temp.upper) ) e = src_list[i].y.upper;
			}
			else if(subset(temp, src_list[i].y) == true) {
				is_end = true;
				b = temp.lower;
				e = temp.upper;
			}
			else if((temp.lower >= src_list[i].y.lower) && (temp.lower <= src_list[i].y.upper)) {
				b = temp.lower;
				e = src_list[i].y.upper;
			}
			else if((temp.upper >= src_list[i].y.lower) && (temp.upper <= src_list[i].y.upper)) {
				b = src_list[i].y.lower;
				e = temp.upper;
			}

			if( (b == temp.lower) && (e == temp.upper) ) is_end = true;
		}
		i++;
	}

	if( b < low ) b = low;
	if( e > hi ) e = hi;

	if( e > b ) *aj = assign_I(b, e);
	else res = false;

	free(res_b);
	free(res_e);
	return(res);
}

int make_exc_list(struct I query, struct ops_list *ops, int cur_loc, int *exc_list)
{
	int i, j = 0;
	struct I dst;

	for( i = 0; i < cur_loc; i++ ) {
		dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
		if( f_loose_subset(dst, query, LOOSE) == true ) {
			exc_list[j] = i;
			j++;
		}
	}

	return(j);
}

void adjust_boundary(int id, struct DotList *temp_algns, struct DotList *algns, int num_algns, FILE *f)
{
	int i, j, num_ortho_algns;
	struct DotList *ortho_algns;
	struct I query, cur_reg;
	int s_loc, e_loc;
	struct slist *sorted;
	int b, e;
	int cur, old;

	j = 0;
	ortho_algns = (struct DotList *) ckalloc(num_algns * sizeof(struct DotList));
	for( i = 0; i < num_algns; i++ ) {
		if( algns[i].sp_id == PAIR ) {
			if( (algns[i].sign != DELETED) && (temp_algns[id].index != algns[i].index) ) {
				assign_algn(ortho_algns, j, algns[i]); 
				j++;
			}
		}
	}
	num_ortho_algns = j;

	sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_ortho_algns);
	for( i = 0; i < num_ortho_algns; i++ ) sorted[i].id = i;

	if( (temp_algns[id].x.upper-temp_algns[id].xr_diff-temp_algns[id].x.lower-temp_algns[id].xl_diff) > 0 )
	{
		sort_init_algns(sorted, ortho_algns, num_ortho_algns, INIT_PAIR); 
		query = assign_I(temp_algns[id].x.lower+temp_algns[id].xl_diff, temp_algns[id].x.upper-temp_algns[id].xr_diff);
		s_loc = search_range_b(sorted, ortho_algns, num_ortho_algns, query.lower, SELF1);
		e_loc = search_range_e(sorted, ortho_algns, num_ortho_algns, query.upper, SELF1);
		for( j = s_loc; j <= e_loc; j++ ) {
			i = sorted[j].id;
			b = ortho_algns[i].x.lower+ortho_algns[i].xl_diff;
			e = ortho_algns[i].x.upper-ortho_algns[i].xr_diff;
			if( (e - b) < DEL_TH ) {}
			else {
				cur_reg = assign_I(b, e);
				if( f_loose_subset(query, cur_reg, STRICT) == true ) {
					temp_algns[id].sign = DELETED;
				}
				else if( subset(cur_reg, query) == true ) {
					temp_algns[id].sign = DELETED;
				}
				else if( proper_overlap(query, cur_reg) == true ) {
					if( query.lower >  cur_reg.lower ) {
						temp_algns[id].xl_diff = cur_reg.upper - temp_algns[id].x.lower;
						if( temp_algns[id].sign == 0 ) {
							if( (temp_algns[id].xl_offset == 0) && (temp_algns[id].xr_offset == 0) ) {
								old = temp_algns[id].y.lower;
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.upper - temp_algns[id].x.lower), NO_GAP_INC);
							}
							else {
								old = find_yloc_one(temp_algns[id], f, temp_algns[id].xl_offset, NO_GAP_INC);
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.upper - temp_algns[id].x.lower + temp_algns[id].xl_offset), NO_GAP_INC);
							}

							if( cur > old ) temp_algns[id].yl_diff = cur - old;
						}
						else if( temp_algns[id].sign == 1 ) {
							if( (temp_algns[id].xl_offset == 0) && (temp_algns[id].xr_offset == 0) ) {
								old = temp_algns[id].y.upper;
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.upper - temp_algns[id].x.lower), NO_GAP_INC);
							}
							else {
								old = find_yloc_one(temp_algns[id], f, temp_algns[id].xl_offset, NO_GAP_INC); // equal to y.upper
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.upper - temp_algns[id].x.lower + temp_algns[id].xl_offset), NO_GAP_INC);
							}

							if( cur < old ) temp_algns[id].yr_diff = old - cur;
						}
					}
					else if( query.upper < cur_reg.upper ) {
						temp_algns[id].xr_diff = temp_algns[id].x.upper - cur_reg.lower;
						if( temp_algns[id].sign == 0 ) {
							if( (temp_algns[id].xl_offset == 0) && (temp_algns[id].xr_offset == 0) ) {
								old = temp_algns[id].y.upper;
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.lower - temp_algns[id].x.lower), NO_GAP_INC);
							}
							else {
								old = find_yloc_one(temp_algns[id], f, temp_algns[id].x.upper - temp_algns[id].x.lower + temp_algns[id].xl_offset, NO_GAP_INC);
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.upper - temp_algns[id].x.lower + temp_algns[id].xl_offset), NO_GAP_INC);
							}

							if( cur < old ) temp_algns[id].yr_diff = old - cur;
						}
						else if( temp_algns[id].sign == 1 ) {
							if( (temp_algns[id].xl_offset == 0) && (temp_algns[id].xr_offset == 0) ) {
								old = temp_algns[id].y.lower;
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.lower - temp_algns[id].x.lower), NO_GAP_INC);
							}
							else {
								old = find_yloc_one(temp_algns[id], f, temp_algns[id].x.upper - temp_algns[id].x.lower + temp_algns[id].xl_offset, NO_GAP_INC);
								cur = find_yloc_one(temp_algns[id], f, abs(cur_reg.lower - temp_algns[id].x.lower + temp_algns[id].xl_offset), NO_GAP_INC);
							}

							if( cur > old ) temp_algns[id].yl_diff = cur - old;
						}
					}
				}
			}
		}
	}				
	else {
		temp_algns[id].sign = DELETED;
	}

	free(sorted);
	free(ortho_algns);	
}

bool check_source_existence(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, int s_new, int e_new, struct DotList *algns, int sp_id, FILE *f, float avg_pid, bool *skip_untagging)
{
	int i = 0, j = 0;
	struct short_alist *src_list, *dst_list;
	struct short_alist *src_new_list, *dst_new_list;
	int num_src_list = 0, num_dst_list = 0;
	int num_src_new = 0, num_dst_new = 0;
	int src_best = 0, dst_best = 0;
	int src_max_width = 0, dst_max_width = 0;
	struct I src_ops_int, dst_ops_int;
	bool res = true;
	int cur_val = 0;

	src_list = (struct short_alist *) ckalloc( (e_loc-s_loc+1) * sizeof(struct short_alist));
	src_new_list = (struct short_alist *) ckalloc( (e_loc-s_loc+1) * sizeof(struct short_alist));
	dst_list = (struct short_alist *) ckalloc( (e_new-s_new+1) * sizeof(struct short_alist));
	dst_new_list = (struct short_alist *) ckalloc( (e_new-s_new+1) * sizeof(struct short_alist));

  for( i = 0; i < (e_loc-s_loc+1); i++ ) {
    src_list[i].id = 0;
    src_list[i].x = assign_I(0, 1);
    src_list[i].y = assign_I(0, 1);
    src_list[i].val = 0;

    src_new_list[i].id = 0;
    src_new_list[i].x = assign_I(0, 1);
    src_new_list[i].y = assign_I(0, 1);
    src_new_list[i].val = 0;
  }

  for( i = 0; i < (e_new-s_new+1); i++ ) {
    dst_list[i].id = 0;
    dst_list[i].x = assign_I(0, 1);
    dst_list[i].y = assign_I(0, 1);
    dst_list[i].val = 0;

    dst_new_list[i].id = 0;
    dst_new_list[i].x = assign_I(0, 1);
    dst_new_list[i].y = assign_I(0, 1);
    dst_new_list[i].val = 0;
  }

	src_ops_int = assign_I(ops[cur_id].srcStart, ops[cur_id].srcEnd);

	num_src_list = chop_match_regions(ops, cur_id, sorted, s_loc, e_loc, algns, sp_id, f, src_list, true, true);
	num_src_new = chop_match_regions(ops, cur_id, sorted, s_loc, e_loc, algns, sp_id, f, src_new_list, true, false);
	for( i = 0; i < num_src_list; i++ ) {
		if( src_max_width < width(src_list[i].y) ) {
			src_max_width = width(src_list[i].y);
		}

		cur_val = (int) ((float)(src_max_width) * ((float)(src_list[i].val)/(float)100));
		if( src_best < cur_val ) {
			src_best = cur_val;
		}
	}
	
	cur_val = 0;
	dst_ops_int = assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd);
	num_dst_list = chop_match_regions(ops, cur_id, sorted, s_new, e_new, algns, sp_id, f, dst_list, false, true);
	num_dst_new = chop_match_regions(ops, cur_id, sorted, s_new, e_new, algns, sp_id, f, dst_new_list, false, false);
	for( i = 0; i < num_dst_list; i++ ) {
		if( dst_max_width < width(dst_list[i].y) ) {
			dst_max_width = width(dst_list[i].y);
		}

		cur_val = (int) ((float)(dst_max_width) * ((float)(dst_list[i].val)/(float)100));
		if( dst_best < cur_val ) {
			dst_best = cur_val;	
		}
	}

	if( src_max_width >= (int)(0.9 * (float)width(src_ops_int)) ) {
		res = true;	
	}
	else {
		if( dst_max_width >= (int)(0.9 * (float)width(dst_ops_int)) ) {
			res = false;
		}
		else if( src_best < dst_best ) {
			res = false;
		}
		else {
			res = true;
			if( debug_mode == TRUE ) {
				if( (src_max_width == 0) && (dst_max_width == 0) ) {
					printf("both source and target do not exist in ops %d-%d and %d-%d\n", ops[cur_id].srcStart, ops[cur_id].srcEnd, ops[cur_id].dstStart, ops[cur_id].dstEnd);
				}	
			}
		}
	}
	
	if( res == true ) {
		dst_best = 0;
		for( i = 0; i < num_dst_new; i++ ) {
			if( (dst_best < dst_new_list[i].val) && (width(dst_new_list[i].y) >= (int)(0.5 * (float)width(dst_ops_int))) ) {
				dst_best = dst_new_list[i].val;
			}
		}
		if( (int)((float)dst_best+0.5) < ((int)(avg_pid+0.5)-PID_DIFF) ) {
			*skip_untagging = true;
		}
	}
	else {
		dst_best = 0;
		for( i = 0; i < num_src_new; i++ ) {
			if( (dst_best < src_new_list[i].val) && (width(src_new_list[i].y) >= (int)(0.5 * (float)width(src_ops_int))) ) {
				dst_best = src_new_list[i].val;
			}
		}
		if( (int)((float)dst_best+0.5) < ((int)(avg_pid+0.5)-PID_DIFF) ) {
			*skip_untagging = true;
		}
	}

	for( i = 0; i < num_src_list; i++ ) {
		for( j = 0; j < num_dst_list; j++ ) {
			if( (proper_overlap(src_list[i].y, dst_list[j].y) == true) && (width(intersect(src_list[i].y, dst_list[j].y)) >= (int)(0.5 * (float)width(src_ops_int))) ) 
			{
				*skip_untagging = true;
			}
		}
	}

	free(src_list);
	free(dst_list);
	free(src_new_list);
	free(dst_new_list);
	return(res);
}

int chop_match_regions(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, struct DotList *algns, int sp_id, FILE *f, struct short_alist *src_list, bool is_source, bool for_untagging)
{
	int i = 0, j = 0;
	struct I cmp, org_cmp;
	int num_list = 0;
	int low = 0, hi = 0;
	struct I ops_reg;
	int *res_b, *res_e;

	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));

	if( (sp_id != SELF1) && (sp_id != SELF2) ) fatalf("species ID should be either SELF1 or SELF2, but %d here\n", sp_id);

	for( i = 0; i < (e_loc-s_loc+1); i++ ) {
		src_list[i].id = 0;
		src_list[i].x = assign_I(0, 1);
		src_list[i].y = assign_I(0, 1);
		src_list[i].val = 0;
	}

	cmp = assign_I(0, 1);
	org_cmp = assign_I(0, 1);
	ops_reg = assign_I(0, 1);
	if( is_source == true ) {
		ops_reg = assign_I(ops[cur_id].srcStart, ops[cur_id].srcEnd);
	}
	else {
		ops_reg = assign_I(ops[cur_id].dstStart, ops[cur_id].dstEnd);
	}

	for( j = s_loc; j <= e_loc; j++ )
	{
		*res_b = -1;
		*res_e = -1;
		i = sorted[j].id;
		if( sp_id == SELF1 ) {
			low = algns[i].x.lower + algns[i].xl_diff;
			hi = algns[i].x.upper - algns[i].xr_diff;
			org_cmp = assign_I(algns[i].x.lower, algns[i].x.upper);
		}
		else {
			low = algns[i].y.lower + algns[i].yl_diff;
			hi = algns[i].y.upper - algns[i].yr_diff;
			org_cmp = assign_I(algns[i].y.lower, algns[i].y.upper);
		}

		if( hi > low ) cmp = assign_I(low, hi);

		if( hi <= low ) {}
		else if( algns[i].sign == DELETED ) {}
		else if( cmp.upper < ops_reg.lower ) {}
		else if( cmp.lower > ops_reg.upper ) {}
		else {
			find_overlapping_ends(ops_reg, cmp, sp_id, algns, i, f, res_b, res_e);
		}

    if( ((*res_b) != -1) && ((*res_e) != -1) && ((*res_b) < (*res_e))) {
			if( for_untagging == true ) {
       	src_list[num_list].id = i;
     	  src_list[num_list].x = assign_I(cmp.lower, cmp.upper);
     	  src_list[num_list].y = assign_I(*res_b, *res_e);
				if( sp_id == SELF1 ) {
					src_list[num_list].val = (int)cal_pid_part_algn(algns, i, abs((*res_b) - algns[i].y.lower), abs(algns[i].y.upper - (*res_e)), f, SELF2);
				}
				else {
					src_list[num_list].val = (int)cal_pid_part_algn(algns, i, abs((*res_b) - algns[i].x.lower), abs(algns[i].x.upper - (*res_e)), f, SELF1);
				}
     	  num_list++;
			}
		}
		else if( for_untagging == false ) {
			*res_b = -2;
			*res_e = -2;
			find_overlapping_ends(ops_reg, org_cmp, sp_id, algns, i, f, res_b, res_e);
    	if( ((*res_b) != -2) && ((*res_e) != -2) && ((*res_b) < (*res_e))) {
     	  src_list[num_list].id = i;
     		src_list[num_list].x = assign_I(cmp.lower, cmp.upper);
     		src_list[num_list].y = assign_I(*res_b, *res_e);
				if( sp_id == SELF1 ) {
					src_list[num_list].val = (int)cal_pid_part_algn(algns, i, abs((*res_b) - algns[i].y.lower), abs(algns[i].y.upper - (*res_e)), f, SELF2);
				}
				else {
					src_list[num_list].val = (int)cal_pid_part_algn(algns, i, abs((*res_b) - algns[i].x.lower), abs(algns[i].x.upper - (*res_e)), f, SELF1);
				}
     		num_list++;
			}
		}
    i++;
  }

	free(res_b);
	free(res_e);
	return(num_list);
}

void find_overlapping_ends(struct I ops_reg, struct I cmp, int sp_id, struct DotList *algns, int i, FILE *f, int *res_b, int *res_e)
{
	int loc = 0;
	int b = -1, e = -1;
	bool for_untagging = true;

	if( ((*res_b) == -2) && ((*res_e) == -2)) {
		b = -2;
		e = -2; // when for_untagging == false
		for_untagging = false;
	}

	if( ((ops_reg.lower >= cmp.lower) && (ops_reg.lower <= cmp.upper)) )
	{
		if( sp_id == SELF1 ) {
			if( algns[i].init_sign == 0 ) b = find_yloc_one(algns[i], f, abs(ops_reg.lower - algns[i].x.lower), NO_GAP_INC);
			else if( algns[i].sign == 1 ) e = find_yloc_one(algns[i], f, abs(ops_reg.lower - algns[i].x.lower), NO_GAP_INC);	
		}
		else {
 	   if( algns[i].init_sign == 0 ) {
 		   loc = find_yloc_one(algns[i], f, abs(ops_reg.lower - algns[i].y.lower), GAP_INC_IN_Y);
        b = find_xloc_one(algns[i], f, loc, GAP_INC);
      }
      else if( algns[i].init_sign == 1 ) {
        loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops_reg.lower), GAP_INC_IN_Y);
        e = find_xloc_one(algns[i], f, loc, GAP_INC);
      }
    }
  }

  if((ops_reg.upper >= cmp.lower) && (ops_reg.upper <= cmp.upper) )
  {
    if( sp_id == SELF1 ) {
      if( algns[i].init_sign == 0 ) e = find_yloc_one(algns[i], f, abs(ops_reg.upper - algns[i].x.lower), NO_GAP_INC);
      else if( algns[i].init_sign == 1 ) b = find_yloc_one(algns[i], f, abs(ops_reg.upper - algns[i].x.lower), NO_GAP_INC);
    }
    else {
      if( algns[i].init_sign == 0 ) {
        loc = find_yloc_one(algns[i], f, abs(ops_reg.upper - algns[i].y.lower), GAP_INC_IN_Y);
        e = find_xloc_one(algns[i], f, loc, GAP_INC);
      }
      else if( algns[i].init_sign == 1 ) {
        loc = find_yloc_one(algns[i], f, abs(algns[i].y.upper - ops_reg.upper), GAP_INC_IN_Y);
        b = find_xloc_one(algns[i], f, loc, GAP_INC);
      }
    }
  }

	if( for_untagging == true ) {
 		if( (b == -1) && (e != -1) ) {
 	  	if( sp_id == SELF1 ) b = algns[i].y.lower + algns[i].yl_diff;
 	  	else b = algns[i].x.lower + algns[i].xl_diff;
 		}
 		else if( (b != -1) && (e == -1) ) {
 	  	if( sp_id == SELF1 ) e = algns[i].y.upper - algns[i].yr_diff;
 	  	else e = algns[i].x.upper - algns[i].xr_diff;
 		}
 		else {
 			if( subset(cmp, ops_reg) == true ) {
 	    	if( sp_id == SELF1 ) {
 	    	   b = algns[i].y.lower + algns[i].yl_diff;
 	     		 e = algns[i].y.upper - algns[i].yr_diff;
 	     	}
 	     	else {
 	      	b = algns[i].x.lower + algns[i].xl_diff;
 	      	e = algns[i].x.upper - algns[i].xr_diff;
 	     	}
 	   	}
	  }
	}
	else {
 		if( (b == -2) && (e != -2) ) {
 	  	if( sp_id == SELF1 ) b = algns[i].y.lower;
 	  	else b = algns[i].x.lower;
 		}
 		else if( (b != -2) && (e == -2) ) {
 	  	if( sp_id == SELF1 ) e = algns[i].y.upper;
 	  	else e = algns[i].x.upper;
 		}
 		else {
 			if( subset(cmp, ops_reg) == true ) {
 	    	if( sp_id == SELF1 ) {
 	    	   b = algns[i].y.lower;
 	     		 e = algns[i].y.upper;
 	     	}
 	     	else {
 	      	b = algns[i].x.lower;
 	      	e = algns[i].x.upper;
 	     	}
 	   	}
	  }
	}

	*res_b = b;
	*res_e = e;
}
