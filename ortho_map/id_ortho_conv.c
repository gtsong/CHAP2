#include "main.h"
#include "id_ortho.h"
#include "id_ortho_conv.h"
#include "regions.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "read_algn.h"
#include "util_algns.h"

extern int debug_mode;

int redo_dups_for_mtom_inc_conv(int num_ops, struct ops_list *ops, int num_init_algns, struct DotList **org_init_algns, FILE *f, int sp_id, int *num_alloc_blocks)
{
  struct DotList *ortho_algns;
  int num_ortho_algns = 0;
  float *temp_val;
	int cur_sp_id;
	int i = 0, j = 0;
	float avg_pid;
	struct slist *sorted;
  int *exc_list;
	int mode = -1;
	int num_added = 0;
	int num_algns = 0;
	int *num_blocks; // the number of allocated blocks for "ortho_algns"
	bool later_dup_exist = false;
	bool *is_untagged;

	num_algns = num_init_algns;

  temp_val = (float *) ckalloc(sizeof(float));
  num_blocks = (int *) ckalloc(sizeof(int));
  is_untagged = (bool *) ckalloc(sizeof(bool));

	if( (*num_alloc_blocks) <= 0 ) {
		*num_alloc_blocks = 1;
	}

  ortho_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * (*num_alloc_blocks));
  initialize_algns(ortho_algns, 0, (*num_alloc_blocks));
	*num_blocks = *num_alloc_blocks;

  if( sp_id == REF_SEQ ) cur_sp_id = SELF1;
  else cur_sp_id = sp_id;

  j = 0;
	*temp_val = (float)0;
  avg_pid = (float)0;
  num_ortho_algns = assign_ortho_algns(ortho_algns, *org_init_algns, num_init_algns, temp_val);
  avg_pid = *temp_val;
  free(temp_val);

  if( (num_ortho_algns <= 0) || (avg_pid < 0) ) {}
  else {
    if( debug_mode == TRUE ) printf("Average pid: %f\n", avg_pid);
    sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_ortho_algns);
		initialize_slist(sorted, 0, num_ortho_algns);
    if( cur_sp_id == SELF1 ) mode = INIT_SELF1;
    else if( cur_sp_id == SELF2 ) mode = INIT_SELF2;
    else {
      fatalf("invalid species id: it must be either SELF1(0) or SELF2(1), but %d here\n", cur_sp_id);
    }

    sort_init_algns(sorted, ortho_algns, num_ortho_algns, mode); // sort by original coordinates in the initial dot plot
    exc_list = (int *) ckalloc(sizeof(int) * num_ops);

    for( i = 0; i < num_ops; i++ ) exc_list[i] = 0;

    for( i = num_ops-1; i >= 0; i-- ) {
			for( j = 0; j < num_ortho_algns; j++ ) sorted[j].id = j;
			sort_init_algns(sorted, ortho_algns, num_ortho_algns, mode);
			num_added = 0;

			*is_untagged = false;
      if(((ops[i].sign == '+') || (ops[i].sign == '-')) && (ops[i].sp_id == cur_sp_id)) {
        add_ortho_intervals_dups(ops, i, sorted, ortho_algns, num_ortho_algns, cur_sp_id, exc_list, f, avg_pid);
 			  for( j = 0; j < num_ortho_algns; j++ ) {
  			 	adjust_boundary_init(*org_init_algns, ortho_algns, j, num_algns, f, SELF2); 
    		}
      }
			else if(((ops[i].sign == 'c') || (ops[i].sign == 'v') ) && (ops[i].sp_id == cur_sp_id) && (ops[i].pid >= avg_pid)) {
				later_dup_exist = false;
				later_dup_exist = check_later_dup_existence(ops, i, num_ops);	
				if( later_dup_exist == false ) {
					num_added = untag_ortho_intervals_conv(ops, i, sorted, &ortho_algns, num_ortho_algns, cur_sp_id, f, org_init_algns, num_algns, num_alloc_blocks, num_blocks, avg_pid, is_untagged);

					if( num_added > 0 ) {
						num_algns = num_algns + num_added;
						num_ortho_algns = num_ortho_algns + num_added;
						sorted = (struct slist *) ckrealloc(sorted, num_ortho_algns * sizeof(struct slist));

						initialize_slist(sorted, 0, num_ortho_algns);
						for( j = 0; j < num_ortho_algns; j++ ) {
							if( (ortho_algns[j].x.lower < 0) || (ortho_algns[j].y.lower < 0 ) || ((ortho_algns[j].x.upper - ortho_algns[j].x.lower) < 0) || ((ortho_algns[j].y.upper - ortho_algns[j].y.lower) < 0)) {
								fatalf("unexpected alignment %d: %d-%d %d-%d\n", j, ortho_algns[j].x.lower, ortho_algns[j].x.upper, ortho_algns[j].y.lower, ortho_algns[j].y.upper);
							} 
						}
						sort_init_algns(sorted, ortho_algns, num_ortho_algns, mode);
					}
					if( (*is_untagged) == true ) {
						ops[i].dir = ops[i].dir * (-1); // mark ops for having untagged alignments
					}
					add_ortho_intervals_dups(ops, i, sorted, ortho_algns, num_ortho_algns, cur_sp_id, exc_list, f, avg_pid);
					if( (*is_untagged) == true ) {
						ops[i].dir = ops[i].dir * (-1); // unmark ops for having untagged alignments
					}

 			   	for( j = 0; j < num_ortho_algns; j++ ) {
  			  	adjust_boundary_init(*org_init_algns, ortho_algns, j, num_algns, f, SELF2); 
    			}
				}
			}
		}

    free(exc_list);
    free(sorted);
  }

  free(ortho_algns);
	free(num_blocks);
	free(is_untagged);

	return(num_algns);
}

int untag_ortho_intervals_conv(struct ops_list *ops, int id, struct slist *sorted, struct DotList **ortho_algns, int num_ortho_algns, int cur_sp_id, FILE *f, struct DotList **org_init_algns, int num_algns, int *num_alloc_blocks, int *num_blocks, float avg_pid, bool *is_untagged)
{
  struct I query, tmp;
  int s_loc = 0, e_loc = 0, s_d = 0, e_d = 0;
	int num_added = 0;
	bool src_exist = true;
	bool *skip_untagging;

	skip_untagging = (bool *) ckalloc(sizeof(bool));
	tmp = assign_I(0, 1);
  query = assign_I(ops[id].srcStart, ops[id].srcEnd);
  s_loc = search_range_b(sorted, *ortho_algns, num_ortho_algns, query.lower, cur_sp_id);
  e_loc = search_range_e(sorted, *ortho_algns, num_ortho_algns, query.upper, cur_sp_id);
  query = assign_I(ops[id].dstStart, ops[id].dstEnd);
  s_d = search_range_b(sorted, *ortho_algns, num_ortho_algns, query.lower, cur_sp_id);
  e_d = search_range_e(sorted, *ortho_algns, num_ortho_algns, query.upper, cur_sp_id);
	
	*skip_untagging = false;
	src_exist = check_source_existence(ops, id, sorted, s_loc, e_loc, s_d, e_d, *ortho_algns, cur_sp_id, f, avg_pid, skip_untagging);
	
	if( src_exist == true ) {
  	query = assign_I(ops[id].dstStart, ops[id].dstEnd);
	}
	else {
  	query = assign_I(ops[id].srcStart, ops[id].srcEnd);
    tmp = assign_I(ops[id].srcStart, ops[id].srcEnd);    
		ops[id].srcStart = ops[id].dstStart;    
		ops[id].srcEnd = ops[id].dstEnd;
    ops[id].dstStart = tmp.lower;
    ops[id].dstEnd = tmp.upper;
    tmp = assign_I(ops[id].src_b, ops[id].src_e);    
		ops[id].src_b = ops[id].dst_b;    
		ops[id].src_e = ops[id].dst_e;
    ops[id].dst_b = tmp.lower;
    ops[id].dst_e = tmp.upper;
	}

	if( (*skip_untagging) == true ) {

	}
	else {
		s_loc = search_range_b(sorted, *ortho_algns, num_ortho_algns, query.lower, cur_sp_id);
		e_loc = search_range_e(sorted, *ortho_algns, num_ortho_algns, query.upper, cur_sp_id);

		num_added = untag_ortho_intervals(ops, id, sorted, s_loc, e_loc, ortho_algns, num_ortho_algns, cur_sp_id, f, org_init_algns, num_algns, num_alloc_blocks, num_blocks, is_untagged); // sorted, ortho_algns, init_algns are updated
	}

	free(skip_untagging);
	return(num_added);
}

int untag_ortho_intervals(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, struct DotList **ortho_algns, int num_ortho, int sp_id, FILE *f, struct DotList **init_algns, int num_init_algns, int *num_alloc_blocks, int *num_blocks, bool *is_untagged)
{
	int i = 0, j = 0;
	struct I cmp, org_cmp;
	int b = 0, e = 1;
	int loc = 0, cur = 0, old = 0;
	int low = 0, hi = 1;
	int num_algns = 0;
	int num_ortho_algns = 0;
	int num_added = 0;
	int cur_diff = 0;
	float cur_pid = (float)0;

	num_algns = num_init_algns;
	num_ortho_algns = num_ortho;
	cmp = assign_I(0, 1);
	org_cmp = assign_I(0, 1);

	*is_untagged = false;
	if( (sp_id != SELF1) && (sp_id != SELF2) ) fatalf("species ID should be either SELF1 or SELF2, but %d here\n", sp_id); 

	if( (ops[cur_id].sign != 'c') && (ops[cur_id].sign != 'v') ) {
		fatalf("op is not conversion: %d\n", ops[cur_id].sign);
	}

	for( j = s_loc; j <= e_loc; j++ ) 
	{
		b = -1;
		e = -1;
		i = sorted[j].id;
		if( sp_id == SELF1 ) {
			low = (*ortho_algns)[i].x.lower + (*ortho_algns)[i].xl_diff;
			hi = (*ortho_algns)[i].x.upper - (*ortho_algns)[i].xr_diff;
			org_cmp = assign_I((*ortho_algns)[i].x.lower, (*ortho_algns)[i].x.upper);
		}
		else {
			low = (*ortho_algns)[i].y.lower + (*ortho_algns)[i].yl_diff;
			hi = (*ortho_algns)[i].y.upper - (*ortho_algns)[i].yr_diff;
			org_cmp = assign_I((*ortho_algns)[i].y.lower, (*ortho_algns)[i].y.upper);
		}

		if( hi > low ) cmp = assign_I(low, hi);

		if( hi <= low ) {}
		else if( (*ortho_algns)[i].sign == DELETED ) {}
		else if( cmp.upper < ops[cur_id].dstStart ) {}
		else if( cmp.lower > ops[cur_id].dstEnd ) {}
		else {
			if( (ops[cur_id].dstStart <= cmp.lower) && (ops[cur_id].dstEnd >= cmp.upper) ) {
				(*ortho_algns)[i].sign = DELETED;
				// delete the entire alignment block
			}
			else if( (cmp.lower <= ops[cur_id].dstStart) && (cmp.upper >= ops[cur_id].dstEnd) ) 
			{
				cur_pid = cal_pid_part_algn((*ortho_algns), i, abs(ops[cur_id].dstStart - org_cmp.lower), abs(org_cmp.upper - ops[cur_id].dstEnd), f, sp_id);
				if( (int)(cur_pid+0.5) >= (int)(ops[cur_id].pid-0.5)) {

				}
				else {
					*is_untagged = true;
					if( (num_ortho_algns + num_added) >= (*num_blocks) ) {
						if( num_ortho_algns < 0 ) {
							fatal("alignments aren't read yet\n");
						}
						else if( (*num_blocks) <= 0 ) {
							*ortho_algns = ckalloc(ALLOC_UNIT * sizeof(struct DotList));
							*num_blocks = ALLOC_UNIT;
						}
						else {
							*ortho_algns = (struct DotList *) ckrealloc(*ortho_algns, ((*num_blocks) + ALLOC_UNIT) * sizeof(struct DotList));
							*num_blocks = (*num_blocks) + ALLOC_UNIT;
						}
						initialize_algns(*ortho_algns, (*num_blocks) - ALLOC_UNIT, *num_blocks);
					}

					if( (num_algns + num_added) >= (*num_alloc_blocks) ) {
						if( num_algns < 0 ) {
							fatal("alignments aren't read yet\n");
						}
						else if( (*num_alloc_blocks) == 0 ) {
							*init_algns = (struct DotList *) ckalloc(ALLOC_UNIT * sizeof(struct DotList));
							*num_alloc_blocks = ALLOC_UNIT;
						}
						else {
							*init_algns = (struct DotList *) ckrealloc(*init_algns, ((*num_alloc_blocks) + ALLOC_UNIT) * sizeof(struct DotList));
							*num_alloc_blocks = (*num_alloc_blocks) + ALLOC_UNIT;
						}
						initialize_algns(*init_algns, (*num_alloc_blocks) - ALLOC_UNIT, *num_alloc_blocks);
					}
//				*ortho_algns = ckrealloc(*ortho_algns, ((num_ortho_algns) + num_added + 1) * sizeof(struct DotList));
//				*init_algns = ckrealloc(*init_algns, (num_algns + num_added + 1) * sizeof(struct DotList));

					assign_algn(*ortho_algns, num_ortho_algns+num_added, (*ortho_algns)[i]);
					assign_algn(*init_algns, num_algns+num_added, (*ortho_algns)[i]);
					if( (*ortho_algns)[i].l_id == -1 ) {
						(*ortho_algns)[num_ortho_algns+num_added].l_id = (*ortho_algns)[i].index;
						(*init_algns)[num_algns+num_added].l_id = (*ortho_algns)[i].index;
					}
					else {
						(*ortho_algns)[num_ortho_algns+num_added].l_id = (*ortho_algns)[i].l_id;
						(*init_algns)[num_algns+num_added].l_id = (*ortho_algns)[i].l_id;
					}
					(*ortho_algns)[num_ortho_algns+num_added].index = num_algns+num_added;
					(*init_algns)[num_algns+num_added].index = num_algns+num_added;
					if( sp_id == SELF1 ) {
						cur_diff = (*ortho_algns)[num_ortho_algns+num_added].xl_diff;
						if( cur_diff < abs(ops[cur_id].dstEnd-(*ortho_algns)[i].x.lower) ) {
							(*ortho_algns)[num_ortho_algns+num_added].xl_diff = abs(ops[cur_id].dstEnd - (*ortho_algns)[i].x.lower);
							cur = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstEnd - (*ortho_algns)[i].x.lower), NO_GAP_INC);
							if( (*ortho_algns)[i].init_sign == 0 ) {
								old = (*ortho_algns)[i].y.lower;
								if( cur > old ) (*ortho_algns)[num_ortho_algns+num_added].yl_diff = cur - old;
							}
							else if( (*ortho_algns)[i].init_sign == 1 ) {
								old = (*ortho_algns)[i].y.upper;
								if( cur < old ) (*ortho_algns)[num_ortho_algns+num_added].yr_diff = old - cur;
							}
						}

						cur_diff = (*ortho_algns)[i].xr_diff;
						if( cur_diff < ((*ortho_algns)[i].x.upper - ops[cur_id].dstStart) ) {
							(*ortho_algns)[i].xr_diff = (*ortho_algns)[i].x.upper - ops[cur_id].dstStart;
							cur = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstStart - (*ortho_algns)[i].x.lower), NO_GAP_INC);
							if( (*ortho_algns)[i].init_sign == 0 ) {
								old = (*ortho_algns)[i].y.upper;
								if( cur < old ) (*ortho_algns)[i].yr_diff = old - cur;
							}
							else if( (*ortho_algns)[i].init_sign == 1 ) {
								old = (*ortho_algns)[i].y.lower;
								if( cur > old ) (*ortho_algns)[i].yl_diff = cur - old;
							}
						}
					}
					else {
						cur_diff = (*ortho_algns)[num_ortho_algns+num_added].yl_diff;
						if( cur_diff < abs(ops[cur_id].dstEnd - (*ortho_algns)[i].y.lower) ) {
							(*ortho_algns)[num_ortho_algns+num_added].yl_diff = abs(ops[cur_id].dstEnd - (*ortho_algns)[i].y.lower);	
							if( (*ortho_algns)[i].init_sign == 0 ) {
								loc = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstEnd - (*ortho_algns)[i].y.lower), GAP_INC_IN_Y);
								cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
								old = (*ortho_algns)[i].x.lower;
								if( cur > old ) (*ortho_algns)[num_ortho_algns+num_added].xl_diff = cur - old;
							}
							else if( (*ortho_algns)[i].init_sign == 1 ) {
								loc = find_yloc_one((*ortho_algns)[i], f, abs((*ortho_algns)[i].y.upper - ops[cur_id].dstEnd), GAP_INC_IN_Y);
								cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
								old = (*ortho_algns)[i].x.upper;
								if( cur < old ) (*ortho_algns)[num_ortho_algns+num_added].xr_diff = old - cur;
							}
						}

						cur_diff = (*ortho_algns)[i].yr_diff;
						if( cur_diff < ((*ortho_algns)[i].y.upper - ops[cur_id].dstStart) ) {
							(*ortho_algns)[i].yr_diff = (*ortho_algns)[i].y.upper - ops[cur_id].dstStart;
							if( (*ortho_algns)[i].init_sign == 0 ) {
								loc = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstStart - (*ortho_algns)[i].y.lower), GAP_INC_IN_Y);
								cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
								old = (*ortho_algns)[i].x.upper;
								if( cur < old ) (*ortho_algns)[i].xr_diff = old - cur;
							}
							else if( (*ortho_algns)[i].init_sign == 1 ) {
								loc = find_yloc_one((*ortho_algns)[i], f, abs((*ortho_algns)[i].y.upper - ops[cur_id].dstStart), GAP_INC_IN_Y);
								cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
								old = (*ortho_algns)[i].x.lower;
								if( cur > old ) (*ortho_algns)[i].xl_diff = cur - old;
							}
						}
					}
					(*init_algns)[num_algns+num_added].xl_diff = (*ortho_algns)[num_ortho_algns+num_added].xl_diff;
					(*init_algns)[num_algns+num_added].xr_diff = (*ortho_algns)[num_ortho_algns+num_added].xr_diff;
					(*init_algns)[num_algns+num_added].yl_diff = (*ortho_algns)[num_ortho_algns+num_added].yl_diff;
					(*init_algns)[num_algns+num_added].yr_diff = (*ortho_algns)[num_ortho_algns+num_added].yr_diff;
					(*init_algns)[(*ortho_algns)[i].index].xl_diff = (*ortho_algns)[i].xl_diff;
					(*init_algns)[(*ortho_algns)[i].index].xr_diff = (*ortho_algns)[i].xr_diff;
					(*init_algns)[(*ortho_algns)[i].index].yl_diff = (*ortho_algns)[i].yl_diff;
					(*init_algns)[(*ortho_algns)[i].index].yr_diff = (*ortho_algns)[i].yr_diff;
					num_added++;
				}
			}
			else {
				if( ((ops[cur_id].dstEnd >= cmp.lower) && (ops[cur_id].dstEnd <= cmp.upper)) ) 
				{
					cur_pid = cal_pid_part_algn((*ortho_algns), i, 0, abs(org_cmp.upper - ops[cur_id].dstEnd), f, sp_id);
					if( (int)(cur_pid+0.5) >= (int)(ops[cur_id].pid-0.5)) {

					}
					else {
						*is_untagged = true;
						if( sp_id == SELF1 ) {
							cur_diff = (*ortho_algns)[i].xl_diff;
							if( cur_diff < abs(ops[cur_id].dstEnd - (*ortho_algns)[i].x.lower) ) {
								(*ortho_algns)[i].xl_diff = abs(ops[cur_id].dstEnd - (*ortho_algns)[i].x.lower);
								cur = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstEnd - (*ortho_algns)[i].x.lower), NO_GAP_INC);
								if( (*ortho_algns)[i].sign == 0 ) {
									old = (*ortho_algns)[i].y.lower;
									if( cur > old ) (*ortho_algns)[i].yl_diff = cur - old;
								}
								else if( (*ortho_algns)[i].sign == 1 ) {
									old = (*ortho_algns)[i].y.upper;
									if( cur < old ) (*ortho_algns)[i].yr_diff = old - cur;
								}
							}
						}
						else {
							cur_diff = (*ortho_algns)[i].xl_diff;
							if( cur_diff < abs(ops[cur_id].dstEnd - (*ortho_algns)[i].y.lower) ) {
								(*ortho_algns)[i].yl_diff = abs(ops[cur_id].dstEnd - (*ortho_algns)[i].y.lower);	
								if( (*ortho_algns)[i].init_sign == 0 ) {
									loc = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstEnd - (*ortho_algns)[i].y.lower), GAP_INC_IN_Y);
									cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
									old = (*ortho_algns)[i].x.lower;
									if( cur > old ) (*ortho_algns)[i].xl_diff = cur - old;
								}
								else if( (*ortho_algns)[i].init_sign == 1 ) {
									loc = find_yloc_one((*ortho_algns)[i], f, abs((*ortho_algns)[i].y.upper - ops[cur_id].dstEnd), GAP_INC_IN_Y);
									cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
									old = (*ortho_algns)[i].x.upper;
									if( cur < old ) (*ortho_algns)[i].xr_diff = old - cur;
								}
							}
						}
					}
				}

				if((ops[cur_id].dstStart >= cmp.lower) && (ops[cur_id].dstStart <= cmp.upper) ) 
				{
					cur_pid = cal_pid_part_algn((*ortho_algns), i, abs(ops[cur_id].dstStart - org_cmp.lower), width(org_cmp), f, sp_id);
					if( (int)(cur_pid+0.5) >= (int)(ops[cur_id].pid-0.5) ) {

					}
					else {
						*is_untagged = true;
						if( sp_id == SELF1 ) {
							cur_diff = (*ortho_algns)[i].xr_diff;
							if( cur_diff < ((*ortho_algns)[i].x.upper - ops[cur_id].dstStart) ) {
								(*ortho_algns)[i].xr_diff = (*ortho_algns)[i].x.upper - ops[cur_id].dstStart;
								cur = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstStart - (*ortho_algns)[i].x.lower), NO_GAP_INC);
								if( (*ortho_algns)[i].init_sign == 0 ) {
									old = (*ortho_algns)[i].y.upper;
									if( cur < old ) (*ortho_algns)[i].yr_diff = old - cur;
								}
								else if( (*ortho_algns)[i].init_sign == 1 ) {
									old = (*ortho_algns)[i].y.lower;
									if( cur > old ) (*ortho_algns)[i].yl_diff = cur - old;
								}
							}
						}
						else {
							cur_diff = (*ortho_algns)[i].yr_diff;
							if( cur_diff < ((*ortho_algns)[i].y.upper - ops[cur_id].dstStart) ) {
								(*ortho_algns)[i].yr_diff = (*ortho_algns)[i].y.upper - ops[cur_id].dstStart;
								if( (*ortho_algns)[i].init_sign == 0 ) {
									loc = find_yloc_one((*ortho_algns)[i], f, abs(ops[cur_id].dstStart - (*ortho_algns)[i].y.lower), GAP_INC_IN_Y);
									cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
									old = (*ortho_algns)[i].x.upper;
									if( cur < old ) (*ortho_algns)[i].xr_diff = old - cur;
								}
								else if( (*ortho_algns)[i].init_sign == 1 ) {
									loc = find_yloc_one((*ortho_algns)[i], f, abs((*ortho_algns)[i].y.upper - ops[cur_id].dstStart), GAP_INC_IN_Y);
									cur = find_xloc_one((*ortho_algns)[i], f, loc, GAP_INC);
									old = (*ortho_algns)[i].x.lower;
									if( cur > old ) (*ortho_algns)[i].xl_diff = cur - old;
								}
							}
						}
					}
				}

				(*init_algns)[(*ortho_algns)[i].index].xl_diff = (*ortho_algns)[i].xl_diff;
				(*init_algns)[(*ortho_algns)[i].index].xr_diff = (*ortho_algns)[i].xr_diff;
				(*init_algns)[(*ortho_algns)[i].index].yl_diff = (*ortho_algns)[i].yl_diff;
				(*init_algns)[(*ortho_algns)[i].index].yr_diff = (*ortho_algns)[i].yr_diff;
			}
		}
		j++;
	}

/*
	if( num_added > 0 ) {
		if( num_ortho_algns < 0 ) {
			fatal("alignments aren't read yet\n");
		}
		else {
//			*ortho_algns = ckrealloc(*ortho_algns, ((num_ortho_algns) + num_added) * sizeof(struct DotList));
		}

		if( num_algns < 0 ) {
			fatal("alignments aren't read yet\n");
		}
		else {
//			*init_algns = ckrealloc(*init_algns, (num_algns + num_added) * sizeof(struct DotList));
		}
	}
*/

//	*num_ortho = num_ortho_algns + num_added;
	return(num_added);
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

bool check_later_dup_existence(struct ops_list *ops, int id, int num_ops)
{
	bool res = false;
	int i = 0;
	struct I src, dst;
	struct I tmp_src, tmp_dst;
	int b = 0, e = 0;

	tmp_src = assign_I(0, 1);
	tmp_dst = assign_I(0, 1);
	src = assign_I(ops[id].srcStart, ops[id].srcEnd);
	dst = assign_I(ops[id].dstStart, ops[id].dstEnd);

	i = 0;
	while( (i < num_ops) && (res == false) ) {
		if( (ops[i].sign == '+') || ((ops[i].sign == 'c') && (i > id)) ) {
			tmp_src = assign_I(ops[i].srcStart, ops[i].srcEnd);
			tmp_dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( (subset(src, tmp_src) == true) && (subset(dst, tmp_dst) == true) )
			{
				b = tmp_dst.lower + (tmp_src.lower - src.lower);
				e = b + width(src);
				tmp_dst = assign_I(b, e);
				if( strict_almost_equal(dst, tmp_dst) == true ) {
					res = true;
				}
			}
			else if( (subset(dst, tmp_src) == true) && (subset(src, tmp_dst) == true) ) 
			{
				b = tmp_dst.lower + (tmp_src.lower - dst.lower);
				e = b + width(dst);
				tmp_dst = assign_I(b, e);
				if( strict_almost_equal(src, tmp_dst) == true ) {
					res = true;
				}
			}
		}
		else if((ops[i].sign == '-') || ((ops[i].sign == 'v') && (i > id))) {
			if( (subset(src, tmp_src) == true) && (subset(dst, tmp_dst) == true) )
			{
				e = tmp_dst.upper - (tmp_src.lower - src.lower);
				b = e - width(src);
				tmp_dst = assign_I(b, e);
				if( strict_almost_equal(dst, tmp_dst) == true ) {
					res = true;
				}
			}
			else if( (subset(dst, tmp_src) == true) && (subset(src, tmp_dst) == true) ) 
			{
				e = tmp_dst.upper - (tmp_src.lower - dst.lower);
				b = e - width(dst);
				tmp_dst = assign_I(b, e);
				if( strict_almost_equal(src, tmp_dst) == true ) {
					res = true;
				}
			}
		}
		i++;
	}

	return(res);
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

float cal_ortho_algn_pid(struct I reg, int sp_id, struct slist *sorted, struct DotList *init_algns, int num_algns, FILE *f)
{
	float res = (float)0;
	float cur_pid = (float)0; 
	int s_loc = 0, e_loc = 0;
	int i = 0, j = 0;
	struct I cmp, org_cmp; 
	int lo = 0, hi = 1;

	cmp = assign_I(0,1);
	org_cmp = assign_I(0,1);
  s_loc = search_range_b(sorted, init_algns, num_algns, reg.lower, sp_id);
  e_loc = search_range_e(sorted, init_algns, num_algns, reg.upper, sp_id);
	
	for( j = s_loc; j <= e_loc; j++ ) {
		i = sorted[j].id;
		cur_pid = (float) 0;
		if( sp_id == SELF1 ) {
      lo = init_algns[i].x.lower + init_algns[i].xl_diff;
      hi = init_algns[i].x.upper - init_algns[i].xr_diff;
      org_cmp = assign_I(init_algns[i].x.lower, init_algns[i].x.upper);
		}
		else if ( sp_id == SELF2 ) {
      lo = init_algns[i].y.lower + init_algns[i].yl_diff;
      hi = init_algns[i].y.upper - init_algns[i].yr_diff;
      org_cmp = assign_I(init_algns[i].y.lower, init_algns[i].y.upper);
		}
		else {
			fatalf("unexpected sp_id %d\n", sp_id);
		}

    if( hi > lo ) cmp = assign_I(lo, hi);

    if( hi <= lo ) {}
    else if( init_algns[i].sign == DELETED ) {}
    else if( cmp.upper < reg.lower ) {}
    else if( cmp.lower > reg.upper ) {}
    else {
      if( (reg.lower <= cmp.lower) && (reg.upper >= cmp.upper) ) {
				cur_pid = init_algns[i].identity;
      }
      else if( (cmp.lower <= reg.lower) && (cmp.upper >= reg.upper) )
      {
        cur_pid = cal_pid_part_algn(init_algns, i, abs(reg.lower - org_cmp.lower), abs(org_cmp.upper - reg.upper), f, sp_id);
			}
			else if( cmp.lower > reg.lower ) {
        cur_pid = cal_pid_part_algn(init_algns, i, abs(cmp.lower - org_cmp.lower), abs(org_cmp.upper - reg.upper), f, sp_id);
			}
			else if( cmp.upper < reg.upper ) {
        cur_pid = cal_pid_part_algn(init_algns, i, abs(reg.lower - org_cmp.lower), abs(org_cmp.upper - cmp.upper), f, sp_id);
			}
		}

		if( cur_pid > res ) {
			res = cur_pid;
		}
	}

	return(res);
}
