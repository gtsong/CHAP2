#include "main.h"
#include "map_algns.h"
#include "regions.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "read_algn.h"
#include "read_maf.h"

extern int debug_mode;

void map_one_to_one(int num_init_algns, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0, k = 0, l = 0;
	struct slist *sorted = NULL;
	struct I cur = {0, 1}, cmp = {0, 1}, ov = {0, 1}, tmp = {0, 1};
	int b = -1, e= -1;
	int flag = TRUE;
	int cur_id = 0, cmp_id = 0, tmp_id = 0;
	struct b_list *a1_info = NULL, *a2_info = NULL;
	int cut_len1 = 0, cut_len2 = 0;
	int end_pos1 = 0, end_pos2 = 0;
	float pid1 = (float) 0, pid2 = (float) 0;
	int beg1 = 0, beg2 = 0;
	int t_b1 = 0, t_b2 = 0;
	int *t_b = NULL;
	int old_pos = 0, cur_pos = 0;
	struct slist h = {0, 0, 0, 0, 0, true};
	char S[BIG] = "", S1[BIG] = "", T[BIG] = "", T1[BIG] = "";

	a1_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	a2_info = (struct b_list *) ckalloc(sizeof(struct b_list));

	if( num_init_algns > 0 ) {
		sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_init_algns);
	}
	t_b = (int *) ckalloc(sizeof(int));

	a1_info->b1 = 0;
	a1_info->e1 = 1;
	a1_info->len1 = 0;
	a1_info->b2 = 0;
	a1_info->e2 = 1;
	a1_info->len2 = 0;
	a1_info->strand = '+';
	a1_info->pid = 0;
	a2_info->b1 = 0;
	a2_info->e1 = 1;
	a2_info->len1 = 0;
	a2_info->b2 = 0;
	a2_info->e2 = 1;
	a2_info->len2 = 0;
	a2_info->strand = '+';
	a2_info->pid = 0;
	*t_b = 0;

	if( num_init_algns > 0 ) {
		for( i = 0; i < num_init_algns; i++ ) {
			sorted[i].id = i;
		}
		sort_init_algns(sorted, init_algns, num_init_algns, SELF1);
	}
	cur = assign_I(0,1);
	cmp = assign_I(0,1);
	ov = assign_I(0,1);

	for( i = 0; i < num_init_algns; i++ ) {
		cur_id = sorted[i].id;
		b = init_algns[cur_id].x.lower + init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;
		e = init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff - init_algns[cur_id].xr_offset;
		if( (init_algns[cur_id].sign == DELETED) || (init_algns[cur_id].sp_id != PAIR) ) {}
		else if( e <= b ) init_algns[cur_id].sign = DELETED;
		else {
			cur = assign_I(b,e);
			beg1 = init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;
			get_nth_algn(S, T, init_algns[cur_id].fid, beg1, f, a1_info, REG);
			if( S[strlen(S)-1] == '\n' ) end_pos1 = strlen(S)-1;
			else end_pos1 = strlen(S);
			k = 0;
			while( (end_pos1 >= 1) && (k < (init_algns[cur_id].xr_diff + init_algns[cur_id].xr_offset)) ) {
				if( S[end_pos1-1] != '-' ) k++;
				end_pos1--;
			}
			S[end_pos1] = '\0';
			flag = TRUE;
			j = i+1;
			while( (flag == TRUE) && (j < num_init_algns) && (init_algns[cur_id].ctg_id1 == init_algns[sorted[j].id].ctg_id1) ) {
				cmp_id = sorted[j].id;
				b = init_algns[cmp_id].x.lower + init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset;
				e = init_algns[cmp_id].x.upper - init_algns[cmp_id].xr_diff - init_algns[cmp_id].xr_offset;
				if( (init_algns[cmp_id].sign == DELETED) || (init_algns[cmp_id].sp_id != PAIR) ) {}
				else if( e <= b ) init_algns[cmp_id].sign = DELETED;
				else {
					cmp = assign_I(b,e);
					if( proper_overlap(cur, cmp) == true ) {
						ov = intersect(cur, cmp);
						beg2 = init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset;
						get_nth_algn(S1, T1, init_algns[cmp_id].fid, beg2, f, a2_info, REG);
						if( S1[strlen(S1)-1] == '\n' ) end_pos2 = strlen(S1)-1;
						else end_pos2 = strlen(S1);
						k = 0;
						while( (end_pos2 >= 1) && (k < (init_algns[cmp_id].xr_diff + init_algns[cmp_id].xr_offset)) ) {
							if( S1[end_pos2-1] != '-' ) k++;
							end_pos2--;
						}
						S1[end_pos2] = '\0';

						cut_len1 = count_ncol(cur, ov, S, end_pos1, t_b);
						t_b1 = *t_b;
						pid1 = cal_pid_maf_beg(S, T, t_b1, cut_len1);
						cut_len2 = count_ncol(cmp, ov, S1, end_pos2, t_b);
						t_b2 = *t_b;
						pid2 = cal_pid_maf_beg(S1, T1, t_b2, cut_len2);

						if( pid1 >= pid2 ) {
							if( (strict_almost_equal(cur, cmp) == true) || (f_loose_subset(cmp, cur, STRICT) == true) ) init_algns[cmp_id].sign = DELETED;						
							else if(subset(cur, cmp) == true) {
								init_algns[cur_id].sign = DELETED;
/* Need to revise */
							}
							else {
								if( cmp.upper > ov.upper ) {
									init_algns[cmp_id].xl_diff = init_algns[cmp_id].xl_diff + width(ov);	
									if( init_algns[cmp_id].sign == 0 ) {
										old_pos = init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_offset;
										cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset, NO_GAP_INC);
										if( cur_pos > old_pos ) init_algns[cmp_id].yl_diff = cur_pos - old_pos;

									}
									else if( init_algns[cmp_id].sign == 1 ) {
				            old_pos = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_offset;
       					    cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff+init_algns[cmp_id].xl_offset, NO_GAP_INC);
            				if( cur_pos < old_pos ) init_algns[cmp_id].yr_diff = old_pos - cur_pos;
									}

									cmp = assign_I(cmp.lower+width(ov), cmp.upper);
									l = j;
									if( (l+1) < num_init_algns ) {
										tmp_id = sorted[l+1].id;
										b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff + init_algns[tmp_id].xl_offset;
										e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff - init_algns[tmp_id].xr_offset;
										tmp = assign_I(b, e);
									}

									while( ((l+1) < num_init_algns) && (cmp.lower >= tmp.lower) && (width(cmp) < width(tmp)) ) {
										h = assign_slist(sorted[l]);
										sorted[l] = assign_slist(sorted[l+1]);
										sorted[l+1] = assign_slist(h);
										l++;
										if ( (l+1) < num_init_algns ) {
											tmp_id = sorted[l+1].id;
											b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff + init_algns[tmp_id].xl_offset;
											e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff - init_algns[tmp_id].xr_offset;
											tmp = assign_I(b, e);
										}
									}

									if( l != j ) j--;
								}
							}
						}
						else if( pid1 < pid2 ) {
							if( (strict_almost_equal(cur, cmp) == true) || (f_loose_subset(cur, cmp, STRICT) == true) ) {
								init_algns[cur_id].sign = DELETED;						
								flag = FALSE;
							}
							else if(subset(cmp, cur) == true) {
/* Need to revise */
								init_algns[cmp_id].sign = DELETED;						
							}
							else {
								if( cur.lower < ov.lower ) {
									init_algns[cur_id].xr_diff = init_algns[cur_id].xr_diff + width(ov);	
				          if( init_algns[cur_id].sign == 0 ) {
    				        old_pos = init_algns[cur_id].y.upper - init_algns[cur_id].yr_offset;
      				      cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);

        				    if( cur_pos < old_pos ) init_algns[cur_id].yr_diff = old_pos - cur_pos;
            			}
          				else if( init_algns[cur_id].sign == 1 ) {
            				old_pos = init_algns[cur_id].y.lower + init_algns[cur_id].yl_offset;
            				cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);
            				if( cur_pos > old_pos ) init_algns[cur_id].yl_diff = cur_pos - old_pos;
          				}
								}
							}
						}
					}
					else {
						flag = FALSE;
					}
				}
				j++;
			}
		}
	}

	free(t_b);
	free(a1_info);
	free(a2_info);
	if( num_init_algns > 0 ) {
		free(sorted);
	}
}

void map_one_to_one_for_sec(int num_init_algns, struct DotList *init_algns, FILE *f)
{
	int i, j, k, l;
	struct slist *sorted;
	struct I cur, cmp, ov, tmp;
	int b, e;
	int flag = TRUE;
	int cur_id, cmp_id, tmp_id;
	char S[BIG], S1[BIG], T[BIG], T1[BIG];
	struct b_list *a1_info, *a2_info;
	int cut_len1 = 0, cut_len2 = 0;
	int end_pos1, end_pos2;
	float pid1, pid2;
	int beg1 = 0, beg2 = 0;
	int t_b1, t_b2;
	int *t_b;
	int old_pos, cur_pos;
	struct slist h;

	a1_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	a2_info = (struct b_list *) ckalloc(sizeof(struct b_list));

	sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_init_algns);
	t_b = (int *) ckalloc(sizeof(int));

	a1_info->b1 = 0;
	a1_info->e1 = 1;
	a1_info->len1 = 0;
	a1_info->b2 = 0;
	a1_info->e2 = 1;
	a1_info->len2 = 0;
	a1_info->strand = '+';
	a1_info->pid = 0;
	a2_info->b1 = 0;
	a2_info->e1 = 1;
	a2_info->len1 = 0;
	a2_info->b2 = 0;
	a2_info->e2 = 1;
	a2_info->len2 = 0;
	a2_info->strand = '+';
	a2_info->pid = 0;
	*t_b = 0;

	sort_init_algns(sorted, init_algns, num_init_algns, SELF1);

	for( i = 0; i < num_init_algns; i++ ) {
		cur_id = sorted[i].id;
		b = init_algns[cur_id].x.lower + init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;
		e = init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff - init_algns[cur_id].xr_offset;
		if( (init_algns[cur_id].sign == DELETED) || (init_algns[cur_id].sp_id != PAIR) ) {}
		else if( e <= b ) init_algns[cur_id].sign = DELETED;
		else {
			cur = assign_I(b,e);
			beg1 = init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;
			get_nth_algn(S, T, init_algns[cur_id].fid, beg1, f, a1_info, REG);
			if( S[strlen(S)-1] == '\n' ) end_pos1 = strlen(S)-1;
			else end_pos1 = strlen(S);
			k = 0;
			while( (end_pos1 >= 1) && (k < (init_algns[cur_id].xr_diff + init_algns[cur_id].xr_offset)) ) {
				if( S[end_pos1-1] != '-' ) k++;
				end_pos1--;
			}
			S[end_pos1] = '\0';
			flag = TRUE;
			j = i+1;
			cmp_id = sorted[j].id;
			while( (flag == TRUE) && (j < num_init_algns) ) {
				b = init_algns[cmp_id].x.lower + init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset;
				e = init_algns[cmp_id].x.upper - init_algns[cmp_id].xr_diff - init_algns[cmp_id].xr_offset;
				if( (init_algns[cmp_id].sign == DELETED) || (init_algns[cmp_id].sp_id != PAIR) ) {}
				else if( e <= b ) init_algns[cmp_id].sign = DELETED;
				else {
					cmp = assign_I(b,e);
					if( proper_overlap(cur, cmp) == true ) {
						ov = intersect(cur, cmp);
						beg2 = init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset;
						get_nth_algn(S1, T1, init_algns[cmp_id].fid, beg2, f, a2_info, REG);
						if( S1[strlen(S1)-1] == '\n' ) end_pos2 = strlen(S1)-1;
						else end_pos2 = strlen(S1);
						k = 0;
						while( (end_pos2 >= 1) && (k < (init_algns[cmp_id].xr_diff + init_algns[cmp_id].xr_offset)) ) {
							if( S1[end_pos2-1] != '-' ) k++;
							end_pos2--;
						}
						S1[end_pos2] = '\0';

						cut_len1 = count_ncol(cur, ov, S, end_pos1, t_b);
						t_b1 = *t_b;
						pid1 = cal_pid_maf_beg(S, T, t_b1, cut_len1);
						cut_len2 = count_ncol(cmp, ov, S1, end_pos2, t_b);
						t_b2 = *t_b;
						pid2 = cal_pid_maf_beg(S1, T1, t_b2, cut_len2);

						if( pid1 >= pid2 ) {
							if( (strict_almost_equal(cur, cmp) == true) || (f_loose_subset(cmp, cur, STRICT) == true) ) init_algns[cmp_id].sign = DELETED;						
							else if(subset(cur, cmp) == true) {
								init_algns[cur_id].sign = DELETED;
/* Need to revise */
							}
							else {
								if( cmp.upper > ov.upper ) {
									init_algns[cmp_id].xl_diff = init_algns[cmp_id].xl_diff + width(ov);	
									if( init_algns[cmp_id].sign == 0 ) {
										old_pos = init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_offset;
										cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset, NO_GAP_INC);
										if( cur_pos > old_pos ) init_algns[cmp_id].yl_diff = cur_pos - old_pos;

									}
									else if( init_algns[cmp_id].sign == 1 ) {
				            old_pos = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_offset;
       					    cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff+init_algns[cmp_id].xl_offset, NO_GAP_INC);
            				if( cur_pos < old_pos ) init_algns[cmp_id].yr_diff = old_pos - cur_pos;
									}

									cmp = assign_I(cmp.lower+width(ov), cmp.upper);
									l = j;
									if( (l+1) < num_init_algns ) {
										tmp_id = sorted[l+1].id;
										b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff + init_algns[tmp_id].xl_offset;
										e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff - init_algns[tmp_id].xr_offset;
										tmp = assign_I(b, e);
									}

									while( ((l+1) < num_init_algns) && (cmp.lower >= tmp.lower) && (width(cmp) < width(tmp)) ) {
										h = assign_slist(sorted[l]);
										sorted[l] = assign_slist(sorted[l+1]);
										sorted[l+1] = assign_slist(h);
										l++;
										if ( (l+1) < num_init_algns ) {
											tmp_id = sorted[l+1].id;
											b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff + init_algns[tmp_id].xl_offset;
											e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff - init_algns[tmp_id].xr_offset;
											tmp = assign_I(b, e);
										}
									}

									if( l != j ) j--;
								}
							}
						}
						else if( pid1 < pid2 ) {
							if( (strict_almost_equal(cur, cmp) == true) || (f_loose_subset(cur, cmp, STRICT) == true) ) {
								init_algns[cur_id].sign = DELETED;						
								flag = FALSE;
							}
							else if(subset(cmp, cur) == true) {
/* Need to revise */
								init_algns[cmp_id].sign = DELETED;						
							}
							else {
								if( cur.lower < ov.lower ) {
									init_algns[cur_id].xr_diff = init_algns[cur_id].xr_diff + width(ov);	
				          if( init_algns[cur_id].sign == 0 ) {
    				        old_pos = init_algns[cur_id].y.upper - init_algns[cur_id].yr_offset;
      				      cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);

        				    if( cur_pos < old_pos ) init_algns[cur_id].yr_diff = old_pos - cur_pos;
            			}
          				else if( init_algns[cur_id].sign == 1 ) {
            				old_pos = init_algns[cur_id].y.lower + init_algns[cur_id].yl_offset;
            				cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);
            				if( cur_pos > old_pos ) init_algns[cur_id].yl_diff = cur_pos - old_pos;
          				}
								}
							}
						}
					}
					else {
						flag = FALSE;
					}
				}
				j++;
				if( j < num_init_algns) cmp_id = sorted[j].id;
			}
		}
	}

	free(t_b);
	free(a1_info);
	free(a2_info);
	free(sorted);
}
