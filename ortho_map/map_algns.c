#include "main.h"
#include "map_algns.h"
#include "regions.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "read_algn.h"
#include "read_maf.h"

#define LONG_ALGN_TH 20000

extern int debug_mode;
extern char S1[BIG], T1[BIG];

void map_one_to_one(int num_init_algns, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0, k = 0, l = 0;
	struct slist *sorted = NULL;
	struct I cur, cmp, ov, tmp;
	struct I cur_y, cmp_y;
	int b = -1, e = -1;
	int y_b = -1, y_e = -1;
	int flag = TRUE;
	int cur_id = -1, cmp_id = -1, tmp_id = -1;
	char S2[BIG], T2[BIG];
	struct b_list *a1_info, *a2_info;
	int cut_len1 = 0, cut_len2 = 0;
	int end_pos1 = 0, end_pos2 = 0;
	float pid1, pid2;
	int beg1 = 0, beg2 = 0;
	int t_b1 = 0, t_b2 = 0;
	int *t_b;
	struct slist h;

	a1_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	a2_info = (struct b_list *) ckalloc(sizeof(struct b_list));

	if( num_init_algns > 0 ) {
		sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_init_algns);
		initialize_slist(sorted, 0, num_init_algns);
		sort_init_algns(sorted, init_algns, num_init_algns, INIT_WIDTH);
	}

	for( i = 0; i < num_init_algns; i++ ) {
		cur_id = sorted[i].id;
		b = init_algns[cur_id].x.lower + init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;
		e = init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff - init_algns[cur_id].xr_offset;
		y_b = init_algns[cur_id].y.lower + init_algns[cur_id].yl_diff + init_algns[cur_id].yl_offset;
		y_e = init_algns[cur_id].y.upper - init_algns[cur_id].yr_diff - init_algns[cur_id].yr_offset;
		if( (init_algns[cur_id].sign == DELETED) || (init_algns[cur_id].sp_id != PAIR) ) {}
		else if( (e <= b) || (y_e <= y_b) ) init_algns[cur_id].sign = DELETED;
		else {
			cur = assign_I(b,e);
			cur_y = assign_I(y_b,y_e);
			for( j = (i+1); j < num_init_algns; j++ ) {
				cmp_id = sorted[j].id;
				b = init_algns[cmp_id].x.lower + init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset;
				e = init_algns[cmp_id].x.upper - init_algns[cmp_id].xr_diff - init_algns[cmp_id].xr_offset;
				y_b = init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_diff + init_algns[cmp_id].yl_offset;
				y_e = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_diff - init_algns[cmp_id].yr_offset;
				if( (init_algns[cmp_id].sign == DELETED) || (init_algns[cmp_id].sp_id != PAIR) ) {}
				else if( (e <= b) || (y_e <= y_b) ) init_algns[cmp_id].sign = DELETED;
				else {
					cmp = assign_I(b,e);
					cmp_y = assign_I(y_b,y_e);
					if( ((init_algns[cur_id].identity - init_algns[cmp_id].identity) >= -1) && ((f_loose_subset(cmp, cur, LOOSE) == true) || (f_loose_subset(cmp_y, cur_y, LOOSE) == true)) && (width(cur) >= LONG_ALGN_TH) && (width(cmp) <= (int)(0.70 * (float)(width(cur)))) ) {
						init_algns[cmp_id].sign = DELETED;
					}
				}
			}
		}	
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

	cur = assign_I(0,1);
	cmp = assign_I(0,1);
	ov = assign_I(0,1);
	tmp = assign_I(0,1);
	pid1 = (float)0;
	pid2 = (float)0;
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
			get_nth_algn(S1, T1, init_algns[cur_id].fid, beg1, f, a1_info, REG);
			if( S1[strlen(S1)-1] == '\n' ) end_pos1 = strlen(S1)-1;
			else end_pos1 = strlen(S1);
			k = 0;
			while( (end_pos1 >= 1) && (k < (init_algns[cur_id].xr_diff + abs(init_algns[cur_id].xr_offset))) ) {
				if( S1[end_pos1-1] != '-' ) k++;
				end_pos1--;
			}
			S1[end_pos1] = '\0';
			flag = TRUE;
			if( (i+1) < num_init_algns ) {
				j = i+1;
				cmp_id = sorted[j].id;
			}
			else {
				j = num_init_algns;
			}

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
						get_nth_algn(S2, T2, init_algns[cmp_id].fid, beg2, f, a2_info, REG);
						if( S2[strlen(S2)-1] == '\n' ) end_pos2 = strlen(S2)-1;
						else end_pos2 = strlen(S2);
						k = 0;
						while( (end_pos2 >= 1) && (k < (init_algns[cmp_id].xr_diff + abs(init_algns[cmp_id].xr_offset))) ) {
							if( S2[end_pos2-1] != '-' ) k++;
							end_pos2--;
						}
						S2[end_pos2] = '\0';

						cut_len1 = count_ncol(cur, ov, S1, end_pos1, t_b);
						t_b1 = *t_b;
						pid1 = cal_pid_maf_beg(S1, T1, t_b1, cut_len1);
						cut_len2 = count_ncol(cmp, ov, S2, end_pos2, t_b);
						t_b2 = *t_b;
						pid2 = cal_pid_maf_beg(S2, T2, t_b2, cut_len2);

						if( ((init_algns[cur_id].identity - init_algns[cmp_id].identity) >= -2) && (f_loose_subset(cmp, cur, LOOSE) == true) && (width(cur) >= LONG_ALGN_TH) && (width(cmp) <= (int)(0.85 * (float)(width(cur)))) ) {
							init_algns[cmp_id].sign = DELETED;
						}
						else if( ((init_algns[cmp_id].identity - init_algns[cur_id].identity) >= -2) && (f_loose_subset(cur, cmp, LOOSE) == true) && (width(cmp) >= LONG_ALGN_TH) && (width(cur) <= (int)(0.85 * (float)(width(cmp)))) ) {
							init_algns[cur_id].sign = DELETED;
						}
						else if( pid1 >= pid2 ) {
/*
							if( ((init_algns[cur_id].identity - init_algns[cmp_id].identity) >= -1) && (width(cur) >= LONG_ALGN_TH) && (width(cmp) <= (int)(0.20 * (float)(width(cur)))) ) 
							{
								if( debug_mode == TRUE ) {
									printf("%d-%d,%d-%d remains\n", init_algns[cur_id].x.lower, init_algns[cur_id].x.upper, init_algns[cur_id].y.lower, init_algns[cur_id].y.upper);
								}
							}
*/

							if( (strict_almost_equal(cur, cmp) == true) || (f_loose_subset(cmp, cur, STRICT) == true) ) init_algns[cmp_id].sign = DELETED;						
							else if(subset(cur, cmp) == true) {
								init_algns[cur_id].sign = DELETED;
							}
							else {
								if( cmp.upper > ov.upper ) {
									adjust_algn_left_diff(init_algns, cmp_id, ov, f);
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
									adjust_algn_right_diff(init_algns, cur_id, ov, f);
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
	if( sorted != NULL ) {
		free(sorted);
	}
}

void adjust_algn_left_diff(struct DotList *init_algns, int cmp_id, struct I ov, FILE *f)
{
	int old_pos = 0, cur_pos = 0;

	init_algns[cmp_id].xl_diff = init_algns[cmp_id].xl_diff + width(ov);
	if( init_algns[cmp_id].sign == 0 ) {
		if( init_algns[cmp_id].xl_offset == 0 ) {
			old_pos = init_algns[cmp_id].y.lower;
			cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff, NO_GAP_INC);
		}
		else {
			old_pos = init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_offset;
			cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset, NO_GAP_INC);
		}

		if( cur_pos > old_pos ) init_algns[cmp_id].yl_diff = cur_pos - old_pos;
	}
	else if( init_algns[cmp_id].sign == 1 ) {
		if( init_algns[cmp_id].xl_offset == 0 ) {
			old_pos = init_algns[cmp_id].y.upper;
			cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff, NO_GAP_INC);
		}
		else {
			old_pos = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_offset;
			cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset, NO_GAP_INC);
		}

		if( cur_pos < old_pos ) init_algns[cmp_id].yr_diff = old_pos - cur_pos;
	}
}

void adjust_algn_right_diff(struct DotList *init_algns, int cur_id, struct I ov, FILE *f)
{
	int old_pos = 0, cur_pos = 0;

	init_algns[cur_id].xr_diff = init_algns[cur_id].xr_diff + width(ov);
	if( init_algns[cur_id].sign == 0 ) {
		if( (init_algns[cur_id].xl_offset == 0) && (init_algns[cur_id].xr_offset == 0) ) {
			old_pos = init_algns[cur_id].y.upper;
			cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);
		}
		else {
			old_pos = init_algns[cur_id].y.upper - init_algns[cur_id].yr_offset;
			cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower + init_algns[cur_id].xl_offset, NO_GAP_INC);
		}

		if( cur_pos < old_pos ) init_algns[cur_id].yr_diff = old_pos - cur_pos;
	}
	else if( init_algns[cur_id].sign == 1 ) {
		if( (init_algns[cur_id].xl_offset == 0) && (init_algns[cur_id].xr_offset == 0) ) {
			old_pos = init_algns[cur_id].y.lower;
			cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);
		}
		else {
			old_pos = init_algns[cur_id].y.lower + init_algns[cur_id].yl_offset;
			cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower + init_algns[cur_id].xl_offset, NO_GAP_INC);
		}

		if( cur_pos > old_pos ) init_algns[cur_id].yl_diff = cur_pos - old_pos;
	}
}

void cut_part_algn(struct DotList *init_algns, int cur_id, int cmp_id, int mode, FILE *f)
{
	int k = 0;
	struct I cur, cmp, ov;
	int b = 0, e = 0;
	int cut_len1 = 0, cut_len2 = 0;
	int end_pos1 = 0, end_pos2 = 0;
	float pid1, pid2;
	int beg1 = 0, beg2 = 0;
	int t_b1 = 0, t_b2 = 0;
	int old_pos = 0, cur_pos = 0;
	struct b_list *a1_info, *a2_info;
	int *t_b;
	char S2[BIG], T2[BIG];
	int len = 0;

	a1_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	a2_info = (struct b_list *) ckalloc(sizeof(struct b_list));
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
  cur = assign_I(0,1);
  cmp = assign_I(0,1);
  ov = assign_I(0,1);
  pid1 = (float)0;
  pid2 = (float)0;

	t_b = (int *) ckalloc(sizeof(int));
	*t_b = 0;

	b = init_algns[cur_id].y.lower + init_algns[cur_id].yl_diff;
	e = init_algns[cur_id].y.upper - init_algns[cur_id].yr_diff;
	if( (init_algns[cur_id].sign == DELETED) || (init_algns[cur_id].sp_id != PAIR) ) {}
	else if( e <= b ) init_algns[cur_id].sign = DELETED;
	else {
		cur = assign_I(b,e);
		beg1 = init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;
		get_nth_algn(S1, T1, init_algns[cur_id].fid, beg1, f, a1_info, REG);
		if( S1[strlen(S1)-1] == '\n' ) end_pos1 = strlen(S1)-1;
		else end_pos1 = strlen(S1);
		k = 0;
		while( (end_pos1 >= 1) && (k < (init_algns[cur_id].yr_diff + abs(init_algns[cur_id].yr_offset))) ) {
			if( S1[end_pos1-1] != '-' ) k++;
			end_pos1--;
		}
		S1[end_pos1] = '\0';

		b = init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_diff;
		e = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_diff;
		if( (init_algns[cmp_id].sign == DELETED) || (init_algns[cmp_id].sp_id != PAIR) ) {}
		else if( e <= b ) init_algns[cmp_id].sign = DELETED;
		else {
			cmp = assign_I(b,e);
			if( (subset(cur, cmp) == true) || (subset(cmp, cur) == true) ) {}
			else if( proper_overlap(cur, cmp) == true ) {
				ov = intersect(cur, cmp);
				beg2 = init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset;
				get_nth_algn(S2, T2, init_algns[cmp_id].fid, beg2, f, a2_info, REG);
				if( S2[strlen(S2)-1] == '\n' ) end_pos2 = strlen(S2)-1;
				else end_pos2 = strlen(S2);
				k = 0;
				while( (end_pos2 >= 1) && (k < (init_algns[cmp_id].yr_diff + abs(init_algns[cmp_id].yr_offset))) ) {
					if( S2[end_pos2-1] != '-' ) k++;
					end_pos2--;
				}
				S2[end_pos2] = '\0';

				cut_len1 = count_ncol(cur, ov, S1, end_pos1, t_b);
				t_b1 = *t_b;
				pid1 = cal_pid_maf_beg(S1, T1, t_b1, cut_len1);
				cut_len2 = count_ncol(cmp, ov, S2, end_pos2, t_b);
				t_b2 = *t_b;
				pid2 = cal_pid_maf_beg(S2, T2, t_b2, cut_len2);

				if( (mode == LEFT_SIDE) || ((mode == TIE) && (pid1 >= pid2)) ) {
					if( (strict_almost_equal(cur, cmp) == true) || (f_loose_subset(cmp, cur, STRICT) == true) ) init_algns[cmp_id].sign = DELETED;						
					else if(subset(cur, cmp) == true) {
						init_algns[cur_id].sign = DELETED;
					}
					else {
						if( cmp.upper > ov.upper ) {
							init_algns[cmp_id].yl_diff = init_algns[cmp_id].yl_diff + width(ov);	
							if( init_algns[cmp_id].sign == 0 ) {
								if( init_algns[cmp_id].yl_offset == 0 ) {
									old_pos = init_algns[cmp_id].x.lower;
									cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].yl_diff, GAP_INC_IN_Y);
									cur_pos = find_xloc_one(init_algns[cmp_id], f, cur_pos, GAP_INC);
								}
								else {
									old_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].yl_offset, GAP_INC_IN_Y);
									old_pos = find_xloc_one(init_algns[cmp_id], f, old_pos, GAP_INC);
									cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].yl_diff + init_algns[cmp_id].yl_offset, GAP_INC_IN_Y);
									cur_pos = find_xloc_one(init_algns[cmp_id], f, cur_pos, GAP_INC);
								}

								if( cur_pos > old_pos ) init_algns[cmp_id].xl_diff = cur_pos - old_pos;
							}
							else if( init_algns[cmp_id].sign == 1 ) {
								if( (init_algns[cmp_id].yl_offset == 0) && (init_algns[cmp_id].yr_offset) ) {
				        	old_pos = init_algns[cmp_id].x.upper;
									len = width(init_algns[cmp_id].y);
       			    	cur_pos = find_yloc_one(init_algns[cmp_id], f, len - init_algns[cmp_id].yl_diff, GAP_INC_IN_Y);
									cur_pos = find_xloc_one(init_algns[cmp_id], f, cur_pos, GAP_INC);
								}
								else {
									len = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_offset - init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_offset;
       			    	old_pos = find_yloc_one(init_algns[cmp_id], f, len - init_algns[cmp_id].yl_offset, GAP_INC_IN_Y);
									old_pos = find_xloc_one(init_algns[cmp_id], f, old_pos, GAP_INC);
       			    	cur_pos = find_yloc_one(init_algns[cmp_id], f, len - init_algns[cmp_id].yl_diff - init_algns[cmp_id].yl_offset, GAP_INC_IN_Y);
									cur_pos = find_xloc_one(init_algns[cmp_id], f, cur_pos, GAP_INC);
								}

            		if( cur_pos < old_pos ) init_algns[cmp_id].xr_diff = old_pos - cur_pos;
							}
						}
						else if( (mode == RIGHT_SIDE) || ( (mode == TIE) && (pid1 < pid2) ) ) 
						{
							if( cur.lower < ov.lower ) {
								init_algns[cur_id].yr_diff = init_algns[cur_id].yr_diff + width(ov);	
				        if( init_algns[cur_id].sign == 0 ) {
									if( (init_algns[cmp_id].yl_offset == 0) && (init_algns[cmp_id].yr_offset == 0) ) {
										len = width(init_algns[cmp_id].y);
				        		old_pos = init_algns[cur_id].x.upper;
       			    		cur_pos = find_yloc_one(init_algns[cur_id], f, len-init_algns[cur_id].yr_diff, GAP_INC_IN_Y);
										cur_pos = find_xloc_one(init_algns[cur_id], f, cur_pos, GAP_INC);
									}
									else {
										len = init_algns[cmp_id].y.upper - init_algns[cmp_id].yr_offset - init_algns[cmp_id].y.lower + init_algns[cmp_id].yl_offset;
       			    		old_pos = find_yloc_one(init_algns[cur_id], f, len-abs(init_algns[cur_id].yr_offset), GAP_INC_IN_Y);
										old_pos = find_xloc_one(init_algns[cur_id], f, old_pos, GAP_INC);
       			    		cur_pos = find_yloc_one(init_algns[cur_id], f, len-(init_algns[cur_id].yr_diff+abs(init_algns[cur_id].yr_offset)), GAP_INC_IN_Y);
										cur_pos = find_xloc_one(init_algns[cur_id], f, cur_pos, GAP_INC);
									}

            			if( cur_pos < old_pos ) init_algns[cur_id].xr_diff = old_pos - cur_pos;
								}
          			else if( init_algns[cur_id].sign == 1 ) {
									if( (init_algns[cmp_id].yl_offset == 0) && (init_algns[cmp_id].yr_offset == 0) ) {
										old_pos = init_algns[cur_id].x.lower;
										cur_pos = find_yloc_one(init_algns[cur_id], f, init_algns[cur_id].yr_diff, GAP_INC_IN_Y);
										cur_pos = find_xloc_one(init_algns[cur_id], f, cur_pos, GAP_INC);
									}
									else {
										old_pos = find_yloc_one(init_algns[cur_id], f, abs(init_algns[cur_id].yr_offset), GAP_INC_IN_Y);
										old_pos = find_xloc_one(init_algns[cur_id], f, old_pos, GAP_INC);
										cur_pos = find_yloc_one(init_algns[cur_id], f, init_algns[cur_id].yr_diff + abs(init_algns[cur_id].yr_offset), GAP_INC_IN_Y);
										cur_pos = find_xloc_one(init_algns[cur_id], f, cur_pos, GAP_INC);
									}

									if( cur_pos > old_pos ) init_algns[cur_id].xl_diff = cur_pos - old_pos;
								}
							}
						}
					}
				}
			}
		}
	}

	free(t_b);
	free(a1_info);
	free(a2_info);
}
