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
	int i, j, k, l;
	struct slist *sorted;
	struct I cur, cmp, ov, tmp;
	int b, e;
	int flag = TRUE;
	int cur_id, cmp_id, tmp_id;
	char S1[BIG], T1[BIG], S2[BIG], T2[BIG];
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
	sort_init_algns(sorted, init_algns, num_init_algns, SELF1);

	for( i = 0; i < num_init_algns; i++ ) {
		cur_id = sorted[i].id;
		b = init_algns[cur_id].x.lower + init_algns[cur_id].xl_diff;
		e = init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff;
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
				b = init_algns[cmp_id].x.lower + init_algns[cmp_id].xl_diff;
				e = init_algns[cmp_id].x.upper - init_algns[cmp_id].xr_diff;
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
										if( init_algns[cmp_id].xl_offset == 0 ) {
											old_pos = init_algns[cmp_id].y.lower;
											cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff, NO_GAP_INC);
										}
										else {
											old_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_offset, NO_GAP_INC);
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
											old_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_offset, NO_GAP_INC);
											cur_pos = find_yloc_one(init_algns[cmp_id], f, init_algns[cmp_id].xl_diff + init_algns[cmp_id].xl_offset, NO_GAP_INC);
										}

            				if( cur_pos < old_pos ) init_algns[cmp_id].yr_diff = old_pos - cur_pos;
									}

									cmp = assign_I(cmp.lower+width(ov), cmp.upper);
									l = j;
									if( (l+1) < num_init_algns ) {
										tmp_id = sorted[l+1].id;
										b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff;
										e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff; 
										tmp = assign_I(b, e);
									}

									while( ((l+1) < num_init_algns) && (cmp.lower >= tmp.lower) && (width(cmp) < width(tmp)) ) {
										h = assign_slist(sorted[l]);
										sorted[l] = assign_slist(sorted[l+1]);
										sorted[l+1] = assign_slist(h);
										l++;
										if ( (l+1) < num_init_algns ) {
											tmp_id = sorted[l+1].id;
											b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff;
											e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff;
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
										if( (init_algns[cur_id].xl_offset == 0) && (init_algns[cur_id].xr_offset == 0) ) {
    				        	old_pos = init_algns[cur_id].y.upper;
      				      	cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower, NO_GAP_INC);
										}
										else {
      				      	old_pos = find_yloc_one(init_algns[cur_id], f, init_algns[cur_id].x.upper - init_algns[cur_id].x.lower + init_algns[cur_id].xl_offset, NO_GAP_INC);
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
      				      	old_pos = find_yloc_one(init_algns[cur_id], f, init_algns[cur_id].x.upper - init_algns[cur_id].x.lower + init_algns[cur_id].xl_offset, NO_GAP_INC);
      				      	cur_pos = find_yloc_one(init_algns[cur_id], f, ov.lower - init_algns[cur_id].x.lower + init_algns[cur_id].xl_offset, NO_GAP_INC);
										}

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

void cut_part_algn(struct DotList *init_algns, int cur_id, int cmp_id, int mode, FILE *f)
{
	int k;
	struct I cur, cmp, ov;
	int b, e;
	int cut_len1 = 0, cut_len2 = 0;
	int end_pos1, end_pos2;
	float pid1, pid2;
	int beg1 = 0, beg2 = 0;
	int t_b1, t_b2;
	int old_pos, cur_pos;
	struct b_list *a1_info, *a2_info;
	int *t_b;
	char S1[BIG], T1[BIG], S2[BIG], T2[BIG];
	int len = 0;

	a1_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	a2_info = (struct b_list *) ckalloc(sizeof(struct b_list));

	t_b = (int *) ckalloc(sizeof(int));

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
