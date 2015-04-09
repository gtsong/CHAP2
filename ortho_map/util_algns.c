#include "main.h"
#include "regions.h"
#include "read_algn.h"
#include "read_maf.h"
#include "util_algns.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"

extern int debug_mode;
extern char S1[BIG], T1[BIG];
extern char S[BIG], T[BIG];

int search_candi_algns(struct DotList *algns, int num_algns, struct I src, struct I dst, struct DotList *candi_algns )
{
	struct DotList *temp_algns; // the sequence starts at 1 and each range of a segment is a form of [a, b), i.e. the nucleotide at 'a' is included and one at 'b' is not
	struct slist *st;
	int i = 0, j = 0;
	int b = 0, e = 0;
	int count = 0;

	temp_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);
	st = (struct slist *) ckalloc(sizeof(struct slist) * num_algns);

	initialize_slist(st, 0, num_algns);
	sort_init_algns(st, algns, num_algns, INIT_SELF1);
	b = search_range_b(st, algns, num_algns, src.lower, SELF1);
	e = search_range_e(st, algns, num_algns, src.upper, SELF1);

	j = 0;
	for( i = b; i <= e; i++ ) {
		assign_algn(temp_algns, j, algns[st[i].id]);	
		j++;
	}
	count = j;
	initialize_slist(st, 0, count);
	sort_init_algns(st, temp_algns, count, INIT_SELF2);
	b = search_range_b(st, temp_algns, count, dst.lower, SELF2);
	e = search_range_e(st, temp_algns, count, dst.upper, SELF2);

	j = 0;
	for( i = b; i <= e; i++ ) {
		assign_algn(candi_algns, j, temp_algns[st[i].id]);	
		j++;
	}
	count = j;

	if( count == 0 ) {
		if(debug_mode == TRUE) {
			printf("candi alignments not found\n");
		}
	}

	free(temp_algns);
	free(st);
	return(count);
}

float pick_sim_level(struct I src, struct I dst, struct DotList *algns, int num_algns, FILE *f)
{
	int i = 0;
	int cur_id = -1;
	struct slist *st;
	int max_val = -1;
	int cur_val = 0;
	int max_id = -1;
	int sign = -1;
	int lo = 0, hi = 1;
	int beg = 0, end = 1, len = 0;
	struct b_list *a_info;
	float res = (float)(-1);
	struct I yval;

	yval = assign_I(0,1);
	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	st = (struct slist *) ckalloc(sizeof(struct slist) * num_algns);
	
	sign = algns[0].sign;
	initialize_slist(st, 0, num_algns);
	sort_init_algns(st, algns, num_algns, INIT_SELF2);
	for( i = 0; i < num_algns; i++ ) {
		cur_id = st[i].id;
		cur_val = 0;
		if( algns[cur_id].sign != sign ) {
			fatalf("different sign: %d vs. %d\n", sign, algns[cur_id].sign);
		}
		else if( src.lower > algns[cur_id].x.upper ) {} 
		else if( src.upper < algns[cur_id].x.lower ) {} 
		else {
			if( (src.lower >= algns[cur_id].x.lower) && (src.lower <= algns[cur_id].x.upper) ) {
				if( sign == 0 ) {
					lo = find_yloc_one(algns[cur_id], f, abs(src.lower - algns[cur_id].x.lower), NO_GAP_INC);
				}
				else { // sign == 1
					hi = find_yloc_one(algns[cur_id], f, abs(src.lower - algns[cur_id].x.lower), NO_GAP_INC);
				}
			}
			else {
				if( sign == 0 ) {
					lo = algns[cur_id].y.lower;
				}
				else { // sign == 1
					hi = algns[cur_id].y.upper;
				}
			}

			if( (src.upper >= algns[cur_id].x.lower) && (src.upper <= algns[cur_id].x.upper) ) {
				if( sign == 0 ) {
					hi = find_yloc_one(algns[cur_id], f, abs(src.upper - algns[cur_id].x.lower), NO_GAP_INC);
				}
				else { // sign == 1
					lo = find_yloc_one(algns[cur_id], f, abs(src.upper - algns[cur_id].x.lower), NO_GAP_INC);
				}
			}
			else {
				if( sign == 0 ) {
					hi = algns[cur_id].y.upper;
				}
				else { // sign == 1
					lo = algns[cur_id].y.lower;
				}
			}

			yval = assign_I(lo, hi);
			if( proper_overlap(yval, dst) == true ) {
				cur_val = width(intersect(yval, dst));
			}

			if( cur_val <= 0 ) {}
			else if( max_val == -1 ) {
				max_id = cur_id;
				max_val = cur_val;
				if( debug_mode == TRUE ) {
					printf("candi: %d-%d\n", yval.lower, yval.upper);
				}
			}
			else if( cur_val > max_val ) {
				if( debug_mode == TRUE ) {
					printf("candi: %d-%d\n", yval.lower, yval.upper);
				}
				max_id = cur_id;
				max_val = cur_val;
			}
		}
	}

	if( max_id == -1 ) {
		if( debug_mode == TRUE ) {
			printf("Warning: alignment for %d-%d %d-%d not found in util_algns.c\n", src.lower, src.upper, dst.lower, dst.upper);
		}
		res = (float)(-1);
	}
	else {
		beg = find_xloc_one(algns[max_id], f, abs(src.lower - algns[max_id].x.lower), NO_GAP_INC);
		end = find_xloc_one(algns[max_id], f, abs(src.upper - algns[max_id].x.lower), NO_GAP_INC);

		strcpy(S1, "");
		strcpy(T1, "");
		get_nth_algn(S1, T1, algns[max_id].fid, 0, f, a_info, REG);
		len = strlen(S1)-1;
		if( end <= beg ) {
			if( (end + 50) < len ) end = end + 50;
			else end = len;

			if( (beg - 50) > 0 ) beg = beg - 50;
			else beg = 0;
		}

		res = cal_pid_maf_beg(S1, T1, beg, end);
		if( res == -1 ) {
			printf("warning! pid is %f\n", res);
		}
	}

	free(st);
	free(a_info);
	return(res);
}

int pick_sim_level_algn(struct DotList *algns, int algn_id, FILE *f, struct exons_list *list1, int num_list1, struct exons_list *list2, int num_list2, int sp_id)
{
	int beg = 0, end = 0, len = 0;
	struct b_list *a_info;
	int res = -1;
	int n = 0;
	int b1 = 0, b2 = 0;
	int i = 0, j = 0, k = 0, cur1 = 0, cur2 = 0;
	int sign = algns[algn_id].init_sign;
	int cur_pid = 0;
	struct I *cur_list_x1, *cur_list_x2;
	struct I *cur_list_y1, *cur_list_y2;
	int num_cur_list_x1 = 0, num_cur_list_x2 = 0;
	int num_cur_list_y1 = 0, num_cur_list_y2 = 0;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	cur_list_x1 = (struct I *) ckalloc(sizeof(struct I) * num_list1);
	cur_list_x2 = (struct I *) ckalloc(sizeof(struct I) * num_list1);
	cur_list_y1 = (struct I *) ckalloc(sizeof(struct I) * num_list2);
	cur_list_y2 = (struct I *) ckalloc(sizeof(struct I) * num_list2);
	
	len = algns[algn_id].x.upper - algns[algn_id].xr_diff - algns[algn_id].xr_offset - algns[algn_id].x.lower - algns[algn_id].xl_offset;
	beg = find_xloc_one(algns[algn_id], f, algns[algn_id].xl_offset+algns[algn_id].xl_diff, NO_GAP_INC);
	end = find_xloc_one(algns[algn_id], f, len, NO_GAP_INC);

	strcpy(S1, "");
	strcpy(T1, "");
	get_nth_algn(S1, T1, algns[algn_id].fid, 0, f, a_info, REG);
	n = strlen(S1);
	cur_pid = cal_pid_maf_beg(S1, T1, beg, end);

	j = 0;
	for( i = 0; i < num_list1; i++ ) {
		if( (list1[i].ctg_id == algns[algn_id].ctg_id1) && (list1[i].val >= (cur_pid-2)) ) {
			cur_list_x1[j] = assign_I(list1[i].reg.lower, list1[i].reg.upper);
			j++;
		}
		
		if( (list1[i].ctg_id == algns[algn_id].ctg_id2) && (list1[i].val >= (cur_pid-2)) ) {
			cur_list_x2[k] = assign_I(list1[i].reg.lower, list1[i].reg.upper);
			k++;
		}
	}
	num_cur_list_x1 = j;
	num_cur_list_x2 = k;
	
	j = 0;
	k = 0;
	for( i = 0; i < num_list2; i++ ) {
		if( (list2[i].ctg_id == algns[algn_id].ctg_id1) && (list2[i].val >= (cur_pid-2)) ) {
			cur_list_y1[j] = assign_I(list2[i].reg.lower, list2[i].reg.upper);
			j++;
		}
		
		if( (list2[i].ctg_id == algns[algn_id].ctg_id2) && (list2[i].val >= (cur_pid-2)) ) {
			cur_list_y2[k] = assign_I(list2[i].reg.lower, list2[i].reg.upper);
			k++;
		}
	}
	num_cur_list_y1 = j;
	num_cur_list_y2 = k;

	b1 = algns[algn_id].x.lower;
	if( sign == 0 ) b2 = algns[algn_id].y.lower;
	else if( sign == 1 ) b2 = algns[algn_id].y.upper;
	else {
		fatalf("unexpected sign for the alignment: %d, %d-%d,%d-%d\n", sign, algns[algn_id].x.lower, algns[algn_id].x.upper, algns[algn_id].y.lower, algns[algn_id].y.upper);
	}

	cur1 = 0;
	cur2 = 0;
	for( i = 0; i < n; i++ ) {
		if( S1[i] != '-') {
			if( sp_id == SELF2 ) {
				if((i < ERR_SM_TH) || ( i > (n-ERR_SM_TH)) ) {
					S1[i] = 'N';
				}
				else if( is_in_skip_int(b1+cur1, cur_list_y1, num_cur_list_y1) == true) {
					S1[i] = 'N';
				}
			}
			else {
				if((i < ERR_SM_TH) || ( i > (n-ERR_SM_TH)) ) {
					S1[i] = 'N';
				}
				else if(is_in_skip_int(b1+cur1, cur_list_x1, num_cur_list_x1) == true) {
					S1[i] = 'N';
				}
			}
			cur1++;
		}

		if( T1[i] != '-') {
			if( sp_id == SELF1 ) {
				if( sign == 0 ) {
					if((i < ERR_SM_TH) || ( i > (n-ERR_SM_TH)) ) {
						T1[i] = 'N';
					}
					else if(is_in_skip_int(b2+cur2, cur_list_x2, num_cur_list_x2) == true) {
						T1[i] = 'N';
					}
				}	
				else {
					if((i < ERR_SM_TH) || ( i > (n-ERR_SM_TH)) ) {
						T1[i] = 'N';
					}
					else if(is_in_skip_int(b2-cur2, cur_list_x2, num_cur_list_x2) == true) {
						T1[i] = 'N';
					}
				}
			}
			else {
				if( sign == 0 ) {
					if((i < ERR_SM_TH) || ( i > (n-ERR_SM_TH)) ) {
						T1[i] = 'N';
					}
					else if(is_in_skip_int(b2+cur2, cur_list_y2, num_cur_list_y2) == true) {
						T1[i] = 'N';
					}
				}	
				else {
					if((i < ERR_SM_TH) || ( i > (n-ERR_TH)) ) {
						T1[i] = 'N';
					}
					else if(is_in_skip_int(b2-cur2, cur_list_y2, num_cur_list_y2) == true) {
						T1[i] = 'N';
					}
				}
			}
			cur2++;
		}
	}
	res = cal_pid_maf_beg(S1, T1, beg, end);
	if( res == 0 ) res = -1;

	free(cur_list_x1);
	free(cur_list_x2);
	free(cur_list_y1);
	free(cur_list_y2);
	free(a_info);
	return(res);
}

bool is_in_skip_int(int loc, struct I *list, int num_list)
{
	int i;
	bool res = false;
	bool is_over = false;

	i = 0;
	while( (i < num_list) && (res == false) && (is_over == false) ) {
		if( in(loc, list[i]) == true ) res = true;
		else if( list[i].lower > loc ) is_over = true;
		i++;
	}

	return res;
}

float cal_pid_part_algn(struct DotList *algns, int algn_id, int left_diff, int right_diff, FILE *f, int sp_id)
{
	int beg = 0, end = 0, len = 0;
	struct b_list *a_info;
	float cur_pid;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	if( sp_id == SELF1 ) {
		len = algns[algn_id].x.upper - algns[algn_id].x.lower - right_diff;
		beg = find_xloc_one(algns[algn_id], f, left_diff, NO_GAP_INC);
		end = find_xloc_one(algns[algn_id], f, len, NO_GAP_INC);
	}
	else if( sp_id == SELF2 ) {
		if( algns[algn_id].init_sign == 0 ) {
			len = algns[algn_id].y.upper - algns[algn_id].y.lower - right_diff;
			beg = find_yloc_one(algns[algn_id], f, left_diff, GAP_INC_IN_Y);
			end = find_yloc_one(algns[algn_id], f, len, GAP_INC_IN_Y);
		}
		else if( algns[algn_id].init_sign == 1 ) {
			len = algns[algn_id].y.upper - algns[algn_id].y.lower - left_diff;
			beg = find_yloc_one(algns[algn_id], f, right_diff, GAP_INC_IN_Y);
			end = find_yloc_one(algns[algn_id], f, len, GAP_INC_IN_Y);
		}
		else {
			fatalf("invalid alignment sign %d\n", algns[algn_id].init_sign);
		}
	}
	else {
		fatalf("wrong sp_id %d\n", sp_id);
	}

	strcpy(S1, "");
	strcpy(T1, "");
	get_nth_algn(S1, T1, algns[algn_id].fid, 0, f, a_info, REG);
	cur_pid = cal_pid_maf_beg(S1, T1, beg, end);

	free(a_info);
	return(cur_pid);
}

void cal_algn_beg_end(struct DotList *algns, int id, struct I ov, int sp_id, FILE *f, int *beg, int *end)
{
	int len = 0;

  if( sp_id == SELF1 ) {
    len = ov.upper - algns[id].x.lower;
    *beg = find_xloc_one(algns[id], f, ov.lower-algns[id].x.lower, NO_GAP_INC);
    *end = find_xloc_one(algns[id], f, len, NO_GAP_INC);
  }
  else if( sp_id == SELF2 ) {
    if( algns[id].init_sign == 0 ) {
      len = ov.upper - algns[id].y.lower;
      *beg = find_yloc_one(algns[id], f, ov.lower-algns[id].y.lower, GAP_INC_IN_Y);
      *end = find_yloc_one(algns[id], f, len, GAP_INC_IN_Y);
    }
    else if( algns[id].init_sign == 1 ) {
      len = algns[id].y.upper - algns[id].y.lower - ov.lower;
      *beg = find_yloc_one(algns[id], f, algns[id].y.upper - ov.upper, GAP_INC_IN_Y);
      *end = find_yloc_one(algns[id], f, len, GAP_INC_IN_Y);
    }
    else {
      fatalf("invalid alignment sign %d\n", algns[id].init_sign);
    }
  }
  else {
    fatalf("wrong sp_id %d\n", sp_id);
  }
}
