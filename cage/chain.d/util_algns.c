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

	sort_init_algns(st, algns, num_algns, INIT_SELF1);
	b = search_range_b(st, algns, num_algns, src.lower, SELF1);
	e = search_range_e(st, algns, num_algns, src.upper, SELF1);

	j = 0;
	for( i = b; i <= e; i++ ) {
		assign_algn(temp_algns, j, algns[st[i].id]);	
		j++;
	}
	count = j;
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
