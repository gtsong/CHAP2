#include "main.h"
#include "regions.h"
#include "get_algn.h"
#include "util.h"
#include "util_i.h"
#include "read_maf.h"
#include "write_maf.h"
#include "read_algn.h"

extern int debug_mode;

void get_nth_algn_ch(char *S1, char *T1, int index, struct DotList *init_algns, FILE *fp)
{
	int old_id, cur_id;
	int *sec_beg, *fst_end;
	int f_len; // f_len is the length of nucleotides already saved in S1 and T1 
	int b;
	struct I temp;
	int y_old, y_cur;
	struct b_list *a_info;
	int i;
	int beg1 = 0, beg2 = 0;
	int rp1_id, rp2_id;
	int gap_len1, gap_len2;
	int t_x1, t_x2, t_y1, t_y2;
	int e1, e2;
	int cid;
	int beg;
	struct I cur_x, cur_y, old_x, old_y;
	struct r_list *rp1, *rp2;
	int fid;

	fid = init_algns[index].fid;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	sec_beg = (int *) ckalloc(sizeof(int));
	fst_end = (int *) ckalloc(sizeof(int));
	rp1 = (struct r_list *) ckalloc(sizeof(struct r_list));
	rp2 = (struct r_list *) ckalloc(sizeof(struct r_list));
	
  rp1->start = 0;
  rp1->end = 1;
  rp1->len = 1;
  rp1->d_rate = 100;
  strcpy(rp1->rn, "none");
  rp2->start = 0;
  rp2->end = 1;
  rp2->len = 1;
  rp2->d_rate = 100;
  strcpy(rp2->rn, "none");

  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;
	*sec_beg = 0;
	*fst_end = 0;
	b = init_algns[index].x.lower + init_algns[index].xl_diff; // 1-base
	beg = find_xloc_one(init_algns[index], fp, init_algns[index].xl_diff, GAP_INC);
	beg1 = beg-init_algns[index].x.lower;
	cid = init_algns[index].c_id;
	if( cid == -1 ) {
		get_nth_algn(S1, T1, fid, beg - init_algns[index].x.lower, fp, a_info, REG);
		cur_id = -1;
	}
	else {
		f_len = 0;
		get_nth_algn(S1, T1, fid, beg - init_algns[index].x.lower, fp, a_info, REG); 
		old_id = index;
		cur_id = init_algns[index].c_id;
		beg = find_xloc_one(init_algns[cur_id], fp, init_algns[cur_id].xl_diff, GAP_INC);
		beg2 = beg-init_algns[cur_id].x.lower;
	}

	while(cur_id != -1) {
		*sec_beg = beg2;

		rp1_id = -1;
		rp2_id = -1;
		if( debug_mode == true ) {
			if( rp1_id != -1 ) {
				printf("A repeat in %d-%d inserted in seq1\n", rp1[rp1_id].start, rp1[rp1_id].end);
			}
			else if( rp2_id != -1 ) {
				printf("A repeat in %d-%d inserted in seq2\n", rp2[rp2_id].start, rp2[rp2_id].end);
			}
			else printf("No repeat insertead but chained\n"); 
		}

		cur_x = assign_I(init_algns[cur_id].x.lower + init_algns[cur_id].xl_diff, init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff);
		cur_y = assign_I(init_algns[cur_id].y.lower + init_algns[cur_id].yl_diff, init_algns[cur_id].y.upper - init_algns[cur_id].yr_diff);
		old_x = assign_I(init_algns[old_id].x.lower + init_algns[old_id].xl_diff, init_algns[old_id].x.upper - init_algns[old_id].xr_diff);
		old_y = assign_I(init_algns[old_id].y.lower + init_algns[old_id].yl_diff, init_algns[old_id].y.upper - init_algns[old_id].yr_diff);

		if( cur_x.lower < old_x.lower ) {
			printf("%d-%d, %d-%d and %d-%d, %d-%d\n", cur_x.lower, cur_x.upper, cur_y.lower, cur_y.upper, old_x.lower, old_x.upper, old_y.lower, old_y.upper);
			fatalf("chaining info error 1\n");
		}
		else if( (proper_overlap(old_x, cur_x) == true) && (proper_overlap(old_y, cur_y) == true) ) 
		{
			temp = intersect(old_x, cur_x);
			if( init_algns[cur_id].sign == 0 ) {
				y_cur = cur_y.lower;
				y_old = find_yloc_one(init_algns[old_id], fp, temp.lower-old_x.lower, NO_GAP_INC);
			}
			else if( init_algns[cur_id].sign == 1 ) {
				y_cur = find_yloc_one(init_algns[old_id], fp, temp.lower-old_x.lower, NO_GAP_INC);
				y_old = cur_y.upper;
			}

			if(y_old <= y_cur) {
				adjust_two_ch_algns(init_algns, init_algns[old_id].index, init_algns[cur_id].index, beg1, beg2, fp, fst_end, sec_beg, true, rp1, rp2, rp1_id, rp2_id);
				S1[f_len + (*fst_end)] = '\0';
				T1[f_len + (*fst_end)] = '\0';
				t_x1 = find_xloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
				t_x2 = find_xloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
				t_y1 = find_yloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
				t_y2 = find_yloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
				gap_len1 = abs(t_x2 - t_x1 - 1);
				gap_len2 = abs(t_y2 - t_y1 - 1);
				fill_N(gap_len1, gap_len2, S1, T1);
			}
			else if( y_old > y_cur ) {
				adjust_two_ch_algns(init_algns, init_algns[old_id].index, init_algns[cur_id].index, beg1, beg2, fp, fst_end, sec_beg, false, rp1, rp2, rp1_id, rp2_id);											
				S1[f_len + (*fst_end)] = '\0';
				T1[f_len + (*fst_end)] = '\0';
				t_x1 = find_xloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
				t_x2 = find_xloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
				t_y1 = find_yloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
				t_y2 = find_yloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
				gap_len1 = abs(t_x2 - t_x1 - 1);
				gap_len2 = abs(t_y2 - t_y1 - 1);
				fill_N(gap_len1, gap_len2, S1, T1);
			}
		}
		else if( proper_overlap(old_x, cur_x) == true ) 
		{ 
			adjust_two_ch_algns(init_algns, init_algns[old_id].index, init_algns[cur_id].index, beg1, beg2, fp, fst_end, sec_beg, true, rp1, rp2, rp1_id, rp2_id);											
			S1[f_len + (*fst_end)] = '\0';
			T1[f_len + (*fst_end)] = '\0';
			t_x1 = find_xloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
			t_x2 = find_xloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
			t_y1 = find_yloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
			t_y2 = find_yloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
			gap_len1 = abs(t_x2 - t_x1 - 1);
			gap_len2 = abs(t_y2 - t_y1 - 1);
			fill_N(gap_len1, gap_len2, S1, T1);
		}
		else if( proper_overlap(old_y, cur_y) == true ) 
		{
			adjust_two_ch_algns(init_algns, init_algns[old_id].index, init_algns[cur_id].index, beg1, beg2, fp, fst_end, sec_beg, false, rp1, rp2, rp1_id, rp2_id);											
			S1[f_len + (*fst_end)] = '\0';
			T1[f_len + (*fst_end)] = '\0';
			t_x1 = find_xloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
			t_x2 = find_xloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
			t_y1 = find_yloc_one(init_algns[old_id], fp, (*fst_end), GAP_INC) - 1; // adjusted endpoint
			t_y2 = find_yloc_one(init_algns[cur_id], fp, (*sec_beg), GAP_INC); // adjusted start point
			gap_len1 = abs(t_x2 - t_x1 - 1);
			gap_len2 = abs(t_y2 - t_y1 - 1);
			fill_N(gap_len1, gap_len2, S1, T1);
		}
		else if( (init_algns[cur_id].sign == 0) && (cur_x.lower >= old_x.upper) && (cur_y.lower >= old_y.upper) ) {
		*sec_beg = beg2;
			gap_len1 = abs(cur_x.lower - old_x.upper);
			gap_len2 = abs(cur_y.lower - old_y.upper);
			fill_N(gap_len1, gap_len2, S1, T1);
		}
		else if( (init_algns[cur_id].sign == 1) && (cur_x.lower >= old_x.upper) && (old_y.lower >= cur_y.upper) ) {
		*sec_beg = beg2;
			gap_len1 = abs(cur_x.lower - old_x.upper);
			gap_len2 = abs(old_y.lower - cur_y.upper);
			fill_N(gap_len1, gap_len2, S1, T1);
		}
		else {
			printf("%d-%d, %d-%d and %d-%d, %d-%d\n", cur_x.lower, cur_x.upper, cur_y.lower, cur_y.upper, old_x.lower, old_x.upper, old_y.lower, old_y.upper);
			fatalf("chaining info error 2\n");
		}

		if( S1[strlen(S1)-1] == '\n' ) f_len = strlen(S1) - 1 - (*sec_beg);
		else f_len = strlen(S1) - (*sec_beg);
		get_nth_algn(S1, T1, init_algns[cur_id].index, (*sec_beg), fp, a_info, EXT);
		old_id = cur_id;
		cur_id = init_algns[old_id].c_id;
		cur_x = assign_I(init_algns[cur_id].x.lower + init_algns[cur_id].xl_diff, init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff);
		cur_y = assign_I(init_algns[cur_id].y.lower + init_algns[cur_id].yl_diff, init_algns[cur_id].y.upper - init_algns[cur_id].yr_diff);
		old_x = assign_I(init_algns[old_id].x.lower + init_algns[old_id].xl_diff, init_algns[old_id].x.upper - init_algns[old_id].xr_diff);
		old_y = assign_I(init_algns[old_id].y.lower + init_algns[old_id].yl_diff, init_algns[old_id].y.upper - init_algns[old_id].yr_diff);

		beg1 = find_xloc_one(init_algns[old_id], fp, *sec_beg, GAP_INC) - old_x.lower;
		if( cur_x.lower < (old_x.lower + beg1) ) {
			beg2 = old_x.lower + beg1 - cur_x.lower;
		}
		else {
			beg = find_xloc_one(init_algns[cur_id], fp, init_algns[cur_id].xl_diff, GAP_INC);
			beg2 = beg-init_algns[cur_id].x.lower;
		}
	}

	if( S1[strlen(S1)-1] == '\n' ) f_len = strlen(S1)-1;
	else f_len = strlen(S1);

	e1 = 0;
	e2 = 0;
	for( i = 0; i < f_len; i++ ) {
		if( strchr("ACGTN", toupper(S1[i])) ) e1++;
		if( strchr("ACGTN", toupper(T1[i])) ) e2++;
	}

	free(rp1);
	free(rp2);
	free(fst_end);
	free(sec_beg);
	free(a_info);
}
