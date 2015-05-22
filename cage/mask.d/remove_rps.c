#include "main.h"
#include "regions.h"
#include "remove_rps.h"
#include "find_merging.h"
#include "util_gen.h"
#include "util_i.h"
#include "util.h"
#include "kd_op.h"

extern int debug_mode;

void remove_overlapped(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, struct kdnode *pair_tree, struct perm_pt *pair_pts, int size_seq1, int size_seq2, struct DotList *pair_alg, int *num_pair, int sp_flag)
{
	struct slist *sorted;
	int *clist;
	int num_list = 0;
	int i = 0, j;
	int num_lines;
	int *longest;
	struct IntList *self_alg;
	int del_id = -1;
	int n_list = 0;
	int s_id = -1;
	struct IntList *self_list;
	int num_self = 0;
	int size;

	if( sp_flag == SP_1 ) size = size_seq1;
	else if( sp_flag == SP_2 ) size = size_seq2;

	longest = (int *) ckalloc(sizeof(int));
	self_alg = (struct IntList *) ckalloc(sizeof(struct IntList));
	num_lines = *num;
	sorted = (struct slist *) ckalloc(num_lines * sizeof(struct slist));
	clist = (int *) ckalloc(num_lines * sizeof(int));
	self_list = (struct IntList *) ckalloc(num_lines * sizeof(struct IntList));

	(*self_alg).reg = assign_I(0, 1);
	strcpy((*self_alg).name1, "");
	strcpy((*self_alg).name2, "");
	for( i = 0; i < num_lines; i++ ) {
		self_list[i].reg = assign_I(0, 1);
		strcpy(self_list[i].name1, "");
		strcpy(self_list[i].name2, "");
	}

	sort_by_yintercept(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		if( dots[sorted[i].id].sign == 2 )
		{
		}
		else if( width(dots[sorted[i].id].x) < S_EF_VALUE) {}
		else
		{
			num_list = find_overlapped(dots, num_lines, sorted, i, clist, longest, self_alg);
			if( num_list > 1 )
			{
				if( debug_mode == TRUE ) {
					printf("Warning 1: masked out a region of %d-%d\n", (*self_alg).reg.lower, (*self_alg).reg.upper);
				}

				for( j = 0; j < num_list; j++ )
				{
					dots[clist[j]].sign = 2;
				}
				strcpy((*self_alg).name1, dots[clist[0]].name1);
				strcpy((*self_alg).name2, dots[clist[0]].name2);
				self_list[num_self].reg = assign_I((*self_alg).reg.lower, (*self_alg).reg.upper);
				strcpy(self_list[num_self].name1, dots[clist[0]].name1);
				strcpy(self_list[num_self].name2, dots[clist[0]].name2);
				num_self++;
				throw_away_rps(dots, num_lines, *self_alg, tree, p_pts, size, X_SIDE);
				throw_away_rps_pair(pair_alg, (*num_pair), *self_alg, pair_tree, pair_pts, size_seq1, size_seq2, sp_flag); 
			}
			else if((abs(dots[sorted[i].id].x.lower - dots[sorted[i].id].y.lower) <= NEW_DIS_TH) || (abs(dots[sorted[i].id].x.upper - dots[sorted[i].id].y.upper) <= NEW_DIS_TH))
			{
				s_id = find_self_alg(dots, num_lines, sorted, sorted[i].id, self_alg, tree, p_pts, size);
				 if( s_id != -1 )
				 {
				 		if( debug_mode == TRUE ) {
							printf("Warning 2: masked out a region of %d-%d\n", (*self_alg).reg.lower, (*self_alg).reg.upper);
						}
						
						strcpy((*self_alg).name1, dots[s_id].name1);
						strcpy((*self_alg).name2, dots[s_id].name2);
						strcpy(self_list[num_self].name1, dots[s_id].name1);
						strcpy(self_list[num_self].name2, dots[s_id].name2);
						self_list[num_self].reg = assign_I((*self_alg).reg.lower, (*self_alg).reg.upper);
						num_self++;
						throw_away_rps(dots, num_lines, *self_alg, tree, p_pts, size, X_SIDE); 
						throw_away_rps_pair(pair_alg, (*num_pair), *self_alg, pair_tree, pair_pts, size_seq1, size_seq2, sp_flag); 
				 }
			}
		}

		if( (strcmp(dots[sorted[i].id].name1, dots[sorted[i].id].name2) == 0 ) && (proper_overlap(dots[sorted[i].id].x, dots[sorted[i].id].y)) )
		{
			if( (width(intersect(dots[sorted[i].id].x, dots[sorted[i].id].y)) >= O_TH) && (compute_closeness(dots, sorted[i].id, SELF_ID) <= T_DIS_TH))
			{
				recal_range(dots, sorted[i].id);
			}
		}

		num_list = 0;
		del_id = -1;
		n_list = 0;
	}

	for( i = 0; i < num_self; i++ )
	{
		throw_away_rps(dots, num_lines, self_list[i], tree, p_pts, size, Y_SIDE); 
	}

	free(longest);
	free(self_alg);
	free(clist);
	free(self_list);
	free(sorted);
}

struct I add_intervals(struct DotList *dots, int id, struct I temp)
{
	struct I res;
	int x1, y1;

	if( dots[id].x.lower < temp.lower )
	{
		x1 = dots[id].x.lower;
	}
	else x1 = temp.lower;

	if( dots[id].y.lower < x1 )
	{
		x1 = dots[id].y.lower;
	}			

	if( dots[id].x.upper > temp.upper )
	{
		y1 = dots[id].x.upper;
	}
	else y1 = temp.upper;

	if( dots[id].y.upper > y1 )
	{
		y1 = dots[id].y.upper;
	}			

	res = assign_I(x1, y1);
	return(res);
}

void recal_range(struct DotList *dots, int id)
{
	int c_prime, b_prime;
	int n;

	n = (int)((float)(width(dots[id].x)) / (float)(dots[id].y.lower - dots[id].x.lower)) + 1;
	b_prime =	dots[id].x.lower + n * (dots[id].y.lower - dots[id].x.lower);  
	c_prime = dots[id].y.lower + width(dots[id].x) - b_prime;	
	dots[id].x = assign_I(dots[id].x.lower, dots[id].x.lower + c_prime);
	dots[id].y = assign_I(b_prime, b_prime + c_prime);
}

int find_self_alg(struct DotList *dots, int num_lines, struct slist *sorted, int id, struct IntList *self_alg, struct kdnode *tree, struct perm_pt *p_pts, int size)
{
	int res = -1;
	int start, end;
	int x1, x2, y1, y2;
	int i, j = 0, n_list = 0;
	int *clist;
	struct perm_pt *temp_plist;
	int *longest;
	int tmp_id = 0;
	char name1[NAME_LEN], name2[NAME_LEN];

	longest = (int *) ckalloc(sizeof(int));
	temp_plist = (struct perm_pt *) ckalloc(sizeof(struct perm_pt) * num_lines);
	clist = (int *) ckalloc(sizeof(int) * num_lines);

	if( dots[id].x.lower < dots[id].y.lower ) x1 = dots[id].x.lower;
	else x1 = dots[id].y.lower;

	if( dots[id].x.upper < dots[id].y.upper ) x2 = dots[id].y.upper;
	else x2 = dots[id].x.upper;

	strcpy(name1, dots[id].name1);
	strcpy(name2, dots[id].name2);

	x1 = x1 - (2*RANGE_TH);
	y1 = x1;
	x2 = x2 + (2*RANGE_TH);
	y2 = size;

	start = find_pred_blk(tree, x1, y1);
	end = find_successor(tree, x2, y2);
	
	for( i = start; i <= end; i++ )
	{
		tmp_id = p_pts[i].id;

		if( (strcmp(name1, dots[tmp_id].name1) == 0) && (strcmp(name2, dots[tmp_id].name2) == 0) ) {
			temp_plist[j] = assign_pm_val(p_pts[i]);
			j++;
		}
	}

	quick_sort_plist_y(temp_plist, 0, j-1);

	i = 0;
	while( (i < j) && ( n_list <= 1) )
	{
		n_list = 0;

		if( (strcmp(name1, dots[temp_plist[i].id].name1) != 0) || (strcmp(name2, dots[temp_plist[i].id].name2) != 0) ) {
			fatalf("%s, %s do not match %s, %s\n", dots[temp_plist[i].id].name1, dots[temp_plist[i].id].name2, name1, name2);	
		}

		if( (dots[temp_plist[i].id].sign != 2) && (overlap(dots[temp_plist[i].id].x, dots[temp_plist[i].id].y) == true) && ((tight_subset(dots[id].x, dots[temp_plist[i].id].x) == true) || (tight_subset(dots[id].y, dots[temp_plist[i].id].x) == true)))
		{
			(*self_alg).reg = assign_I(dots[temp_plist[i].id].x.lower, dots[temp_plist[i].id].x.upper);
			n_list = find_clist(dots, num_lines, sorted, temp_plist[i].id, clist, longest, self_alg, NO_OL);
		} 
		i++;
	}

	if( n_list > 1 ) res = temp_plist[i-1].id; 

	free(longest);
	free(temp_plist);
	free(clist);
	return(res);
}

void throw_away_rps(struct DotList *dots, int num_lines, struct IntList self_alg, struct kdnode *tree, struct perm_pt *p_pts, int size, int flag) 
{
	int start = 0, end = 1;
	int x1 = 0, x2 = 0, y1 = 1, y2 = 1;
	int i = 0, j = 0;
	int *longest = NULL;
	struct perm_pt *temp_plist = NULL;

	longest = (int *) ckalloc(sizeof(int));
	temp_plist = (struct perm_pt *) ckalloc(num_lines * sizeof(struct perm_pt));

	if( flag == X_SIDE )
	{
		x1 = self_alg.reg.lower - (2*RANGE_TH);
		y1 = self_alg.reg.lower - (2*RANGE_TH);
		x2 = self_alg.reg.upper + (3*RANGE_TH);
		y2 = size;
	}
//	else if( flag == Y_SIDE )
	else // flag == Y_SIDE )
	{
		x1 = 1;
		y1 = self_alg.reg.lower - (3*RANGE_TH);
		x2 = self_alg.reg.upper + (2*RANGE_TH);
		y2 = self_alg.reg.upper + (2*RANGE_TH);
	}

	start = find_pred_blk(tree, x1, y1);
	end = find_successor(tree, x2, y2);
	
	for( i = start; i <= end; i++ )
	{
		if( ((flag == X_SIDE) && (strcmp(dots[p_pts[i].id].name1, self_alg.name1) == 0) ) || ((flag == Y_SIDE) && (strcmp(dots[p_pts[i].id].name2, self_alg.name2) == 0)) ) {
			temp_plist[j] = assign_pm_val(p_pts[i]);
			j++;
		}
	}

	for( i = 0; i < j; i++ )
	{
		if(flag == X_SIDE)
		{
			if((dots[temp_plist[i].id].sign != 2) && (strict_subset(dots[temp_plist[i].id].x, self_alg.reg) == true))
			{
				if( strcmp(dots[temp_plist[i].id].name1, self_alg.name1) != 0 ) {
					fatalf("Error: checking on different contigs %s, %s\n", dots[temp_plist[i].id].name1, self_alg.name1);
				}

				if(check_neighbors(dots, i, temp_plist, j, self_alg.reg, flag) == true) dots[temp_plist[i].id].sign = 2;
			} 
		}
		else if(flag == Y_SIDE)
		{
			if((dots[temp_plist[i].id].sign != 2) && (strict_subset(dots[temp_plist[i].id].y, self_alg.reg) == true))
			{
				if( strcmp(dots[temp_plist[i].id].name2, self_alg.name2) != 0 ) {
					fatalf("Error: checking on different contigs %s, %s\n", dots[temp_plist[i].id].name2, self_alg.name2);
				}

				if(check_neighbors(dots, i, temp_plist, j, self_alg.reg, flag) == true) dots[temp_plist[i].id].sign = 2;
			} 
		}
	}

	free(longest);
	free(temp_plist);
}

void throw_away_rps_pair(struct DotList *dots, int num_lines, struct IntList self_alg, struct kdnode *tree, struct perm_pt *p_pts, int size_seq1, int size_seq2, int sp_flag) 
{
	int start = 0, end = 1;
	int x1 = 0, x2 = 0, y1 = 1, y2 = 1;
	int i = 0, j = 0;
	int *longest = NULL;
	struct perm_pt *temp_plist = NULL;

	longest = (int *) ckalloc(sizeof(int));
	temp_plist = (struct perm_pt *) ckalloc(num_lines * sizeof(struct perm_pt));
	
	if( sp_flag == SP_1 )
	{
		x1 = self_alg.reg.lower - (2*RANGE_TH);
		y1 = 1;
		x2 = self_alg.reg.upper + (3*RANGE_TH);
		y2 = size_seq2;
	}
	else // sp_flag == SP_2 
//	else if( sp_flag == SP_2 )
	{
		x1 = 1;
		y1 = self_alg.reg.lower - (2*RANGE_TH);
		x2 = size_seq1;
		y2 = self_alg.reg.upper + (3*RANGE_TH);
	}

	start = find_pred_blk(tree, x1, y1);
	end = find_successor(tree, x2, y2);
	
	for( i = start; i <= end; i++ )
	{
		if( ((sp_flag == SP_1) && (strcmp(dots[p_pts[i].id].name1, self_alg.name1) == 0)) || ((sp_flag == SP_2) && (strcmp(dots[p_pts[i].id].name2, self_alg.name2) == 0)) ) {
			temp_plist[j] = assign_pm_val(p_pts[i]);
			j++;
		}
	}

	for( i = 0; i < j; i++ )
	{
		if(sp_flag == SP_1)
		{
			if((dots[temp_plist[i].id].sign != 2) && (strict_subset(dots[temp_plist[i].id].x, self_alg.reg) == true))
			{
				if( strcmp(dots[temp_plist[i].id].name1, self_alg.name1) != 0 ) {
					fatalf("Error: checking on different contigs %s, %s\n", dots[temp_plist[i].id].name1, self_alg.name1);
				}

				if(check_neighbors(dots, i, temp_plist, j, self_alg.reg, X_SIDE) == true) dots[temp_plist[i].id].sign = 2;
			} 
		}
		else if(sp_flag == SP_2)
		{
			if((dots[temp_plist[i].id].sign != 2) && (strict_subset(dots[temp_plist[i].id].y, self_alg.reg) == true))
			{
				if( strcmp(dots[temp_plist[i].id].name2, self_alg.name2) != 0 ) {
					fatalf("Error: checking on different contigs %s, %s\n", dots[temp_plist[i].id].name2, self_alg.name2);
				}

				if(check_neighbors(dots, i, temp_plist, j, self_alg.reg, Y_SIDE) == true) dots[temp_plist[i].id].sign = 2;
			} 
		}
	}

	free(longest);
	free(temp_plist);
}

void reduce_rps(struct DotList *dots, int num_lines, struct slist *sorted, int cur, struct kdnode *tree, struct perm_pt *p_pts, int size)
{
	int start, end;
	int x1, x2, y1, y2;
	int id;
	int i, j = 0, n_list = 0;
	int *clist;
	int *longest;
	struct perm_pt *temp_plist;

	longest = (int *) ckalloc(sizeof(int));
	clist = (int *) ckalloc(num_lines * sizeof(int));
	temp_plist = (struct perm_pt *) ckalloc(num_lines * sizeof(struct perm_pt));

	id = sorted[cur].id;
	x1 = dots[id].x.lower - RANGE_TH;
	x2 = dots[id].x.upper + RANGE_TH;
	y1 = dots[id].y.lower - RANGE_TH;
	y2 = size;

	start = find_pred_blk(tree, x1, y1);
	end = find_successor(tree, x2, y2);
	
	for( i = start; i <= end; i++ )
	{
		temp_plist[j] = assign_pm_val(p_pts[i]);
		j++;
	}
	quick_sort_plist_x(temp_plist, 0, j-1);

	for( i = 0; i < j; i++ )
	{
		if( (temp_plist[i].id != id) && (dots[temp_plist[i].id].sign != 2) && (proper_overlap(dots[temp_plist[i].id].y, dots[id].y) == true))
		{
			n_list = find_rps(dots, temp_plist, i, 0, j-1, clist, longest);
			if( n_list > 1 )
			{
				for( j = 0; j < n_list; j++ )
				{
					if( clist[j] != (*longest) )
					{
						dots[clist[j]].sign = 2;
					}
				}
			}
		} 
	}
	free(longest);
	free(clist);
	free(temp_plist); 
}

int find_overlapped(struct DotList *dots, int num, struct slist *sorted, int cur, int *clist, int *longest, struct IntList *self_alg)
{
	int n_list = 0;
	int cur_id;

	cur_id = sorted[cur].id;

	if( (proper_overlap(dots[cur_id].x, dots[cur_id].y) == true) && (strcmp(dots[cur_id].name1, dots[cur_id].name2) == 0) )
	{
		if( (width(intersect(dots[cur_id].x, dots[cur_id].y)) >= O_TH) && (compute_closeness(dots, cur_id, SELF_ID) <= T_DIS_TH) ) 
		{
			n_list = find_clist(dots, num, sorted, cur_id, clist, longest, self_alg, SELF); 
		}
	}
	else if( (dots[cur_id].sign != 2) && (width(dots[cur_id].x) >= MID2_SELF) ) {
		remove_odd_algns(dots, num, sorted, cur_id);
	}

	return(n_list);
}

void remove_odd_algns(struct DotList *dots, int num, struct slist *sorted, int cur_id)
{
	int i = 0, j = 0;
	float r = (float)0;
	int a1 = 0, b1 = 0, a2 = 1, b2 = 1;
	float c1 = (float)0, c2 = (float)0;
	int x1 = 0, y1 = 0;

	if( dots[cur_id].sign == 0 ) {
		x1 = dots[cur_id].x.lower;
		y1 = dots[cur_id].y.lower;
	}
	else if( dots[cur_id].sign == 1 ) {
		x1 = dots[cur_id].x.lower;
		y1 = dots[cur_id].y.upper;
	}
	r = (float)((float)width(dots[cur_id].x)/(float)width(dots[cur_id].y));
	if( dots[cur_id].sign == 1 ) {
		r = (float)(-1) * r;
	}

	for( i = 0; i < num; i++ )
	{
		j = sorted[i].id;
		if( dots[j].sign == 0 ) {
			a1 = dots[j].x.lower;
			b1 = dots[j].y.lower;
			a2 = dots[j].x.upper;
			b2 = dots[j].y.upper;
		}
		else if( dots[j].sign == 1 ) {
			a1 = dots[j].x.lower;
			b1 = dots[j].y.upper;
			a2 = dots[j].x.upper;
			b2 = dots[j].y.lower;
		}

		if( (j == cur_id) || (dots[j].sign == 2) || (dots[cur_id].sign == 2) ) {}
		else if((width(dots[j].x) >= MID2_SELF) && (abs(dots[cur_id].identity - dots[j].identity) > 10) && (strcmp(dots[cur_id].name1, dots[j].name1) == 0) && (strcmp(dots[cur_id].name2, dots[j].name2) == 0) && (((subset(dots[cur_id].x, dots[j].x) == true) && (subset(dots[cur_id].y, dots[j].y) == true)) || ((subset(dots[j].x, dots[cur_id].x) == true) && (subset(dots[j].y, dots[cur_id].y) == true)))) {
			c1 = (float)(a1-x1)-r*((float)(b1-y1));
			c2 = (float)(a2-x1)-r*((float)(b2-y1));

			if( (c1*c2) < ((float)0) ) {
				if(dots[cur_id].identity < dots[j].identity) {
					dots[cur_id].sign = 2;
				}
				else {
					dots[j].sign = 2;
				}
			}
		} 
	}
}

void remove_rps(struct DotList *dots, int *num)
{
	struct slist *sorted;
	int *clist;
	int num_list = 0;
	int i = 0, j;
	int num_lines;
	int *longest;
	struct IntList *self_alg;

	longest = (int *) ckalloc(sizeof(int));
	self_alg = (struct IntList *) ckalloc(sizeof(struct IntList));
	num_lines = *num;
	sorted = (struct slist *) ckalloc(num_lines * sizeof(struct slist));
	clist = (int *) ckalloc(num_lines * sizeof(int));

	sort_by_yintercept(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		if( dots[sorted[i].id].sign == 2 )
		{
		}
		else
		{
			num_list = find_clist(dots, num_lines, sorted, sorted[i].id, clist, longest, self_alg, NO_OL);
			if( num_list > 1 )
			{
				for( j = 0; j < num_list; j++ )
				{
					if( clist[j] != (*longest) )
					{
						dots[clist[j]].sign = 2;
					}
				}
			}
			num_list = 0;
		}
	}

	overwrite_dots(num, dots);

	free(longest);
	free(self_alg);
	free(sorted);
	free(clist);
}

int find_rps(struct DotList *dots, struct perm_pt *p_pts, int cur, int start, int end, int *clist, int *longest)
{
	int n_list = 0;
	int cur_n_list;
	int i;
	int j;
	int closeness;
	bool d_flag; 
	int s_len;

	s_len = width(dots[p_pts[cur].id].y);
	*longest = p_pts[cur].id;
	clist[n_list] = p_pts[cur].id;
	n_list++;

	for( i = start; i <= end; i++ )
	{
		j = 0;
		cur_n_list = n_list;
		d_flag = false;
		if( check_contain(p_pts[i].id, clist, cur_n_list) == false )
		{
			while((d_flag == false) && (j < cur_n_list))
			{
				if((dots[p_pts[i].id].sign != 2) && (strcmp(dots[p_pts[i].id].name1, dots[clist[j]].name1) == 0) && (strcmp(dots[p_pts[i].id].name2, dots[clist[j]].name2) == 0) && (dots[p_pts[i].id].sign == dots[clist[j]].sign) && (proper_overlap(dots[p_pts[i].id].y, dots[clist[j]].y) == true))
				{
					if( width(dots[p_pts[i].id].y) < s_len ) 
					{
						s_len = width(dots[p_pts[i].id].y);
					}

					if( width(intersect(dots[p_pts[i].id].y, dots[clist[j]].y)) >= s_len )
					{
						closeness = compute_closeness(dots, clist[j], p_pts[i].id);
						if( closeness <= T_DIS_TH )
						{
							if( width(dots[*longest].x) < width(dots[p_pts[i].id].x) )
							{
								if( ( dots[*longest].identity - 2 ) <= dots[p_pts[i].id].identity )
								{
									*longest = p_pts[i].id;
								}
							}
							clist[n_list] = p_pts[i].id;
							n_list++;
							d_flag = true;
						}
					}
				}
				j++;
			}
		}
	}
	return(n_list);
}

int find_clist(struct DotList *dots, int num, struct slist *sorted, int cur, int *clist, int *longest, struct IntList *self_alg, int flag)
{
	int n_list = 0;
	int cur_n_list = 0;
	int i = 0, j = 0;
	int d = 0, closeness = 0;
	bool d_flag = false;
	int temp_l = 0, temp_u = 0;
	int dis_th = 0;
	float w_ratio = (float)0;

	(*self_alg).reg = assign_I(dots[cur].x.lower, dots[cur].x.upper);
	strcpy((*self_alg).name1, dots[cur].name1);
	strcpy((*self_alg).name2, dots[cur].name2);
	*longest = cur;
	clist[n_list] = cur;
	n_list++;

	if( width((*self_alg).reg) >= LONG_SELF ) dis_th = T_DIS_TH;
	else if( width((*self_alg).reg) >= MID_SELF ) dis_th = MID_SELF;
	else if( width((*self_alg).reg) >= MID2_SELF ) dis_th = MID2_SELF;
	else dis_th = SELF_T_DIS_TH;

	for( i = 0; i < num; i++ )
	{
		j = 0;
		cur_n_list = n_list;
		d_flag = false;
		if( check_contain(sorted[i].id, clist, cur_n_list) == false )
		{
			while((d_flag == false) && (j < cur_n_list))
			{
				if((dots[sorted[i].id].sign != 2) && (strcmp(dots[sorted[i].id].name1, dots[clist[j]].name1) == 0) && (strcmp(dots[sorted[i].id].name2, dots[clist[j]].name2) == 0) && (dots[sorted[i].id].sign == dots[clist[j]].sign))
				{
					d = overlap_len(dots, sorted[i].id, clist[j], (*self_alg).reg);
					w_ratio = ((float)d)/((float)(width(dots[clist[j]].x)));
					closeness = compute_closeness(dots, clist[j], sorted[i].id);

					if( (closeness <= dis_th) && ((d >= (2*RP_BD)) || (w_ratio >= 0.4)))
					{
						if((flag == OVERLAP) || (flag == SELF))
						{
							if( proper_overlap(dots[sorted[i].id].x, dots[sorted[i].id].y) == false )
							{
								if( (*longest) == cur )
								{
									*longest = sorted[i].id;
								}
								else 
								{
									if( width(dots[*longest].x) < width(dots[sorted[i].id].x) )
									{
										*longest = sorted[i].id;
									}
								}
							}

							if( flag == SELF )
							{
								temp_l = (*self_alg).reg.lower;
								temp_u = (*self_alg).reg.upper;
								if( (dots[sorted[i].id].x.lower < temp_l) && (proper_overlap(dots[sorted[i].id].x, (*self_alg).reg) == true) )
								{
									temp_l = dots[sorted[i].id].x.lower;	
								}

								if( (dots[sorted[i].id].y.lower < temp_l) && (proper_overlap(dots[sorted[i].id].y, (*self_alg).reg) == true) )
								{
									temp_l = dots[sorted[i].id].y.lower;	
								}

								if( (dots[sorted[i].id].x.upper > temp_u) && ( proper_overlap(dots[sorted[i].id].x, (*self_alg).reg) == true) )
								{
									temp_u = dots[sorted[i].id].x.upper;	
								}

								if( (dots[sorted[i].id].y.upper > temp_u) && (proper_overlap(dots[sorted[i].id].y, (*self_alg).reg) == true))
								{
									temp_u = dots[sorted[i].id].y.upper;	
								}

								if( width(dots[*longest].x) < width(dots[sorted[i].id].x) )
								{
									*longest = sorted[i].id;
								}
								(*self_alg).reg = assign_I(temp_l, temp_u);
							}
						}
						else if( flag == NO_OL )
						{
							if( width(dots[*longest].x) < width(dots[sorted[i].id].x) )
							{
								if( ( dots[*longest].identity - 2 ) <= dots[sorted[i].id].identity )
								{
									*longest = sorted[i].id;
								}
							}
							temp_l = (*self_alg).reg.lower;
							temp_u = (*self_alg).reg.upper;
							if( (dots[sorted[i].id].x.lower < temp_l) && (proper_overlap(dots[sorted[i].id].x, (*self_alg).reg) == true) )
							{
								temp_l = dots[sorted[i].id].x.lower;	
							}

							if( (dots[sorted[i].id].x.upper > temp_u) && ( proper_overlap(dots[sorted[i].id].x, (*self_alg).reg) == true) )
							{
								temp_u = dots[sorted[i].id].x.upper;	
							}
							(*self_alg).reg = assign_I(temp_l, temp_u);
						}
						clist[n_list] = sorted[i].id;
						n_list++;
						d_flag = true;
					}
				}
				j++;
			}
		}
	}

	n_list = check_longest_one(dots, clist, self_alg, flag, n_list);

	return(n_list);
}


int check_longest_one(struct DotList *dots, int *clist, struct IntList *self_alg, int flag, int num_list)
{
	int num = num_list;
	int l_id;
	int sec_id;
	int max_len, sec_len;
	bool remove_longest = false;
	int i, j;
	int temp_l, temp_u;

	max_len = 0;
	sec_len = 0;	
	for( i = 0; i < num; i++ )
	{
		if( max_len < width(dots[clist[i]].x) )
		{
			l_id = i;
			max_len = width(dots[clist[i]].x);
		}
	}

	for( i = 0; i < num; i++ )
	{
		if( i != l_id )
		{
			if( sec_len < width(dots[clist[i]].x) )
			{
				sec_id = i;
				sec_len = width(dots[clist[i]].x);
			}
		}
	}

	if( max_len > (2*sec_len) )
	{
		remove_longest = true;	

		j = 0;
		for( i = 0; i < num; i++ )
		{
			if( i == l_id ){}
			else
			{
				clist[j] = clist[i];
			}
		}
		num = num - 1;
	}

	if( remove_longest == true )
	{
		temp_l = dots[clist[0]].x.lower;
		temp_u = dots[clist[0]].x.upper;

		for( i = 0; i < num; i++ )
		{
			if( flag == SELF )
			{
				if( dots[clist[i]].x.lower < temp_l )
				{
					temp_l = dots[clist[i]].x.lower;	
				}

				if( dots[clist[i]].y.lower < temp_l )
				{
					temp_l = dots[clist[i]].y.lower;	
				}

				if( dots[clist[i]].x.upper > temp_u )
				{
					temp_u = dots[clist[i]].x.upper;	
				}

				if( dots[clist[i]].y.upper > temp_u )
				{
					temp_u = dots[clist[i]].y.upper;	
				}

				(*self_alg).reg = assign_I(temp_l, temp_u);
			}
			else if( flag == NO_OL )
			{
				if( dots[clist[i]].x.lower < temp_l )
				{
					temp_l = dots[clist[i]].x.lower;	
				}

				if( dots[clist[i]].x.upper > temp_u )
				{
					temp_u = dots[clist[i]].x.upper;	
				}
				(*self_alg).reg = assign_I(temp_l, temp_u);
			}
		}
	}
	return(num);
}

bool check_contain(int cur, int *clist, int num_list)
{
	int i;
	bool res = false;

	for( i = 0; i < num_list; i++ )
	{
		if( cur == clist[i] ) res = true;
	}

	return(res);
}

int overlap_len(struct DotList *dots, int i, int id, struct I self)
{
	int len = 0;
	int temp_len;

	if( proper_overlap(dots[i].x, dots[id].x) == true )
	{
		len = width(intersect(dots[i].x, dots[id].x));	
	}

	if( proper_overlap(dots[i].y, dots[id].y) == true )
	{
		temp_len = width(intersect(dots[i].y, dots[id].y));
		if( len < temp_len )
		{
			len = temp_len;
		}
	}

	if( proper_overlap(dots[i].x, self) == true )
	{
		temp_len = width(intersect(dots[i].x, self));	
		if( len < temp_len )
		{
			len = temp_len;
		}
	}

	if( proper_overlap(dots[i].y, self) == true )
	{
		temp_len = width(intersect(dots[i].y, self));	
		if( len < temp_len )
		{
			len = temp_len;
		}
	}

	return(len);
}

bool check_neighbors(struct DotList *dots, int id, struct perm_pt *list, int num, struct I self_alg, int flag)
{
	int i;
	bool res = false;
	int dis_th;
	int closeness;
	int len;
	int cur_id;
	struct I cur_reg;

	len = width(self_alg);
	cur_id = list[id].id;
	
	if( len >= LONG_SELF ) dis_th = T_DIS_TH;
	else if( len >= MID_SELF ) dis_th = MID_SELF;
	else if( len >= MID2_SELF ) dis_th = MID2_SELF;
	else dis_th = SELF_T_DIS_TH;
	
	i = 0;
	while((i < num) && (res == false))
	{
		if( (cur_id != list[i].id) && (is_same_contigs(dots, cur_id, list[i].id) == true) ) {
			closeness = compute_closeness(dots, cur_id, list[i].id);
			if( flag == X_SIDE ) cur_reg = assign_I(dots[list[i].id].x.lower, dots[list[i].id].x.upper);
			else if( flag == Y_SIDE ) cur_reg = assign_I(dots[list[i].id].y.lower, dots[list[i].id].y.upper);

			if((d_tight_subset(cur_reg, self_alg) == true) && ( closeness <= dis_th ) ) res = true;
		}
		i++;
	}

	return res;
}

bool is_same_contigs(struct DotList *dots, int id1, int id2)
{
	bool res = false;

	if( (strcmp(dots[id1].name1, dots[id2].name1) == 0) && (strcmp(dots[id1].name2, dots[id2].name2) == 0) ) {
		res = true;
	}

	return(res);
}
