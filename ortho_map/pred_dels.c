#include "main.h"
#include "regions.h"
#include "pred_dels.h"
#include "deal_gaps.h"
#include "find_merging.h"
#include "util_gen.h"
#include "find_dup_copy.h"
#include "check_copy.h"
#include "util.h"
#include "util_i.h"

extern int debug_mode;

int pred_dels(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct gap_list *del_list, int pid)
{
	struct slist *sorted;
	int i = 0, j;
	int num_lines;
	int *sd;
	bool *is_x;
	int xval, yval;
	int w_sid, w_fid, h_sid, h_fid;
	int opt_id = -1;
	int num_del = 0;
	struct gap_list *gp;
	int last_gid = 0;
	int *cid;
	bool *x_ins;
	int temp_score;
	struct I del_reg;

	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));
	gp = (struct gap_list *) ckalloc(sizeof(struct gap_list));
	cid = (int *) ckalloc(sizeof(int));
	x_ins = (bool *) ckalloc(sizeof(bool));
	sorted = (struct slist *) ckalloc(sizeof(struct slist) * (*num));

	num_lines = *num;

	initialize_slist(sorted, 0, num_lines);
	sort_list(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		opt_id = -1;
		if( (dots[sorted[i].id].sign == 2) || (dots[sorted[i].id].sign == 20) || (dots[sorted[i].id].sign == 21) ) {}
		else
		{
			xval = dots[sorted[i].id].x.lower;
			if( dots[sorted[i].id].sign == 0 )
			{
				yval = dots[sorted[i].id].y.lower;
			}
			else 
			{
				yval = dots[sorted[i].id].y.upper;
			}	
			
			w_sid = find_id(dots, tree, size, sorted[i].id, xval, yval, W_SID);
			w_fid = find_id(dots, tree, size, sorted[i].id, xval, yval, W_FID);
			h_sid = find_id(dots, tree, size, sorted[i].id, xval, yval, H_SID);
			h_fid = find_id(dots, tree, size, sorted[i].id, xval, yval, H_FID);

			opt_id = find_opt_del(dots, sorted[i].id, p_pts, w_sid, w_fid, h_sid, h_fid, gp);
		}

		if( opt_id != -1 ) 
		{
			temp_score = get_score_copy(dots, num_lines, *gp, cid, x_ins);

			if( (*gp).type == 1 ) 
			{
				*x_ins = true;
			}
			else if( (*gp).type == 2 )
			{
				*x_ins = false;
			}
			else 
			{
			}

			if( temp_score == -1 ) {
				if( (*gp).x1 >= (*gp).x2 ) del_reg = assign_I((*gp).x2-DEL_TH-1, (*gp).x1+DEL_TH+1);
				else del_reg = assign_I((*gp).x1-DEL_TH-1, (*gp).x2+DEL_TH+1);
			}

			if( temp_score != -1 )
			{
				opt_id = -1;
			}
			else if( check_whole_del_region_inclusion(dots, num_lines, del_reg, (*gp).id1, (*gp).id2) == true )
			{
				opt_id = -1;
			}
			else
			{	
				last_gid = assign_gap_list(del_list, num_del, last_gid, *gp);
				num_del++;
			}
		}	
	}

	if( num_del > 0 )
	{
		check_overlapped_del(dots, del_list, num_del);
		remove_dif_del(del_list, num_del, last_gid, pid, dots);
	}

	j = 0;
	for( i = 0; i < num_del; i++ )
	{
		if( del_list[i].gid == -1 ){}
		else
		{
			del_list[j] = assign_glist(del_list[i]);
			j++;
		}
	}
	num_del = j;

	free(sorted);
	free(x_ins);
	free(cid);
	free(gp);
	free(sd);
	free(is_x);
	return(num_del);
}

void remove_dif_del(struct gap_list *del_list, int num_del, int last_gid, int pid, struct DotList *dots)
{
	int i, j;
	int min_width = -1; // the most common width should be the minimum
	int max_pid = -1;
	int num_count = 0;

	for( i = 1; i <= last_gid; i++ )
	{
		num_count = 0;
		min_width = -1;
		for( j = 0; j < num_del; j++ )
		{
			if( del_list[j].gid == i )
			{
				if( min_width == -1 ) 
				{
					min_width = abs(del_list[j].y2 - del_list[j].y1);
				}
				else if( min_width > abs(del_list[j].y2 - del_list[j].y1) )
				{
					min_width = abs(del_list[j].y2 - del_list[j].y1);
				}

				if( max_pid < dots[del_list[j].id1].identity ) {
					max_pid = dots[del_list[j].id1].identity;
				}

				if( max_pid < dots[del_list[j].id2].identity ) {
					max_pid = dots[del_list[j].id2].identity;
				}
			}
		}

		for( j = 0; j < num_del; j++ )
		{
			if( del_list[j].gid == i )
			{
				if( max_pid < pid) {
					del_list[j].gid = -1;
				}
				else if( ((min_width >= 1000) && (abs(del_list[j].y2-del_list[j].y1) - min_width <= ERR_TH)) || ((((float)abs(del_list[j].y2 - del_list[j].y1) - (float)min_width) / ((float)min_width)) <= 0.05) ) 
				{
					num_count++;
/*
					if( dots[del_list[j].id1].pair_self == SELF ) {
						if( pair_self == -1 ) pair_self = SELF;
						else if( pair_self == PAIR ) pair_self = SELF_PAIR;
					}
					else if( dots[del_list[j].id1].pair_self == PAIR ) {
						if( pair_self == -1 ) pair_self = PAIR;
						else if( pair_self == SELF ) pair_self = SELF_PAIR;
					}
*/
				}
				else 
				{
//					del_list[j].gid = -1;
				}
			}
		}
		
//		if( (num_count <= 1) || (pair_self != SELF_PAIR)) {
		if( num_count <= 1) {
			for( j = 0; j < num_del; j++ ) {
				if( del_list[j].gid == i ) {
					del_list[j].gid = -1;
				}
			}
		}
	}
}

void check_overlapped_del(struct DotList *dots, struct gap_list *del_list, int num_del)
{
	int j = 0, k = 0;
	bool is_res = false;
	struct I left, right; // left: the left side region in deleted point
								 // rigth: the right side region in deleted point
	struct gap_list gp;

	left = assign_I(0, 1);
	right = assign_I(0, 1);

	for( j = 0; j < num_del; j++ )
	{
		gp = assign_glist(del_list[j]);

		if( gp.type == 1 ) // deleted in x region
		{
			left = assign_I(dots[gp.id1].x.lower, dots[gp.id1].x.upper);
			right = assign_I(dots[gp.id2].x.lower, dots[gp.id2].x.upper);
		}
		else if( gp.type == 2 ) // deleted in y region
		{
			if( dots[gp.id1].sign == 0 )
			{
				left = assign_I(dots[gp.id1].y.lower, dots[gp.id1].y.upper);
				right = assign_I(dots[gp.id2].y.lower, dots[gp.id2].y.upper);
			}
			else if( dots[gp.id1].sign == 1 )
			{
				left = assign_I(dots[gp.id2].y.lower, dots[gp.id2].y.upper);
				right = assign_I(dots[gp.id1].y.lower, dots[gp.id1].y.upper);
			}
		}

/*
		if( (dots[gp.id1].pair_self == PAIR) || (dots[gp.id2].pair_self == PAIR) ) 
		{
			is_res = true;
		}
*/
		if( (width(dots[gp.id1].x) >= EFFEC_DEL_LEN_TH) && (width(dots[gp.id2].x) >= EFFEC_DEL_LEN_TH) && (abs(dots[gp.id1].identity - dots[gp.id2].identity) >= DEL_ID_TH) ) 
		{
			is_res = true;
		}
		if( (dots[gp.id1].lock == BOTH_LOCK) || (dots[gp.id1].lock == LEFT_LOCK) || (dots[gp.id2].lock == RIGHT_LOCK) || (dots[gp.id2].lock == BOTH_LOCK)) 
		{
			is_res = true;
		}
		else
		{
			is_res = false;
		}

		if( is_res == true ) 
		{
			for( k = 0; k < num_del; k++ )
			{
				if( del_list[k].gid == del_list[j].gid ) del_list[k].gid = -1;
			}

			del_list[j].gid = -1;
		}
		else {
		}
	}
}

bool check_inv_same_del(struct gap_list *del_list, int num_del, int cur, int del_id)
{
	int i;
	bool res = false;

	for( i = 0; i < num_del; i++ )
	{
		if( (del_list[i].id1 == cur) || (del_list[i].id2 == cur) )
		{
			if( del_list[i].gid == del_id )
			{
				res = true;
			}
		}
	}

	return(res);
}

int assign_gap_list(struct gap_list *del_list, int loc, int last_gid, struct gap_list gp)
{
	int j = 0;
	bool is_end = false;
	struct I reg_temp, reg_comp;
	struct I temp, comp;
	int wid_x, wid_xt;
	int mx_wid;
	int res = last_gid;

	comp = assign_I(gp.x1, gp.x2);

	wid_x = width(comp);

	if( loc == 0 )
	{
		res = last_gid + 1;
		del_list[loc].gid = res;
	}
	else
	{
		while( (j < loc) && (is_end == false) )
		{
			temp = assign_I(del_list[j].x1, del_list[j].x2);

			wid_xt = width(temp);

			if( wid_x > wid_xt ) mx_wid = wid_x;
			else mx_wid = wid_xt;

			if( overlap(comp, temp) == true )
			{
				if( (comp.lower == temp.lower) || (comp.lower == temp.upper) || (comp.upper == temp.lower) || (comp.upper == temp.upper ) )
				{
					del_list[loc].gid = del_list[j].gid;
					is_end = true;
				}
				else {
					reg_comp = assign_I(comp.lower-DEL_TH, comp.upper+DEL_TH);
					reg_temp = assign_I(temp.lower-DEL_TH, temp.upper+DEL_TH);
					if( (subset(comp, reg_temp) == true) || (subset(temp, reg_comp) == true) ) 
					{
						del_list[loc].gid = del_list[j].gid;
						is_end = true;
					}
				}
			}
			j++;

			if( ( is_end == false ) && ( j == loc ) )
			{
				res = last_gid + 1;
				del_list[loc].gid = res;
			}
		}
	}

	del_list[loc].type = gp.type;
	del_list[loc].id1 = gp.id1;
	del_list[loc].id2 = gp.id2;
	del_list[loc].x1 = gp.x1;
	del_list[loc].x2 = gp.x2;
	del_list[loc].y1 = gp.y1;
	del_list[loc].y2 = gp.y2;

	return(res);
}

int find_opt_del(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct gap_list *gp)
{
	int i = 0;
	int max_id = -1;
	bool *is_x;
	int *sd;
	int d = 0;
	struct gap_list gps;
	int start = 0, mid1 = -1, mid2 = -1, end = 0;
	int temp_score = -1;
	int min_score = 100;
	int closeness = MDIS_THRESHOLD;
	int min_closeness = MDIS_THRESHOLD;
	int len_gap = MDIS_THRESHOLD;
	int min_gap = MDIS_THRESHOLD;

	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));

	*sd = 0;
	*is_x = true;

	if( w_sid < h_sid )
	{
		start = w_sid;
		if( w_fid < h_sid )
		{
			mid1 = w_fid;
			mid2 = h_sid;
			end = h_fid;
		}
		else 
		{
			if( w_fid < h_fid ) end = h_fid;
			else end = w_fid;
		}
	}
	else 
	{
		start = h_sid;
		if( h_fid < w_sid )
		{
			mid1 = h_fid;
			mid2 = w_sid;
			end = w_fid;
		}
		else
		{
			if( h_fid < w_fid ) end = w_fid;
			else end = h_fid;
		}
	}

	for( i = start; i <= end; i++ )
	{
		if( st[i].id == id ) {}	
		else if( dots[st[i].id].sp_id != dots[id].sp_id ) {}
		else if( dots[st[i].id].ctg_id1 != dots[id].ctg_id1 ) {}
		else if( dots[st[i].id].ctg_id2 != dots[id].ctg_id2 ) {}
		else if( abs(dots[st[i].id].identity - dots[id].identity) > PID_DIFF ) {}
		else if( (width(dots[st[i].id].x) >= EFFEC_DEL_LEN_TH) && (width(dots[id].x) >= EFFEC_DEL_LEN_TH) && (abs(dots[st[i].id].identity - dots[id].identity) >= DEL_ID_TH) ) {}
		else if( (dots[st[i].id].sign == 2) || (dots[st[i].id].sign == 20) || (dots[st[i].id].sign == 21) ) {}
		else if( dots[st[i].id].x.lower > dots[id].x.lower ) {}
		else if( (dots[st[i].id].sign == 0) && (dots[st[i].id].y.lower > dots[id].y.lower )) {}
		else if( (dots[st[i].id].sign == 1) && (dots[st[i].id].y.lower < dots[id].y.lower )) {}
		else if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
			if((dots[st[i].id].sign != 2) && (dots[st[i].id].sign == dots[id].sign) && ((d = distance(dots, st[i].id, id, is_x, sd)) <= MDIS_THRESHOLD))
			{
				if( (d <= 0) || (d >= MAX_DEL) ){}
				else
				{
					if( (*sd) <= DEL_TH )
					{
						if( check_candi(dots, id, st[i].id, CHECK_DEL) == false ) {}
						else
						{
							gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
							if((gps.type == -1) || (gps.type == 0) || (gps.type == 3))
							{
							}
							else
							{
								gps.gid = 0;
								closeness = *sd;
								len_gap = d;
								temp_score = abs(dots[id].identity - dots[st[i].id].identity);
								if( min_closeness > closeness )
								{
									min_closeness = closeness;
									min_score = temp_score;
									max_id = st[i].id;
									min_gap = d;
									*gp = gps;
								}
								else if( min_closeness == closeness )
								{
									if( min_gap > d )
									{
										min_closeness = closeness;
										min_score = temp_score;
										max_id = st[i].id;
										min_gap = d;
										*gp = gps;
									}
									else if( min_gap == d )
									{
										if( min_score > temp_score ) 
										{
											min_closeness = closeness;
											min_score = temp_score;
											max_id = st[i].id;
											min_gap = d;
											*gp = gps;
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

	free(sd);
	free(is_x); 

	if( max_id == -1 ) return(-1);
	else 
	{
		return(max_id);
	}
}

int remove_false_dels(struct DotList *dots, struct gap_list *del_list, int num_del, struct ID_List *dlist, int num_dup)
{
	int i = 0, j = 0, k = 0;
	struct I temp;
	struct I ins_reg;

	for( i = 0; i < num_del; i++ )
	{
		temp = assign_I(del_list[i].y1, del_list[i].y2);

		for( j = 0; j < num_dup; j++ )
		{
			if( dlist[j].is_x == true )
			{
				ins_reg = assign_I(dots[dlist[j].m_id].x.lower, dots[dlist[j].m_id].x.upper);
			}
			else
			{
				ins_reg = assign_I(dots[dlist[j].m_id].y.lower, dots[dlist[j].m_id].y.upper);
			}

			if( f_loose_overlap(temp, ins_reg, STRICT) == true )
			{
				for( k = 0; k < num_del; k++ )
				{
					if( k != i )
					{
						if( del_list[k].gid == del_list[i].gid ) del_list[k].gid = -1;
					}
				}
				del_list[i].gid = -1;
			}
		}
	}

	j = 0;
	for( i = 0; i < num_del; i++ )
	{
		if( del_list[i].gid != -1 )
		{
			del_list[j] = assign_glist(del_list[i]);
			j++;
		}
	}

	return(j);
}

int find_redo_del_list(struct I src, struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct gap_list *del_list)
{ 
	int num_del = 0;
	struct perm_pt *m_pts;
	struct kdnode *m_tree;
	int num_lines = *num;
	int i, j;
	int xval, yval;
	int w_sid, w_fid, h_sid, h_fid;
	int start, end;
	int cur_id, tmp_id;
	struct I cmp_reg;
	bool is_skip;
	int b, e;
	struct gap_list cur_del;
	int min_width;
	int min_id;
	struct I tmp_reg;
	int mid1, mid2;

	m_pts = (struct perm_pt *) ckalloc(sizeof(struct perm_pt) * num_lines);	
	m_tree = (struct kdnode *) ckalloc(sizeof(struct kdnode));

	assign_perm(m_pts, num_lines, dots, RIGHT);
	m_tree = build_kd(m_pts, 0, num_lines-1);

	cmp_reg = assign_I(src.lower-DEL_TH, src.lower+DEL_TH);

	w_sid = find_id(dots, m_tree, size, -1, src.lower, src.lower, W_SID);
	w_fid = find_id(dots, m_tree, size, -1, src.lower, src.lower, W_FID);
	h_sid = find_id(dots, m_tree, size, -1, src.lower, src.lower, H_SID);
	h_fid = find_id(dots, m_tree, size, -1, src.lower, src.lower, H_FID);

  if( w_sid < h_sid )
  {
    start = w_sid;
    if( w_fid < h_sid )
    {
      mid1 = w_fid;
      mid2 = h_sid;
      end = h_fid;
    }
    else
    {
      if( w_fid < h_fid ) end = h_fid;
      else end = w_fid;
    }
  }
  else
  {
    start = h_sid;
    if( h_fid < w_sid )
    {
      mid1 = h_fid;
      mid2 = w_sid;
      end = w_fid;
    }
    else
    {
      if( h_fid < w_fid ) end = w_fid;
      else end = h_fid;
    }
  }

	for( i = start; i <= end; i++ ) {
		cur_id = m_pts[i].id;
		is_skip = true;
		min_id = -1;
		min_width = -1;
		if( dots[cur_id].sign == 0 ) {
			xval = dots[cur_id].x.lower;
			yval = dots[cur_id].y.lower;
			if( in(dots[cur_id].x.lower, cmp_reg) == true ) {
				is_skip = false;
				b = find_id(dots, tree, size, cur_id, xval, yval, H_SID);
				e = find_id(dots, tree, size, cur_id, xval, yval, H_FID);
			}
			else if( in(dots[cur_id].y.lower, cmp_reg) == true ) {
				is_skip = false;
				b = find_id(dots, tree, size, cur_id, xval, yval, W_SID);
				e = find_id(dots, tree, size, cur_id, xval, yval, W_FID);
			}
		}
		else if( dots[cur_id].sign == 1 ) {
			xval = dots[cur_id].x.lower;
			yval = dots[cur_id].y.upper;
			if( in(dots[cur_id].x.lower, cmp_reg) == true ) {
				is_skip = false;
				b = find_id(dots, tree, size, cur_id, xval, yval, H_SID);
				e = find_id(dots, tree, size, cur_id, xval, yval, H_FID);
			}
			else if( in(dots[cur_id].y.upper, cmp_reg) == true ) {
				is_skip = false;
				b = find_id(dots, tree, size, cur_id, xval, yval, W_SID);
				e = find_id(dots, tree, size, cur_id, xval, yval, W_FID);
			}
		}

		if( is_skip == false ) {
			for( j = b; j <= e; j++ ) {
				tmp_id = p_pts[j].id;	
				if( (tmp_id != cur_id) && ( dots[cur_id].sign == dots[tmp_id].sign ) && (dots[cur_id].sign == 0) ) { 
					if( (in(dots[cur_id].x.lower, cmp_reg) == true) && (in(dots[tmp_id].x.upper, cmp_reg) == true) && (dots[cur_id].y.lower > dots[tmp_id].y.upper)) 
					{
						tmp_reg = assign_I(dots[tmp_id].y.upper, dots[cur_id].y.lower);
						if( (min_id == -1) || ((min_width != -1) && (min_width > width(tmp_reg)))) 
						{
							min_id = tmp_id;
							min_width = width(tmp_reg);
							cur_del.type = 1; // x region is deleted
							cur_del.y1 = tmp_reg.lower;
							cur_del.y2 = tmp_reg.upper;
							cur_del.id1 = tmp_id;
							cur_del.id2 = cur_id;

							if( dots[cur_id].x.lower < dots[tmp_id].x.upper) {
								cur_del.x1 = dots[cur_id].x.lower;
								cur_del.x2 = dots[tmp_id].x.upper;
							}
							else if( dots[cur_id].x.lower > dots[tmp_id].x.upper) {
								cur_del.x1 = dots[tmp_id].x.upper;
								cur_del.x2 = dots[cur_id].x.lower;
							}
							else {
								cur_del.x1 = dots[cur_id].x.lower;
								cur_del.x2 = cur_del.x1 + 1;
							}
						}
					}
					else if( (in(dots[cur_id].y.lower, cmp_reg) == true) && (in(dots[tmp_id].y.upper, cmp_reg) == true) && (dots[cur_id].x.lower > dots[tmp_id].x.upper)) 
					{
						tmp_reg = assign_I(dots[tmp_id].x.upper, dots[cur_id].x.lower);
						if( (min_id == -1) || ((min_width != -1) && (min_width > width(tmp_reg)))) 
						{
							min_id = tmp_id;
							min_width = width(tmp_reg);
							cur_del.type = 2; // y region is deleted
							cur_del.y1 = tmp_reg.lower;
							cur_del.y2 = tmp_reg.upper;
							cur_del.id1 = tmp_id;
							cur_del.id2 = cur_id;

							if( dots[cur_id].y.lower < dots[tmp_id].y.upper) {
								cur_del.x1 = dots[cur_id].y.lower;
								cur_del.x2 = dots[tmp_id].y.upper;
							}
							else if( dots[cur_id].y.lower > dots[tmp_id].y.upper) {
								cur_del.x1 = dots[tmp_id].y.upper;
								cur_del.x2 = dots[cur_id].y.lower;
							}
							else {
								cur_del.x1 = dots[cur_id].y.lower;
								cur_del.x2 = cur_del.x1 + 1;
							}
						}
					}
				}
				else if( (tmp_id != cur_id) && (dots[cur_id].sign == dots[tmp_id].sign ) && (dots[cur_id].sign == 1) ) 
				{ 
          if( (in(dots[cur_id].x.lower, cmp_reg) == true) && (in(dots[tmp_id].x.upper, cmp_reg) == true) && (dots[tmp_id].y.lower > dots[cur_id].y.upper)) 
					{              	
						tmp_reg = assign_I(dots[cur_id].y.upper, dots[tmp_id].y.lower);     		        
						if( (min_id == -1) || ((min_width != -1) && (min_width > width(tmp_reg))))
           	{
             	min_id = tmp_id;
             	min_width = width(tmp_reg);
             	cur_del.type = 1; // x region is deleted
             	cur_del.y1 = tmp_reg.lower;
             	cur_del.y2 = tmp_reg.upper;
             	cur_del.id1 = tmp_id;
             	cur_del.id2 = cur_id;

             	if( dots[cur_id].x.lower < dots[tmp_id].x.upper) {
               	cur_del.x1 = dots[cur_id].x.lower;
               	cur_del.x2 = dots[tmp_id].x.upper;
             	}
             	else if( dots[cur_id].x.lower > dots[tmp_id].x.upper) {
               	cur_del.x1 = dots[tmp_id].x.upper;
               	cur_del.x2 = dots[cur_id].x.lower;
             	}
             	else {
               	cur_del.x1 = dots[cur_id].x.lower;
               	cur_del.x2 = cur_del.x1 + 1;
             	}
            }
					}
          else if( (in(dots[cur_id].y.upper, cmp_reg) == true) && (in(dots[tmp_id].y.lower, cmp_reg) == true) && (dots[cur_id].x.lower > dots[tmp_id].x.upper))             
					{              
						tmp_reg = assign_I(dots[tmp_id].x.upper, dots[cur_id].x.lower);
						if( (min_id == -1) || ((min_width != -1) && (min_width > width(tmp_reg))))
						{
               min_id = tmp_id;
               min_width = width(tmp_reg);
               cur_del.type = 2; // y region is deleted
               cur_del.y1 = tmp_reg.lower;
               cur_del.y2 = tmp_reg.upper;
               cur_del.id1 = tmp_id;
               cur_del.id2 = cur_id;

               if( dots[cur_id].y.upper < dots[tmp_id].y.lower) {
                 cur_del.x1 = dots[cur_id].y.upper;
                 cur_del.x2 = dots[tmp_id].y.lower;
               }
               else if( dots[cur_id].y.upper > dots[tmp_id].y.lower) {
                 cur_del.x1 = dots[tmp_id].y.lower;
                 cur_del.x2 = dots[cur_id].y.upper;
               }
               else {
                 cur_del.x1 = dots[cur_id].y.upper;
                 cur_del.x2 = cur_del.x1 + 1;
               }
             }
					}
				}
			}
		}

		if( min_width != -1 ) {
			del_list[num_del] = assign_glist(cur_del);
			del_list[num_del].gid = 1;
			num_del++;
		}
	}

	if( num_del > 0 ) {
		check_overlapped_del(dots, del_list, num_del);
		remove_dif_del(del_list, num_del, 1, 0, dots); // del_list[i].gid is always 1 and percentage identity is not necessary to be considered in redoing steps, so 0 is passed.
	}

	j = 0;
	for( i = 0; i < num_del; i++ )
	{
		if( del_list[i].gid == -1 ){}
		else
		{
			del_list[j] = assign_glist(del_list[i]);
			j++;
		}
	}
	num_del = j;

	free(m_pts);
	free_kd(m_tree);
	return(num_del);
}
