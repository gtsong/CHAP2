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

int pred_dels(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct gap_list *del_list)
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

	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));
	gp = (struct gap_list *) ckalloc(sizeof(struct gap_list));
	cid = (int *) ckalloc(sizeof(int));
	x_ins = (bool *) ckalloc(sizeof(bool));
	sorted = (struct slist *) ckalloc(sizeof(struct slist) * (*num));

	num_lines = *num;

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

			if( temp_score != -1 )
			{
				opt_id = -1;
			}
			else if( check_whole_regions_inclusion(dots, num_lines, (*gp).id1, (*gp).id1, (*gp).id2, *x_ins) == true )
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
		check_overlapped_del(dots, num_lines, del_list, num_del);
		remove_dif_del(del_list, num_del, last_gid);
	}

	j = 0;
	for( i = 0; i < num_del; i++ )
	{
		if( del_list[i].gid == -1 ){}
		else
		{
			del_list[j] = del_list[i];
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

void remove_dif_del(struct gap_list *del_list, int num_del, int last_gid)
{
	int i, j;
	int min_width = -1;

	for( i = 1; i <= last_gid; i++ )
	{
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
			}
		}

		for( j = 0; j < num_del; j++ )
		{
			if( del_list[j].gid == i )
			{
				if( ((float)abs(del_list[j].y2 - del_list[j].y1) - (float)min_width) / ((float)min_width) <= 0.05 ) {}
				else 
				{
					del_list[j].gid = -1;
				}
			}
		}
	}
}

void check_overlapped_del(struct DotList *dots, int num, struct gap_list *del_list, int num_del)
{
	int i, j, k;
	bool is_res;
	struct I left, right; // left: the left side region in deleted point
								 // rigth: the right side region in deleted point
	struct gap_list gp;
	struct I temp_x, temp_y;

	for( j = 0; j < num_del; j++ )
	{
		gp = del_list[j];

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

		if( (dots[gp.id1].pair_self == PAIR) || (dots[gp.id2].pair_self == PAIR) ) 
		{
			is_res = true;
		}
		else if( abs(dots[gp.id1].identity - dots[gp.id2].identity) >= DEL_ID_TH ) 
		{
			is_res = true;
		}
		else if( (dots[gp.id1].lock == BOTH_LOCK) || (dots[gp.id1].lock == LEFT_LOCK) || (dots[gp.id2].lock == RIGHT_LOCK) || (dots[gp.id2].lock == BOTH_LOCK)) 
		{
			is_res = true;
		}
		else
		{
			is_res = false;
		}

		i = 0;
		while((is_res == false) && (i < num))
		{
			if( dots[i].sign == 2 ) {}
			else if( (i == gp.id1) || (i == gp.id2) ) {}
			else if( check_inv_same_del(del_list, num_del, i, gp.gid) == true ) {}
			else
			{
				if( in(left.lower, dots[i].x) == true )
				{
					temp_x = assign_I(left.lower, dots[i].x.upper);
				}
				else if( in(right.upper, dots[i].x) == true )
				{
					temp_x = assign_I(dots[i].x.lower, right.upper);
				}
				else temp_x = assign_I(dots[i].x.lower, dots[i].x.upper);

				if( in(left.lower, dots[i].y) == true )
				{
					temp_y = assign_I(left.lower, dots[i].y.upper);
				}
				else if( in(right.upper, dots[i].y) == true )
				{
					temp_y = assign_I(dots[i].y.lower, right.upper);
				}
				else temp_y = assign_I(dots[i].y.lower, dots[i].y.upper);
			
				if( (subset(temp_x, left) == true) || (subset(temp_y, left) == true) )
				{
				}
				else if( (proper_overlap(temp_x, left) == true) && (temp_x.upper > (left.upper + DEL_OVERLAP_TH))) 
				{				
					if(width(intersect(temp_x, left)) > DEL_OVERLAP_TH)
					{
						is_res = true;
					}
				}
				else if( (proper_overlap(temp_y, left) == true) && (temp_y.upper > (left.upper + DEL_OVERLAP_TH)))
				{
					if( width(intersect(temp_y, left)) > DEL_OVERLAP_TH)
					{
						is_res = true;
					}
				}
				else{}
				
				if( (subset(temp_x, right) == true) || (subset(temp_y, right) == true) )
				{
				}
				else if( (proper_overlap(temp_x, right) == true) && (temp_x.lower < (right.lower - DEL_OVERLAP_TH)))
				{
					if(width(intersect(temp_x, right)) > DEL_OVERLAP_TH)
					{
						is_res = true;
					}
				}
				else if( (proper_overlap(temp_y, right) == true) && (temp_y.lower < (right.lower - DEL_OVERLAP_TH)))
				{
					if(width(intersect(temp_y, right)) > DEL_OVERLAP_TH)
					{
						is_res = true;
					}
				}
			}
			i++;
		}

		if( is_res == true ) 
		{
			for( k = 0; k < num_del; k++ )
			{
				if( del_list[k].gid == del_list[j].gid ) del_list[k].gid = -1;
			}

			del_list[j].gid = -1;
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
	struct I temp_x, comp_x;
	struct I temp_y, comp_y;
	int wid_x, wid_y, wid_xt, wid_yt;
	int mx_wid, my_wid;
	int res = last_gid;

	comp_x = assign_I(gp.x1, gp.x2);
	comp_y = assign_I(gp.y1, gp.y2);

	wid_x = width(comp_x);
	wid_y = width(comp_y);

	if( loc == 0 )
	{
		res = last_gid + 1;
		del_list[loc].gid = res;
	}
	else
	{
		while( (j < loc) && (is_end == false) )
		{
			temp_x = assign_I(del_list[j].x1, del_list[j].x2);
			temp_y = assign_I(del_list[j].y1, del_list[j].y2);

			wid_xt = width(temp_x);
			wid_yt = width(temp_y);

			if( wid_x > wid_xt ) mx_wid = wid_x;
			else mx_wid = wid_xt;

			if( wid_y > wid_yt ) my_wid = wid_x;
			else my_wid = wid_yt;

			if( proper_overlap(comp_x, temp_x) == true )
			{
				if( (comp_x.lower == temp_x.lower) || (comp_x.lower == temp_x.upper) || (comp_x.upper == temp_x.lower) || (comp_x.upper == temp_x.upper ) )
				{
					del_list[loc].gid = del_list[j].gid;
					is_end = true;
				}
				else if( (((float)width(intersect(comp_x, temp_x)))/((float)mx_wid))	>= 0.95 )
				{
					del_list[loc].gid = del_list[j].gid;
					is_end = true;
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
	int i;
	int max_id = -1;
	bool *is_x;
	int *sd;
	int d;
	struct gap_list gps;
	int start, mid1 = -1, mid2 = -1, end;
	int temp_score = -1;
	int min_score = 100;
	int closeness;
	int min_closeness = MDIS_THRESHOLD;
	int len_gap;
	int min_gap = MDIS_THRESHOLD;

	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));

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

int remove_false_dels(struct DotList *dots, struct gap_list *del_list, int num_del, struct ID_List *dlist, int num_dup )
{
	int i, j;
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
				for( j = 0; j < num_del; j++ )
				{
					if( j != i )
					{
						if( del_list[j].gid == del_list[i].gid ) del_list[j].gid = -1;
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
			del_list[j].gid = del_list[i].gid;
			del_list[j].type = del_list[i].type;
			del_list[j].id1 = del_list[i].id1;
			del_list[j].id2 = del_list[i].id2;
			del_list[j].x1 = del_list[i].x1;
			del_list[j].x2 = del_list[i].x2;
			del_list[j].y1 = del_list[i].y1;
			del_list[j].y2 = del_list[i].y2;
			j++;
		}
	}

	return(j);
}
