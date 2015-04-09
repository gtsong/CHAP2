#include "main.h"
#include "regions.h"
#include "find_gene_loss.h"
#include "util_gen.h"
#include "find_dup_copy.h"
#include "deal_gaps.h"
#include "find_merging.h"
#include "chain_pair_alg.h"
#include "util.h"
#include "util_i.h"

void check_gene_loss(int *num_list, struct DotList *dots, int sp_code, int rm_sp, int left_sp, int *cur_num, struct ops_list *ops)
{
  int i, j = 0, l = 0, k = 0;
	int loc;
  int num_pair = 0, num_rm_sp = 0, num_left_sp = 0;
  struct kdnode *tree;
  struct perm_pt *p_pts;
  struct kdnode *ch_tree;
  struct perm_pt *ch_p_pts;
	int size = 0;
  struct DotList *pair_alg;
  struct DotList *rm_sp_alg;
  struct DotList *left_sp_alg;
	int max_x = 0, max_y = 0;
  struct slist *sorted;
  int opt_id;
  int xval, yval;
  int w_sid, w_fid, h_sid, h_fid;
  struct gap_list *gp;
	int *num_ch_pair;
	int num_clist;
	struct chain_list *clist;

  for( i = 0; i < (*num_list); i++ )
  {
    if(dots[i].sp_id == sp_code)
    {
      num_pair++;
			if( max_x < dots[i].x.upper ) max_x = dots[i].x.upper;
			if( max_y < dots[i].y.upper ) max_y = dots[i].y.upper;
    }

    if( dots[i].sp_id == rm_sp )
    {
      num_rm_sp++;
    }

    if( dots[i].sp_id == left_sp )
    {
      num_left_sp++;
    }
  }

	pair_alg = (struct DotList *) ckalloc(num_pair * sizeof(struct DotList));
	sorted = (struct slist *) ckalloc(num_pair * sizeof(struct slist));
	clist = (struct chain_list *) ckalloc(num_pair * sizeof(struct chain_list));
	rm_sp_alg = (struct DotList *) ckalloc(num_rm_sp * sizeof(struct DotList));
	left_sp_alg = (struct DotList *) ckalloc(num_left_sp * sizeof(struct DotList));

  for( i = 0; i < (*num_list); i++ )
  {
    if( dots[i].sp_id == sp_code )
    {
			if( size < dots[i].y.upper ) size = dots[i].y.upper;

      assign_algn(pair_alg, j, dots[i]);
      j++;
    }

    if( dots[i].sp_id == rm_sp )
    {
      assign_algn(rm_sp_alg, l, dots[i]);
      l++;
    }

    if( dots[i].sp_id == left_sp )
    {
      assign_algn(left_sp_alg, k, dots[i]);
      k++;
    }
  }

  p_pts = (struct perm_pt *) ckalloc(num_pair * sizeof(struct perm_pt));
  ch_p_pts = (struct perm_pt *) ckalloc(num_pair * sizeof(struct perm_pt));
	num_ch_pair = (int *) ckalloc(sizeof(int));
  gp = (struct gap_list *) ckalloc(sizeof(struct gap_list));

	*num_ch_pair = num_pair;
  assign_perm(p_pts, num_pair, pair_alg, LEFT);
	tree = build_kd(p_pts, 0, num_pair);
	num_clist = chain_pair_alg(pair_alg, num_ch_pair, tree, p_pts, max_x+1, max_y+1, clist);

	assign_perm(ch_p_pts, *num_ch_pair, pair_alg, LEFT);
	ch_tree = build_kd(ch_p_pts, 0, *num_ch_pair);

  sort_list(sorted, pair_alg, *num_ch_pair);
  for( i = 0; i < (*num_ch_pair) ; i++ )
  {
    opt_id = -1;
    loc = sorted[i].id;

    if( pair_alg[loc].sign == 2 ) {}
    else
    {
      xval = pair_alg[loc].x.lower;
      if( pair_alg[loc].sign == 0 ) yval = pair_alg[loc].y.lower;
      else yval = pair_alg[loc].y.upper;

      w_sid = find_id(dots, tree, size, sorted[i].id, xval, yval, W_SID);
	    w_fid = find_id(dots, tree, size, sorted[i].id, xval, yval, W_FID);
	    h_sid = find_id(dots, tree, size, sorted[i].id, xval, yval, H_SID);
	    h_fid = find_id(dots, tree, size, sorted[i].id, xval, yval, H_FID);

			if( w_sid >= (*num_ch_pair) ) w_sid = (*num_ch_pair) - 1;
			if( w_fid >= (*num_ch_pair) ) w_fid = (*num_ch_pair) - 1;
			if( h_sid >= (*num_ch_pair) ) h_sid = (*num_ch_pair) - 1;
			if( h_fid >= (*num_ch_pair) ) h_fid = (*num_ch_pair) - 1;

	    opt_id = find_opt_gene_loss(pair_alg, loc, p_pts, w_sid, w_fid, h_sid, h_fid, gp, rm_sp_alg, num_rm_sp);

			if( opt_id != -1 )
			{
  			ops[*cur_num].sign = 'l';
  			ops[*cur_num].src_b = (*gp).y1;
  			ops[*cur_num].src_e = (*gp).y2;
  			ops[*cur_num].dst_b = 0;
  			ops[*cur_num].dst_e = 0;
  			ops[*cur_num].sp_id = sp_code;
				*cur_num = (*cur_num) + 1;
			}
	  }
  }

	free(clist);
	free(sorted);
	free(pair_alg); 
	free(rm_sp_alg);
	free(left_sp_alg);
	free(gp);
	free(p_pts);
	free_kd(tree);
	free(num_ch_pair);
  free(ch_p_pts);
  free_kd(ch_tree);
}

int find_opt_gene_loss(struct DotList *pair_alg, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct gap_list *gp, struct DotList *rm_sp_alg, int num_rm_sp)
{
	struct I temp;
	int num_algns;
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
	int len, len1, len2;
	int op_len, op_len_x, op_len_y;
	int m_th;

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
		else if( (pair_alg[st[i].id].sign == 2) || (pair_alg[st[i].id].sign == 20) || (pair_alg[st[i].id].sign == 21) ) {}
		else if( pair_alg[st[i].id].x.lower > pair_alg[id].x.lower ) {}
		else if( (pair_alg[st[i].id].sign == 0) && (pair_alg[st[i].id].y.lower > pair_alg[id].y.lower )) {}
		else if( (pair_alg[st[i].id].sign == 1) && (pair_alg[st[i].id].y.lower < pair_alg[id].y.lower )) {}
		else if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
			if((pair_alg[st[i].id].sign != 2) && (pair_alg[st[i].id].sign == pair_alg[id].sign) && ((d = distance(pair_alg, st[i].id, id, is_x, sd)) <= MDIS_THRESHOLD))
			{
				if( (d <= 0) || (d >= MAX_DEL) ) {}
				else
				{
					len1 = width(pair_alg[st[i].id].x);
					len2 = width(pair_alg[id].x);

					if( len1 > len2 )
					{
						len = len2;
					}
					else len = len1;

					if( (len1 >= LG_TH) && (len2 >= LG_TH) ) m_th = L_M_TH;
					else m_th = M_TH;

					if( ((*sd) <= len) || ( ((*sd) <= L_M_TH) && (len1 >= LG_TH) && (len2 >= LG_TH)) )
					{
						op_len = 0;
						op_len_x = 0;
						op_len_y = 0;

						if( proper_overlap(pair_alg[st[i].id].x, pair_alg[id].x) == true )
						{
							op_len_x = width(intersect(pair_alg[st[i].id].x, pair_alg[id].x));
							op_len = op_len_x;
						}

						if( proper_overlap(pair_alg[st[i].id].y, pair_alg[id].y) == true )
						{
							op_len_y = width(intersect(pair_alg[st[i].id].y, pair_alg[id].y));
							if( op_len < op_len_y )
							{
								op_len = op_len_y;
							}
						}

						if( ((*sd) > m_th) || (op_len > m_th) )
						{
							if( (strict_almost_equal(pair_alg[st[i].id].x, pair_alg[id].x) == true) || (strict_almost_equal(pair_alg[st[i].id].y, pair_alg[id].y) == true ))
							{
								gps.type = -1;
							}
							else if( ((*is_x) == true) && (proper_overlap(pair_alg[st[i].id].x, pair_alg[id].x) == true )) 
							{
								gps = define_gap_new_type(pair_alg, st[i].id, id, (*is_x));
							}
							else if( ((*is_x) == false) && (proper_overlap(pair_alg[st[i].id].y, pair_alg[id].y) == true ) )
							{
								gps = define_gap_new_type(pair_alg, st[i].id, id, (*is_x));
							}
					    else if( (op_len_x >= op_len_y) && ((strict_subset(pair_alg[st[i].id].x, pair_alg[id].x) == true) || (strict_subset(pair_alg[id].x, pair_alg[st[i].id].x) == true)) )
			        {
			          gps = define_gap_new_type_inc(pair_alg, st[i].id, id, true);
							}
              else if( (op_len_x <= op_len_y) && ((strict_subset(pair_alg[st[i].id].y, pair_alg[id].y) == true) || (strict_subset(pair_alg[id].y, pair_alg[st[i].id].y) == true)) )
							{
								gps = define_gap_new_type_inc(pair_alg, st[i].id, id, false);
							}
							else
							{
								gps.type = -1;
							}	
						}
						else if( (check_candi(pair_alg, id, st[i].id, CHECK_INS_DUP) == false) && ( (*sd) > TD_TH) )
						{
							gps.type = -1;
						}
						else
						{
							gps = define_gap(pair_alg, st[i].id, id, d, *sd, *is_x);
						}
						
						if( (gps.type == 2) || (gps.type == 12) || (gps.type == 22) ) 
						{
							gps = redefine_for_del(pair_alg, gps);
							temp = assign_I(gps.x1, gps.x2);
							num_algns = check_inc_algns( temp, rm_sp_alg, num_rm_sp );

							if( num_algns <= 0 ) gps.type = -1;
						}
						
						if( (gps.type == -1) || (gps.type == 3) || (gps.type == 0) ) {}
						else if( (gps.type == 1) || (gps.type == 11) || (gps.type == 21) ) {}
						else
						{
							gps.gid = 0;
							closeness = *sd;
							len_gap = d;
							temp_score = abs(pair_alg[id].identity - pair_alg[st[i].id].identity);
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

	free(sd);
	free(is_x); 

	if( max_id == -1 ) return(-1);
	else 
	{
		return(max_id);
	}
}

int check_inc_algns(struct I reg, struct DotList *dots, int num)
{
	int i;
	int count = 0;

	for( i = 0; i < num; i++ )
	{
		if( dots[i].pair_self == PAIR ) {}
		else if( (loose_subset(dots[i].x, reg) == true ) || (loose_subset(dots[i].y, reg) == true ) ) count++;
		else if( ( proper_overlap(dots[i].x, reg) == true ) && ( width(intersect(dots[i].x, reg)) >= RP_BD ) ) count++;
		else if( ( proper_overlap(dots[i].y, reg) == true ) && ( width(intersect(dots[i].y, reg)) >= RP_BD ) ) count++;
	}

	return(count);
}
