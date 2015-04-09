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
#include "util_ops.h"
#include "map_algns.h"
#include "update_init_algns.h"

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

	initialize_slist(sorted, 0, (*num_ch_pair));
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

void iden_gene_loss(int num_init_algns, struct DotList *init_algns, int num_ops, struct ops_list *ops, FILE *f, struct ops_list **del_ops, int *num_del_ops, int *num_alloc)
{
	struct DotList *pair_algns;
	int num_pair_algns = 0;
	int i = 0, j = 0, cur_len = 0;
	struct I src, dst;
	int index = -1;
	int index1 = -1, index2 = -1;
	int side = TIE;
  struct slist *sorted;
	int cur_id = -1, cmp_id = -1;
	int num_del = 0;
	int xl_diff1 = 0, xr_diff1 = 0, xl_diff2 = 0, xr_diff2 = 0;
	int lo = 0, hi = 1; 

	src = assign_I(0, 1);
	dst = assign_I(0, 1);
	j = 0;
	num_del = *num_del_ops;
	for( i = 0; i < num_init_algns; i++ ) {
		if( (init_algns[i].sign != DELETED) && (init_algns[i].sp_id == PAIR) ) 
		{
			cur_len = (init_algns[i].x.upper - init_algns[i].xr_diff) - (init_algns[i].x.lower + init_algns[i].xl_diff);
			if( cur_len > 0 ) {
				j++;
			}
		}
	}

	num_pair_algns = j;
	pair_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_pair_algns);
	sorted = (struct slist *) ckalloc(sizeof(struct slist) * num_pair_algns);

	initialize_algns(pair_algns, 0, num_pair_algns);

	j = 0;
	for( i = 0; i < num_init_algns; i++ ) {
		if( (init_algns[i].sign != DELETED) && (init_algns[i].sp_id == PAIR) ) 
		{
			cur_len = (init_algns[i].x.upper - init_algns[i].xr_diff) - (init_algns[i].x.lower + init_algns[i].xl_diff);
			if( cur_len > 0 ) {
				assign_algn(pair_algns, j, init_algns[i]);
				pair_algns[j].x = assign_I(init_algns[i].x.lower+init_algns[i].xl_diff, init_algns[i].x.upper-init_algns[i].xr_diff);
				pair_algns[j].y = assign_I(init_algns[i].y.lower+init_algns[i].yl_diff, init_algns[i].y.upper-init_algns[i].yr_diff);
				j++;
			}
		}
	}

	initialize_slist(sorted, 0, num_pair_algns);
  sort_init_algns(sorted, pair_algns, num_pair_algns, SELF2);

	num_pair_algns = j;

	for( i = 0; i < num_pair_algns; i++ ) {
		cur_id = sorted[i].id;
		for( j = (i+1); j < num_pair_algns; j++ ) {
			cmp_id = sorted[j].id;
			lo = 0;
			hi = 1;
			xl_diff1 = -1;
			xr_diff1 = -1;
			xl_diff2 = -1;
			xr_diff2 = -1;
			side = TIE;
			if( i == j ) {}
			else if( (pair_algns[cur_id].sign == DELETED) || (pair_algns[cur_id].sign == DELETED)) {}
			else if( proper_overlap(pair_algns[cur_id].y, pair_algns[cmp_id].y) == true ) 
			{
				lo = 0;
				hi = 1;
				index1 = pair_algns[cur_id].index;
				index2 = pair_algns[cmp_id].index;
				xl_diff1 = init_algns[index1].xl_diff;
				xr_diff1 = init_algns[index1].xr_diff;
				xl_diff2 = init_algns[index2].xl_diff;
				xr_diff2 = init_algns[index2].xr_diff;

				if( (strict_almost_equal(pair_algns[cur_id].y, pair_algns[cmp_id].y) == true) || (strict_almost_equal(pair_algns[cmp_id].y, pair_algns[cur_id].y) == true) )
				{
					src = assign_I(pair_algns[cur_id].x.lower, pair_algns[cur_id].x.upper);
					dst = assign_I(pair_algns[cmp_id].y.lower, pair_algns[cmp_id].y.upper);

					side = which_side_in_ops_list(pair_algns[cur_id], pair_algns[cmp_id], num_ops, ops);
					if( side == LEFT_SIDE ) {
						index = pair_algns[cmp_id].index;
						init_algns[index].sign = DELETED;
					}
					else if( side == RIGHT_SIDE ) {
						index = pair_algns[cur_id].index;
						init_algns[index].sign = DELETED;
					}
					else {
						if( pair_algns[cur_id].identity >= pair_algns[cmp_id].identity ) 
						{
							index = pair_algns[cmp_id].index;
							init_algns[index].sign = DELETED;
						}
						else {
							index = pair_algns[cur_id].index;
							init_algns[index].sign = DELETED;
						}
					}
				}
				else if( f_loose_subset(pair_algns[cur_id].y, pair_algns[cmp_id].y, STRICT) == true ) 
				{
					src = assign_I(pair_algns[cur_id].x.lower, pair_algns[cur_id].x.upper);
					index = pair_algns[cur_id].index;
					init_algns[index].sign = DELETED;
				}
				else if( f_loose_subset(pair_algns[cmp_id].y, pair_algns[cur_id].y, STRICT) == true ) 
				{
					src = assign_I(pair_algns[cur_id].x.lower, pair_algns[cur_id].x.upper);
					index = pair_algns[cmp_id].index;
					init_algns[index].sign = DELETED;
				}
				else { // overlapping case
					side = which_side_in_ops_list(pair_algns[cur_id], pair_algns[cmp_id], num_ops, ops);
					if( side == LEFT_SIDE ) {
						cut_part_algn(init_algns, pair_algns[cur_id].index, pair_algns[cmp_id].index,  LEFT_SIDE, f);	// LEFT_SIDE should remain and the other side is removed
					} 
					else if( side == RIGHT_SIDE ) {
						cut_part_algn(init_algns, pair_algns[cur_id].index, pair_algns[cmp_id].index, RIGHT_SIDE, f);
					}
					else {
						cut_part_algn(init_algns, pair_algns[cur_id].index, pair_algns[cmp_id].index, TIE, f);
					}
				}					

				if( (xl_diff1 == -1) || (xr_diff1 == -1) || (xl_diff2 == -1) || (xr_diff2 == -1) ) {}
				else if( init_algns[index1].sign == DELETED ) {
					lo = init_algns[index1].x.lower + init_algns[index1].xl_diff;
					hi = init_algns[index1].x.upper - init_algns[index1].xr_diff;
				}
				else if( init_algns[index2].sign == DELETED ) {
					lo = init_algns[index2].x.lower + init_algns[index2].xl_diff;
					hi = init_algns[index2].x.upper - init_algns[index2].xr_diff;
				}
				else if( xl_diff1 != init_algns[index1].xl_diff ) {
					lo = init_algns[index1].x.lower + xl_diff1;
					hi = init_algns[index1].x.lower + init_algns[index1].xl_diff;
				}
				else if( xr_diff1 != init_algns[index1].xr_diff ) {
					lo = init_algns[index1].x.upper - init_algns[index1].xr_diff;
					hi = init_algns[index1].x.upper - xr_diff1;
				}
				else if( xl_diff2 != init_algns[index2].xl_diff ) {
					lo = init_algns[index2].x.lower + xl_diff2;
					hi = init_algns[index2].x.lower + init_algns[index2].xl_diff;
				}
				else if( xl_diff2 != init_algns[index2].xl_diff ) {
					lo = init_algns[index2].x.upper - xr_diff2;
					hi = init_algns[index2].x.upper - init_algns[index2].xr_diff;
				}

				if( ((lo == 0) && (hi == 1)) || ( (hi - lo) <= DEL_TH ) ) {}
				else {
					if( num_del >= (*num_alloc) ) {
						*del_ops = (struct ops_list *) ckrealloc(*del_ops, ((*num_alloc)+ALLOC_UNIT) * sizeof(struct ops_list));
						init_ops(*del_ops, (*num_alloc), (*num_alloc)+ALLOC_UNIT);
						*num_alloc = *num_alloc + ALLOC_UNIT;
					}
					(*del_ops)[num_del].dstStart = lo;
					(*del_ops)[num_del].dstEnd = hi;
					(*del_ops)[num_del].sign = 'l';
					(*del_ops)[num_del].sp_id = SELF1;
					num_del++;
				}
			}
		}
	}

	for( i = (*num_del_ops); i < num_del; i++ ) {
		src = assign_I((*del_ops)[i].dstStart, (*del_ops)[i].dstEnd);
		for( j = (i+1); j < num_del; j++ ) {
			dst = assign_I((*del_ops)[j].dstStart, (*del_ops)[j].dstEnd);
			if( f_loose_subset(src, dst, STRICT) == true ) {
				(*del_ops)[i].sp_id = -1;
			}
		}
	}

	j = 0;
	for( i = 0; i < num_del; i++ ) {
		if( (*del_ops)[i].sp_id != -1 ) {
			(*del_ops)[j].dstStart = (*del_ops)[i].dstStart;
			(*del_ops)[j].dstEnd = (*del_ops)[i].dstEnd;
			(*del_ops)[j].sign = 'l';
			(*del_ops)[j].sp_id = SELF1;
			j++;
		}
	}

	*num_del_ops = j;

	free(sorted);
	free(pair_algns);
}

int which_side_in_ops_list(struct DotList algn1, struct DotList algn2, int num_ops, struct ops_list *ops)
{
	int res = TIE;
	int i = 0;
	struct I src, dst;

	src = assign_I(0,1);
	dst = assign_I(0,1);
	i = 0;
	while( (i < num_ops) && (res == TIE) ) {
		if(ops[i].sp_id == SELF1)
		{
			src = assign_I(ops[i].srcStart, ops[i].srcEnd);
			dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( (strict_almost_equal(algn1.x, src) == true) || (f_loose_subset(src, algn1.x, STRICT) == true) ) {
				res = LEFT_SIDE;
			}
			else if( (strict_almost_equal(algn2.x, src) == true) || (f_loose_subset(src, algn2.x, STRICT) == true) ) {
				res = RIGHT_SIDE;
			}
		}
		i++;
	}

	i = 0;
	while( (i < num_ops) && (res == TIE) ) {
		if(ops[i].sp_id == SELF1)
		{
			src = assign_I(ops[i].srcStart, ops[i].srcEnd);
			dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( (strict_almost_equal(algn1.x, src) == true) || (f_loose_subset(src, algn1.x, STRICT) == true) || (strict_almost_equal(algn1.x, dst) == true) || (f_loose_subset(dst, algn1.x, STRICT) == true) ) {
				res = LEFT_SIDE;
			}
			else if( (strict_almost_equal(algn2.x, src) == true) || (f_loose_subset(src, algn2.x, STRICT) == true) || (strict_almost_equal(algn2.x, dst) == true) || (f_loose_subset(dst, algn2.x, STRICT) == true) ) {
				res = RIGHT_SIDE;
			}
		}
		i++;
	}

	i = 0;
	while( (i < num_ops) && (res == TIE) ) {
		if(ops[i].sp_id == SELF1)
		{
			src = assign_I(ops[i].srcStart, ops[i].srcEnd);
			dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( (proper_overlap(algn1.x, src) == true) && (width(intersect(algn1.x, src)) >= DEL_TH) ) {
				res = LEFT_SIDE;	
			}
			else if( (proper_overlap(algn2.x, src) == true) && (width(intersect(algn2.x, src)) >= DEL_TH) ) {
				res = RIGHT_SIDE;	
			}
		}
		i++;
	}

	i = 0;
	while( (i < num_ops) && (res == TIE) ) {
		if(ops[i].sp_id == SELF1)
		{
			src = assign_I(ops[i].srcStart, ops[i].srcEnd);
			dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( ((proper_overlap(algn1.x, src) == true) && (width(intersect(algn1.x, src)) >= DEL_TH)) || ((proper_overlap(algn1.x, dst) == true) && (width(intersect(algn1.x, dst)) >= DEL_TH )) ) {
				res = LEFT_SIDE;	
			}
			else if( ((proper_overlap(algn2.x, src) == true) && (width(intersect(algn2.x, src)) >= DEL_TH)) || ((proper_overlap(algn2.x, dst) == true) && (width(intersect(algn2.x, dst)) >= DEL_TH)) ) {
				res = RIGHT_SIDE;	
			}
		}
		i++;
	}

	return(res);
}

void update_init_algns_for_loss(struct ops_list *del_ops, int num_del_ops, int num_init_algns, struct DotList *init_algns, int sp_id, FILE *f)
{
	int i = 0, j = 0;
	struct I del_reg, cur_reg;

	del_reg = assign_I(-1,0);
	cur_reg = assign_I(-1,0);
	for( i = 0; i < num_del_ops; i++ ) {
		del_reg = assign_I(del_ops[i].dstStart, del_ops[i].dstEnd);	
		for( j = 0; j < num_init_algns; j++ ) {
			if( sp_id == SELF1 ) {
				cur_reg = assign_I(init_algns[j].x.lower+init_algns[j].xl_diff, init_algns[j].x.upper-init_algns[j].xr_diff);
			}
			else {
				cur_reg = assign_I(init_algns[j].y.lower+init_algns[j].yl_diff, init_algns[j].y.upper-init_algns[j].yr_diff);
			}

			if( (init_algns[j].sign != DELETED) && (init_algns[j].sp_id == PAIR) && (proper_overlap(del_reg, cur_reg) == true)) {
				if( f_loose_subset(cur_reg, del_reg, LOOSE) == true ) {
					init_algns[j].sign = DELETED;
				}
				else if(subset(del_reg, assign_I(cur_reg.lower+DEL_TH, cur_reg.upper-DEL_TH)) == true) {}
				else {
					if( sp_id == SELF1 ) {
						adjust_init_algn(init_algns, j, del_reg, true, f, DUP, true);
					}
					else if( sp_id == SELF2 ) {
						adjust_init_algn(init_algns, j, del_reg, false, f, DUP, false);
					}
				}
			}
		}
	}
}
