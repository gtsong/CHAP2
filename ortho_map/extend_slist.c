#include "main.h"
#include "regions.h"
#include "extend_slist.h"
#include "find_dup_copy.h"
#include "find_merging.h"
#include "deal_gaps.h"
#include "util_gen.h"
#include "check_copy.h"
#include "kd_op.h"
#include "read_algn.h"
#include "util.h"
#include "util_i.h"

void extend_slist(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct ID_List *slist, int num_dup, int res_id, bool is_x, FILE *fp, struct DotList *init_dots)
{
	int i = 0;
	int num_lines;
	bool *l_f_is_x, *r_f_is_x;
	int l_opt_id = -1, r_opt_id = -1;
	struct I ins_reg, temp;

	l_f_is_x = (bool *) ckalloc(sizeof(bool));
	r_f_is_x = (bool *) ckalloc(sizeof(bool));

	*l_f_is_x = true;
	*r_f_is_x = true;

	num_lines = *num;
	
	ins_reg = assign_I(0, 1);
	temp = assign_I(0, 1);
	
	if( is_x == true )
	{
		ins_reg = assign_I(dots[res_id].x.lower, dots[res_id].x.upper);
	}
	else if( is_x == false )
	{
		ins_reg = assign_I(dots[res_id].y.lower, dots[res_id].y.upper);
	}

	for( i = 0; i < num_dup; i++ )
	{
		l_opt_id = -1;
		r_opt_id = -1;
		if( slist[i].is_x == true ) 
		{
			temp = assign_I(dots[slist[i].m_id].x.lower, dots[slist[i].m_id].x.upper);
		}
		else
		{
			temp = assign_I(dots[slist[i].m_id].y.lower, dots[slist[i].m_id].y.upper);
		}

		if( strict_almost_equal(ins_reg, temp) == true ) {}
		else
		{
			l_opt_id =  find_alt_ins_id(slist[i].left_id, dots, tree, p_pts, res_id, is_x, l_f_is_x, size, fp, init_dots);

			if( l_opt_id != -1 )
			{
				if( check_whole_regions_inclusion(dots, num_lines, res_id, l_opt_id, slist[i].left_id, is_x) == true )
				{
					l_opt_id = -1;	
				}			
				else if( in_slist(slist[i].left_id, l_opt_id, i, num_dup, slist, dots) == true )
				{
					l_opt_id = -1;
				}
			}

			r_opt_id =  find_alt_ins_id(slist[i].right_id, dots, tree, p_pts, res_id, is_x, r_f_is_x, size, fp, init_dots);

			if( r_opt_id != -1 )
			{
				if( check_whole_regions_inclusion(dots, num_lines, res_id, r_opt_id, slist[i].left_id, is_x) == true )
				{
					r_opt_id = -1;	
				}			
				else if( in_slist(slist[i].right_id, r_opt_id, i, num_dup, slist, dots) == true )
				{
					r_opt_id = -1;
				}
			}
		}

		if( r_opt_id != -1 )
		{
			if( dots[slist[i].right_id].x.lower < dots[r_opt_id].x.lower )
			{
				slist[i].left_id = slist[i].right_id;
				slist[i].right_id = r_opt_id;
				slist[i].f_is_x = *r_f_is_x;
				slist[i].m_id = res_id;
				slist[i].is_x = is_x;
			}
			else if( dots[slist[i].right_id].x.lower > dots[r_opt_id].x.lower )
			{
				slist[i].left_id = r_opt_id;
				slist[i].f_is_x = *r_f_is_x;
				slist[i].m_id = res_id;
				slist[i].is_x = is_x;
			}
			else
			{	
				if( ((dots[r_opt_id].sign == 0) && ( dots[slist[i].right_id].y.lower < dots[r_opt_id].y.lower )) || ( (dots[r_opt_id].sign == 1) && ( dots[slist[i].right_id].y.lower > dots[r_opt_id].y.lower ) ) )
				{
					slist[i].left_id = slist[i].right_id;
					slist[i].right_id = r_opt_id;
					slist[i].f_is_x = *r_f_is_x;
					slist[i].m_id = res_id;
					slist[i].is_x = is_x;
				}
				else if( ((dots[r_opt_id].sign == 0) && ( dots[slist[i].right_id].y.lower > dots[r_opt_id].y.lower )) || ( (dots[r_opt_id].sign == 1) && ( dots[slist[i].right_id].y.lower < dots[r_opt_id].y.lower ) ) )
				{
					slist[i].left_id = r_opt_id;
					slist[i].f_is_x = *r_f_is_x;
					slist[i].m_id = res_id;
					slist[i].is_x = is_x;
				}	
				else {} 
			}
		}
		else if( l_opt_id != -1 )
		{
			if( dots[slist[i].left_id].x.lower < dots[l_opt_id].x.lower )
			{
				slist[i].right_id = l_opt_id;
				slist[i].f_is_x = *l_f_is_x;
				slist[i].m_id = res_id;
				slist[i].is_x = is_x;
			}	
			else if( dots[slist[i].left_id].x.lower > dots[l_opt_id].x.lower )
			{
				slist[i].right_id = slist[i].left_id;
				slist[i].left_id = l_opt_id;
				slist[i].f_is_x = *l_f_is_x;
				slist[i].m_id = res_id;
				slist[i].is_x = is_x;
			}
			else
			{	
				if( ((dots[l_opt_id].sign == 0) && ( dots[slist[i].left_id].y.lower < dots[l_opt_id].y.lower )) || ( (dots[l_opt_id].sign == 1) && ( dots[slist[i].left_id].y.lower > dots[l_opt_id].y.lower ) ) )
				{
					slist[i].right_id = l_opt_id;
					slist[i].f_is_x = *l_f_is_x;
					slist[i].m_id = res_id;
					slist[i].is_x = is_x;
				}
				else if( ((dots[l_opt_id].sign == 0) && ( dots[slist[i].left_id].y.lower > dots[l_opt_id].y.lower )) || ( (dots[l_opt_id].sign == 1) && ( dots[slist[i].left_id].y.lower < dots[l_opt_id].y.lower ) ) )
				{
					slist[i].right_id = slist[i].left_id;
					slist[i].left_id = l_opt_id;
					slist[i].f_is_x = *l_f_is_x;
					slist[i].m_id = res_id;
					slist[i].is_x = is_x;
				}	
				else {} 
			}
		}
	}

	free(l_f_is_x);
	free(r_f_is_x);
}

int find_gap_for_ins(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid , int h_sid, int h_fid, int res_id, bool res_is_x, bool *f_is_x, FILE *fp, struct DotList *init_dots)
{
	int i = 0;
	int min_score = 1000;
	int max_id = -1;
	bool *is_x;
	int *sd;
	int d = 0;
	struct gap_list gps;
	int temp_score = -1;
	int start = 0, mid1 = -1, mid2 = -1, end = 0;
	int len1 = 0, len2 = 0, len = 0;
	int op_len = 0, op_len_x = 0, op_len_y = 0;
	int m_th = 0;
	int y_old = 0, y_cur = 0;
	struct I temp;
	int closeness = 0;

	sd = (int *) ckalloc(sizeof(int));
	is_x = (bool *) ckalloc(sizeof(bool));

	*is_x = true;
	*sd = 0;

	gps.type = -1;
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
		if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
			if( st[i].id == id ) {}
			else if( dots[st[i].id].sign == 2) {}
			else if( dots[st[i].id].x.lower > dots[id].x.lower ) {}
			else if( (dots[st[i].id].sign == 0) && (dots[st[i].id].y.lower > dots[id].y.lower )) {}
			else if( (dots[st[i].id].sign == 1) && (dots[st[i].id].y.lower < dots[id].y.lower )) {}
      else if( subset(dots[st[i].id].m_x, dots[id].x) || subset(dots[st[i].id].m_y, dots[id].y) || subset(dots[id].m_x, dots[st[i].id].x) || subset(dots[id].m_y, dots[st[i].id].y) ) {}
			else 
			{
				if((dots[st[i].id].sign != 2) && (dots[st[i].id].sign == dots[id].sign) && (dots[st[i].id].sp_id == dots[id].sp_id) && ((d = distance(dots, st[i].id, id, is_x, sd)) <= MDIS_THRESHOLD))
				{
					if( d <= M_THRESHOLD ) 
					{
					}
					else if( d > M_THRESHOLD )
					{
						len1 = width(dots[st[i].id].x);
						len2 = width(dots[id].x);

						if( len1 > len2 ) 
						{
							len = len2;
						}
						else len = len1;

						if( (len1 >= LG_TH) && (len2 >= LG_TH)) m_th = L_M_TH;
						else m_th = M_TH;

						if(((*sd) <= len) || (((*sd) <= L_M_TH) && (len1 >= LG_TH) && (len2 >= LG_TH)) )
						{
							op_len = 0;
							op_len_x = 0;
							op_len_y = 0;

							if( proper_overlap(dots[st[i].id].x, dots[id].x) == true )
							{
								op_len_x = width(intersect(dots[st[i].id].x, dots[id].x));
								op_len = op_len_x;
							}

							if( proper_overlap(dots[st[i].id].y, dots[id].y) == true )
							{
								op_len_y = width(intersect(dots[st[i].id].y, dots[id].y));
								if( op_len < op_len_y )
								{
									op_len = op_len_y;
								}
							}

							if( ((*sd) > m_th) || (op_len > m_th) )
							{
								if( (strict_almost_equal(dots[st[i].id].x, dots[id].x) == true) || (strict_almost_equal(dots[st[i].id].y, dots[id].y) == true ) ) 
								{
									gps.type = -1;
								}
	              else if( ((closeness = compute_closeness(dots, st[i].id, id)) <= C_OP_TH) && ((*sd) > m_th) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true )) gps.type = -1;
 		            else if( ((*sd) > m_th) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true )) {
                	gps.type = -1;
              	}
								else if( (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true) )
              	{
                	temp = intersect(dots[st[i].id].x, dots[id].x);
	                if( dots[id].sign == 0 ) {
	                  y_cur = init_dots[dots[id].index].y.lower + init_dots[dots[id].index].yl_diff;
	 	                y_old = find_yloc_one_ch(init_dots, dots[st[i].id], fp, width(temp), NO_GAP_INC);
                  	if( y_old == -1 ) {
                    	y_cur = dots[id].y.lower;
                    	y_old = dots[st[i].id].y.upper - width(temp);
                  	}
                	}
                	else if( dots[id].sign == 1 ) {
                  	y_old = init_dots[dots[id].index].y.upper - init_dots[dots[id].index].yr_diff;
                  	y_cur = find_yloc_one_ch(init_dots, dots[st[i].id], fp, width(temp), NO_GAP_INC);
                  	if( y_cur == -1 ) {
                    	y_old = dots[id].y.upper;
                    	y_cur = dots[st[i].id].y.lower + width(temp);
                  	}
                	}

                	if( y_old >= y_cur ) {
										gps = define_gap_new_type(dots, st[i].id, id, false);
										if( gps.type == -1 ) gps = define_gap_new_type(dots, st[i].id, id, true);
									}
                	else {
										 gps = define_gap_new_type(dots, st[i].id, id, true);
										 if( gps.type == -1 ) gps = define_gap_new_type(dots, st[i].id, id, false);
									}
              	}
								else if(proper_overlap(dots[st[i].id].x, dots[id].x) == true ) 
								{
									gps = define_gap_new_type(dots, st[i].id, id, ((*is_x)));
								}
								else if( proper_overlap(dots[st[i].id].y, dots[id].y) == true ) 
								{
									gps = define_gap_new_type(dots, st[i].id, id, ((*is_x)));
								}
     		      	else if( ((closeness = compute_closeness(dots, st[i].id, id)) > C_OP_TH) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true ))
								{
	              	if( ((subset(dots[st[i].id].x, dots[id].x) == true) || (subset(dots[id].x, dots[st[i].id].x) == true)) || ((subset(dots[st[i].id].y, dots[id].y) == true) || (subset(dots[id].y, dots[st[i].id].y) == true)) ) {
                		gps.type = -1;
              		}
              		else {
                		temp = intersect(dots[id].x, dots[st[i].id].x);
		                if( dots[id].sign == 0 ) {
            		      y_cur = init_dots[dots[id].index].y.lower + init_dots[dots[id].index].yl_diff;
             		    	y_old = find_yloc_one_ch(init_dots, dots[st[i].id], fp, width(temp), NO_GAP_INC);            

	  	                if( y_old == -1 ) {
                    		y_cur = dots[id].y.lower;
                    		y_old = dots[st[i].id].y.upper - width(temp);
                  		}
                		}
                		else if( dots[id].sign == 1 ) {
                  		y_cur = find_yloc_one_ch(init_dots, dots[st[i].id], fp, width(temp), NO_GAP_INC);              
                  		y_old = init_dots[dots[id].index].y.upper - init_dots[dots[id].index].yr_diff;
                  		if( y_cur == -1 ) {
                    		y_old = dots[id].y.upper;
                    		y_cur = dots[st[i].id].y.lower + width(temp);
                  		}
                		}
                		else gps.type = -1;
                
										if( (dots[id].sign == 0) || (dots[id].sign == 1) ) {
		                  if( y_old >= y_cur ) gps = define_gap_new_type(dots, st[i].id, id, false);
   		              	else gps = define_gap_new_type(dots, st[i].id, id, true);      
										}
									}
								}
								else if( (check_candi(dots, id, st[i].id, CHECK_INS_DUP) == false) && ( (*sd) > TD_TH) ) 
								{
									gps.type = -1;
								}
								else
								{
              		if( (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true ) )              
									{
      						 	gps.type = -1;              
									}              
									else if( proper_overlap(dots[st[i].id].x, dots[id].x) == true)
									{                
										gps = define_gap_new_type(dots, st[i].id, id, true);              						}              
									else if( proper_overlap(dots[st[i].id].y, dots[id].y) == true )
									{                
										gps = define_gap_new_type(dots, st[i].id, id, false);      
									}
									else gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
								}
							}

							if((gps.type == -1) || (gps.type == 3) || (gps.type == 0))
							{
							}
							else
							{
								gps.gid = 0;
								temp_score = check_match_gap_aln(dots, gps, res_id, res_is_x);
								if( temp_score == -1 ) {
						      if( (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true) ) {
                		if( (gps.type == 21) || (gps.type == 22) ){
                  		if( gps.type == 21 ) {
                    		gps = define_gap_new_type(dots, st[i].id, id, false);
                  		}
                  		else if( gps.type == 22 ) {
                    		gps = define_gap_new_type(dots, st[i].id, id, true);
                  		}

                  		if((gps.type == -1) || (gps.type == 3) || (gps.type == 0)) {}
                  		else {
                    		gps.gid = 0;
												temp_score = check_match_gap_aln(dots, gps, res_id, res_is_x);
											}
										}
									}
								}

								if( temp_score != -1 )
								{
									if( min_score >= temp_score )
									{
										min_score = temp_score;
										max_id = st[i].id;
										
										if( (gps.type == 1) || (gps.type == 11) || (gps.type == 21) )
										{
											*f_is_x = false;
										}
										else if( (gps.type == 2) || (gps.type == 12) || (gps.type == 22) )
										{
											*f_is_x = true;
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

	free(is_x);
	free(sd);

	if( max_id == -1 ) return(-1);
	else 
	{
		return(max_id);
	}
}

bool in_slist(int l_id, int r_id, int cur_index, int num_dup, struct ID_List *slist, struct DotList *dots)
{
	int i = 0;
	struct I temp;
	bool res = false;

	temp = assign_I(0, 1);

	for( i = 0; i < num_dup; i++ )
	{
		if( i != cur_index )
		{
			if( slist[i].is_x == true ) temp = assign_I(dots[slist[i].m_id].x.lower, dots[slist[i].m_id].x.upper);
			else temp = assign_I(dots[slist[i].m_id].y.lower, dots[slist[i].m_id].y.upper);

			if( (slist[i].left_id == l_id) && (slist[i].right_id == r_id) )
			{
				res = true;
			}
			else if( (slist[i].left_id == r_id) && (slist[i].right_id == l_id) )
			{
				res = true;
			}
		}
	}

	return(res);
}

int find_alt_ins_id(int id, struct DotList *dots, struct kdnode *tree, struct perm_pt *p_pts, int res_id, bool is_x, bool *f_is_x, int size, FILE *fp, struct DotList *init_dots)
{
	int xval = 0, yval = 0;
	int opt_id = -1;
	int w_sid = 0, w_fid = 0, h_sid = 0, h_fid = 0;

	if( id == res_id ) {}
	else
	{
		xval = dots[id].x.lower;
		if( dots[id].sign == 0 )
		{
			yval = dots[id].y.lower;
		}
		else 
		{
			yval = dots[id].y.upper;
		}	
	
		w_sid = find_id(dots, tree, size, id, xval, yval, W_SID);
		w_fid = find_id(dots, tree, size, id, xval, yval, W_FID);
		h_sid = find_id(dots, tree, size, id, xval, yval, H_SID);
		h_fid = find_id(dots, tree, size, id, xval, yval, H_FID);
		opt_id = find_gap_for_ins(dots, id, p_pts, w_sid, w_fid, h_sid, h_fid, res_id, is_x, f_is_x, fp, init_dots);
	}

	return(opt_id);
}
