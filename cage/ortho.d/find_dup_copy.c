#include "main.h"
#include "regions.h"
#include "find_dup_copy.h"
#include "find_merging.h"
#include "deal_gaps.h"
#include "util_gen.h"
#include "util_i.h"
#include "ins_dup_copy.h"
#include "check_copy.h"
#include "kd_op.h"
#include "read_algn.h"
#include "util.h"

extern int debug_mode;

int find_dup_copy(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct ID_List *dlist, FILE *fp, struct DotList *init_dots)
{
	struct slist *sorted;
	int i = 0;
	int num_lines;
	bool *x_ins;
	bool *f_is_x;
	int *t_ins;
	bool *is_t_ins;
	int xval, yval;
	int w_sid, w_fid, h_sid, h_fid;
	int opt_id;
	int *cid;
	int num_ins_copy = 0;
	int new_id;

	sorted = (struct slist *) ckalloc(sizeof(struct slist) * (*num));
	x_ins = (bool *) ckalloc(sizeof(bool));
	cid = (int *) ckalloc(sizeof(int));
	f_is_x = (bool *) ckalloc(sizeof(bool));
	t_ins = (int *) ckalloc(sizeof(int));
	is_t_ins = (bool *) ckalloc(sizeof(bool));

	num_lines = *num;
	sort_list(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		opt_id = -1;
		if( dots[sorted[i].id].sign == 2) {}
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

			opt_id = find_opt_du_copy(dots, num_lines, sorted[i].id, p_pts, tree, size, w_sid, w_fid, h_sid, h_fid, cid, x_ins, f_is_x, t_ins, fp, init_dots);
		}

		if( opt_id != -1 ) 
		{
			if( ((*t_ins) == 11) || ((*t_ins) == 12) || ((*t_ins) == 21) || ((*t_ins) == 22) )
			{
				*is_t_ins = true;	
			}
			else *is_t_ins = false;

			dlist[num_ins_copy].is_x = *x_ins;
			dlist[num_ins_copy].m_id = *cid;
			dlist[num_ins_copy].left_id = opt_id;
			dlist[num_ins_copy].right_id = sorted[i].id;
			dlist[num_ins_copy].f_is_x = *f_is_x;
			dlist[num_ins_copy].is_t_ins = *is_t_ins;

			if( (*is_t_ins) == false ) 
			{
				new_id = check_inclusion_close_dup(dlist[num_ins_copy].m_id, dots, num_lines, x_ins, is_t_ins);
				if( new_id != -1 )
				{
					dlist[num_ins_copy].m_id = new_id;
					dlist[num_ins_copy].is_x = *x_ins;
					dlist[num_ins_copy].is_t_ins = *is_t_ins;
				}
			}

			if( dots[dlist[num_ins_copy].left_id].lock == -1 )
			{
				dots[dlist[num_ins_copy].left_id].lock = RIGHT_LOCK;
			}
			else if( dots[dlist[num_ins_copy].left_id].lock == LEFT_LOCK )
			{
				dots[dlist[num_ins_copy].left_id].lock = BOTH_LOCK;
			}

			if( dots[dlist[num_ins_copy].right_id].lock == -1 )
			{
				dots[dlist[num_ins_copy].right_id].lock = LEFT_LOCK;
			}
			else if( dots[dlist[num_ins_copy].right_id].lock == RIGHT_LOCK )
			{
				dots[dlist[num_ins_copy].right_id].lock = BOTH_LOCK;
			}

			num_ins_copy++;

			if( ((*t_ins) == 21) || ((*t_ins) == 22 ) )
			{
				if( is_tandem(dots[*cid]) == true ) {
					if( (*x_ins) == true )
					{
						dlist[num_ins_copy].is_x = false;
					}
					else dlist[num_ins_copy].is_x = true;

					dlist[num_ins_copy].m_id = dlist[num_ins_copy-1].m_id;
					dlist[num_ins_copy].left_id = dlist[num_ins_copy-1].left_id;
					dlist[num_ins_copy].right_id = dlist[num_ins_copy-1].right_id;
					dlist[num_ins_copy].f_is_x = dlist[num_ins_copy-1].f_is_x;
					dlist[num_ins_copy].is_t_ins = dlist[num_ins_copy-1].is_t_ins;
					num_ins_copy++;

					if( debug_mode == TRUE ) printf("Double ==> ");
				}
			}

			if( debug_mode == TRUE ) {
				printf("Insertion of a duplicated copy: %d-%d %d-%d, %d\n" , dots[dlist[num_ins_copy-1].m_id].x.lower, dots[dlist[num_ins_copy-1].m_id].x.upper, dots[dlist[num_ins_copy-1].m_id].y.lower, dots[dlist[num_ins_copy-1].m_id].y.upper, dlist[num_ins_copy-1].m_id);
				printf("%d-%d %d-%d\n", dots[opt_id].x.lower, dots[opt_id].x.upper, dots[opt_id].y.lower, dots[opt_id].y.upper);
				printf("%d-%d %d-%d, %d, %d, %d\n", dots[sorted[i].id].x.lower, dots[sorted[i].id].x.upper, dots[sorted[i].id].y.lower, dots[sorted[i].id].y.upper, dots[*cid].identity, dots[opt_id].identity, dots[sorted[i].id].identity);
			}
	
			if( (debug_mode == TRUE) && ( dlist[num_ins_copy-1].is_t_ins == true) )
			{
				printf("Tandem Insertion\n");
			}
		}	
	}

	free(sorted);
	free(is_t_ins);
	free(t_ins);
	free(f_is_x);
	free(cid);
	free(x_ins);

	return(num_ins_copy);
}

int find_opt_du_copy(struct DotList *dots, int num_lines, int id, struct perm_pt *st, struct kdnode *tree, int size, int w_sid, int w_fid, int h_sid, int h_fid, int *cid, bool *x_ins, bool *f_is_x, int *t_ins, FILE *fp, struct DotList *init_dots)
{
	int i;
	int min_score = 1000;
	int max_id = -1;
	bool *is_x;
	int *sd;
	int d;
	struct gap_list gps;
	int temp_score;
	int y_cur, y_old;
	struct I temp;
	int closeness;
	int start, mid1 = -1, mid2 = -1, end;
	int len1, len2, len;
	int opt_cid;
	int op_len = 0, op_len_x, op_len_y;
	int m_th;
	int from = 0, to = 1;

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
		if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
			if( st[i].id == id ) {}
			else if( dots[st[i].id].sign == 2) {}
			else if( dots[st[i].id].x.lower > dots[id].x.lower ) {}
			else if( dots[st[i].id].ctg_id1 != dots[id].ctg_id1 ) {}
			else if( dots[st[i].id].ctg_id2 != dots[id].ctg_id2 ) {}
			else if( (dots[st[i].id].sign == 0) && (dots[st[i].id].y.lower > dots[id].y.lower )) {}
			else if( (dots[st[i].id].sign == 1) && (dots[st[i].id].y.lower < dots[id].y.lower )) {}
			else if( subset(dots[st[i].id].m_x, dots[id].x) || subset(dots[st[i].id].m_y, dots[id].y) || subset(dots[id].m_x, dots[st[i].id].x) || subset(dots[id].m_y, dots[st[i].id].y) ) {}
			else if( ((dots[st[i].id].pair_self == SELF) && (is_tandem(dots[st[i].id]) == true)) && ((dots[id].pair_self == SELF) && (is_tandem(dots[id]) == true))) {}
			else 
			{
				if((dots[st[i].id].sign != 2) && (dots[st[i].id].sign == dots[id].sign) && (dots[st[i].id].sp_id == dots[id].sp_id) && ((d = distance(dots, st[i].id, id, is_x, sd)) <= MDIS_THRESHOLD))
				{
					len1 = width(dots[st[i].id].x);
					len2 = width(dots[id].x);

					if( len1 > len2 ) len = len2;
					else len = len1;

					if( (len1 >= LG_TH) && (len2 >= LG_TH)) m_th = L_M_TH;
					else m_th = M_TH;

					if((*sd) <= len) 
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
							if( (strict_almost_equal(dots[st[i].id].x, dots[id].x) == true) || (strict_almost_equal(dots[st[i].id].y, dots[id].y) == true ) ) gps.type = -1;
//            	else if( ((*sd) > m_th) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true )) {
            	else if( (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true ) && (tandem_exist(dots, st, tree, size, st[i].id, id) == false)) { 
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
							else if( proper_overlap(dots[st[i].id].x, dots[id].x) == true ) 
							{
								gps = define_gap_new_type(dots, st[i].id, id, true);
							}
							else if( proper_overlap(dots[st[i].id].y, dots[id].y) == true ) 
							{
								gps = define_gap_new_type(dots, st[i].id, id, false);
							}
							else gps.type = -1;
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
 							
									if( debug_mode == TRUE ) printf("Gap: %d-%d, %d-%d\n", dots[st[i].id].x.lower, dots[st[i].id].x.upper, dots[st[i].id].y.lower, dots[st[i].id].y.upper);
								}
							}
						}          
						else if( (check_candi(dots, id, st[i].id, CHECK_INS_DUP) == false) && ( (*sd) > TD_TH)) gps.type = -1;          
						else          
						{            
							if( (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true ) ) 
							{              
								gps.type = -1;	
							}            
							else if( proper_overlap(dots[st[i].id].x, dots[id].x) == true )  
							{              
								gps = define_gap_new_type(dots, st[i].id, id, true);            
							}            
							else if( proper_overlap(dots[st[i].id].y, dots[id].y) == true )            
							{              
								gps = define_gap_new_type(dots, st[i].id, id, false);            
							}
							else gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
						}
					}
					else gps.type = -1;

					if((gps.type == -1) || (gps.type == 3) || (gps.type == 0))
					{
					}
					else
					{
						gps.gid = 0;

						temp_score = get_score_copy(dots, num_lines, gps, cid, x_ins);
						if( temp_score != -1 )
						{
							if( check_whole_regions_inclusion(dots, num_lines, *cid, st[i].id, id, *x_ins) == true )
							{
								temp_score = -1;	
							}			
						}

						if( (temp_score == -1) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true) && ((gps.type == 21) || (gps.type == 22) ) ) {
							from = gps.y1;
							to = gps.y2;
							gps.y1 = from - abs(to - from);
							gps.y2 = from;
							temp_score = get_score_copy(dots, num_lines, gps, cid, x_ins);
							if( temp_score != -1 )
							{
								if( check_whole_regions_inclusion(dots, num_lines, *cid, st[i].id, id, *x_ins) == true )
								{
									temp_score = -1;	
								}			
							}

							if( temp_score == -1 ) {
								gps.y1 = from - abs(to-from)/2;
								gps.y2 = from + abs(to-from)/2;
								temp_score = get_score_copy(dots, num_lines, gps, cid, x_ins);
								if( temp_score != -1 )
								{
									if( check_whole_regions_inclusion(dots, num_lines, *cid, st[i].id, id, *x_ins) == true )
									{
										temp_score = -1;	
									}			
								}
							}
						}

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
										temp_score = get_score_copy(dots, num_lines, gps, cid, x_ins);
										if( temp_score != -1 )
										{
											if( check_whole_regions_inclusion(dots, num_lines, *cid, st[i].id, id, *x_ins) == true )
											{
												temp_score = -1;	
											}			
										}
									}
								}
							}
						}

						if( temp_score != -1 )
						{
							if( min_score > temp_score )
							{
								min_score = temp_score;
								max_id = st[i].id;
								opt_cid = *cid;
								
								if( gps.type == 1 )
								{
									*f_is_x = false;
								}
								else if( gps.type == 2 )
								{
									*f_is_x = true;
								}

								if( (gps.type == 11) || (gps.type == 21) )
								{
									*f_is_x = false;
									*t_ins = gps.type;
								}
								else if( (gps.type == 12) || (gps.type == 22) )
								{
									*f_is_x = true;
									*t_ins = gps.type;
								}
								else
								{
									*t_ins = -1;
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
		*cid = opt_cid;
		return(max_id);
	}
}

int find_id(struct DotList *dots, struct kdnode *tree, int size, int id, int xval, int yval, int option)
{
	int res_id = 0;
	int x = 0, y = 1;
	int len = 0;
	int sign = -1;

	if( id == -1 ) { // to redo 'del' events beforing ancestral steps
		len = RANGE_TH;
		sign = -1;
	}
	else {
		if( width(dots[id].x) > width(dots[id].y) ) len = width(dots[id].x);
		else len = width(dots[id].y);

		sign = dots[id].sign;
	}

	len = len + 100;

	if( option == W_SID )
	{
		x = xval - MDIS_THRESHOLD;
		if( len > RANGE_TH )
		{
			y = yval - len;
		}
		else y = yval - RANGE_TH;

		if( x <= 0 ) x = 1;
		if( y <= 0 ) y = 1;
		res_id = find_pred_blk(tree, x, y);
	}
	else if( option == W_FID )
	{
		if( len > RANGE_TH )
		{
			x = xval + len;
			y = yval + len;
		}
		else 
		{
			x = xval + RANGE_TH;
			y = yval + RANGE_TH;
		}
		if( x >= size ) x = size;
		if( y >= size ) y = size;
		res_id = find_successor(tree, x, y);
	}
	else if( option == H_SID )
	{
		if( len > RANGE_TH )
		{
			x = xval - len;
		}
		else x = xval - RANGE_TH;

		if( sign == 0 ) y = yval - MDIS_THRESHOLD;
		else if( (sign == 1) || (sign == -1) ) 
		{	
			if( len > RANGE_TH )
			{
				y = yval - len;
			}
			else y = yval - RANGE_TH;
		}
		if( x <= 0 ) x = 1;
		if( y <= 0 ) y = 1;
		res_id = find_pred_blk(tree, x, y);
	}
	else if( option == H_FID )
	{
		if( len > RANGE_TH )
		{
			x = xval + len;
		}
		else x = xval + RANGE_TH;

		if( sign == 0 ) 
		{
			if( len > RANGE_TH )
			{
				y = yval + len;
			}
			else y = yval + RANGE_TH;
		}
		else if( (sign == 1) || (sign == -1) ) y = yval + MDIS_THRESHOLD;
		if( x >= size ) x = size;
		if( y >= size ) y = size;
		res_id = find_successor(tree, x, y);
	}
	else if( option == OW_SID )
	{
		x = xval - (2 * RANGE_TH);
		if( sign == 0)	y = xval - (2 * RANGE_TH);
		else if( sign == 1 ) y = yval - MDIS_THRESHOLD;
		if( x <= 0 ) x = 1;
		if( y <= 0 ) y = 1;
		res_id = find_pred_blk(tree, x, y);
	}
//	else if( option == OH_FID )
	else // option == OH_FID
	{
		x = xval + (2 * RANGE_TH);
		y = yval + (2 * RANGE_TH);
		if( x >= size ) x = size;
		if( y >= size ) y = size;
		res_id = find_successor(tree, x, y);
	}

	return(res_id);
}

int check_inclusion_close_dup(int id, struct DotList *dots, int num_lines, bool *x_ins, bool *t_ins)
{
	int res = -1;
	int i = 0;
	int temp_res = id;
	struct I temp;

	if( (*x_ins) == true ) 
	{
		temp = assign_I(dots[id].x.lower, dots[id].x.upper);
	}
	else temp = assign_I(dots[id].y.lower, dots[id].y.upper);

	while( ( i < num_lines) && (res == -1) )
	{
		if( dots[i].pair_self == PAIR ) 
		{
		}
		else if( dots[i].sign == 2 ) {}
		else if( (i != id) && (dots[i].sign == 0) && ((dots[i].y.lower - dots[i].x.upper) <= THRESHOLD) )
		{
			if( strict_almost_equal(temp, dots[i].x) == true)
			{
				temp_res = i;
				res = i;
				*x_ins = true;
			}

			if( strict_almost_equal(temp, dots[i].y) == true)
			{
				temp_res = i;
				res = i;
				*x_ins = false;
			}
		}
		i++;
	}
	
	if( res != -1 )
	{
		if( is_tandem(dots[res]) == true )
		{
			*t_ins = true;
		}
		else 
		{
			*t_ins = false;
		}
	}
	else
	{
		i = 0;
		while( (i < num_lines) && (res == -1) )
		{
			if( dots[i].pair_self == PAIR ) {}
			else if( dots[i].sign == 2 ) {}
			else if( dots[temp_res].pair_self == PAIR )
			{
				if( strict_almost_equal(temp, dots[i].x) == true )
				{
					res = i;
					*x_ins = true;
				}
				else if( strict_almost_equal(temp, dots[i].y) == true )
				{
					res = i;
					*x_ins = false;
				}
			}
			i++;
		}

		if(res != -1) {
			if( is_tandem(dots[res]) == true ) *t_ins = true;
			else *t_ins = false;
		}
	}
	return( res);
}

bool check_whole_regions_inclusion(struct DotList *dots, int num_lines, int mid, int left_id, int right_id, bool is_x)			
{
	int i;
	bool res = false;
	struct I temp;

	if( is_x == true ) temp = assign_I(dots[mid].x.lower, dots[mid].x.upper);
	else temp = assign_I(dots[mid].y.lower, dots[mid].y.upper);

	for( i = 0; i < num_lines; i++ )
	{
		if( is_x == true )
		{
			if( strict_overlap(dots[left_id].y, dots[right_id].y, (M_TH/2)+1) == true )
			{
				if( (loose_subset(dots[mid].x, dots[i].x) == true) && (loose_subset(dots[left_id].x, dots[i].x) == true) && (loose_subset(dots[right_id].x, dots[i].x) == true))
				{
					res = true;	
				}
				else if( (loose_subset(dots[mid].x, dots[i].y) == true) && (loose_subset(dots[left_id].x, dots[i].y) == true) && (loose_subset(dots[right_id].x, dots[i].y) == true))
				{
					res = true;
				}
			}
			else if( strict_overlap(dots[left_id].x, dots[right_id].x, (M_TH/2)+1) == true )
			{
				if( (loose_subset(dots[mid].x, dots[i].x) == true) && (loose_subset(dots[left_id].y, dots[i].x) == true) && (loose_subset(dots[right_id].y, dots[i].x) == true))
				{
					res = true;	
				}
				else if( (loose_subset(dots[mid].x, dots[i].y) == true) && (loose_subset(dots[left_id].y, dots[i].y) == true) && (loose_subset(dots[right_id].y, dots[i].y) == true))
				{
					res = true;
				}
			}
			else
			{
			}
		}
		else
		{
			if( strict_overlap(dots[left_id].x, dots[right_id].x, (M_TH/2)+1) == true )
			{
				if( (loose_subset(dots[mid].y, dots[i].x) == true) && (loose_subset(dots[left_id].y, dots[i].x) == true) && (loose_subset(dots[right_id].y, dots[i].x) == true))
				{
					res = true;	
				}
				else if( (loose_subset(dots[mid].y, dots[i].y) == true) && (loose_subset(dots[left_id].y, dots[i].y) == true) && (loose_subset(dots[right_id].y, dots[i].y) == true))
				{
					res = true;
				}
			}
			else if( strict_overlap(dots[left_id].y, dots[right_id].y, (M_TH/2)+1) == true )
			{
				if( (loose_subset(dots[mid].y, dots[i].x) == true) && (loose_subset(dots[left_id].x, dots[i].x) == true) && (loose_subset(dots[right_id].x, dots[i].x) == true))
				{
					res = true;	
				}
				else if( (loose_subset(dots[mid].y, dots[i].y) == true) && (loose_subset(dots[left_id].x, dots[i].y) == true) && (loose_subset(dots[right_id].x, dots[i].y) == true))
				{
					res = true;
				}
			}
			else
			{
			}
		}
	}

	return(res);
}

int inserted_dup_copy(int id, struct ID_List *dlist, int *num_dup, struct DotList *dots, int num_lines, int mode, int *add_info, struct kdnode *tree, struct perm_pt *p_pts, int size, FILE *fp, struct DotList *init_dots) 
{
	bool *x_ins;
	bool *f_is_x;
	int *t_ins;
	bool *is_t_ins;
	int xval, yval;
	int w_sid, w_fid, h_sid, h_fid;
	int opt_id;
	int *cid;
	int num_ins_copy = 0;
	int new_id;
	int sp_state;

	x_ins = (bool *) ckalloc(sizeof(bool));
	cid = (int *) ckalloc(sizeof(int));
	f_is_x = (bool *) ckalloc(sizeof(bool));
	t_ins = (int *) ckalloc(sizeof(int));
	is_t_ins = (bool *) ckalloc(sizeof(bool));

	opt_id = -1;
	if( dots[id].sign == 2) {}
	else if( dots[id].pair_self == PAIR ) {}
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

		opt_id = find_opt_du_copy(dots, num_lines, id, p_pts, tree, size, w_sid, w_fid, h_sid, h_fid, cid, x_ins, f_is_x, t_ins, fp, init_dots);
	}

	if( opt_id != -1 ) 
	{
		if( ((*t_ins) == 11) || ((*t_ins) == 12) || ((*t_ins) == 21) || ((*t_ins) == 22) )
		{
			*is_t_ins = true;	
		}
		else *is_t_ins = false;

		dlist[num_ins_copy].is_x = *x_ins;
		dlist[num_ins_copy].m_id = *cid;
		dlist[num_ins_copy].left_id = opt_id;
		dlist[num_ins_copy].right_id = id;
		dlist[num_ins_copy].f_is_x = *f_is_x;
		dlist[num_ins_copy].is_t_ins = *is_t_ins;

		if( (*is_t_ins) == false ) 
		{
			new_id = check_inclusion_close_dup(dlist[num_ins_copy].m_id, dots, num_lines, x_ins, is_t_ins);
			if( new_id != -1 )
			{
				dlist[num_ins_copy].m_id = new_id;
				dlist[num_ins_copy].is_x = *x_ins;
				dlist[num_ins_copy].is_t_ins = *is_t_ins;
			}
		}

		if( dots[dlist[num_ins_copy].left_id].lock == -1 )
		{
			dots[dlist[num_ins_copy].left_id].lock = RIGHT_LOCK;
		}
		else if( dots[dlist[num_ins_copy].left_id].lock == LEFT_LOCK )
		{
			dots[dlist[num_ins_copy].left_id].lock = BOTH_LOCK;
		}

		if( dots[dlist[num_ins_copy].right_id].lock == -1 )
		{
			dots[dlist[num_ins_copy].right_id].lock = LEFT_LOCK;
		}
		else if( dots[dlist[num_ins_copy].right_id].lock == RIGHT_LOCK )
		{
			dots[dlist[num_ins_copy].right_id].lock = BOTH_LOCK;
		}

		num_ins_copy++;

		if( ((*t_ins) == 21) || ((*t_ins) == 22 ) )
		{
			if( (*x_ins) == true )
			{
				dlist[num_ins_copy].is_x = false;
			}
			else dlist[num_ins_copy].is_x = true;

			dlist[num_ins_copy].m_id = dlist[num_ins_copy-1].m_id;
			dlist[num_ins_copy].left_id = dlist[num_ins_copy-1].left_id;
			dlist[num_ins_copy].right_id = dlist[num_ins_copy-1].right_id;
			dlist[num_ins_copy].f_is_x = dlist[num_ins_copy-1].f_is_x;
			dlist[num_ins_copy].is_t_ins = dlist[num_ins_copy-1].is_t_ins;
			num_ins_copy++;

 			if( debug_mode == true ) printf("Double ==> ");
		}

		if( debug_mode == true ) {
			printf("Insertion of a duplicated copy: %d-%d %d-%d, %d\n", dots[dlist[num_ins_copy-1].m_id].x.lower, dots[dlist[num_ins_copy-1].m_id].x.upper, dots[dlist[num_ins_copy-1].m_id].y.lower, dots[dlist[num_ins_copy-1].m_id].y.upper, dlist[num_ins_copy-1].m_id);
			printf("%d-%d %d-%d\n", dots[opt_id].x.lower, dots[opt_id].x.upper, dots[opt_id].y.lower, dots[opt_id].y.upper);
			printf("%d-%d %d-%d, %d, %d, %d\n", dots[id].x.lower, dots[id].x.upper, dots[id].y.lower, dots[id].y.upper, dots[*cid].identity, dots[opt_id].identity, dots[id].identity);
		}

		if( (debug_mode == TRUE) && (dlist[num_ins_copy-1].is_t_ins == true) )
		{
			printf("Tandem Insertion\n");
		}
	}

	sp_state = check_ins_dup_copy(id, dlist, num_ins_copy, dots, mode, add_info);	
	
	*num_dup = num_ins_copy;
	free(is_t_ins);
	free(t_ins);
	free(f_is_x);
	free(cid);
	free(x_ins);

	return(sp_state);
}

bool check_whole_del_region_inclusion(struct DotList *dots, int num_lines, struct I del_reg)
{
	int i;
	bool res = false;

	for( i = 0; i < num_lines; i++ )
	{
		if( (subset(del_reg, dots[i].x) == true) || (subset(del_reg, dots[i].y) == true ) ) res = true;
	}

	return(res);
}

bool tandem_exist(struct DotList *dots, struct perm_pt *p_pts, struct kdnode *tree, int size, int id1, int id2)
{
	bool res = false;
	struct I reg1, reg2;
	int sid = 0, eid = 0;
	int i = 0;
	int cur_id = 0;

	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);
	
	if( (dots[id1].sign == dots[id2].sign) && (proper_overlap(dots[id1].x, dots[id2].x) == true ) && (proper_overlap(dots[id1].y, dots[id2].y) == true) ) {
		reg1 = intersect(dots[id1].x, dots[id2].x);			
		sid = find_id_len(tree, size, width(reg1), reg1.lower, reg1.lower, W_SID);
		eid = find_id_len(tree, size, width(reg1), reg1.upper, reg1.upper, W_FID);

		i = sid;
		while( (i <= eid) && (res == false) ) {
			cur_id = p_pts[i].id;
			if( is_tandem(dots[cur_id]) == true ) res = true;
			i++;
		}

		if( res == false ) {
			reg2 = intersect(dots[id1].y, dots[id2].y);
			sid = find_id_len(tree, size, width(reg2), reg2.lower, reg2.lower, W_SID);
			eid = find_id_len(tree, size, width(reg2), reg2.upper, reg2.upper, W_FID);

			i = sid;
			while( (i <= eid) && (res == false)) {
				cur_id = p_pts[i].id;
				if( is_tandem(dots[cur_id]) == true ) res = true;
				i++;
			}
		}	
	}

	return(res);
}

int find_id_len(struct kdnode *tree, int size, int len, int xval, int yval, int option)
{
	int res_id = 0;
	int x = 0, y = 1;

	if( option == W_SID )
	{
		if( len > RANGE_TH )
		{
			x = xval - len;
			y = yval - len;
		}
		else {
			x = xval - RANGE_TH;
			y = yval - RANGE_TH;
		}

		if( x <= 0 ) x = 1;
		if( y <= 0 ) y = 1;
		res_id = find_pred_blk(tree, x, y);
	}
	else if( option == W_FID )
	{
		if( len > RANGE_TH )
		{
			x = xval + len;
			y = yval + len;
		}
		else 
		{
			x = xval + RANGE_TH;
			y = yval + RANGE_TH;
		}
		if( x >= size ) x = size;
		if( y >= size ) y = size;
		res_id = find_successor(tree, x, y);
	}

	return(res_id);
}
