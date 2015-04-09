#include "regions.h"
#include "check_copy.h"
#include "find_merging.h"
#include "util_gen.h"
#include "deal_gaps.h"
#include "adjust_algn.h"
#include "util.h"
#include "util_i.h"

int get_score_copy(struct DotList *dots, int num_lines, struct gap_list gps, int *cid, bool *x_ins)
{
	int score = -1;

	if( (gps.type == 0) || (gps.type == -1)) 
	{
		return(-1);
	}
	else
	{
		if( is_back_to_td(dots, gps) ) return(-1);
		else if( (score = check_insertion_copy(dots, num_lines, gps, cid, x_ins)) != -1 ) return(score);
		else return(-1);
	}
}

bool is_back_to_td(struct DotList *dots, struct gap_list gp) 
{
	struct I x1, x2, y1, y2;
	struct I gap;
	struct I x, y;
	int low, up;
	int pid_diff = 0;

	pid_diff = abs(dots[gp.id1].identity - dots[gp.id2].identity);

	x1 = assign_I(dots[gp.id1].x.lower, dots[gp.id1].x.upper);
	x2 = assign_I(dots[gp.id2].x.lower, dots[gp.id2].x.upper);
	y1 = assign_I(dots[gp.id1].y.lower, dots[gp.id1].y.upper);
	y2 = assign_I(dots[gp.id2].y.lower, dots[gp.id2].y.upper);
	gap = assign_I(gp.y1, gp.y2);
	
	if( (gp.type == 1) || (gp.type == 11) || (gp.type == 21) ) {
		if( proper_overlap(x1, gap) && (gap.lower < x1.lower) ) x1 = assign_I(gap.lower, x1.upper - width(gap));
		if( proper_overlap(x2, gap) && (gap.lower < x2.lower) ) x2 = assign_I(gap.lower, x2.upper - width(gap));
		y1 = assign_I(y1.lower - width(gap), y1.upper - width(gap));
		y2 = assign_I(y2.lower - width(gap), y2.upper - width(gap));
	}
	else if( (gp.type == 2) || (gp.type == 12) || (gp.type == 22) ) {
		if( proper_overlap(y1, gap) && (gap.lower < y1.lower) ) y1 = assign_I(gap.lower, y1.upper - width(gap));
		if( proper_overlap(y2, gap) && (gap.lower < y2.lower) ) y2 = assign_I(gap.lower, y2.upper - width(gap));
	}
	else return(false);

	if( x1.lower <= x2.lower ) low = x1.lower;
	else low = x2.lower;

	if( x1.upper <= x2.upper ) up = x2.upper;
	else up = x1.upper;

	x = assign_I(low, up);

	if( y1.lower <= y2.lower ) low = y1.lower;
	else low = y2.lower;

	if( y1.upper <= y2.upper ) up = y2.upper;
	else up = y1.upper;

	y = assign_I(low, up);

	if( proper_overlap(x, y) && (width(intersect(x,y)) >= ALT_EFFEC_VALUE) && (pid_diff > PID_DIFF_TH_STRICT) ) {
		return(true);
	}
	else return(false);

}

int check_insertion_copy(struct DotList *dots, int num_lines, struct gap_list gp, int *cid, bool *x_ins)
{
	int pid_algns;
	int score;
	int i;
	struct I x1, x2, y1, y2;
	int min_score = 1000;
	int max_pid = 0;
	int pid;
	struct I temp;
	int type, f_type;

	x1 = assign_I(dots[gp.id1].x.lower, dots[gp.id1].x.upper);
	x2 = assign_I(dots[gp.id2].x.lower, dots[gp.id2].x.upper);
	y1 = assign_I(dots[gp.id1].y.lower, dots[gp.id1].y.upper);
	y2 = assign_I(dots[gp.id2].y.lower, dots[gp.id2].y.upper);

	if( abs(dots[gp.id1].identity - dots[gp.id2].identity) > 5 ) pid_algns = 116;
	else pid_algns = ((dots[gp.id1].identity * width(x1)) + (dots[gp.id2].identity * width(x2))) / (width(x1) + width(x2));

	temp = assign_I(gp.y1, gp.y2);
	f_type = gp.type;
	for(i = 0; i < num_lines; i++)
	{
		type = gp.type;
		if( dots[i].pair_self == PAIR ) {}
		else if( (i == gp.id1) || ( i == gp.id2) ) {}
		else if( dots[i].sign == 2 ) {} 
		else if( proper_overlap(dots[i].x, assign_I(gp.y1, gp.y2)) == true) 
		{
			if( check_again_equal(dots[i].x, gp, dots) == true )
			{
				if( dots[i].identity >= (pid_algns-3) ) 
				{
					pid = dots[i].identity;
					score = (int)((((float)(abs(dots[i].x.lower - gp.y1) + abs(dots[i].x.upper - gp.y2)))/((float)(gp.y2 - gp.y1)))*100); 

					if( ((type == 11) || (type == 12) || (type == 21) || (type == 22)) && (is_tandem(dots[i]) == false) ) {
						if( (type == 11) || (type == 21) ) type = 1;
						else type = 2;
					}

					if( (type == 11) || (type == 12) || (type == 21) || (type == 22) )
					{
						max_pid = 101;
						*cid = i;
						*x_ins = true;
						min_score = -2;
						f_type = type;
					}
					else if( pid > max_pid )
					{
						max_pid = pid;
						min_score = score;
						*cid = i;
						*x_ins = true;
						f_type = type;
					}
					else if( pid == max_pid )
					{
						if( score < min_score )
						{
							min_score = score;
							*cid = i;
							*x_ins = true;
							f_type = type;
						}
					}
				}
			}
		}

		type = gp.type;
		if( dots[i].pair_self == PAIR ) {}
		else if( (i == gp.id1) || ( i == gp.id2) ) {}
		else if( dots[i].sign == 2 ) {}
		else if( proper_overlap(dots[i].y, assign_I(gp.y1, gp.y2)) == true )
		{
			if(check_again_equal(dots[i].y, gp, dots) == true)
			{
				if( dots[i].identity >= (pid_algns-3) ) 
				{
					pid = dots[i].identity;
					score = (int)((((float)(abs(dots[i].y.lower - gp.y1) + abs(dots[i].y.upper - gp.y2)))/((float)(gp.y2 - gp.y1)))*100); 
					if( ((type == 11) || (type == 12) || (type == 21) || (type == 22)) && (is_tandem(dots[i]) == false) ) {
						if( (type == 11) || (type == 21) ) type = 1;
						else type = 2;
					}

					if( (type == 11) || (type == 12) || (type == 21) || (type == 22) )
					{
						max_pid = 101;
						min_score = -2;
						*cid = i;
						*x_ins = false;
						f_type = type;
					}
					else if( pid > max_pid )
					{
						max_pid = pid;
						min_score = score;
						*cid = i;
						*x_ins = false;
						f_type = type;
					}
					else if( pid == max_pid )
					{
						if( score < min_score )
						{
							*cid = i;
							*x_ins = false;
							min_score = score;
							f_type = type;
						}
					}
				}
			}
		}
	}

	if( min_score == 1000 ) score = -1;
	else score = min_score;

	gp.type = f_type;
	return(score);
}

bool check_boundary(struct I x, int from, int to)
{
	bool res;
	int len;

	len = to - from;

	if( width(x) <= S_RP_BD ) 
	{
		res = false;
	}
	else if( width(x) < (5*RP_BD) )
	{
		if(in( to, x ))
		{
			if(in( from, x ))
			{
				if((float)(((float)len) / ((float)width(x))) >= 0.8)
				{
					res = true;			
				}
				else res = false;
			}
			else
			{
				if((float)((float)(to - x.lower)) / ((float)width(x)) >= 0.8)
				{
					if((float)((float)(x.lower - from)) / ((float)width(x)) <= 0.2)
					{
						res = true;	
					}
					else res = false;
				}
				else res = false;
			}
		}
		else 
		{
			if((float)((float)(to - x.upper) / ((float)width(x))) <= 0.2)
			{
				if(in( from, x ))
				{
					if((float)((float)(x.upper - from)) / ((float)width(x)) >= 0.8)
					{
						if((float)((float)(from - x.lower)) / ((float)width(x)) <= 0.2)
						{
							res = true;	
						}
						else res = false;
					}
					else res = false;
				}
				else
				{
					if((float)((float)(x.lower - from)) / ((float)width(x)) <= 0.2)
					{
						res = true;			
					}
					else res = false;
				}
			}
			else res = false;
		}
	}
	else
	{
		if(in( to, x ))
		{
			if(in( from, x ))
			{
				if(((x.upper - to) <= RP_BD) && (from - x.lower <= RP_BD))
				{
					res = true;			
				}
				else res = false;
			}
			else
			{
				if((x.upper - to) <= RP_BD)
				{
					if((x.lower - from) <= RP_BD)
					{
						res = true;	
					}
					else res = false;
				}
				else res = false;
			}
		}
		else 
		{
			if((to - x.upper) <= RP_BD)
			{
				if(in( from, x ))
				{
					if( (from - x.lower) <= RP_BD )
					{
						res = true;	
					}
					else res = false;
				}
				else
				{
					if((x.lower - from) <= RP_BD)
					{
						res = true;			
					}
					else res = false;
				}
			}
			else res = false;
		}
	}
	return(res);
}

bool check_again_equal(struct I cur, struct gap_list gp, struct DotList *dots)
{
	int res;
	int op_len;
	struct DotList aln_1, aln_2; 
	struct I *temp_x, *temp_y;
	int *from;
	int *to;
	int new_len;
	struct I new_gap;

	temp_x = (struct I *) ckalloc(sizeof(struct I));
	temp_y = (struct I *) ckalloc(sizeof(struct I));
	from = (int *) ckalloc(sizeof(int));
	to = (int *) ckalloc(sizeof(int));

	(*temp_x) = assign_I(0,1);
	(*temp_y) = assign_I(0,1);
	*from = 1;
	*to = 0;

	aln_1.x = assign_I(dots[gp.id1].x.lower, dots[gp.id1].x.upper);
	aln_1.y = assign_I(dots[gp.id1].y.lower, dots[gp.id1].y.upper);
	aln_1.sign = dots[gp.id1].sign;
	aln_2.x = assign_I(dots[gp.id2].x.lower, dots[gp.id2].x.upper);
	aln_2.y = assign_I(dots[gp.id2].y.lower, dots[gp.id2].y.upper);
	aln_2.sign = dots[gp.id2].sign;

	if( gp.type == 1 )
	{
		op_len = abs(gp.x1-gp.x2);
		if( proper_overlap(dots[gp.id1].x, dots[gp.id2].x) == true )
		{
			op_len = width(intersect(dots[gp.id1].x, dots[gp.id2].x));
			if( (proper_overlap(dots[gp.id1].y, cur) == true) || (proper_overlap(dots[gp.id2].y, cur) == true) )
			{
				*temp_x = assign_I(aln_1.x.lower, aln_1.x.upper);
				*temp_y = assign_I(aln_1.y.lower, aln_1.y.upper);
				if( proper_overlap(dots[gp.id1].y, cur) == true ) 
				{
					res = adjust_algn_reg(temp_x, temp_y, cur, aln_1.sign);

					if( res == 0 )
					{
						aln_1.x = assign_I((*temp_x).lower, (*temp_x).upper);
						aln_1.y = assign_I((*temp_y).lower, (*temp_y).upper);
					}
				}

				*temp_x = assign_I(aln_2.x.lower, aln_2.x.upper);
				*temp_y = assign_I(aln_2.y.lower, aln_2.y.upper);
				if( proper_overlap(dots[gp.id2].y, cur) == true ) 
				{
					res = adjust_algn_reg(temp_x, temp_y, cur, aln_2.sign);

					if( res == 0 )
					{
						aln_2.x = assign_I((*temp_x).lower, (*temp_x).upper);
						aln_2.y = assign_I((*temp_y).lower, (*temp_y).upper);
					}
				}

				new_len = get_starting_loc_reg(aln_1, aln_2, false, from, to);

				if( (*from) > (*to) )
				{
					free(to);
					free(from);
					free(temp_x);
					free(temp_y);
					return false;
				}
				else
				{
					if( new_len < op_len )
					{
						
						if( (new_len <= M_TH) || ((new_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
						{
							new_gap = assign_I((*from), (*to));
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return( strict_almost_equal(cur, new_gap) );
						}
						else
						{
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return false;
						}
					}
					else
					{
						if( (op_len <= M_TH) || ((op_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
						{
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return( strict_almost_equal(cur, assign_I(gp.y1, gp.y2) ) );
						}
						else 
						{
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return false;
						}
					}
				}
			}
			else
			{
				if( (op_len <= M_TH) || ((op_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
				{
					free(to);
					free(from);
					free(temp_x);
					free(temp_y);
					return( strict_almost_equal(cur, assign_I(gp.y1, gp.y2) ) );
				}
				else 
				{
					free(to);
					free(from);
					free(temp_x);
					free(temp_y);
					return false;
				}
			}
		}
		else
		{
			if( (op_len <= M_TH) || ((op_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
			{
				free(to);
				free(from);
				free(temp_x);
				free(temp_y);
				return( strict_almost_equal(cur, assign_I(gp.y1, gp.y2) ) );
			}
			else 
			{
				free(to);
				free(from);
				free(temp_x);
				free(temp_y);
				return false;
			}
		}
	}
	else if( gp.type == 2 )
	{
		op_len = abs(gp.x1-gp.x2);
		if( proper_overlap(dots[gp.id1].y, dots[gp.id2].y) == true )
		{
			op_len = width(intersect(dots[gp.id1].y, dots[gp.id2].y));

			if( (proper_overlap(dots[gp.id1].x, cur) == true) || (proper_overlap(dots[gp.id2].x, cur) == true) )
			{
				*temp_x = assign_I(aln_1.x.lower, aln_1.x.upper);
				*temp_y = assign_I(aln_1.y.lower, aln_1.y.upper);
				if( proper_overlap(dots[gp.id1].x, cur) == true ) 
				{
					res = adjust_algn_reg(temp_x, temp_y, cur, aln_1.sign);

					if( res == 0 )
					{
						aln_1.x = assign_I((*temp_x).lower, (*temp_x).upper);
						aln_1.y = assign_I((*temp_y).lower, (*temp_y).upper);
					}
				}

				*temp_x = assign_I(aln_2.x.lower, aln_2.x.upper);
				*temp_y = assign_I(aln_2.y.lower, aln_2.y.upper);

				if( proper_overlap(dots[gp.id2].x, cur) == true ) 
				{
					res = adjust_algn_reg(temp_x, temp_y, cur, aln_2.sign);

					if( res == 0 )
					{
						aln_2.x = assign_I((*temp_x).lower, (*temp_x).upper);
						aln_2.y = assign_I((*temp_y).lower, (*temp_y).upper);
					}
				}

				new_len = get_starting_loc_reg(aln_1, aln_2, true, from, to);

				if( (*from) > (*to) )
				{
					free(to);
					free(from);
					free(temp_x);
					free(temp_y);
					return false;
				}
				else
				{
					if( new_len <= op_len )
					{
						if( (new_len <= M_TH) || ((new_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
						{
							new_gap = assign_I((*from), (*to));
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return( strict_almost_equal(cur, new_gap) );
						}
						else 
						{
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return false;
						}
					}
					else
					{
						if( (op_len <= M_TH) || ((op_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
						{
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return( strict_almost_equal(cur, assign_I(gp.y1, gp.y2) ) );
						}
						else
						{
							free(to);
							free(from);
							free(temp_x);
							free(temp_y);
							return false;
						}
					}
				}
			}
			else
			{
				if( (op_len <= M_TH) || ((op_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
				{
					free(to);
					free(from);
					free(temp_x);
					free(temp_y);
					return( strict_almost_equal(cur, assign_I(gp.y1, gp.y2) ) );
				}
				else
				{
					free(to);
					free(from);
					free(temp_x);
					free(temp_y);
					return false;
				}
			}
		}
		else
		{
			if( (op_len <= M_TH) || ((op_len <= L_M_TH) && ( ( width(dots[gp.id1].x) >= LG_TH) && (width(dots[gp.id2].x) >= LG_TH ) )) )
			{
				free(to);
				free(from);
				free(temp_x);
				free(temp_y);
				return( strict_almost_equal(cur, assign_I(gp.y1, gp.y2) ) );
			}
			else 
			{
				free(to);
				free(from);
				free(temp_x);
				free(temp_y);
				return false;
			}
		}
	}
	else 
	{
 		if( ( gp.type == 11 ) || ( gp.type == 12) || ( gp.type == 21 ) || (gp.type == 22))
		{
			free(to);
			free(from);
			free(temp_x);
			free(temp_y);
			if( strict_almost_equal(cur, assign_I(gp.y1, gp.y2)) ) return true;
			else if( ((gp.y2 + gp.offset) > gp.y1) && (strict_almost_equal(cur, assign_I(gp.y1, gp.y2 + gp.offset))) ) return true;
			else return false;
		}
		else 
		{
			free(to);
			free(from);
			free(temp_x);
			free(temp_y);
			return false;
		}
	}

	free(to);
	free(from);
	free(temp_x);
	free(temp_y);
	return false;
}

int check_match_gap_aln(struct DotList *dots, struct gap_list gp, int res_id, bool is_x)
{
	int pid_algns;
	struct I x1, x2, y1, y2;
	int score = -1;
	struct I temp, ins_reg;

	x1 = assign_I(dots[gp.id1].x.lower, dots[gp.id1].x.upper);
	x2 = assign_I(dots[gp.id2].x.lower, dots[gp.id2].x.upper);
	y1 = assign_I(dots[gp.id1].y.lower, dots[gp.id1].y.upper);
	y2 = assign_I(dots[gp.id2].y.lower, dots[gp.id2].y.upper);

	pid_algns = ((dots[gp.id1].identity * width(x1)) + (dots[gp.id2].identity * width(x2))) / (width(x1) + width(x2));

	temp = assign_I(gp.y1, gp.y2);

	if( is_x == true ) ins_reg = dots[res_id].x;
	else ins_reg = dots[res_id].y;

	if( ((gp.type == 11) || (gp.type == 12) || (gp.type == 21) || (gp.type == 22)) && (is_tandem(dots[res_id]) == false) ) {}
	else if( (dots[res_id].identity >= (pid_algns-10)) && (check_again_equal(ins_reg, gp, dots) == true) )
	{
		score = (int)((((float)(abs(ins_reg.lower - gp.y1) + abs(ins_reg.upper - gp.y2)))/((float)(gp.y2 - gp.y1)))*100); 
	}

	return(score);
}
