#include "main.h"
#include "regions.h"
#include "deal_gaps.h"
#include "util.h"
#include "util_i.h"

struct gap_list define_gap_new_type_inc(struct DotList *dots, int loc_id, int comp_id, bool is_x)
{
	struct gap_list gp;
	int len_dif;
	int sm_id, lg_id;
	int diff;

	gp.type = -1;
	if( is_x == true ) // the overlap is in x region
	{
		if( strict_subset(dots[loc_id].x, dots[comp_id].x) == true )
		{
			sm_id = loc_id;
			lg_id = comp_id;
		}
		else if( strict_subset(dots[comp_id].x, dots[loc_id].x) == true )
		{
			sm_id = comp_id;
			lg_id = loc_id;
		}	

		if( (dots[sm_id].y.lower < dots[lg_id].y.lower) && (abs(dots[sm_id].x.lower - dots[lg_id].x.lower) <= ERR_SM_TH) )
		{
			if( dots[sm_id].sign == 0 ) 
			{
				diff = dots[lg_id].x.lower - dots[sm_id].x.lower;
				len_dif = abs(dots[lg_id].y.lower - dots[sm_id].y.lower) + diff;
				gp.y1 = dots[sm_id].y.upper;
				gp.y2 = gp.y1 + len_dif;
				gp.id1 = sm_id;
				gp.id2 = lg_id;
				gp.type = 11;
			}
		}
		else if( (dots[sm_id].y.lower > dots[lg_id].y.lower) && (abs(dots[sm_id].x.lower - dots[lg_id].x.lower) <= ERR_SM_TH) )
		{
			if( dots[sm_id].sign == 1 ) 
			{
				diff = dots[lg_id].x.lower - dots[sm_id].x.lower;
				len_dif = abs(dots[lg_id].y.upper - dots[sm_id].y.upper) + diff;
				gp.y2 = dots[sm_id].y.lower;
				gp.y1 = gp.y2 - len_dif;
				gp.id1 = sm_id;
				gp.id2 = lg_id;
				gp.type = 11;
			}
		}
		else if( (dots[sm_id].y.lower < dots[lg_id].y.lower) && (abs(dots[sm_id].x.upper - dots[lg_id].x.upper) <= ERR_SM_TH) )
		{
			if( dots[sm_id].sign == 1 ) 
			{
				diff = dots[lg_id].x.upper - dots[sm_id].x.upper;
				len_dif = abs(dots[lg_id].y.lower - dots[sm_id].y.lower) + diff;
				gp.y1 = dots[sm_id].y.upper;
				gp.y2 = gp.y1 + len_dif;
				gp.id1 = lg_id;
				gp.id2 = sm_id;
				gp.type = 11;
			}
		}
		else if( (dots[sm_id].y.lower > dots[lg_id].y.lower) && (abs(dots[sm_id].x.upper - dots[lg_id].x.upper) <= ERR_SM_TH) )
		{	
			if( dots[sm_id].sign == 0 )
			{
				diff = dots[lg_id].x.upper - dots[sm_id].x.upper;
				len_dif = abs(dots[lg_id].y.upper - dots[sm_id].y.upper) + diff;
				gp.y2 = dots[sm_id].y.lower;
				gp.y1 = gp.y2 - len_dif;
				gp.id1 = lg_id;
				gp.id2 = sm_id;
				gp.type = 11;
			}	
		}
	}
	else if( is_x == false ) // the overlap is in y region
	{
		if( strict_subset(dots[loc_id].y, dots[comp_id].y) == true )
		{
			sm_id = loc_id;
			lg_id = comp_id;
		}
		else if( strict_subset(dots[comp_id].y, dots[loc_id].y) == true )
		{
			sm_id = comp_id;
			lg_id = loc_id;
		}	

		if( (dots[sm_id].x.lower < dots[lg_id].x.lower) && (abs(dots[sm_id].y.lower - dots[lg_id].y.lower) <= ERR_SM_TH) )
		{
			if( dots[sm_id].sign == 0 ) 
			{
				diff = dots[lg_id].y.lower - dots[sm_id].y.lower;
				len_dif = abs(dots[lg_id].x.lower - dots[sm_id].x.lower) + diff;
				gp.y1 = dots[sm_id].x.upper;
				gp.y2 = gp.y1 + len_dif;
				gp.id1 = sm_id;
				gp.id2 = lg_id;
				gp.type = 12;
			}
		}
		else if( (dots[sm_id].x.lower > dots[lg_id].x.lower) && (abs(dots[sm_id].y.lower - dots[lg_id].y.lower) <= ERR_SM_TH) )
		{
			if( dots[sm_id].sign == 1 ) 
			{
				diff = dots[lg_id].y.lower - dots[sm_id].y.lower;
				len_dif = abs(dots[lg_id].x.upper - dots[sm_id].x.upper) + diff;
				gp.y2 = dots[sm_id].x.lower;
				gp.y1 = gp.y2 - len_dif;
				gp.id1 = lg_id;
				gp.id2 = sm_id;
				gp.type = 12;
			}
		}
		if( (dots[sm_id].x.lower < dots[lg_id].x.lower) && (abs(dots[sm_id].y.upper - dots[lg_id].y.upper) <= ERR_SM_TH) )
		{
			if( dots[sm_id].sign == 1 ) 
			{
				diff = dots[lg_id].y.upper - dots[sm_id].y.upper;
				len_dif = abs(dots[lg_id].x.lower - dots[sm_id].x.lower) + diff;
				gp.y1 = dots[sm_id].x.upper;
				gp.y2 = gp.y1 + len_dif;
				gp.id1 = sm_id;
				gp.id2 = lg_id;
				gp.type = 12;
			}
		}
		else if( (dots[sm_id].x.lower > dots[lg_id].x.lower) && (abs(dots[sm_id].y.upper - dots[lg_id].y.upper) <= ERR_SM_TH) )
		{	
			if( dots[sm_id].sign == 0 )
			{
				diff = dots[lg_id].y.upper - dots[sm_id].y.upper;
				len_dif = abs(dots[lg_id].x.upper - dots[sm_id].x.upper) + diff;
				gp.y2 = dots[sm_id].x.lower;
				gp.y1 = gp.y2 - len_dif;
				gp.id1 = lg_id;
				gp.id2 = sm_id;
				gp.type = 12;
			}	
		}
	}
	
	return(gp);
}

/* when two alignments have an overlapped region */
struct gap_list define_gap_new_type(struct DotList *dots, int loc_id, int comp_id, bool is_x)
{
	struct gap_list gp;
	struct I temp;
	int len_x, len_y;

	gp.id1 = loc_id;
	gp.id2 = comp_id;

	gp.type = -1;
	gp.x1 = 0;
	gp.x2 = 1;
	gp.y1 = 0;
	gp.y2 = 1;
	gp.offset = 0;
	if( is_x == true ) // the overlap of x region is larger than y's
	{
		if( proper_overlap(dots[loc_id].x, dots[comp_id].x) == true )
		{
			temp = intersect(dots[loc_id].x, dots[comp_id].x);
			gp.type = 21; // the gap is in y side
			
			if( dots[loc_id].y.lower <= dots[comp_id].y.lower )
			{
				gp.y1 = dots[loc_id].y.upper;
				gp.y2 = dots[comp_id].y.lower + width(temp);
				len_x = width(dots[comp_id].x);
				len_y = width(dots[comp_id].y);
				gp.offset = len_y - len_x;

				if( dots[loc_id].sign == 0 )
				{
					gp.x1 = dots[loc_id].x.upper;
					gp.x2 = gp.x1 + 1;
				}
				else if( dots[loc_id].sign == 1 )
				{
					gp.x1 = dots[loc_id].x.lower;
					gp.x2 = gp.x1 + 1;
				}
				else gp.type = -1;
			}
			else
			{
				gp.y1 = dots[comp_id].y.upper;
				gp.y2 = dots[loc_id].y.lower + width(temp);
				len_x = width(dots[loc_id].x);
				len_y = width(dots[loc_id].y);
				gp.offset = len_y - len_x;

				if( dots[comp_id].sign == 0 )
				{
					gp.x1 = dots[comp_id].x.upper;
					gp.x2 = gp.x1 + 1;
				}
				else if( dots[comp_id].sign == 1 )
				{
					gp.x1 = dots[comp_id].x.lower;
					gp.x2 = gp.x1 + 1;
				}
				else gp.type = -1;
			}
		}
		else
		{
			gp.type = -1;
		}
	}
	else 
	{
		if( proper_overlap(dots[loc_id].y, dots[comp_id].y) == true )
		{
			temp = intersect(dots[loc_id].y, dots[comp_id].y);
			gp.type = 22; // the gap is in x side

			if( dots[loc_id].x.lower <= dots[comp_id].x.lower )
			{
				gp.y1 = dots[loc_id].x.upper;
				gp.y2 = dots[comp_id].x.lower + width(temp);
				len_x = width(dots[comp_id].x);
				len_y = width(dots[comp_id].y);
				gp.offset = len_x - len_y;

				if( dots[loc_id].sign == 0 )
				{
					gp.x1 = dots[loc_id].y.upper;
					gp.x2 = gp.x1 + 1;
				}
				else if( dots[loc_id].sign == 1 )
				{
					gp.x1 = dots[loc_id].y.lower;
					gp.x2 = gp.x1 + 1;
				}
				else gp.type = -1;
			}
			else 
			{
				gp.y1 = dots[comp_id].x.upper;
				gp.y2 = dots[loc_id].x.lower + width(temp);
				len_x = width(dots[loc_id].x);
				len_y = width(dots[loc_id].y);
				gp.offset = len_x - len_y;
			}

		}
		else
		{
			gp.type = -1;
		}
	}

	if( (gp.type != -1) && (gp.y2 <= gp.y1) ) {
		gp.type = -1;
	}

	if( gp.type != -1 ) {
		temp = assign_I(gp.y1, gp.y2);
		if( ( strict_almost_equal(temp, dots[comp_id].x) == true ) || ( strict_almost_equal(temp, dots[comp_id].y) == true ) || ( strict_almost_equal(temp, dots[loc_id].x) == true ) || (strict_almost_equal(temp, dots[loc_id].y) == true ))
		{
			gp.type = -1;
		}
	}	

	return(gp);
}

struct gap_list define_gap(struct DotList *dots, int loc_id, int comp_id, int d, int sd, bool is_x)
{
	int from, to;
	int from_1, to_1;
	struct gap_list gp;
	struct I temp;
	int m_th;

	gp.type = 3;
	gp.id1 = loc_id;
	gp.id2 = comp_id;

	if( ( width(dots[loc_id].x) >= LG_TH ) && ( width(dots[comp_id].x) >= LG_TH ) )
	{
		m_th = L_M_TH;
	}
	else m_th = M_TH;
	
	if( overlap(dots[loc_id].x, dots[comp_id].x) == true )
	{
		if(width(intersect(dots[loc_id].x, dots[comp_id].x)) == 0)
		{
			if( dots[loc_id].x.upper == dots[comp_id].x.lower )
			{
				gp.x1 = dots[loc_id].x.upper;
				gp.x2 = gp.x1+1;
			}
			else
			{
			}
		}
		else
		{
			temp = intersect(dots[loc_id].x, dots[comp_id].x);
			gp.x1 = temp.lower;
			gp.x2 = temp.upper;
		}

		if(gp.type == 3) gp.type = 1;
	}

	if( overlap(dots[loc_id].y, dots[comp_id].y) == true)
	{
		if(gp.type == 1) gp.type = 0;
		else if(gp.type == 3) gp.type = 2;

		if(width(intersect(dots[loc_id].y, dots[comp_id].y)) == 0)
		{
			if(dots[loc_id].y.upper == dots[comp_id].y.lower)
			{
				if( gp.type == 0 )
				{
					gp.y1 = dots[loc_id].y.upper;
					gp.y2 = gp.y1+1;
				}
				else if( gp.type == 2 )
				{
					gp.x1 = dots[loc_id].y.upper;
					gp.x2 = gp.x1+1;
				}
			}
			else if(dots[loc_id].y.lower == dots[comp_id].y.upper)
			{
				if( gp.type == 0 )
				{
					gp.y1 = dots[loc_id].y.lower;
					gp.y2 = gp.y1+1;
				}
				else if( gp.type == 2)
				{
					gp.x1 = dots[loc_id].y.lower;
					gp.x2 = gp.x1+1;
				}
			}
			else
			{
			}
		}
		else
		{
			temp = intersect(dots[loc_id].y, dots[comp_id].y);
			if( gp.type == 0 )
			{
				gp.y1 = temp.lower;
				gp.y2 = temp.upper;
			}
			else if( gp.type == 2 )
			{
				gp.x1 = temp.lower;
				gp.x2 = temp.upper;
			}
		}
	}

	if((gp.type == 1) || (gp.type == 2) || (gp.type == 3)) 
	{
		from = get_starting_loc(dots, loc_id, comp_id, is_x);
		to = from + d;
		if( (gp.type == 1) || (gp.type == 2) || ((gp.type == 3) && (sd <= m_th)))
		{
			gp.y1 = from;
			gp.y2 = to;

			if( gp.type == 3 )
			{
				if( is_x == true ) gp.type = 2;
				else gp.type = 1;

				from_1 = get_starting_loc(dots, loc_id, comp_id, (!is_x));
				to_1 = from_1 + sd;
				gp.x1 = from_1;
				gp.x2 = to_1;
			}

			if( d <= MIN_GAP ) gp.type = 0;
		}
		else if((gp.type == 3) && (sd > m_th))
		{
			gp.type = -1;
			gp.x1 = 0;
			gp.x2 = 0;
			gp.y1 = 0;
			gp.y2 = 0;
		}
	}

	return(gp);
}

bool check_inclusion_alignments(struct gap_list gp, struct DotList *dots, int num)
{
	struct I x, y;
	bool res = false;
	int i;

	if( gp.x1 >= gp.x2 )
	{
	}
	else
	{
		x = assign_I(gp.x1, gp.x2);
	}

	if( gp.y1 >= gp.y2 )
	{
	}
	else
	{
		y = assign_I(gp.y1, gp.y2);
	}

	for( i = 0; i < num; i++ )
	{
		if( dots[i].sign != 2 )
		{
			if( (subset(dots[i].x, x) == true) && (dots[i].identity > dots[gp.id1].identity) )
			{
				res = true;
			}
			else if( (subset(dots[i].y, y) == true) && (dots[i].identity > dots[gp.id2].identity))
			{
				res = true;
			}
		}
	}
	
	return(res);
}

int check_into_own(struct DotList *dots, int loc_id, int comp_id)
{
	struct I temp;
	int res;

	if( (dots[loc_id].sign == 0) && (dots[comp_id].sign == 0) )
	{
		if( (strict_overlap(dots[loc_id].x, dots[loc_id].y, 10*T_OP_TH) == true) && (strict_overlap(dots[comp_id].x, dots[comp_id].y, 10*T_OP_TH) == true ) && (dots[loc_id].x.upper < dots[comp_id].x.lower) )
		{
			temp = assign_I(dots[loc_id].x.upper, dots[comp_id].x.lower);
			if( strict_almost_equal(dots[loc_id].y, temp) == true )
			{
				res = COPY_OWN;	
			}
			else res = NON_COPY;
		}
		else if( (strict_overlap(dots[loc_id].x, dots[loc_id].y, 10*T_OP_TH) == true) && (strict_overlap(dots[comp_id].x, dots[comp_id].y, 10*T_OP_TH) == true) && (dots[comp_id].x.upper < dots[loc_id].x.lower ) )
		{
			temp = assign_I(dots[comp_id].x.upper, dots[loc_id].x.lower);
			if( strict_almost_equal(dots[comp_id].y, temp) == true )
			{
				res = COPY_OWN;
			}
			else res = NON_COPY;
		}
		else 
		{
			res = NON_COPY;
		}
	}
	else if( (dots[loc_id].sign == 1) && (dots[comp_id].sign == 1) )
	{
		if( (strict_overlap(dots[loc_id].x, dots[comp_id].x, 10*T_OP_TH) == true ) && (strict_overlap(dots[loc_id].y, dots[comp_id].y, 10*T_OP_TH) == true ) && (dots[comp_id].x.upper < dots[comp_id].y.lower))
		{
			temp = assign_I(dots[comp_id].x.upper, dots[comp_id].y.lower);
			if( strict_almost_equal(dots[loc_id].y, temp) == true )
			{
				res = COPY_OWN_INV;
			}
			else res = NON_COPY;
		}
		else if( (strict_overlap(dots[loc_id].x, dots[comp_id].x, 10*T_OP_TH) == true ) &&(strict_overlap(dots[loc_id].y, dots[comp_id].y, 10*T_OP_TH) == true ) && (dots[loc_id].x.upper < dots[loc_id].y.lower))
		{
			temp = assign_I(dots[loc_id].x.upper, dots[loc_id].y.lower);
			if( strict_almost_equal(dots[comp_id].y, temp) == true )
			{
				res = COPY_OWN_INV;
			}
			else res = NON_COPY;
		}
		else
		{
			res = NON_COPY;
		}
	}
	else res = NON_COPY;

	return(res);
}

struct gap_list redefine_for_del(struct DotList *dots, struct gap_list gps) // when the type of a gap is 2, 12, or 22
{
	struct I temp;
	int loc_id, comp_id;
	int len;
	int x;

	loc_id = gps.id1;
	comp_id = gps.id2;
	if( (gps.type == 12) || (gps.type == 22) )
	{
		if( proper_overlap(dots[loc_id].y, dots[comp_id].y) == true )
		{
			temp = intersect(dots[loc_id].y, dots[comp_id].y);
			gps.x1 = gps.y1;
			gps.x2 = gps.y2;
			gps.y1 = temp.lower;
			gps.y2 = temp.upper;
		}
	}
	else if( gps.type == 2 )
	{
		len = abs(gps.y2 - gps.y1);
		x = gps.y1;
		gps.y1 = gps.x1;
		gps.y2 = gps.x1 + len;
		gps.x1 = x;
		gps.x2 = x + len;
	}
	return(gps);
}

int get_starting_loc(struct DotList *dots, int loc, int comp, bool is_x)
{
	int from;

	if( is_x == true )
	{
		if( dots[loc].x.lower < dots[comp].x.lower )
		{
			from = dots[loc].x.upper;
		}
		else
		{
			from = dots[comp].x.upper;
		}
	}
	else
	{
		if( dots[loc].sign == 0 )
		{
			if( dots[loc].x.lower < dots[comp].x.lower )
			{
				from = dots[loc].y.upper;
			}
			else if( dots[loc].x.lower > dots[comp].x.lower )
			{
				from = dots[comp].y.upper;
			}
			else 
			{
				if( dots[loc].x.upper <= dots[comp].x.upper )
				{
					from = dots[loc].y.upper;
				}
				else
				{
					from = dots[comp].y.upper;
				}
			}
		}		
		else if( dots[loc].sign == 1 )
		{
			if( dots[loc].x.lower < dots[comp].x.lower )
			{
				from = dots[comp].y.upper;
			}
			else if( dots[loc].x.lower > dots[comp].x.lower )
			{
				from = dots[loc].y.upper;
			}
			else 
			{
				if( dots[loc].x.upper <= dots[comp].x.upper )
				{
					from = dots[comp].y.upper;
				}
				else
				{
					from = dots[loc].y.upper;
				}
			}
		}
	}
	return(from);
}

int get_starting_loc_reg(struct DotList aln_1, struct DotList aln_2, bool is_x, int *from, int *to)
{
	int len;

	if( is_x == true )
	{
		if( aln_1.x.lower < aln_2.x.lower )
		{
			*from = aln_1.x.upper;
			*to = aln_2.x.lower;
		}
		else
		{
			*from = aln_2.x.upper;
			*to = aln_1.x.lower;
		}

		if( abs(aln_1.y.upper - aln_2.y.lower) > abs(aln_2.y.upper - aln_1.y.lower) )
		{
			len = abs(aln_2.y.upper - aln_1.y.lower);
		}
		else
		{
			len = abs(aln_1.y.upper - aln_2.y.lower);
		}
	}
	else
	{
		if( aln_1.sign == 0 )
		{
			if( aln_1.x.lower < aln_2.x.lower )
			{
				*from = aln_1.y.upper;
				*to = aln_2.y.lower;
			}
			else if( aln_1.x.lower > aln_2.x.lower )
			{
				*from = aln_2.y.upper;
				*to = aln_1.y.lower;
			}
			else 
			{
				if( aln_1.x.upper <= aln_2.x.upper )
				{
					*from = aln_1.y.upper;
					*to = aln_2.y.lower;
				}
				else
				{
					*from = aln_2.y.upper;
					*to = aln_1.y.lower;
				}
			}
		}		
		else if( aln_1.sign == 1 )
		{
			if( aln_1.x.lower < aln_2.x.lower )
			{
				*from = aln_2.y.upper;
				*to = aln_1.y.lower;
			}
			else if( aln_1.x.lower > aln_2.x.lower )
			{
				*from = aln_1.y.upper;
				*to = aln_2.y.lower;
			}
			else 
			{
				if( aln_1.x.upper <= aln_2.x.upper )
				{
					*from = aln_2.y.upper;
					*to = aln_1.y.lower;
				}
				else
				{
					*from = aln_1.y.upper;
					*to = aln_2.y.lower;
				}
			}
		}

		if( abs(aln_1.x.upper - aln_2.x.lower) > abs(aln_2.x.upper - aln_1.x.lower) )
		{
			len = abs(aln_2.x.upper - aln_1.x.lower);
		}
		else
		{
			len = abs(aln_1.x.upper - aln_2.x.lower);
		}
	}
	return(len);
}
