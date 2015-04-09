#include "main.h"
#include "regions.h"
#include "find_merging.h"
#include "util_i.h"

int distance(struct DotList *dots, int loc_id, int comp_id, bool *is_x, int *sd, int flag)
{
	int x, y;
	int t_x1, t_x2, t_y1, t_y2;
	int d;
	int res;
	int left_id, right_id;

	if( dots[loc_id].x.lower < dots[comp_id].x.lower )
	{
		left_id = loc_id;
		right_id = comp_id;
	}
	else if( dots[loc_id].x.lower > dots[comp_id].x.lower )
	{
		left_id = comp_id;
		right_id = loc_id;
	}
	else
	{
		if( dots[loc_id].x.upper <= dots[comp_id].x.upper )
		{
			left_id = loc_id;
			right_id = comp_id;
		}
		else
		{
			left_id = comp_id;
			right_id = loc_id;
		}
	}

	t_x1 = dots[left_id].x.upper; 
	t_x2 = dots[right_id].x.lower;

	if( dots[left_id].sign == 0 )
	{
		t_y1 = dots[left_id].y.upper;
		t_y2 = dots[right_id].y.lower;
	}	
	else if( dots[left_id].sign == 1 )
	{
		t_y1 = dots[left_id].y.lower;
		t_y2 = dots[right_id].y.upper;
	}
	
	d = compute_distance(dots[left_id].x, dots[left_id].y, dots[right_id].x, dots[right_id].y, dots[left_id].sign);

	if( overlap(dots[left_id].x, dots[right_id].x) == true )
	{
		if( (d > 2*DIS_THRESHOLD) && (flag == TRANSITIVE) ) {
			return (MDIS_THRESHOLD+1);
		}
		else x = 0;
	}
 	else x = abs(t_x1 - t_x2);

	if( overlap(dots[left_id].y, dots[right_id].y) == true)
	{
		if( (d > 2*DIS_THRESHOLD) && (flag == TRANSITIVE) ) {
			return (MDIS_THRESHOLD+1);
		}
		else y = 0;
	}
	else y = abs(t_y1 - t_y2);

	if( x >= y ) 
	{
		*is_x = true;
		*sd = y;
		res = x;
	}	
	else 
	{
		*is_x = false;
		*sd = x;
		res = y;
	}

	return(res);
}

bool check_candi(struct DotList *dots, int id, int comp_id, bool is_x)
{
	struct I t1, t2;
	
	if( is_x == true )
	{
		t1 = assign_I(dots[id].y.lower, dots[id].y.upper);
		t2 = assign_I(dots[comp_id].y.lower, dots[id].y.upper);
	}
	else
	{
		t1 = assign_I(dots[id].x.lower, dots[id].x.upper);
		t2 = assign_I(dots[comp_id].x.lower, dots[comp_id].x.upper);
	}

	if( proper_overlap(t1, t2) == true)
	{
		if(	width(intersect(t1, t2)) >= (3*RP_BD)) return(false);
		else return(true);
	}
	else return(true);
}

void merging_step(struct DotList *dots, int loc, int comp)
{
	struct I t_x, t_y;
	struct I x, y;
	int from, to;
	
	x = assign_I(dots[loc].x.lower, dots[loc].x.upper);
	y = assign_I(dots[loc].y.lower, dots[loc].y.upper);
	t_x = assign_I(dots[comp].x.lower, dots[comp].x.upper);
	t_y = assign_I(dots[comp].y.lower, dots[comp].y.upper);

	dots[loc].identity = ((dots[loc].identity * width(x)) + (dots[comp].identity * width(t_x))) / (width(x) + width(t_x));

	if( x.lower <= t_x.lower ) from = x.lower;
	else from = t_x.lower;

	if( t_x.upper >= x.upper ) to = t_x.upper;
	else to = x.upper;

	dots[loc].x = assign_I(from, to);

	if( y.lower <= t_y.lower ) from = y.lower;
	else from = t_y.lower;

	if( t_y.upper >= y.upper ) to = t_y.upper;
	else to = y.upper;
	dots[loc].y = assign_I(from, to);

	dots[comp].sign = 2;
}

int m_val(struct slist *st, struct DotList *dots, int i)
{
	if(st[i].is_x) return(dots[st[i].id].x.lower);
	else return(dots[st[i].id].y.lower);
}

int compute_closeness(struct DotList *dots, int loc_id, int comp_id)
{
	int res;
	struct I self_alg;

	self_alg = assign_I(1, 100000);

	if( comp_id == SELF_ID )
	{
		res = compute_distance(dots[loc_id].x, dots[loc_id].y, self_alg, self_alg, dots[loc_id].sign);
	}
	else
	{
		res = compute_distance(dots[loc_id].x, dots[loc_id].y, dots[comp_id].x, dots[comp_id].y, dots[loc_id].sign);
	}

	return(res);
}
