#include "main.h"
#include "regions.h"
#include "adjust_algn.h"
#include "util_gen.h"
#include "util_i.h"
#include "find_merging.h"

void adjust_alignment(struct DotList *dots, int i, struct I to)
{
	int wide = 0;
	int flag = 0;
	int temp_val = 0;
	bool is_td = false;
	int len_x = 0, len_y = 0;

	is_td = is_tandem(dots[i]);
	len_x = width(dots[i].x);
	len_y = width(dots[i].y);

	if(proper_overlap(dots[i].x, to))
	{
		wide = width(intersect(dots[i].x, to));
		if( (wide + 5) > width(dots[i].x))
		{
			dots[i].sign = 2;
		}
		else if( in(dots[i].x.lower, intersect(dots[i].x, to)))
		{
			dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
			dots[i].xl_diff = dots[i].xl_diff + wide;
			flag = 1;
		}
		else if( in(dots[i].x.upper, intersect(dots[i].x, to)))
		{
			dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
			dots[i].xr_diff = dots[i].xr_diff + wide;
			flag = 2;
		}
		else
		{
			if( abs(dots[i].x.lower - to.lower) < abs(dots[i].x.upper - to.upper) )
			{
				temp_val = dots[i].x.lower;
				dots[i].x = assign_I(to.lower+wide, dots[i].x.upper);
				wide = wide + abs(temp_val - to.lower);
				dots[i].xl_diff = dots[i].xl_diff + wide;
				flag = 3;
			}
			else
			{
				temp_val = dots[i].x.upper;
				dots[i].x = assign_I(dots[i].x.lower, to.upper-wide);
				wide = wide + abs(temp_val - to.upper);
				dots[i].xr_diff = dots[i].xr_diff + wide;
				flag = 4;
			}
		}

		if( dots[i].sign == 0 )
		{
			if( (wide + 5) > width(dots[i].y))
			{
				dots[i].sign = 2;
			}
			else if( flag == 1 )
			{
				dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				dots[i].yl_diff = dots[i].yl_diff + wide;
			}
			else if( flag == 2 )
			{
				dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				dots[i].yr_diff = dots[i].yr_diff + wide;
			}
			else
			{
				if( flag == 3 )
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
					dots[i].yl_diff = dots[i].yl_diff + wide;
				}
				else if( flag == 4 )
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
					dots[i].yr_diff = dots[i].yr_diff + wide;
				}
			}
		}
		else if( dots[i].sign == 1)
		{
			if( (wide + 5) > width(dots[i].y))
			{
				dots[i].sign = 2;
			}
			else if( flag == 1 )
			{
				dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				dots[i].yr_diff = dots[i].yr_diff + wide;
			}
			else if( flag == 2 )
			{
				dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				dots[i].yl_diff = dots[i].yl_diff + wide;
			}
			else
			{
				if( flag == 3 )
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
					dots[i].yr_diff = dots[i].yr_diff + wide;
				}
				else if( flag == 4 )
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
					dots[i].yl_diff = dots[i].yl_diff + wide;
				}
			}
		}

		if( (width(dots[i].x) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= EFFECTIVE_VALUE))
		{
			dots[i].sign = 2;
		}

		if( (is_td) && (dots[i].sign != DELETED)) {
			dots[i].x = assign_I(dots[i].y.lower - ((width(dots[i].y)*len_x)/len_y), dots[i].y.lower);
		}
	}
	else if(proper_overlap(dots[i].y, to))
	{
		wide = width(intersect(dots[i].y, to));
		if( (wide + 5) > width(dots[i].y))
		{
			dots[i].sign = 2;
		}
		else if( in(dots[i].y.lower, intersect(dots[i].y, to)))
		{
			dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
			dots[i].yl_diff = dots[i].yl_diff + wide;
			flag = 1;
		}
		else if( in(dots[i].y.upper, intersect(dots[i].y, to)))
		{
			dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
			dots[i].yr_diff = dots[i].yr_diff + wide;
			flag = 2;
		}
		else
		{
			if( abs(dots[i].y.lower - to.lower) < abs(dots[i].y.upper - to.upper) )
			{
				temp_val = dots[i].y.lower;
				dots[i].y = assign_I(to.lower+wide, dots[i].y.upper);
				dots[i].yl_diff = dots[i].yl_diff + wide;
				wide = wide + abs(temp_val - to.lower);
				flag = 3;
			}
			else
			{
				temp_val = dots[i].y.upper;
				dots[i].y = assign_I(dots[i].y.lower, to.upper-wide);
				dots[i].yr_diff = dots[i].yr_diff + wide;
				wide = wide + abs(temp_val - to.upper);
				flag = 4;
			}
		}

		if( dots[i].sign == 0 )
		{
			if( (wide + 5) > width(dots[i].x))
			{
				dots[i].sign = 2;
			}
			else if( flag == 1 )
			{
				dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				dots[i].xl_diff = dots[i].xl_diff + wide;
			}
			else if( flag == 2 )
			{
				dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				dots[i].xr_diff = dots[i].xr_diff + wide;
			}
			else
			{
				if( flag == 3 )
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
					dots[i].xl_diff = dots[i].xl_diff + wide;
				}
				else if( flag == 4 )
				{
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
					dots[i].xr_diff = dots[i].xr_diff + wide;
				}
			}
		}
		else if( dots[i].sign == 1)
		{
			if( (wide + 5) > width(dots[i].x))
			{
				dots[i].sign = 2;
			}
			else if( flag == 1 )
			{	
				dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				dots[i].xr_diff = dots[i].xr_diff + wide;
			}
			else if( flag == 2 )
			{
				dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				dots[i].xl_diff = dots[i].xl_diff + wide;
			}
			else
			{
				if( flag == 3 )
				{
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
					dots[i].xr_diff = dots[i].xr_diff + wide;
				}
				else if( flag == 4 )
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
					dots[i].xl_diff = dots[i].xl_diff + wide;
				}
			}
		}

		if( (width(dots[i].y) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= EFFECTIVE_VALUE)) {
			dots[i].sign = 2;
		}

		if( (is_td) && (dots[i].sign != DELETED)) {
			dots[i].y = assign_I(dots[i].x.upper, dots[i].x.upper + ((width(dots[i].x)*len_y)/len_x));
		}
	}
}

int adjust_algn_reg(struct I *reg_x, struct I *reg_y, struct I to, int sign)
{
	int wide = 0;
	struct I temp_x, temp_y;
	int res = 0;
	int flag = 0;
	int temp_val = 0;

	temp_x = (*reg_x);
	temp_y = (*reg_y);

	if(proper_overlap(temp_x, to))
	{
		wide = width(intersect(temp_x, to));
		if( (wide + 5) > width(temp_x))
		{
			res = -1;
		}
		else if( in(temp_x.lower, intersect(temp_x, to)))
		{
			temp_x = assign_I(temp_x.lower+wide, temp_x.upper);
			flag = 1;
		}
		else if( in(temp_x.upper, intersect(temp_x, to)))
		{
			temp_x = assign_I(temp_x.lower, temp_x.upper-wide);
			flag = 2;
		}
		else
		{
			if( abs(temp_x.lower - to.lower) < abs(temp_x.upper - to.upper) )
			{
				temp_val = temp_x.lower;
				temp_x = assign_I(to.lower+wide, temp_x.upper);
				wide = wide + abs(temp_val - to.lower);
				flag = 3;
			}
			else
			{
				temp_val = temp_x.upper;
				temp_x = assign_I(temp_x.lower, to.upper-wide);
				wide = wide + abs(temp_val - to.upper);
				flag = 4;
			}
		}

		if( sign == 0 )
		{
			if( (wide + 5) > width(temp_y))
			{
				res = -1;
			}
			else if( flag == 1 )
			{
				temp_y = assign_I(temp_y.lower+wide, temp_y.upper);
			}
			else if( flag == 2 )
			{
				temp_y = assign_I(temp_y.lower, temp_y.upper-wide);
			}
			else
			{
				if( flag == 3 )
				{
					temp_y = assign_I(temp_y.lower+wide, temp_y.upper);
				}
				else if( flag == 4 )
				{
					temp_y = assign_I(temp_y.lower, temp_y.upper-wide);
				}
			}
		}
		else if( sign == 1)
		{
			if( (wide + 5) > width(temp_y))
			{
				res = -1;
			}
			else if( flag == 1 )
			{
				temp_y = assign_I(temp_y.lower, temp_y.upper-wide);
			}
			else if( flag == 2 )
			{
				temp_y = assign_I(temp_y.lower+wide, temp_y.upper);
			}
			else
			{
				if( flag == 3 )
				{
					temp_y = assign_I(temp_y.lower, temp_y.upper-wide);
				}
				else if( flag == 4 )
				{
					temp_y = assign_I(temp_y.lower+wide, temp_y.upper);
				}
			}
		}

		if( (width(temp_x) <= EFFECTIVE_VALUE) || (abs(temp_x.lower - temp_y.lower) <= EFFECTIVE_VALUE))
		{
			res = -1;
		}
	}
	else if(proper_overlap(temp_y, to))
	{
		wide = width(intersect(temp_y, to));
		if( (wide + 5) > width(temp_y))
		{
			res = -1;
		}
		else if( in(temp_y.lower, intersect(temp_y, to)))
		{
			temp_y = assign_I(temp_y.lower+wide, temp_y.upper);
			flag = 1;
		}
		else if( in(temp_y.upper, intersect(temp_y, to)))
		{
			temp_y = assign_I(temp_y.lower, temp_y.upper-wide);
			flag = 2;
		}
		else
		{
			if( abs(temp_y.lower - to.lower) < abs(temp_y.upper - to.upper) )
			{
				temp_val = temp_y.lower;
				temp_y = assign_I(to.lower+wide, temp_y.upper);
				wide = wide + abs(temp_val - to.lower);
				flag = 3;
			}
			else
			{
				temp_val = temp_y.upper;
				temp_y = assign_I(temp_y.lower, to.upper-wide);
				wide = wide + abs(temp_val - to.upper);
				flag = 4;
			}
		}

		if( sign == 0 )
		{
			if( (wide + 5) > width(temp_x))
			{
				res = -1;
			}
			else if( flag == 1 )
			{
				temp_x = assign_I(temp_x.lower+wide, temp_x.upper);
			}
			else if( flag == 2 )
			{
				temp_x = assign_I(temp_x.lower, temp_x.upper-wide);
			}
			else
			{
				if( flag == 3 )
				{
					temp_x = assign_I(temp_x.lower+wide, temp_x.upper);
				}
				else if( flag == 4 )
				{
					temp_x = assign_I(temp_x.lower, temp_x.upper-wide);
				}
			}
		}
		else if( sign == 1 )
		{
			if( (wide + 5) > width(temp_x))
			{
				res = -1;
			}
			else if( flag == 1 )
			{	
				temp_x = assign_I(temp_x.lower, temp_x.upper-wide);
			}
			else if( flag == 2 )
			{
				temp_x = assign_I(temp_x.lower+wide, temp_x.upper);
			}
			else
			{
				if( flag == 3 )
				{
					temp_x = assign_I(temp_x.lower, temp_x.upper-wide);
				}
				else if( flag == 4 )
				{
					temp_x = assign_I(temp_x.lower+wide, temp_x.upper);
				}
			}
		}

		if( (width(temp_y) <= EFFECTIVE_VALUE) || (abs(temp_x.lower - temp_y.lower) <= EFFECTIVE_VALUE)) {
			 res = -1;
		}
	}

	*reg_x = temp_x;
	*reg_y = temp_y;
	return(res);
}
