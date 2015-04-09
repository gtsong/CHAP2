#include "main.h"
#include "regions.h"
#include "check_overlap.h"
#include "util_i.h"

struct I change_value(struct I temp, int distance, int sign, bool is_x)
{
	struct I res;

	if( sign == 0)
	{
		if( is_x == true)
		{
			res = assign_I(temp.lower + distance, temp.upper + distance);		
		}
		else 
		{
			res = assign_I(temp.lower - distance, temp.upper - distance);		
		}
	}
	else 
	{
		res = assign_I((distance - temp.upper), (distance - temp.lower));
	}
	return res;
}

bool find_overlapped_region(struct I *temp_ptr, struct I *pair_ptr, struct I cur_reg, int id, struct DotList *dots, int sign, bool is_x, int distance) 
{
	int difference;
	struct I temp, pair;	

	temp = assign_I((*temp_ptr).lower, (*temp_ptr).upper);
	pair = assign_I((*pair_ptr).lower, (*pair_ptr).upper);

	if( in(cur_reg.upper, temp) )
	{
		temp = assign_I(temp.lower, cur_reg.upper);
		difference = width(temp);
		if( dots[id].sign == 0)
		{
			if(sign == 0)
			{
				pair = assign_I(pair.lower, pair.lower + difference);
			}
			else 
			{
				pair = assign_I((pair.upper - difference), pair.upper);
			}
		}
		else
		{
			if(sign == 0)
			{
				pair = assign_I((pair.upper - difference), pair.upper);
			}
			else 
			{
				pair = assign_I(pair.lower, pair.lower + difference);
			}
		}
	}
	else if( in(cur_reg.lower, temp))
	{
		temp = assign_I(cur_reg.lower, temp.upper);
		difference = width(temp);
		if( dots[id].sign == 0)
		{
			if(sign == 0)
			{
				pair = assign_I((pair.upper - difference), pair.upper);
			}
			else 
			{
				pair = assign_I(pair.lower, pair.lower + difference);
			}
		}
		else
		{
			if(sign == 0)
			{
				pair = assign_I(pair.lower, pair.lower + difference);
			}
			else 
			{
				pair = assign_I((pair.upper - difference), pair.upper);
			}
		}
	}
	else
	{
		return false;
	}

	temp = change_value(temp, distance, dots[id].sign, is_x);

	*temp_ptr = assign_I(temp.lower, temp.upper);
	*pair_ptr = assign_I(pair.lower, pair.upper);
	return true;
}
