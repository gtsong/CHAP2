#include "main.h"
#include "regions.h"
#include "rollback.h"
#include "pred_dels.h"
#include "pred_ops.h"
#include "adjust_algn.h"
#include "find_merging.h"
#include "util_i.h"

extern int debug_mode;

void shift_regions(int con, struct I young, int *num, struct DotList *dots)
{
	int i;
	int num_loop;

	num_loop = *num;

	for(i = 0; i < num_loop; i++)
	{
		if( dots[i].sign == DELETED ) {}
		else if( dots[i].sign == 10 )
		{
			dots[i].sign = 0;
		}
		else if( dots[i].sign == 11 )
		{
			dots[i].sign = 1;
		}
		else if( con == 0 )
		{
			if( dots[i].x.lower >= young.lower )
			{
				dots[i].x = assign_I(dots[i].x.lower - width(young), dots[i].x.upper - width(young));
				if( is_tandem(dots[i]) ) dots[i].y = assign_I(dots[i].y.lower - width(young), dots[i].y.upper - width(young));
			}

			if( (!is_tandem(dots[i])) && (dots[i].y.lower >= young.lower) )
			{
				dots[i].y = assign_I(dots[i].y.lower - width(young), dots[i].y.upper - width(young));
			}
			
			if( abs(dots[i].x.lower - dots[i].y.lower) <= EFFECTIVE_VALUE )
			{
				dots[i].sign = DELETED;
			}
		}
	}
}

void shift_regions_pair(bool is_x, int con, struct I young, int *num, struct DotList *dots)
{
	int i;
	int num_loop;

	num_loop = *num;

	for(i = 0; i < num_loop; i++)
	{
		if( dots[i].sign == DELETED ) {}
		else if( dots[i].sign == 10 )
		{
			dots[i].sign = 0;
		}
		else if( dots[i].sign == 11 )
		{
			dots[i].sign = 1;
		}
		else if( con == 0 )
		{
			if( is_x == true )
			{
				if( dots[i].x.lower >= young.lower )
				{
					dots[i].x = assign_I(dots[i].x.lower - width(young), dots[i].x.upper - width(young));
				}
			}
			else
			{
				if( dots[i].y.lower >= young.lower )
				{
					dots[i].y = assign_I(dots[i].y.lower - width(young), dots[i].y.upper - width(young));
				}
			}
		}
	}
}

int rollback_step_dup_no_overlap(bool is_x, int id, int *num, struct DotList *dots)
{
	int i;
	struct I to;
	int num_lines;
	bool is_td;

	is_td = is_tandem(dots[id]);

	if(is_x) 
	{
		to = assign_I(dots[id].y.lower, dots[id].y.upper);
	}
	else
	{
		to = assign_I(dots[id].x.lower, dots[id].x.upper);
	}
	
	num_lines = *num;
	if( debug_mode == TRUE ) printf("Removed: %d-%d: %d\n", to.lower, to.upper, dots[id].identity);
	for( i = 0 ; i < num_lines; i++ ) // Excluding the test for split region newly
	{
		if( (dots[i].sign == 20) || (dots[i].sign == 21) )
		{
			if( dots[i].sign == 20 ) dots[i].sign = 0;
			else if( dots[i].sign == 21 ) dots[i].sign = 1;
		}
		else if( (i != id)  && (is_td) && (is_tandem(dots[i]) == true)) {
			if( f_loose_subset(dots[i].x, to, STRICT) == true ) {
			}
			else if( f_loose_subset(dots[i].y, to, STRICT) == true ) {
			}
			else if( (proper_overlap(dots[i].x, to) == true ) || (proper_overlap(dots[i].y, to) == true )) {
				adjust_alignment(dots, i, to);
			}
		}
		else if(f_loose_subset(dots[i].x, to, STRICT) || f_loose_subset(dots[i].y, to, STRICT))
		{
			dots[i].sign = DELETED; // The alignment is deleted
		}
		else if( (proper_overlap(dots[i].x, to) == true ) || (proper_overlap(dots[i].y, to) == true ))
		{
			adjust_alignment(dots, i, to);
		}
		else
		{
		}
	}

	shift_regions(0, to, num, dots);
	overwrite_dots(num, dots);
	return(0);
}

void chop_off_self(bool x_or_y, int id, int num, struct DotList *dots)
{
	int i;
	struct I to;
	int num_lines;
	int wide;

	if(x_or_y)  // if x_or_y is true, the region of y is chopped off
	{
		to = assign_I(dots[id].y.lower, dots[id].y.upper);
	}
	else
	{
		to = assign_I(dots[id].x.lower, dots[id].x.upper);
	}
	
	num_lines = num;
	for( i = 0 ; i < num_lines; i++ ) // Excluding the test for split region newly
	{
		if(loose_subset(dots[i].x, to) || loose_subset(dots[i].y, to))
		{
			dots[i].sign = DELETED; // The line is deleted
		}
		else if(proper_overlap(dots[i].x, to))
		{
			wide = width(intersect(dots[i].x, to));
			if( in(dots[i].x.lower, intersect(dots[i].x, to)))
			{
				dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
			}
			else
			{
				dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
			}
			if( dots[i].sign == 0 )
			{
				if( in(dots[i].x.lower, intersect(dots[i].x, to)))
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				}
				else
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				}
			}
			else if( dots[i].sign == 1)
			{
				if( in(dots[i].x.lower, intersect(dots[i].x, to)))
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				}
				else
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				}
			}

			if( (width(dots[i].x) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= 0))
			{
				dots[i].sign = DELETED;
			}
		}
		else if(proper_overlap(dots[i].y, to))
		{
			wide = width(intersect(dots[i].y, to));
			if( in(dots[i].y.lower, intersect(dots[i].y, to)))
			{
				dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
			}
			else
			{
				dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
			}
			if( dots[i].sign == 0 )
			{
				if( in(dots[i].y.lower, intersect(dots[i].y, to)))
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				}
				else
				{
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				}
			}
			else if( dots[i].sign == 1)
			{
				if( in(dots[i].y.lower, intersect(dots[i].y, to)))
				{	
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				}
				else
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				}
			}

			if( (width(dots[i].y) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= 0)) {
				dots[i].sign = DELETED;
			}
		}
		else
		{
		}
	}
}

int rollback_step_conversion(bool is_x, int id, int *num, struct DotList *dots)
{
	struct I to;
	int num_lines;

	if(is_x) 
	{
		to = assign_I(dots[id].y.lower, dots[id].y.upper);
	}
	else
	{
		to = assign_I(dots[id].x.lower, dots[id].x.upper);
	}
	
	int i;
	num_lines = *num;
	for( i = 0 ; i < num_lines; i++ ) // Excluding the test for split region newly
	{
		if(subset(dots[i].x, to) || subset(dots[i].y, to))
		{
			dots[i].sign = DELETED; // The line is deleted
		}
	}

	shift_regions(width(to), to, num, dots);
	overwrite_dots(num, dots);
	return(width(to));
}

struct I redefine_del_list(int *extra, int *cur_count, int *del_loc, int num_del, struct gap_list *del_list, int count, struct DotList *dots)
{
	int i;
	int min_len = MDIS_THRESHOLD;
	struct I res;
	int cur_num = 0;
	int mid_loc;
	int temp_extra;

	for( i = 0; i < num_del; i++ )
	{
		if( del_list[i].gid == count )
		{
			if( abs(del_list[i].y2 - del_list[i].y1) != 0 )
			{
				if( min_len > abs(del_list[i].y2 - del_list[i].y1) )
				{
					min_len = abs(del_list[i].y2 - del_list[i].y1);
					if( del_list[i].y1 > del_list[i].y2)
					{		
						res = assign_I(del_list[i].y2, del_list[i].y1);
						if( del_list[i].x1 <= del_list[i].x2 )
						{
							*del_loc = del_list[i].x1;
						}
						else	
						{
							*del_loc = del_list[i].x2;
						}
						mid_loc = (del_list[i].x1 + del_list[i].x2)/2;
					}
					else
					{
						res = assign_I(del_list[i].y1, del_list[i].y2);
						if( del_list[i].x1 <= del_list[i].x2 )
						{
							*del_loc = del_list[i].x1;
						}
						else	
						{
							*del_loc = del_list[i].x2;
						}
						mid_loc = (del_list[i].x1 + del_list[i].x2)/2;
					}

					if( del_list[i].type == 1 )
					{
						if( proper_overlap(dots[del_list[i].id1].x, dots[del_list[i].id2].x))
						{
							temp_extra = width(intersect(dots[del_list[i].id1].x, dots[del_list[i].id2].x));	
						}
						else temp_extra = 0;
					}
					else if( del_list[i].type == 2)
					{
						if( proper_overlap(dots[del_list[i].id1].y, dots[del_list[i].id2].y))
						{
							temp_extra = width(intersect(dots[del_list[i].id1].y, dots[del_list[i].id2].y));	
						}
						else temp_extra = 0;
					}
				}
			}
			else del_list[i].gid = -2;
		}
	}

	for( i = 0; i < num_del; i++ )
	{
		if( del_list[i].gid == count )
		{
			if( abs(del_list[i].y2 - del_list[i].y1) > (min_len + DEL_DIF_TH) )
			{
				del_list[i].gid = -2; // remove from the deletion list
			}
			else if( in(mid_loc, assign_I(del_list[i].x1, del_list[i].x2)) != true )
			{
				del_list[i].gid = -2; // remove from the deletion list
			}
			else cur_num++;
		}
	}

	*extra = temp_extra;
	*cur_count = cur_num;
	return(res);
}

void chop_off(int mode, struct I seg, int num, struct DotList *dots)
{
	int i;
	struct I temp;
	int num_lines;
	int wide;

	num_lines = num;
	for( i = 0 ; i < num_lines; i++ ) // Excluding the test for split region newly
	{
		if( mode == PAIR_1 )
		{
			temp = assign_I(dots[i].x.lower, dots[i].x.upper);
		}
		else if( mode == PAIR_2 )
		{
			temp = assign_I(dots[i].y.lower, dots[i].y.upper);
		}

		if(loose_subset(temp, seg) == true)
		{
			dots[i].sign = 2; // The line is deleted
		}
		else if((mode == PAIR_1) && (proper_overlap(dots[i].x, seg)))
		{
			wide = width(intersect(dots[i].x, seg));
			if( in(dots[i].x.lower, intersect(dots[i].x, seg)))
			{
				dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
			}
			else
			{
				dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
			}
			if( dots[i].sign == 0 )
			{
				if( in(dots[i].x.lower, intersect(dots[i].x, seg)))
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				}
				else
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				}
			}
			else if( dots[i].sign == 1)
			{
				if( in(dots[i].x.lower, intersect(dots[i].x, seg)))
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				}
				else
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				}
			}

			if( (width(dots[i].x) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= 0) )
			{
				dots[i].sign = DELETED;
			}
		}
		else if((mode == PAIR_2) && (proper_overlap(dots[i].y, seg)))
		{
			wide = width(intersect(dots[i].y, seg));
			if( in(dots[i].y.lower, intersect(dots[i].y, seg)))
			{
				dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
			}
			else
			{
				dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
			}
			if( dots[i].sign == 0 )
			{
				if( in(dots[i].y.lower, intersect(dots[i].y, seg)))
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				}
				else
				{
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				}
			}
			else if( dots[i].sign == 1)
			{
				if( in(dots[i].y.lower, intersect(dots[i].y, seg)))
				{	
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				}
				else
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				}
			}

			if( (width(dots[i].x) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= 0) )
			{
				dots[i].sign = DELETED;
			}
		}
	}
}

int rollback_step_pairwise(bool is_x, struct I reg, int *num, struct DotList *dots)
{
	int i;
	int num_lines;
	int wide;
	struct I cur;

	num_lines = *num;
	for( i = 0 ; i < num_lines; i++ ) // Excluding the test for split region newly
	{
		if( is_x == true )
		{
			cur = assign_I(dots[i].x.lower, dots[i].x.upper);
		}
		else
		{
			cur = assign_I(dots[i].y.lower, dots[i].y.upper);
		}

		if(loose_subset(cur, reg))
		{
			dots[i].sign = DELETED; // The line is deleted
		}
		else if((is_x == true) && (proper_overlap(dots[i].x, reg) == true))
		{
			wide = width(intersect(dots[i].x, reg));

			if( (wide+5) >= width(dots[i].x) )
			{
				dots[i].sign = DELETED;
			}
			else if( in(dots[i].x.lower, intersect(dots[i].x, reg)))
			{
				dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
				dots[i].left_diff = dots[i].left_diff + wide;
			}
			else if( in(dots[i].x.upper, intersect(dots[i].x, reg)))
			{
				dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
				dots[i].right_diff = dots[i].right_diff + wide;
			}
			else
			{
				if( abs(dots[i].x.lower - reg.lower) < abs(dots[i].x.upper - reg.upper) )
				{
					dots[i].x = assign_I(reg.lower+wide, dots[i].x.upper);
					wide = wide + abs(dots[i].x.lower - reg.lower);
					dots[i].left_diff = dots[i].left_diff + wide;
				}
				else
				{
					dots[i].x = assign_I(dots[i].x.lower, reg.upper-wide);
					wide = wide + abs(dots[i].x.upper - reg.upper);
					dots[i].right_diff = dots[i].right_diff + wide;
				}	
			}

			if( dots[i].sign == 0 )
			{
				if( (wide+5) >= width(dots[i].y) )
				{
					dots[i].sign = DELETED;
				}
				else if( in(dots[i].x.lower, intersect(dots[i].x, reg)))
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				}
				else if( in(dots[i].x.upper, intersect(dots[i].x, reg)))
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				}
				else 
				{
					if( abs(dots[i].x.lower - reg.lower) < abs(dots[i].x.upper - reg.upper) )
					{
						dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
					}
					else
					{
						dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
					}	
				}
			}
			else if( dots[i].sign == 1)
			{
				if( (wide+5) >= width(dots[i].y) )
				{
					dots[i].sign = DELETED;
				}
				else if( in(dots[i].x.lower, intersect(dots[i].x, reg)))
				{
					dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
				}
				else if( in(dots[i].x.upper, intersect(dots[i].x, reg)))
				{
					dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
				}
				else 
				{
					if( abs(dots[i].x.lower - reg.lower) < abs(dots[i].x.upper - reg.upper) )
					{
						dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
					}
					else
					{
						dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
					}	
				}
			}

			if( width(dots[i].x) <= EFFECTIVE_VALUE )
			{
				dots[i].sign = DELETED;
			}
		}
		else if((is_x == false) && (proper_overlap(dots[i].y, reg) == true))
		{
			wide = width(intersect(dots[i].y, reg));
			if( (wide+5) >= width(dots[i].y) )
			{
				dots[i].sign = DELETED;
			}
			else if( in(dots[i].y.lower, intersect(dots[i].y, reg)))
			{
				dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
			}
			else if( in(dots[i].y.upper, intersect(dots[i].y, reg)))
			{
				dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
			}
			else 
			{
				if( abs(dots[i].y.lower - reg.lower) < abs(dots[i].y.upper - reg.upper) )
				{
					dots[i].y = assign_I(reg.lower+wide, dots[i].y.lower);
					wide = wide + abs(dots[i].y.lower - reg.lower);
				}
				else
				{
					dots[i].y = assign_I(dots[i].y.lower, reg.upper-wide);
					wide = wide + abs(dots[i].y.upper - reg.upper);
				}	
			}

			if( dots[i].sign == 0 )
			{
				if( (wide+5) >= width(dots[i].x))
				{
					dots[i].sign = DELETED;
				}
				else if( in(dots[i].y.lower, intersect(dots[i].y, reg)))
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
					dots[i].left_diff = dots[i].left_diff + wide;
				}
				else if( in(dots[i].y.upper, intersect(dots[i].y, reg)))
				{
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
					dots[i].right_diff = dots[i].right_diff + wide;
				}
				else 
				{
					if( abs(dots[i].y.lower - reg.lower) < abs(dots[i].y.upper - reg.upper) )
					{
						dots[i].y = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
						dots[i].left_diff = dots[i].left_diff + wide;
					}
					else
					{
						dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
						dots[i].right_diff = dots[i].right_diff + wide;
					}	
				}
			}
			else if( dots[i].sign == 1)
			{
				if( (wide+5) >= width(dots[i].x) )
				{
					dots[i].sign = DELETED;
				}
				else if( in(dots[i].y.lower, intersect(dots[i].y, reg)))
				{	
					dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
					dots[i].right_diff = dots[i].right_diff + wide;
				}
				else if( in(dots[i].y.upper, intersect(dots[i].y, reg)))
				{
					dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
					dots[i].left_diff = dots[i].left_diff + wide;
				}
				else 
				{
					if( abs(dots[i].y.lower - reg.lower) < abs(dots[i].y.upper - reg.upper) )
					{
						dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
						dots[i].right_diff = dots[i].right_diff + wide;
					}
					else
					{
						dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
						dots[i].left_diff = dots[i].left_diff + wide;
					}	
				}
			}

			if( width(dots[i].y) <= EFFECTIVE_VALUE)
			{
				dots[i].sign = DELETED;
			}
		}
		else
		{
		}
	}

	shift_regions_pair(is_x, 0, reg, num, dots);
	overwrite_dots(num, dots);
	return(0);
}

int rollback_step_dup_overlap(bool is_x, int id, int *num, struct DotList *dots)
{  
	int i;  
	struct I to;  
	int num_lines;  
	int wide;  
	
	if(is_x)  
	{    
		to = assign_I(dots[id].y.lower, dots[id].y.upper);
	}  
	else  
	{    
		to = assign_I(dots[id].x.lower, dots[id].x.upper);
	}

  num_lines = *num;
	for( i = 0 ; i < num_lines; i++ ) // Excluding the test for split region newly  
	{    
		if(f_loose_subset(dots[i].x, to, STRICT) || f_loose_subset(dots[i].y, to, STRICT)) 
		{
			dots[i].sign = DELETED; // The line is deleted    
		} 
    else if(proper_overlap(dots[i].x, to))
    {
      wide = width(intersect(dots[i].x, to));
      if( in(dots[i].x.lower, intersect(dots[i].x, to)))
      {
        dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
      }
      else
      {
        dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
      }

 	    if( dots[i].sign == 0 )
      {
        if( in(dots[i].x.lower, intersect(dots[i].x, to)))
        {
          dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
        }
        else
        {
          dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
        }
      }
 	 		else if( dots[i].sign == 1)
      {
        if( in(dots[i].x.lower, intersect(dots[i].x, to)))
        {
          dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
        }
        else
        {
          dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
        }
      }

      if( (width(dots[i].x) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= EFFECTIVE_VALUE))
      {
        dots[i].sign = DELETED;
      }
    }
    else if(proper_overlap(dots[i].y, to))
   	{
  		wide = width(intersect(dots[i].y, to));
     	if( in(dots[i].y.lower, intersect(dots[i].y, to)))
     	{
       	dots[i].y = assign_I(dots[i].y.lower+wide, dots[i].y.upper);
     	}
     	else
     	{
       	dots[i].y = assign_I(dots[i].y.lower, dots[i].y.upper-wide);
     	}

     	if( dots[i].sign == 0 )
     	{
       	if( in(dots[i].y.lower, intersect(dots[i].y, to)))
       	{
         	dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
      	}
       	else
       	{
         	dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
       	}
     	}
     	else if( dots[i].sign == 1)
     	{
       	if( in(dots[i].y.lower, intersect(dots[i].y, to)))
       	{
         	dots[i].x = assign_I(dots[i].x.lower, dots[i].x.upper-wide);
       	}
       	else
       	{
        	dots[i].x = assign_I(dots[i].x.lower+wide, dots[i].x.upper);
       	}
     	}

      if( (width(dots[i].y) <= EFFECTIVE_VALUE) || (abs(dots[i].x.lower - dots[i].y.lower) <= EFFECTIVE_VALUE)) 
			{
        dots[i].sign = DELETED;
      }
    }
    else{}
	}

	shift_regions(0, to, num, dots);
	overwrite_dots(num, dots);
	return(0);
}

int rollback_step_del(int cur, int num_del, struct gap_list *del_list, int *num, struct DotList *dots, int num_ops, struct ops_list *ops)
{
  int i;
	int max_width = -1;
	int cut_point;
	int b, e;
	int cur_gid;

	max_width = abs(del_list[cur].y2 - del_list[cur].y1);
	cut_point = (del_list[cur].x2 + del_list[cur].x1)/2;

	for( i = 0; i < num_del; i++ )
	{
		if( i == cur ) {}
		else if( del_list[i].gid == del_list[cur].gid )
		{
			if( max_width < abs(del_list[i].y2 - del_list[i].y1) )
			{
				max_width = abs(del_list[i].y2 - del_list[i].y1);
				cut_point = (del_list[i].x2 + del_list[cur].x1)/2;
			}
		}
	}

	for(i = 0; i < *num; i++)
  {
		if( dots[i].sign == DELETED ) {}
		else if( check_inv_same_del(del_list, num_del, i, del_list[cur].gid) == true ) {}
		else if( ( i == del_list[cur].id1 ) || ( i == del_list[cur].id2 ) ) {}
		else
		{
			if( dots[i].x.lower >= cut_point )
 		  {
				dots[i].x = assign_I(dots[i].x.lower + max_width, dots[i].x.upper + max_width);
			}

	   	if( dots[i].y.lower >= cut_point )
 	  	{
				dots[i].y = assign_I(dots[i].y.lower + max_width, dots[i].y.upper + max_width);
			}
		}
	}

	cur_gid = del_list[cur].gid;

	for(i = 0; i < num_del; i++ )
	{
		if( (dots[del_list[i].id1].sign != DELETED) && (dots[del_list[i].id2].sign != DELETED))
		{
			if( del_list[i].gid == cur_gid )
			{
				if( debug_mode == TRUE ) {
					printf("Merged: %d-%d, %d-%d and %d-%d, %d-%d\n", dots[del_list[i].id1].x.lower, dots[del_list[i].id1].x.upper, dots[del_list[i].id1].y.lower, dots[del_list[i].id1].y.upper, dots[del_list[i].id2].x.lower, dots[del_list[i].id2].x.upper, dots[del_list[i].id2].y.lower, dots[del_list[i].id2].y.upper);
				}

				del_list[i].gid = -1; 

				if( dots[del_list[i].id1].x.lower < dots[del_list[i].id2].x.lower )
				{
					b = dots[del_list[i].id1].x.lower;
				}
				else
				{
					b = dots[del_list[i].id2].x.lower;
				}

				if( dots[del_list[i].id1].x.upper < dots[del_list[i].id2].x.upper )
				{
					e = dots[del_list[i].id2].x.upper;
				}
				else
				{
					e = dots[del_list[i].id1].x.upper;
				}
				
				if( del_list[i].type == 1 ) 
				{
					e = e + max_width;
				}
				dots[del_list[i].id1].x = assign_I(b, e);
				
				if( dots[del_list[i].id1].y.lower < dots[del_list[i].id2].y.lower )
				{
					b = dots[del_list[i].id1].y.lower;
				}
				else
				{
					b = dots[del_list[i].id2].y.lower;
				}

				if( dots[del_list[i].id1].y.upper < dots[del_list[i].id2].y.upper )
				{
					e = dots[del_list[i].id2].y.upper;
				}
				else
				{
					e = dots[del_list[i].id1].y.upper;
				}

				if( del_list[i].type == 2 )
				{
					e = e + max_width;
				}
				else if( b >= cut_point ) 
				{
					b = b + max_width;
					e = e + max_width;
				}

				dots[del_list[i].id1].y = assign_I(b, e);

				dots[del_list[i].id2].sign = DELETED;
			}
		}
	}
	del_list[cur].gid = -1; // initialization of the current item of the deletion list

  overwrite_dots(num, dots);
	generate_ops_del('d', cut_point, max_width, num_ops, ops);

	return(max_width);
}
