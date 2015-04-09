#include "main.h"
#include "regions.h"
#include "util_i.h"

bool strict_subset(struct I reg1, struct I reg2)
{
  struct I temp;
  int len;

  if( width(reg1) > width(reg2) ) len = width(reg2);
  else len = width(reg1);

  if( len <= MIN_LEN )
  {
    temp = assign_I(reg2.lower - ERR_SM_TH, reg2.upper + ERR_SM_TH );
  }
  else if( len <= DIS_THRESHOLD )
  {
    temp = assign_I(reg2.lower - ERR_LG_TH, reg2.upper + ERR_LG_TH );
  }
  else
  {
    temp = assign_I(reg2.lower - ERR_TH, reg2.upper + ERR_TH );
  }

  return(subset(reg1, temp));
}

bool proper_in(int r, struct I reg)
{
	if(in(r, reg)) 
	{
		if(r == reg.lower || r == reg.upper) return false;
		else return true;
	}
	else return false;
}

bool proper_overlap(struct I reg1, struct I reg2)
{
	if(overlap(reg1, reg2)) 
	{
		if(width(intersect(reg1, reg2)) == 0) return false;
		else return true;
	}
	else return false;
}

bool almost_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	temp = assign_I((reg2.lower) - THRESHOLD, (reg2.upper) + THRESHOLD);
	if( (width(reg1) > 5*THRESHOLD) && (width(reg2) > 5*THRESHOLD) )
	{
		if(subset(reg1, temp)) return(true);
		else if( proper_overlap(reg1, reg2) )
		{
			if(width(intersect(reg1, reg2)) > ((0.7)*width(reg1))) return(true);
			else return(false);				
		}
		else return(false);
	}
	else return(subset(reg1, temp));
}

bool d_tight_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	if( width(reg1) >= TIGHT_BASE_TH )
	{
		temp = assign_I((reg2.lower - TIGHT_TH), (reg2.upper + TIGHT_TH));
	}
	else if( width(reg1) >= (TIGHT_BASE_TH/4) )
	{
		temp = assign_I((reg2.lower - (TIGHT_TH/2)), (reg2.upper + (TIGHT_TH/2)));
	}
	else
	{
		temp = assign_I((reg2.lower - (TIGHT_TH/6)), (reg2.upper + (TIGHT_TH/6)));
	}

	return(subset(reg1, temp));
}

bool s_tight_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	temp = assign_I((reg2.lower - (TIGHT_TH/20)), (reg2.upper + (TIGHT_TH/20)));

	return(subset(reg1, temp));
}

bool tight_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	temp = assign_I((reg2.lower - TIGHT_TH), (reg2.upper + TIGHT_TH));

	return(subset(reg1, temp));
}

bool loosen_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	if( width(reg1) > (DIS_THRESHOLD) ) 
	{
		temp = assign_I((reg2.lower) - 4*THRESHOLD, (reg2.upper) + 4*THRESHOLD);
	}
	else temp = assign_I((reg2.lower) - 2*THRESHOLD, (reg2.upper) + 2*THRESHOLD);
	return(subset(reg1, temp));
}

bool too_loosen_subset(struct I reg1, struct I reg2)
{
	struct I temp;

	if( width(reg1) > LOOSEN_T ) 
	{
		temp = assign_I((reg2.lower) - (LOOSEN_T/2), (reg2.upper) + (LOOSEN_T/2));
	}
	else
	{
		if( width(reg1) > (LOOSEN_T/2) )
		{	
			temp = assign_I((reg2.lower) - (LOOSEN_T/4), (reg2.upper) + (LOOSEN_T/4));
		}
		else
		{
			temp = assign_I((reg2.lower) - (LOOSEN_T/8), (reg2.upper) + (LOOSEN_T/8));
		}
	}
	return(subset(reg1, temp));
}

void init_array(int *array, int num)
{
	int i;

	for( i = 0 ; i < num; i++ ) array[i] = 0;
}

void cut_off(int *num, struct DotList *dots)
{
  int i;
	int j = 0;

  for(i = 0 ; i < *num; i++)
  {
    if( dots[i].sign == 0 || dots[i].sign == 1 )
    {
			if( dots[i].identity >= P_IDT )
			{
				dots[j].x = assign_I(dots[i].x.lower, dots[i].x.upper);
				dots[j].y = assign_I(dots[i].y.lower, dots[i].y.upper);
				dots[j].sign = dots[i].sign;
				dots[j].identity = dots[i].identity;
				j++;
			}
    }
  }
  *num = j; // Modify the number of lines
}

void overwrite_dots(int *num, struct DotList *dots)
{
  int i;
	int j = 0;

  for(i = 0 ; i < *num; i++)
  {
    if( dots[i].sign == 0 || dots[i].sign == 1 )
    {
			dots[j].x = assign_I(dots[i].x.lower, dots[i].x.upper);
			dots[j].y = assign_I(dots[i].y.lower, dots[i].y.upper);
	    dots[j].sign = dots[i].sign;
	    dots[j].identity = dots[i].identity;
	    dots[j].l_pid = dots[i].l_pid;
      dots[j].l_id = dots[i].l_id;
      dots[j].lock = -1;
      dots[j].c_id = dots[i].c_id;
      dots[j].fid = dots[i].fid;
      dots[j].left_diff = dots[i].left_diff;
      dots[j].right_diff = dots[i].right_diff;
      j++;
    }
  }
  *num = j; // Modify the number of lines
}

int compute_distance(struct I x1, struct I y1, struct I x2, struct I y2, int sign)
{
  int c, d;

  if( sign == 0 )
  {
    c = y2.lower - x2.lower;
    d = abs(x1.lower - y1.lower + c);
  }
  else
  {
    c = (-1)*(x2.lower + y2.upper);
    d = abs(x1.lower + y1.upper + c);
  }

  return d;
}
