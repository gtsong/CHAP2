#include "main.h"
#include "regions.h"
#include "adjust_plot.h"
#include "find_merging.h"
#include "check_repeats.h"
#include "util_i.h"

int count_lower(char *buf, int len)
{
	int i;
	int res = 0;

	for( i = 0; i < len; i++ )
	{
		if(islower(buf[i])) res++;
	}

	return res;
}

float compute_portion_rps(int from, int to, struct r_list *rp, int n_rps)
{
	int start, end;
	int sum_nc_rps;
	int size_gap;
	float f;

	if( from < rp[0].start ) start = from;
	else start = rp[0].start;

	if( to > rp[n_rps-1].end ) end = to;
	else end = rp[n_rps-1].end;
	
	sum_nc_rps = sum_len_rps(rp, n_rps);
	size_gap = end - start + 1;
	f = (float)(((float) sum_nc_rps) / ((float) (size_gap)));

	return(f);
}

int sum_len_rps(struct r_list *rp, int n_rps)
{
	int i;
	int sum = 0;

	for( i = 0; i < n_rps; i++ )
	{
		sum = sum + abs(rp[i].end - rp[i].start + 1); 
	}

	return(sum);
}

/* if there exists no repeat in a gap
	 return -1
	 otherwise return the number of repeats with the list of repeats
	 *cr is a ratio to express how correspondent the boundary of an interspersed repeat and and a gap are
*/
int obtain_repeats_list(int from, int to, struct r_list *rp, struct r_list *rp_list, int num_list)
{
	struct I i_rp;
	int i = 0, j = 0;
	int a = 0, b = 0;
	int id = 0;
	float rate = (float) 0;
	int within_b = -1;
	bool is_over = false;
	int pre_end = from; 

	i_rp = assign_I(0, 1);
	i = 0;
	while((i < num_list) && (is_over == false))
	{
		rate = rp_list[i].d_rate;
		a = rp_list[i].start;
		b = rp_list[i].end;
		id = rp_list[i].id;
		i_rp = assign_I(a,b);

		if((!in(to, i_rp)) && (to >  b))
		{
			if( width(i_rp) <= S_RP_BD ) {}
			else if( in(from, i_rp) ) 
			{
				if( width(i_rp) < (5*RP_BD) )
				{
					if(((float)(((float)(abs(b - from))) / ((float)width(i_rp))) >= RR_C))
					{
						within_b = 1;
					}
					else if((float)(((float)(abs(b - from))) / ((float)width(i_rp))) <= RR_NC)
					{
						within_b = 0;
					}
					else
					{
						within_b = -1;
					}
				}
				else
				{
					if(abs(from - a) <= RP_BD)
					{
						within_b = 1;
					}
					else if(abs(b - from) <= RP_BD)
					{
						within_b = 0;
					}
					else
					{
						within_b = -1;
					}
				}

				if( within_b == 1 )
				{
					if( abs(a - pre_end) > (3*S_RP_BD) )
					{
						return(-1);
					}
					else pre_end = b;

					if( (j >= 1) && (rate == 0 ) ) {
						// Simple_repeat or Low_complexity
					}
					else {
						rp[j].start = a;
						rp[j].end = b;
						rp[j].d_rate = rate;
						rp[j].id = id;
						j++;
					}

					if( j >= Max_Rps ) 
					{
						return(-1);
					}			
				}
				else if( within_b == -1 )
				{
					return(-1);
				}
			}
			else if( overlap(assign_I(from, to), i_rp ))
			{
				if( j == 0 )
				{
					if( width(i_rp) < (5*RP_BD) )
					{
						if((float)(((float)(abs(a - from))) / ((float)width(i_rp))) <= RR_NC)
						{
							within_b = 1;
						}
						else
						{
							within_b = -1;
						}
					}
					else
					{
						if((abs(a - from)) <= RP_BD)
						{
							within_b = 1;
						}
						else
						{
							within_b = -1;
						}
					}
				}

				if( within_b == -1 )
				{
					return(-1);
				}
				else 
				{
					if(abs(a - pre_end) > (3*S_RP_BD))
					{
						return(-1);
					}
					else pre_end = b;

					if( (j >= 1) && (rate == 0 ) ) {
						// Simple_repeat or Low_complexity
					}
					else {
						rp[j].start = a;
						rp[j].end = b;
						rp[j].d_rate = rate;
						rp[j].id = id;
						j++;
					}

					if( j >= Max_Rps ) 
					{
						return(-1);
					}			
				}
			}
		}
		else is_over = true;
		i++;
	}

	if((to > a) && (is_over == true) )
	{
		if( (width(i_rp) < (5*RP_BD)) && (width(i_rp) >= S_RP_BD) )
		{
			if(((float)(((float)(abs(a - from))) / ((float)width(i_rp))) <= RR_NC) && ((float)(((float)(abs(to - b))) / ((float)width(i_rp))) <= RR_NC))
			{
				within_b = 1;
			}
			else if((float)(((float)(to - a)) / ((float)width(i_rp))) <= RR_NC)
			{
				within_b = 0;
			}
			else
			{
				within_b = -1;
			}
		}
		else if( width(i_rp) >= (5*RP_BD))
		{
			if((abs(a - from) <= RP_BD) && (abs(to - b) <= RP_BD))
			{
				within_b = 1;
			}
			else if((to - a) <= RP_BD)
			{
				within_b = 0;
			}
			else
			{
				within_b = -1;
			}
		}
		else within_b = 0;
	}
	else 
	{
		within_b = 0;
	}

	if( within_b == 1 )
	{
		if( abs(a - pre_end) > (3*S_RP_BD) )
		{
			return(-1);
		}

		if( (j >= 1) && (rate == 0 ) ) {
			// Simple_repeat or Low_complexity
		}
		else {
			if( i > 0 ) {
				rp[j].start = a;
				rp[j].end = b;
				rp[j].d_rate = rate;
				rp[j].id = rp_list[i-1].id;
				j++;	
			}
		}

		if( j >= Max_Rps ) 
		{
			return(-1);
		}			
	}
	else if( within_b == 0 )
	{
		if( abs(to - pre_end) > (3*S_RP_BD) )
		{
			return(-1);
		}
	}
	else if( within_b == -1 )
	{
		return(-1);
	}

	return j;
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

