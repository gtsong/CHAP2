#include "main.h"
#include "regions.h"
#include "adjust_pair.h"

int mark_outparalog(int num_pair, struct DotList *pair_alg, int id)
{
	int i;
	int num_checked = 0;

	for( i = 0; i < num_pair; i++ )
	{
		if( (pair_alg[i].sign == 0) || (pair_alg[i].sign == 1) )
		{
			if( loose_overlap(pair_alg[id].x, pair_alg[i].x) == true )
			{
				if( pair_alg[i].sign == 0 )
				{
					pair_alg[i].sign = OUT_PAR;
				}
				else if( pair_alg[i].sign == 1 )
				{
					pair_alg[i].sign = OUT_PAR_REV;
				}
				num_checked++;
			}
			else if( loose_overlap(pair_alg[id].y, pair_alg[i].y) == true )
			{
				if( pair_alg[i].sign == 0 )
				{
					pair_alg[i].sign = OUT_PAR;
				}
				else if( pair_alg[i].sign == 1 )
				{
					pair_alg[i].sign = OUT_PAR_REV;
				}
				num_checked++;
			}
		}
	}
	return(num_checked);
}
