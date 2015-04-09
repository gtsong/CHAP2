#include "main.h"
#include "regions.h"
#include "find_own_dup_copy.h"
#include "find_dup_copy.h"
#include "deal_gaps.h"
#include "util_gen.h"
#include "util_i.h"
#include "util.h"
#include "kd_op.h"

void find_own_dup_copy(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size)
{
	struct slist *sorted;
	int i = 0;
	int num_lines;
	int xval = 0, yval = 0;
	int w_sid = 0, h_fid = 0;

	sorted = (struct slist *) ckalloc(sizeof(struct slist) * (*num));

	num_lines = *num;
	sort_list(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		if( (dots[sorted[i].id].sign == 2) || (dots[sorted[i].id].pair_self == PAIR) || (dots[sorted[i].id].rp1_id != -1))
		{
		}
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
			
			w_sid = find_id(dots, tree, size, sorted[i].id, xval, yval, OW_SID);
			h_fid = find_id(dots, tree, size, sorted[i].id, xval, yval, OH_FID);

			find_opt_own_dup_copy(dots, sorted[i].id, p_pts, w_sid, h_fid);
		}
	}

	free(sorted);
}

void find_opt_own_dup_copy(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int h_fid)
{
	int i = 0;
	bool *is_x;
	int *sd;
	int start = 0, mid1 = -1, mid2 = -1, end = 0;
	int len1 = 0, len2 = 0;

	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));

	*is_x = true;
	*sd = 0;

	start = w_sid;
	end = h_fid;
	for( i = start; i <= end; i++ )
	{
		if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
			if( st[i].id == id ) {}
			else if( dots[st[i].id].pair_self == PAIR ) {}
			else if( (dots[st[i].id].ctg_id1 != dots[id].ctg_id1) || (dots[st[i].id].ctg_id2 != dots[id].ctg_id2) ) {}
			else if( abs(dots[st[i].id].identity - dots[id].identity) > DIFF_PID ) {}
			else if( dots[st[i].id].l_id != -1 ) 
			{
			}
			else if( dots[id].l_id != -1 ) 
			{
			}
			else if( dots[st[i].id].x.lower > dots[id].x.lower ) 
			{
			}
			else if( (dots[st[i].id].sign == 0) && (dots[st[i].id].y.lower > dots[id].y.lower )) 
			{
			}
			else if( (dots[st[i].id].sign == 1) && (dots[st[i].id].y.lower < dots[id].y.lower )) 
			{
				if(check_into_own(dots, st[i].id, id) != NON_COPY )
				{
					dots[st[i].id].lock = BOTH_LOCK;
					dots[id].lock = BOTH_LOCK;
					dots[st[i].id].l_id = id;
					dots[st[i].id].m_x = assign_I(dots[st[i].id].x.lower, dots[st[i].id].x.upper);
					dots[st[i].id].m_y = assign_I(dots[st[i].id].y.lower, dots[st[i].id].y.upper);
					dots[st[i].id].m_pid = dots[st[i].id].identity;
					dots[st[i].id].identity = ((dots[id].identity * width(dots[id].x)) + (dots[st[i].id].identity * width(dots[st[i].id].x))) / (width(dots[id].x) + width(dots[st[i].id].x));
					if( dots[st[i].id].sign == 0 )
					{
						dots[st[i].id].y = assign_I(dots[st[i].id].m_x.upper, dots[id].x.upper); 
					}
					else if( dots[st[i].id].sign == 1 )
					{
						dots[st[i].id].y = assign_I(dots[st[i].id].m_x.upper, dots[st[i].id].y.upper); 
					}
					dots[st[i].id].x = assign_I(dots[st[i].id].m_x.lower, dots[st[i].id].m_x.lower + width(dots[st[i].id].y)); 
					dots[id].sign = 2;
				}
			}
			else 
			{
				if((dots[st[i].id].sign != 2) && (dots[st[i].id].sign == dots[id].sign))
				{
					len1 = width(dots[st[i].id].x);
					len2 = width(dots[id].x);
					if(check_into_own(dots, st[i].id, id) != NON_COPY )
					{
						dots[st[i].id].lock = BOTH_LOCK;
						dots[id].lock = BOTH_LOCK;
						dots[st[i].id].l_id = id;
						dots[st[i].id].m_x = assign_I(dots[st[i].id].x.lower, dots[st[i].id].x.upper);
						dots[st[i].id].m_y = assign_I(dots[st[i].id].y.lower, dots[st[i].id].y.upper);
						dots[st[i].id].m_pid = dots[st[i].id].identity;
						dots[st[i].id].identity = ((dots[id].identity * width(dots[id].x)) + (dots[st[i].id].identity * width(dots[st[i].id].x))) / (width(dots[id].x) + width(dots[st[i].id].x));
						if( dots[st[i].id].sign == 0 )
						{
							dots[st[i].id].y = assign_I(dots[st[i].id].m_x.upper, dots[id].x.upper); 
						}
						else if( dots[st[i].id].sign == 1 )
						{
							dots[st[i].id].y = assign_I(dots[st[i].id].m_x.upper, dots[st[i].id].y.upper); 
						}
						dots[st[i].id].x = assign_I(dots[st[i].id].m_x.lower, dots[st[i].id].m_x.lower + width(dots[st[i].id].y)); 
						dots[id].sign = 2;
					}
				}
			}
		}
	}

	free(sd);
	free(is_x);
}
