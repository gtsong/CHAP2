#include "main.h"
#include "extend_slist.h"
#include "regions.h"
#include "pred_regions.h"
#include "util.h"
#include "util_i.h"

/* Need to consider more about the direction of conversion */
int find_covering_regions(int youngest, bool *is_x, int *regs, int num_list, struct DotList *dots)
{
	int num = 0;
	int i;

	for(i = 0; i < num_list; i++)
	{
		if( i == youngest) {}
		else
		{
			if(subset(dots[youngest].y, dots[i].y))
			{
				regs[num] = i;
				is_x[num] = false;
				num++;
			}
			else if(subset(dots[youngest].y, dots[i].x))
			{
				regs[num] = i;
				is_x[num] = true;
				num++;
			}
		}
	}
	return num;
}

/* 0: no overlap
 * 1: x part has overlap
 * 2: y part has overlap
 * 3: both x and y parts have overlap - need overwriting
 */
int is_left_to_right_count(int *num_x, int *num_y, int id, int num_list, struct DotList *org, int t_val, struct ID_List *dlist, int num_dup, int flag, struct kdnode *tree, struct perm_pt *p_pts, int size, FILE *fp, struct DotList *init_algns)
{
	int i = 0;
	int res = NO_OVERLAP;
	bool x_fully_belong = false;
	bool y_fully_belong = false;

	*num_x = 0;
	*num_y = 0;

	while( i < num_list )
	{
		if( (i == id) || (org[i].sign == 2) ) 
		{
			i++;
		}
		else if( (flag == TANDEM_CHECK) && ((f_loose_subset(org[i].x, org[id].x, t_val) || f_loose_subset(org[i].y, org[id].x, t_val) || f_loose_subset(org[i].x, org[id].x, t_val) || f_loose_subset(org[i].y, org[id].y, t_val)) || (is_s_list(org, id, i, dlist, num_dup, tree, p_pts, size, fp, init_algns) == true)) )  // check whether dots[id] is in the suspend list or not
		{
			i++;
		}
		else
		{
			if( (fully_subset(org[id].x, org[i].x) == true) || (fully_subset(org[id].x, org[i].y) == true) )
			{
				x_fully_belong = true;
			}

			if( ( ((org[i].pair_self==PAIR) && ((org[i].identity-2) >= org[id].identity)) && (((f_loose_subset(org[i].x, org[id].x, t_val) == true) && (f_loose_subset(org[id].x, org[i].x, t_val) == true)) || ((f_loose_subset(org[i].y, org[id].x, t_val) == true) && (f_loose_subset(org[id].x, org[i].y, t_val) == true))) ) || ((f_loose_subset(org[i].x, org[id].x, t_val) == false) && ((width(org[i].x) >= 3*D_OP_TH) && (f_loose_overlap(org[id].x, org[i].x, t_val) == true))) || ((f_loose_subset(org[i].y, org[id].x, t_val) == false) && ((width(org[i].y) >= 3*D_OP_TH) && (f_loose_overlap(org[id].x, org[i].y, t_val) == true)))) 
			{
				(*num_x)++;
				if( (res == OVERLAP_IN_Y) || (res == OVERLAP_BOTH) )
				{
					res = OVERLAP_BOTH;
				}
				else
				{
					res = OVERLAP_IN_X;
				}
			}

			if( (fully_subset(org[id].y, org[i].x) == true) || (fully_subset(org[id].y, org[i].y) == true) )
			{
				y_fully_belong = true;
			}

			if( (((org[i].pair_self==PAIR) && ((org[i].identity-2) >= org[id].identity)) && (((f_loose_subset(org[i].x, org[id].y, t_val) == true) && (f_loose_subset(org[id].y, org[i].x, t_val) == true)) || ((f_loose_subset(org[i].y, org[id].y, t_val) == true) && (f_loose_subset(org[id].y, org[i].y, t_val) == true))) ) || ((f_loose_subset(org[i].x, org[id].y, t_val) == false) && ((width(org[i].x) >= 3*D_OP_TH) && (f_loose_overlap(org[id].y, org[i].x, t_val) == true))) || ((f_loose_subset(org[i].y, org[id].y, t_val) == false) && ((width(org[i].y) >= 3*D_OP_TH) && (f_loose_overlap(org[id].y, org[i].y, t_val) == true)))) 
			{
				(*num_y)++;
				if( (res == OVERLAP_IN_X) || (res == OVERLAP_BOTH) )
				{
					res = OVERLAP_BOTH;
				}
				else 
				{
					res = OVERLAP_IN_Y;
				}
			}

			if( org[i].l_id != -1 )
			{
				if( ((f_loose_subset(org[org[i].l_id].y, org[id].x, t_val) == false) && (f_loose_overlap(org[id].x, org[org[i].l_id].y, t_val) == true)) ) 
				{
					(*num_x)++;
					if( (res == OVERLAP_IN_Y) || (res == OVERLAP_BOTH) )
					{
						res = OVERLAP_BOTH;
					}
					else
					{
						res = OVERLAP_IN_X;
					}
				}

				if( ((f_loose_subset(org[org[i].l_id].y, org[id].y, t_val) == false) && (f_loose_overlap(org[id].y, org[org[i].l_id].y, t_val) == true)) ) 
				{
					(*num_y)++;
					if( (res == 1) || (res == 3) )
					{
						res = OVERLAP_BOTH;
					}
					else 
					{
						res = OVERLAP_IN_Y;
					}
				}
			}

			i++;
		}
	}
	
	if( x_fully_belong == true ) 
	{
		*num_x = num_list;
	}

	if( y_fully_belong == true ) 
	{
		*num_y = num_list;
	}

	return(res);
}

int is_left_to_right_dup_count(int *num_x, int *num_y, int id, int num_list, struct DotList *org, int t_val, int r_id, int l_id)
{
	int i = 0;
	int res = 0;
	bool x_fully_belong = false;
	bool y_fully_belong = false;

	*num_x = 0;
	*num_y = 0;
	while( i < num_list )
	{
		if( (i == id) || (i == r_id) || (i == l_id) || (org[i].sign == 2) ) i++;
		else
		{
			if( (fully_subset(org[id].x, org[i].x) == true) || (fully_subset(org[id].x, org[i].y) == true) )
			{
				x_fully_belong = true;
			}

			if( ((f_loose_subset(org[i].x, org[id].x, t_val) == false) && (f_loose_overlap(org[id].x, org[i].x, t_val) == true)) || ((f_loose_subset(org[i].y, org[id].x, t_val) == false) && (f_loose_overlap(org[id].x, org[i].y, t_val) == true))) 
			{
				(*num_x)++;
				if( (res == 2) || (res == 3) )
				{
					res = 3;
				}
				else
				{
					res = 1;
				}
			}

			if( (fully_subset(org[id].y, org[i].x) == true) || (fully_subset(org[id].y, org[i].y) == true) )
			{
				y_fully_belong = true;
			}

			if( ((f_loose_subset(org[i].x, org[id].y, t_val) == false) && (f_loose_overlap(org[id].y, org[i].x, t_val) == true)) || ((f_loose_subset(org[i].y, org[id].y, t_val) == false) && (f_loose_overlap(org[id].y, org[i].y, t_val) == true))) 
			{
				(*num_y)++;
				if( (res == 1) || (res == 3) )
				{
					res = 3;
				}
				else 
				{
					res = 2;
				}
			}
			i++;
		}
	}
	
	if( x_fully_belong == true ) 
	{
		*num_x = num_list;
	}

	if( y_fully_belong == true ) 
	{
		*num_y = num_list;
	}

	return(res);
}

void get_elm_list(int *num_id, int *elm_id, int id, bool is_x, int num_list, struct DotList *dots)
{
	int i = 0;
	struct I temp;
	int t_val = LOOSE;

	*num_id = 0;
	if( is_x == true )
	{
		temp = assign_I(dots[id].x.lower, dots[id].x.upper);
	}
	else temp = assign_I(dots[id].y.lower, dots[id].y.upper);

	while( i < num_list )
	{
		if( (i == id) || (dots[i].sign == 2) ) i++;
		else
		{
			if( ((f_loose_subset(dots[i].x, temp, t_val) == false) && (f_loose_overlap(temp, dots[i].x, t_val) == true)) || ((f_loose_subset(dots[i].y, temp, t_val) == false) && (f_loose_overlap(temp, dots[i].y, t_val) == true))) 
			{
				elm_id[*num_id] = i;
				(*num_id)++;
			}
			i++;
		}
	}
}

int is_left_to_right(int id, int num_list, struct DotList *org)
{
	int i = 0;
	int res = 0;

	while( i < num_list )
	{
		if( (i == id ) || (org[i].sign == 2)) i++;
		else
		{
			if( ((loose_subset(org[i].x, org[id].x) == false) && (loose_overlap(org[id].x, org[i].x) == true)) || ((loose_subset(org[i].y, org[id].x) == false) && (loose_overlap(org[id].x, org[i].y) == true))) 
			{
				if( res == 2 )
				{
					res = 3;
					return(res);
				}
				else
				{
					res = 1;
				}
			}
			if( ((loose_subset(org[i].x, org[id].y) == false) && (loose_overlap(org[id].y, org[i].x) == true)) || ((loose_subset(org[i].y, org[id].y) == false) && (loose_overlap(org[id].y, org[i].y) == true))) 
			{
				if( res == 1 )
				{
					res = 3;
					return(res);
				}
				else 
				{
					res = 2;
				}
			}
			i++;
		}
	}
	return(res);
}

bool check_negligible_line(int id, int num_list, struct DotList *dots)
{
	int i = 0;

	while( i < num_list )
	{
		if( i == id ) i++;
		else
		{
			if(proper_subset(dots[id].x, dots[i].x) && proper_subset(dots[id].y, dots[i].x) )
			{
				return true;
			}
			else if(proper_subset(dots[id].x, dots[i].y) && proper_subset(dots[id].y, dots[i].y) )
			{
				return true;
			}
			i++;
		}	
	}
	return false;
}

int is_left_to_right_again(int id, int num_list, struct DotList *dots)
{
	int i = 0;
	bool is_conv = false;
	bool is_overlap_x = false;
	bool is_overlap_y = false;
	bool is_subset_x = false;
	bool is_subset_y = false;

	while( (i < num_list) && (is_conv == false) )
	{
		if( (i == id) || (dots[i].sign == 2) ) i++;
		else
		{
			if((proper_overlap(dots[id].x, dots[i].x) && (!subset(dots[id].x, dots[i].x))) || (proper_overlap(dots[id].x, dots[i].y) && (!subset(dots[id].x, dots[i].y))))
			{
				if(is_overlap_x == false) is_overlap_x = true;
				if(is_subset_y || is_overlap_y) 
				{
					is_conv = true;
				}
			}

			if(proper_subset(dots[id].x, dots[i].x) || proper_subset(dots[id].x, dots[i].y))
			{
				if(is_subset_x == false) is_subset_x = true;
				if(is_overlap_y == true) 
				{
					is_conv = true;
				}
			}

			if((proper_overlap(dots[id].y, dots[i].x) && (!subset(dots[id].y, dots[i].x))) || (proper_overlap(dots[id].y, dots[i].y) && (!subset(dots[id].y, dots[i].y))))
			{
				if(is_overlap_y == false) is_overlap_y = true;
				if(is_subset_x || is_overlap_x)
				{
					is_conv = true;
				}
			}

			if(proper_subset(dots[id].y, dots[i].x) || proper_subset(dots[id].y, dots[i].y))
			{
				if(is_subset_y == false) is_subset_y = true;
				if(is_overlap_x == true) 
				{
					is_conv = true;
				}
			}
			i++;
		}
	}

	i = 0;
	while( i < num_list )
	{
		if( i == id ) i++;
		else
		{
			if(equal(dots[id].y, dots[i].x) || equal(dots[id].y, dots[i].y))
			{
				if(is_conv == true) return(CON_X_TO_Y);
				else return(DUP_X_TO_Y);
			}
			else if(equal(dots[id].x, dots[i].x) || equal(dots[id].x, dots[i].y) )
			{
				if(is_conv == true) return(CON_Y_TO_X);
				else return(DUP_Y_TO_X);
			}
			i++;
		}
	}

	i = 0;
	while( i < num_list )
	{
		if( i == id ) i++;
		else
		{
			if(subset(dots[id].x, dots[i].x) || subset(dots[id].x, dots[i].y) )
			{
				if(is_conv == true) return(CON_X_TO_Y);
				else return(DUP_X_TO_Y);	
			}
			else if(subset(dots[id].y, dots[i].x) || subset(dots[id].y, dots[i].y))
			{
				if(is_conv == true) return(CON_X_TO_Y);
				else return(DUP_Y_TO_X);
			}
			i++;
		}
	}
	
	if(is_conv == true) return(CON_X_TO_Y);
	else return(DUP_X_TO_Y);
}

int is_left_to_right_count_strict(int *num_x, int *num_y, int id, int num_list, struct DotList *org)
{
	int i = 0;
	int res = 0;

	*num_x = 0;
	*num_y = 0;
	while( i < num_list )
	{
		if( (i == id) || (org[i].sign == 2) ) i++;
		else
		{
			if( ((subset(org[i].x, org[id].x) == false) && (proper_overlap(org[id].x, org[i].x) == true)) || ((subset(org[i].y, org[id].x) == false) && (proper_overlap(org[id].x, org[i].y) == true))) 
			{
				(*num_x)++;
				if( (res == 2) || (res == 3) )
				{
					res = 3;
				}
				else
				{
					res = 1;
				}
			}

			if( ((subset(org[i].x, org[id].y) == false) && (overlap(org[id].y, org[i].x) == true)) || ((subset(org[i].y, org[id].y) == false) && (overlap(org[id].y, org[i].y) == true))) 
			{
				(*num_y)++;
				if( (res == 1) || (res == 3) )
				{
					res = 3;
				}
				else 
				{
					res = 2;
				}
			}
			i++;
		}
	}
	return(res);
}

bool is_s_list(struct DotList *dots, int ins_id, int cur_id, struct ID_List *dlist, int num_dup, struct kdnode *tree, struct perm_pt *p_pts, int size, FILE *fp, struct DotList *init_dots)
{
	bool res = false;
	int i = 0;
	struct I temp;
	int x_opt_id, y_opt_id;
	bool *f_is_x;

	f_is_x = (bool *) ckalloc(sizeof(bool));

	x_opt_id = find_alt_ins_id(cur_id, dots, tree, p_pts, ins_id, true, f_is_x, size, fp, init_dots);
	y_opt_id = find_alt_ins_id(cur_id, dots, tree, p_pts, ins_id, false, f_is_x, size, fp, init_dots);

	if( (x_opt_id != -1) || (y_opt_id != -1) )
	{
		res = true;
	}

	while( (i < num_dup) && (res == false))
	{
		if( dlist[i].is_x == true )
		{
			temp = assign_I(dots[dlist[i].m_id].x.lower, dots[dlist[i].m_id].x.upper);
		}
		else
		{
			temp = assign_I(dots[dlist[i].m_id].y.lower, dots[dlist[i].m_id].y.upper);
		}

		if( (dlist[i].m_id == ins_id) || (strict_almost_equal(dots[ins_id].x, temp) == true) || (strict_almost_equal(dots[ins_id].y, temp) == true) )
		{
			if( (cur_id == dlist[i].left_id) || (cur_id == dlist[i].right_id) )
			{
				res = true;
			}
		}
		i++;
	}

	free(f_is_x);
	return(res);
}
