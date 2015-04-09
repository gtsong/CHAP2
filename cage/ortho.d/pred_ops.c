#include "main.h"
#include "regions.h"
#include "pred_ops.h"
#include "rollback.h"
#include "pred_regions.h"
#include "adjust_for_indels.h"
#include "util_i.h"

void predict_op(bool is_x, int id, int *num_list, struct DotList *dots, int overlap, int num_ops, struct ops_list *ops)
{
	char op_ch;
	int conversion;
	int pred_op;
	bool is_x_to_y;

	conversion = 0;
	if( dots[id].sign == 0 )
	{
		op_ch = '+';
	}
	else if( dots[id].sign == 1)
	{
		op_ch = '-';
	}
	else
	{
	}
	
	if( overlap == -1 ) pred_op = 3;
	else pred_op = 2;

	if( dots[id].l_id != -1 ) is_x_to_y = true;
	else if( is_x == true ) is_x_to_y = false;
	else is_x_to_y = true;
	pred_dup(conversion, op_ch, pred_op, is_x_to_y, id, num_list, dots, num_ops, ops);
}

void pred_dup(int con, char op_ch, int pred_op, bool is_x_to_y, int id, int *num_list, struct DotList *dots, int num_ops, struct ops_list *ops)
{
	int wide;
	struct I from, to;
	int flag = DEL;
	int i;
	int sp_id;

	sp_id = dots[id].sp_id;

	if((dots[id].l_id == -1) && (proper_overlap(dots[id].x, dots[id].y) == true) && (width(intersect(dots[id].x, dots[id].y)) <= THRESHOLD))
	{
		dots[id].y = assign_I(dots[id].x.upper, dots[id].x.upper + width(dots[id].x));
	}

	for( i = 0; i < *num_list; i++ )
	{
		if( i != id )
		{
			if( (dots[i].l_id != -1) && (dots[i].sign != 2) )
			{
				dots[i].x = assign_I(dots[i].m_x.lower, dots[i].m_x.upper);
				dots[i].y = assign_I(dots[i].m_y.lower, dots[i].m_y.upper);
				dots[dots[i].l_id].sign = dots[i].sign;
				dots[i].l_id = -1;
				dots[i].identity = dots[i].m_pid;
				dots[i].m_x = assign_I(0,1);
				dots[i].m_y = assign_I(0,1);
			}
		}
	}

	if( dots[id].l_id != -1 )
	{
		from = assign_I(dots[id].x.lower, dots[id].x.upper);
		to = assign_I(dots[id].y.lower, dots[id].y.upper);
		dots[id].sign = 2;
		flag = NONE;	
	}
	else if( is_x_to_y )
	{
		from = assign_I(dots[id].x.lower, dots[id].x.upper);
		to = assign_I(dots[id].y.lower, dots[id].y.upper);
	}
	else
	{
		from = assign_I(dots[id].y.lower, dots[id].y.upper);
		to = assign_I(dots[id].x.lower, dots[id].x.upper);
	}

	if( pred_op == 0 ) 
	{
		wide = rollback_step_dup_no_overlap(is_x_to_y, id, num_list, dots);
	}
	else if(pred_op == 2)
	{	
		wide = rollback_step_dup_no_overlap(is_x_to_y, id, num_list, dots);
	}
	else if(pred_op == 3)
	{
		wide = rollback_step_dup_overlap(is_x_to_y, id, num_list, dots);
	}
	else if(pred_op == 4)
	{
		wide = rollback_step_conversion(is_x_to_y, id, num_list, dots);
		if( con > 0 ) wide = con;
	}

	generate_ops(op_ch, wide, is_x_to_y, from, to, flag, num_ops, ops, sp_id);
}

void generate_ops_del(char op, int location, int wide, int num_ops, struct ops_list *ops)
{
	ops[num_ops].sign = op;
	ops[num_ops].src_b = location;
	ops[num_ops].src_e = location+wide;
	ops[num_ops].dst_b = 0;
	ops[num_ops].dst_e = 0;
}

void generate_ops(char op, int wide, bool direction, struct I from, struct I to, int flag, int num_ops, struct ops_list *ops, int sp_id)
{
	if( flag == NONE ) {}
	else
	{
	}

	if( wide == 0 )
	{
		if( direction == true )
		{
			ops[num_ops].sign = op;
			ops[num_ops].src_b = from.lower;
			ops[num_ops].src_e = from.upper;
			ops[num_ops].dst_b = to.lower;
			ops[num_ops].dst_e = to.upper;
			ops[num_ops].sp_id = sp_id;
		}
		else
		{
			ops[num_ops].sign = op;
			ops[num_ops].src_b = from.lower - width(to);
			ops[num_ops].src_e = from.upper - width(to);
			ops[num_ops].dst_b = to.lower;
			ops[num_ops].dst_e = to.upper;
			ops[num_ops].sp_id = sp_id;
		}
	}
	else
	{
	}
}
