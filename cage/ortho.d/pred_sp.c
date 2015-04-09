#include "main.h"
#include "pred_sp.h"
#include "regions.h"
#include "pred_ops.h"
#include "rollback.h"
#include "pred_regions.h"
#include "adjust_for_indels.h"
#include "find_gene_loss.h"
#include "util_i.h"

void predict_sp_op(int sp_code, int rm_sp, int left_sp, int *num_list, struct DotList *dots, int *cur_num, struct ops_list *ops)
{
	char op_ch;
	int r_st = -1, r_end = -1; // the range of a removed species
	struct I temp_reg;
	int len;
	int i = 0;

	check_gene_loss(num_list, dots, sp_code, rm_sp, left_sp, cur_num, ops);

	op_ch = 's';
	for( i = 0; i < (*num_list); i++ )
	{
		if( dots[i].sp_id == sp_code )
		{
			if( ( r_st == -1 ) && ( r_end == -1 ) )
			{
				r_st = dots[i].y.lower;
				r_end = dots[i].y.upper;
			}
			else 
			{
				if( dots[i].y.lower < r_st ) r_st = dots[i].y.lower;
				if( dots[i].y.upper > r_end ) r_end = dots[i].y.upper;
			}
		}
	}

	temp_reg = assign_I(r_st, r_end);
	len = r_end - r_st + 1;

	for( i = 0; i < (*num_list); i++ )
	{
		if( (proper_overlap(temp_reg, dots[i].x) == true) || (proper_overlap(temp_reg, dots[i].y) == true) )
		{
			dots[i].sign = 2;	
		}
		else 
		{
			if( dots[i].x.lower > r_st )
			{
				dots[i].x = assign_I(dots[i].x.lower - len, dots[i].x.upper - len);
				dots[i].y = assign_I(dots[i].y.lower - len, dots[i].y.upper - len);
			}
			else if( dots[i].y.lower > r_st )
			{
				dots[i].y = assign_I(dots[i].y.lower - len, dots[i].y.upper - len);
			}
		}
	}

	overwrite_dots(num_list, dots);
	ops[*cur_num].sign = op_ch;
	ops[*cur_num].src_b = r_st;
	ops[*cur_num].src_e = r_end;
	ops[*cur_num].dst_b = 0;
	ops[*cur_num].dst_e = 0;
	ops[*cur_num].sp_id = rm_sp;
}

void rollback_sp(struct DotList *dots, int *num_dots, struct DotList *init_dots, int num_init_dots, int cur_num, struct ops_list *ops, int sp_id)
{
	int i;
	char op_ch;

	op_ch = 's';
	for( i = 0; i < (*num_dots); i++ ) {
		if( (dots[i].sign != DELETED) && ( (dots[i].sp_id == sp_id) || (dots[i].sp_id == PAIR) ) ) dots[i].sign = DELETED;
	}

	for( i = 0; i < num_init_dots; i++ ) {
		if( (init_dots[i].sign != DELETED) && ( (init_dots[i].sp_id == sp_id) || (init_dots[i].sp_id == PAIR) ) ) init_dots[i].sign = DELETED;
	}

	overwrite_dots(num_dots, dots);
	ops[cur_num].sign = op_ch;
	ops[cur_num].src_b = 0;
	ops[cur_num].src_e = 0;
	ops[cur_num].dst_b = 0;
	ops[cur_num].dst_e = 0;
	ops[cur_num].sp_id = sp_id;
}
