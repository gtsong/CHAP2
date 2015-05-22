#include "main.h"
#include "apply_ops.h"
#include "regions.h"
#include "util.h"
#include "util_i.h"
#include "util_ops.h"
#include "util_gen.h"

void print_conv_on_dup(int avg_pid, struct DotList alg, struct cv_list *cv, int num_cv, int *num_ops, struct ops_list *ops, int run_mode)
{
	int i;
	struct I src, dst;
	int sign;
	int ops_id;

	ops_id = *num_ops;

	for( i = 0; i < num_cv; i++ ) {
		src = assign_I(cv[i].a1, cv[i].a2);
		dst = assign_I(cv[i].b1, cv[i].b2);
		if( cv[i].ori == '+' ) sign = 0;
		else sign = 1;

		if( (run_mode == INF_DUP) && (alg.sign != sign) ) {}
		else if( (run_mode == INF_DUP) && (((subset(src, alg.x) == true) && (subset(dst, alg.y) == true) ) || ((subset(src, alg.y) == true) && (subset(dst, alg.x) == true))) ) 
		{
			ops[ops_id].id = alg.index;
			if( cv[i].ori == '+' ) ops[ops_id].sign = 'c';
			if( cv[i].ori == '-' ) ops[ops_id].sign = 'v';
			ops[ops_id].dir = cv[i].dir;
			ops[ops_id].sp_id = alg.sp_id;
			ops[ops_id].pid = cv[i].conv_pid;
			if( cv[i].dir == 1 )  {
				ops[ops_id].src_b = cv[i].b1;
				ops[ops_id].src_e = cv[i].b2;
				ops[ops_id].ctg_id1 = cv[i].ctg_id2;
				ops[ops_id].dst_b = cv[i].a1;
				ops[ops_id].dst_e = cv[i].a2;
				ops[ops_id].ctg_id2 = cv[i].ctg_id1;
			}
			else {
				ops[ops_id].src_b = cv[i].a1;
				ops[ops_id].src_e = cv[i].a2;
				ops[ops_id].ctg_id1 = cv[i].ctg_id1;
				ops[ops_id].dst_b = cv[i].b1;
				ops[ops_id].dst_e = cv[i].b2;
				ops[ops_id].ctg_id2 = cv[i].ctg_id2;
			}
			ops_id++;
		}
		else if( (run_mode == CONTENT_ORTHO) && (cv[i].ori != 'd') && ((cv[i].conv_pid) > alg.identity) && (cv[i].conv_pid > avg_pid)) {
			ops[ops_id].id = alg.index;
			if( cv[i].ori == '+' ) ops[ops_id].sign = 'c';
			if( cv[i].ori == '-' ) ops[ops_id].sign = 'v';
			ops[ops_id].dir = cv[i].dir;
			ops[ops_id].pid = cv[i].conv_pid;
//			ops[ops_id].sp_id = alg.sp_id;
			ops[ops_id].sp_id = cv[i].sp_id;
			if( cv[i].dir == 1 )  {
				ops[ops_id].src_b = cv[i].b1;
				ops[ops_id].src_e = cv[i].b2;
				ops[ops_id].ctg_id1 = cv[i].ctg_id2;
				ops[ops_id].dst_b = cv[i].a1;
				ops[ops_id].dst_e = cv[i].a2;
				ops[ops_id].ctg_id2 = cv[i].ctg_id1;
			}
			else {
				ops[ops_id].src_b = cv[i].a1;
				ops[ops_id].src_e = cv[i].a2;
				ops[ops_id].ctg_id1 = cv[i].ctg_id1;
				ops[ops_id].dst_b = cv[i].b1;
				ops[ops_id].dst_e = cv[i].b2;
				ops[ops_id].ctg_id2 = cv[i].ctg_id2;
			}
			cv[i].ori = 'd';
			ops_id++;
		}
	}
	*num_ops = ops_id;
}

int update_conv(struct DotList *algns, int rm_id, bool is_x, struct cv_list *cv, int num_cv)
{
	struct I rm_reg;
	struct I src, dst;
	int i = 0, j = 0;

	if( is_x == true ) rm_reg = assign_I(algns[rm_id].x.lower, algns[rm_id].x.upper);
	else rm_reg = assign_I(algns[rm_id].y.lower, algns[rm_id].y.upper);
	
	for( i = 0; i < num_cv; i++ ) {
		src = assign_I(cv[i].a1, cv[i].a2);
		dst = assign_I(cv[i].b1, cv[i].b2);
		if( (proper_overlap(src, rm_reg) == true) || (proper_overlap(dst, rm_reg) == true ) ) {
			cv[i].fid = -1;
		}
		else {
			if(src.lower >= rm_reg.lower) {
				cv[i].a1 = cv[i].a1 - width(rm_reg);
				cv[i].a2 = cv[i].a2 - width(rm_reg);
			}

			if(dst.lower >= rm_reg.lower) {
				cv[i].b1 = cv[i].b1 - width(rm_reg);
				cv[i].b2 = cv[i].b2 - width(rm_reg);
			}
		}
	}

	j = 0; 
	for( i = 0; i < num_cv; i++ ) 
	{
		if( cv[i].fid != -1 ) {
			cv[j] = assign_conv(cv[i]);
			j++;
		}
	}

	return(j);
}

int update_conv_del(struct I reg, struct cv_list *cv, int num_cv)
{
	int i, j;
	struct I ins_reg;
	struct I src, dst;

	ins_reg = assign_I(reg.lower-DEL_TH, reg.lower+DEL_TH);

	for( i = 0; i < num_cv; i++ ) {
		src = assign_I(cv[i].a1, cv[i].a2);
		dst = assign_I(cv[i].b1, cv[i].b2);
		if( (proper_overlap(src, ins_reg) == true) || (proper_overlap(dst, ins_reg) == true ) ) {
			cv[i].fid = -1;
		}
		else {
			if(src.lower >= reg.lower) {
				cv[i].a1 = cv[i].a1 + width(reg);
				cv[i].a2 = cv[i].a2 + width(reg);
			}

			if(dst.lower >= reg.lower) {
				cv[i].b1 = cv[i].b1 + width(reg);
				cv[i].b2 = cv[i].b2 + width(reg);
			}
		}
	}

	j = 0; 
	for( i = 0; i < num_cv; i++ ) 
	{
		if( cv[i].fid != -1 ) {
			cv[j] = assign_conv(cv[i]);
			j++;
		}
	}

	return(j);
}

int update_exons(struct DotList *algns, int rm_id, bool is_x, struct exons_list *exons, int num_exons)
{
	struct I rm_reg;
	struct I cur_reg;
	struct I res;
	int i, j;
	int init_len;

	if( is_x == true ) rm_reg = assign_I(algns[rm_id].x.lower, algns[rm_id].x.upper);
	else rm_reg = assign_I(algns[rm_id].y.lower, algns[rm_id].y.upper);
	
	for( i = 0; i < num_exons; i++ ) {
		cur_reg = assign_I(exons[i].reg.lower, exons[i].reg.upper);
		if( (subset(cur_reg, rm_reg) == true ) || (subset(rm_reg, cur_reg) == true) ) {
			exons[i].fid = -1;
		}
		else if( proper_overlap(cur_reg, rm_reg) == true ) {
			if( cur_reg.lower < rm_reg.lower ) {
				res = assign_I(cur_reg.lower, rm_reg.lower);
			}
			else {
				res = assign_I(rm_reg.upper, cur_reg.upper);
			}

			if( width(res) <= EXON_TH ) {
				exons[i].fid = -1;
			}
			else {
				init_len = (int)((((float)100)/((float)exons[i].val)) * ((float)width(exons[i].reg)) );
				exons[i].val = (int)(((float)100)*((float)width(res))/((float)init_len));
				if( exons[i].val <= 5 ) exons[i].fid = -1;
				exons[i].reg = assign_I(res.lower, res.upper);
			}
		}
		else {
			if(cur_reg.lower >= rm_reg.lower) {
				exons[i].reg = assign_I(exons[i].reg.lower - width(rm_reg), exons[i].reg.upper - width(rm_reg));
			}
		}
	}

	for( i = 0; i < num_exons; i++ ) {
		if( exons[i].reg.upper < exons[i].reg.lower ) {
			exons[i].fid = -1;
		}
	}

	j = 0; 
	for( i = 0; i < num_exons; i++ ) 
	{
		if( exons[i].fid != -1 ) {
			exons[j] = assign_exons(exons[i]);
			j++;
		}
	}

	return(j);
}

int update_exons_del(struct I reg, struct exons_list *exons, int num_exons)
{
	int i;
	struct I cur_reg;
	int b = 0, e = 1;
	int j = 0;

	for( i = 0; i < num_exons; i++ ) {
		cur_reg = assign_I(exons[i].reg.lower, exons[i].reg.upper);
		if( cur_reg.lower >= reg.lower ) b = cur_reg.lower + width(reg);
		else b = cur_reg.lower;

		if( cur_reg.upper >= reg.lower ) e = cur_reg.upper + width(reg);
		else e = cur_reg.upper;

		exons[i].reg = assign_I(b, e);
	}

	for( i = 0; i < num_exons; i++ ) {
		if( exons[i].reg.upper < exons[i].reg.lower ) {
			exons[i].fid = -1;
		}
	}

	j = 0; 
	for( i = 0; i < num_exons; i++ ) 
	{
		if( exons[i].fid != -1 ) {
			exons[j] = assign_exons(exons[i]);
			j++;
		}
	}

	return(j);
}

int check_exons_bound(struct DotList algn, struct exons_list *exons, int num_exons, struct exons_list *genes, int num_genes)
{
	int i;
	int res_left = -1;
	int res_right = -1;
	struct I left, right;
	int res;

	left = assign_I(algn.x.lower, algn.x.upper);
	right = assign_I(algn.y.lower, algn.y.upper);

	i = 0;
	while((i < num_exons) && ((res_left == -1) || (res_right == -1))) {
		if( (f_loose_subset(exons[i].reg, left, STRICT) == true) || (f_loose_subset(left, exons[i].reg, STRICT) == true) ) {}
		else if( f_loose_overlap(exons[i].reg, left, STRICT) == true ) res_left = 0;
		
		if( (f_loose_subset(exons[i].reg, right, STRICT) == true) || (f_loose_subset(right, exons[i].reg, STRICT) == true) ) {}
		else if( f_loose_overlap(exons[i].reg, right, STRICT) == true ) res_right = 0;
		i++;
	}

	if( res_left == res_right ) {
		res_left = -1;
		res_right = -1;
		i = 0;
		while((i < num_genes) && ((res_left == -1) || (res_right == -1))) {
			if( (f_loose_subset(genes[i].reg, left, STRICT) == true) || (f_loose_subset(left, genes[i].reg, STRICT) == true) ) {}
			else if( f_loose_overlap(genes[i].reg, left, STRICT) == true ) res_left = 0;
		
			if( (f_loose_subset(genes[i].reg, right, STRICT) == true) || (f_loose_subset(right, genes[i].reg, STRICT) == true) ) {}
			else if( f_loose_overlap(genes[i].reg, right, STRICT) == true ) res_right = 0;
			i++;
		}
	}
	
	if( res_left == res_right ) res = TIE;
	else if( res_left == -1 ) res = LEFT_SIDE;
	else if( res_right == -1 ) res = RIGHT_SIDE;
	else res = TIE;
	
	return(res);
}

int cal_cur_pos_ops(int num_ops, struct ops_list *ops, struct ops_list *new_ops, int sp_id, int seq1_len)
{
	int i, j;
	int num_new_ops;
	struct I dup_reg;
	int cut_point;
	int del_len;
	int dup_len;

	j = 0;
	for( i = 0; i < num_ops; i++ ) {
		if( (ops[i].sign == '+') || (ops[i].sign == '-') || (ops[i].sign == 'd') || (ops[i].sign == 'c') || (ops[i].sign == 'v')) {
			new_ops[j] = assign_ops(ops[i]);
			new_ops[j].srcStart = ops[i].src_b;
			new_ops[j].srcEnd = ops[i].src_e;
			new_ops[j].ctg_id1 = ops[i].ctg_id1;
			new_ops[j].ctg_id2 = ops[i].ctg_id2;
			new_ops[j].id = i;
			if( ops[i].sign != 'd' ) {
				new_ops[j].dstStart = ops[i].dst_b;
				new_ops[j].dstEnd = ops[i].dst_e;
				if( ((ops[i].sign == '+') || (ops[i].sign == '-') ) && (new_ops[j].srcStart >= new_ops[j].dstStart) ) {
					dup_len = new_ops[j].dst_e - new_ops[j].dst_b;
					new_ops[j].srcStart = new_ops[j].srcStart + dup_len;
					new_ops[j].srcEnd = new_ops[j].srcEnd + dup_len;
				}
				else if( ((ops[i].sign == '+') || (ops[i].sign == '-') ) && (strict_almost_equal(assign_I(new_ops[j].srcStart, new_ops[j].srcEnd), assign_I(new_ops[j].dstStart, new_ops[j].dstEnd)) == true) ) { // tandem dup
					dup_len = new_ops[j].dst_e - new_ops[j].dst_b;
					new_ops[j].srcStart = new_ops[j].srcStart + dup_len;
					new_ops[j].srcEnd = new_ops[j].srcEnd + dup_len;
				}
			}
			j++;
		}
	}

	num_new_ops = j;
	for( i = (num_new_ops-1); i >= 0; i-- )
	{
		if( (new_ops[i].sign == '+') || (new_ops[i].sign == '-') ) {
			cut_point = -1;
			del_len = 0;
			dup_reg = assign_I(new_ops[i].dstStart, new_ops[i].dstEnd);
			dup_len = new_ops[i].dst_e - new_ops[i].dst_b;
		}
		else if(new_ops[i].sign == 'd') {
			dup_len = 0;	
			cut_point = new_ops[i].srcStart;
			del_len = abs(new_ops[i].src_e - new_ops[i].src_b);
		}

		j = i;
		for( j = (num_new_ops-1); j > i; j-- )
		{
			if( (new_ops[i].sign == '+') || (new_ops[i].sign == '-') ) 
			{
				if( (new_ops[j].sign == '+') || (new_ops[j].sign == '-') || (new_ops[j].sign == 'c') || (new_ops[j].sign == 'v'))
				{
					if( new_ops[j].srcStart >= dup_reg.lower ) new_ops[j].srcStart = new_ops[j].srcStart + dup_len;
					if( new_ops[j].srcEnd >= dup_reg.lower ) new_ops[j].srcEnd = new_ops[j].srcEnd + dup_len;
					if( new_ops[j].dstStart >= dup_reg.lower ) new_ops[j].dstStart = new_ops[j].dstStart + dup_len;
					if( new_ops[j].dstEnd >= dup_reg.lower ) new_ops[j].dstEnd = new_ops[j].dstEnd + dup_len;
				}
				else if(new_ops[j].sign == 'd')  
				{
					if( new_ops[j].srcStart >= dup_reg.lower ) new_ops[j].srcStart = new_ops[j].srcStart + dup_len;
					if( new_ops[j].srcEnd >= dup_reg.lower ) new_ops[j].srcEnd = new_ops[j].srcEnd + dup_len;
				}
				else {
					fatalf("unallowed operation in the new_ops list: %d\n", new_ops[j].sign);
				}
			}
			else if( new_ops[i].sign == 'd' ) {
				if( (new_ops[j].sign == '+') || (new_ops[j].sign == '-') || (new_ops[j].sign == 'c') || (new_ops[j].sign == 'v' ) )
				{
					if( new_ops[j].srcStart >= (cut_point + del_len)) {
						new_ops[j].srcStart = new_ops[j].srcStart - del_len;
					}
					else if( new_ops[j].srcStart >= cut_point ) {
						new_ops[j].srcStart = cut_point;
					}
			
					if( new_ops[j].srcEnd >= (cut_point + del_len)) {
						new_ops[j].srcEnd = new_ops[j].srcEnd - del_len;
					}
					else if( new_ops[j].srcEnd >= cut_point ) {
						new_ops[j].srcEnd = cut_point;
					}

					if( new_ops[j].dstStart >= (cut_point + del_len)) {
						new_ops[j].dstStart = new_ops[j].dstStart - del_len;
					}
					else if( new_ops[j].dstStart >= cut_point ) {
						new_ops[j].dstStart = cut_point;
					}

					if( new_ops[j].dstEnd >= (cut_point + del_len)) {
						new_ops[j].dstEnd = new_ops[j].dstEnd - del_len;
					}
					else if( new_ops[j].dstEnd >= cut_point ) {
						new_ops[j].dstEnd = cut_point;
					}
				}
				else if(new_ops[j].sign == 'd')  
				{
					if( new_ops[j].srcStart >= (cut_point + del_len)) {
						new_ops[j].srcStart = new_ops[j].srcStart - del_len;
					}
					else if( new_ops[j].srcStart >= cut_point ) {
						new_ops[j].srcStart = cut_point;
					}

					if( new_ops[j].srcEnd >= (cut_point + del_len)) {
						new_ops[j].srcEnd = new_ops[j].srcEnd - del_len;
					}
					else if( new_ops[j].srcEnd >= cut_point ) {
						new_ops[j].srcEnd = cut_point;
					}
				}
				else {
					fatalf("unallowed operation in the new_ops list: %d\n", new_ops[j].sign);
				}
			}
		}
	}

	j = 0;
	for( i = 0; i < num_new_ops; i++ ) {
		if( (new_ops[i].sp_id == sp_id) || (sp_id == SELF_PAIR) ) {
			new_ops[j] = assign_ops(new_ops[i]);
			if( new_ops[i].sp_id == SELF2 ) {
				new_ops[j].srcStart = new_ops[j].srcStart - seq1_len;
				new_ops[j].srcEnd = new_ops[j].srcEnd - seq1_len;
				if( new_ops[j].sign != 'd' ) {
					new_ops[j].dstStart = new_ops[j].dstStart - seq1_len;
					new_ops[j].dstEnd = new_ops[j].dstEnd - seq1_len;
				}
			}
			j++;
		}
	}
	num_new_ops = j;

	return(num_new_ops);
}

void update_cur_pos_ops(struct ops_list *ops, int *num_ops, int size1)
{
	int num_new1 = 0, num_new2 = 0;
	struct ops_list *new_ops1 = NULL, *new_ops2 = NULL;
	int i = 0, j = 0;

	num_new1 = *num_ops;
	num_new2 = *num_ops;
	if( (*num_ops) > 0 ) {
		new_ops1 = (struct ops_list *) ckalloc(num_new1 * sizeof(struct ops_list));
		new_ops2 = (struct ops_list *) ckalloc(num_new2 * sizeof(struct ops_list));

		init_ops(new_ops1, 0, num_new1);
		init_ops(new_ops2, 0, num_new2);

		num_new1 = cal_cur_pos_ops(num_new1, ops, new_ops1, SELF1, 0);
		num_new2 = cal_cur_pos_ops(num_new2, ops, new_ops2, SELF2, size1);

		j = 0;
		for( i = 0; i < num_new1; i++ ) {
			ops[j] = assign_ops(new_ops1[i]);
			j++;
		}

		for( i = 0; i < num_new2; i++ ) {
			ops[j] = assign_ops(new_ops2[i]);
			j++;
		}
		*num_ops = j;
		free(new_ops1);
		free(new_ops2);
	}
}

