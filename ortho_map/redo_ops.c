#include "main.h"
#include "redo_ops.h"
#include "pred_regions.h"
#include "util_gen.h"
#include "regions.h"
#include "kd_tree.h"
#include "find_dup_copy.h"
#include "ins_dup_copy.h"
#include "extend_slist.h"
#include "handle_tandem_dup.h"
#include "util.h"
#include "util_i.h"
#include "find_merging.h"
#include "update_init_algns.h"
#include "pred_ops.h"
#include "pred_dels.h"
#include "rollback.h"
#include "apply_ops.h"

extern int debug_mode;

void redo_ops(struct ops_list *ops, int num_ops, int *num_algns, struct DotList *algns, struct cv_list *cv, int *num_cv, int size, struct DotList *init_algns, int num_init_algns, FILE *fp, int size1)
{
	struct DotList *ops_algns;
	int num_ops_algns;
	struct slist rm_list;
	int num_ins_regs;
	int *num_suspend_pairs;
	struct ID_List *suspend_list;
	int temp_count = 0;
	struct ops_list *new_ops;
	int *num_new_ops;

	struct kdnode *tree;
	struct perm_pt *p_pts;
	struct kdnode *m_tree;
	struct perm_pt *m_p_pts;

	int i = 0;
	struct I src, dst;
	char ori;
	int cur_num_cv;
	struct gap_list *del_list;
	int num_del;

	ops_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_ops);
	new_ops = (struct ops_list *) ckalloc(sizeof(struct ops_list) * num_ops);
	num_new_ops = (int *) ckalloc(sizeof(int));
	num_suspend_pairs = (int *) ckalloc(sizeof(int));
	suspend_list = (struct ID_List *) ckalloc(sizeof(struct ID_List) * (*num_algns));
	del_list = (struct gap_list *) ckalloc(sizeof(struct gap_list) * (*num_algns));

	*num_new_ops = 0;
	cur_num_cv = *num_cv;
	for( i = 0; i < num_ops; i++ ) {
		tree = (struct kdnode *) ckalloc(sizeof(struct kdnode));
		p_pts = (struct perm_pt *) ckalloc((*num_algns) * sizeof(struct perm_pt));
		m_tree = (struct kdnode *) ckalloc(sizeof(struct kdnode));
		m_p_pts = (struct perm_pt *) ckalloc((*num_algns) * sizeof(struct perm_pt));
		temp_count++;

		ori = ops[i].sign;
		src = assign_I(ops[i].src_b, ops[i].src_e);
		if( ori != 'd' ) dst = assign_I(ops[i].dst_b, ops[i].dst_e);
		else dst = assign_I(0,1);
		
		num_ops_algns = make_ops_algns(ops_algns, ops, i, num_ops-1);
		assign_perm(p_pts, (*num_algns), algns, LEFT);

		tree = build_kd(p_pts, 0, (*num_algns)-1);

//		handle_tandem_dup_for_redo(algns, num_algns, init_algns, ops_algns, num_ops_algns);
		handle_tandem_dup(algns, num_algns, init_algns);

		assign_perm(m_p_pts, (*num_algns), algns, LEFT);  
		m_tree = build_kd(m_p_pts, 0, (*num_algns)-1);
		num_ins_regs = find_dup_copy(algns, num_algns, m_tree, m_p_pts, size, suspend_list, fp, init_algns);
		*num_suspend_pairs = num_ins_regs;

		if( ori == 'd' ) {
			num_del = 0;
			num_del = find_redo_del_list(src, algns, num_algns, m_tree, m_p_pts, size, del_list);
			if( num_del > 0 ) {
				rm_list = rollback_step_del(num_del, del_list, num_algns, algns, num_new_ops, new_ops, init_algns, fp, src.lower, width(src), size1);
    		cur_num_cv = update_conv_del(src, cv, cur_num_cv);
			}
			else {
				rm_list = rollback_step_del(num_del, del_list, num_algns, algns, num_new_ops, new_ops, init_algns, fp, src.lower, width(src), size1);
    		cur_num_cv = update_conv_del(src, cv, cur_num_cv);
				if(debug_mode == TRUE) printf("Nothing is redone in %d-%d for del", src.lower, src.upper);
			}

			rm_list.id = DEL_EXIST;
		}
		else {
			rm_list = find_mapping_algn(ori, src, dst, algns, *num_algns, size, suspend_list, num_suspend_pairs, fp, init_algns, tree, p_pts);

	  	if( (rm_list.id != NO_EXIST) &&( rm_list.id != SP_EVENT) )
  		{
   	 		if( num_ins_regs > 0 )
    		{
      		extend_slist(algns, num_algns, m_tree, m_p_pts, size, suspend_list, num_ins_regs, rm_list.id, rm_list.is_x, fp, init_algns);
    		}
  		}
		}

	  if( (rm_list.id != SP_EVENT) && ((rm_list.id != NO_EXIST) && (rm_list.id != DEL_EXIST)) && ((rm_list.sp_state == INSERTED_IN_X) || (rm_list.sp_state == INSERTED_IN_Y) || (rm_list.sp_state == INS_IN_X_TAN) || (rm_list.sp_state == INS_IN_Y_TAN)) )
 	 	{ 
 	  	num_ins_regs = *num_suspend_pairs;
  	  rollback_ins_dup(rm_list.sp_state, rm_list.id, algns, suspend_list, num_ins_regs, init_algns, fp); // when the inferred duplication occurs in an inserted region of the suspend list
		}
    
    if( rm_list.id == NO_EXIST ) {}
    else if( rm_list.id == DEL_EXIST )
    {
      size = size + rm_list.val;
    }
    else
    { 
      rollback_init_dots(algns, rm_list.id, rm_list.is_x, *num_algns, init_algns, num_init_algns, fp, *num_new_ops, new_ops, FIRST_RUN, size1);
			
    	cur_num_cv = update_conv(algns, rm_list.id, rm_list.is_x, cv, cur_num_cv);
		 	predict_op(rm_list.is_x, rm_list.id, num_algns, algns, rm_list.sp_state, *num_new_ops, new_ops);
      *num_new_ops = (*num_new_ops) + 1;
    } 

		free(m_p_pts);
		free_kd(m_tree);
		free(p_pts);
		free_kd(tree);
	}	
	
	*num_cv = cur_num_cv;
	free(ops_algns);
	free(del_list);
	free(suspend_list);
	free(num_suspend_pairs);
	free(num_new_ops);
	free(new_ops);
}

int search_algns_for_op(char ori, struct I src, struct I dst, struct DotList *algns, int num_algns, bool *is_x)
{
  struct DotList *temp_algns; // the sequence starts at 1 and each range of a segment is a form of [a, b), i.e. the nucleotide at 'a' is included and one at 'b' is not
  struct slist *st;
  int i, j;
  int b, e;
  int count;
	int sign;
	int res = -1;
	int min_diff, cur_diff;
	int min_id = -1;
	int cur_id;

	min_diff = width(src) + width(dst) + 1;

	if( ori == '+' ) sign = 0;
	else if( ori == '-' ) sign = 1;

  temp_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);
  st = (struct slist *) ckalloc(sizeof(struct slist) * num_algns);

	initialize_slist(st, 0, num_algns);
  sort_init_algns(st, algns, num_algns, INIT_SELF1);
  b = search_range_b(st, algns, num_algns, src.lower, SELF1);
  e = search_range_e(st, algns, num_algns, src.upper, SELF1);

  j = 0;
  for( i = b; i <= e; i++ ) {
    assign_algn(temp_algns, j, algns[st[i].id]);
		temp_algns[j].l_id = st[i].id;
    j++;
  }
  count = j;

	for( i = 0; i < num_algns; i++ ) st[i].id = i;
  sort_init_algns(st, temp_algns, count, INIT_SELF2);
  b = search_range_b(st, temp_algns, count, dst.lower, SELF2);
  e = search_range_e(st, temp_algns, count, dst.upper, SELF2);

  j = 0;
  for( i = b; i <= e; i++ ) {
		cur_id = st[i].id;
		if( sign == temp_algns[cur_id].sign ) {
			if( (almost_equal(temp_algns[cur_id].x, src) == true) && (almost_equal(temp_algns[cur_id].y, dst) == true) ) {
				cur_diff = abs(temp_algns[cur_id].x.lower-src.lower) + abs(temp_algns[cur_id].x.upper-src.upper) + abs(temp_algns[cur_id].y.lower-dst.lower) + abs(temp_algns[cur_id].y.upper-dst.upper);
				if( (min_diff == -1) || (min_diff > cur_diff) ) {
					min_diff = abs(temp_algns[cur_id].x.lower-src.lower) + abs(temp_algns[cur_id].x.upper-src.upper) + abs(temp_algns[cur_id].y.lower-dst.lower) + abs(temp_algns[cur_id].y.upper-dst.upper);
					min_id = cur_id;
					*is_x = false;
					res = temp_algns[cur_id].l_id;
				}
			}
			else if( (almost_equal(temp_algns[cur_id].x, dst) == true) && (almost_equal(temp_algns[cur_id].y, src) == true) ) {
				cur_diff = abs(temp_algns[cur_id].x.lower-dst.lower) + abs(temp_algns[cur_id].x.upper-dst.upper) + abs(temp_algns[cur_id].y.lower-src.lower) + abs(temp_algns[cur_id].y.upper-src.upper);
				if( (min_diff == -1) || (min_diff > cur_diff) ) {
					min_diff = abs(temp_algns[cur_id].x.lower-dst.lower) + abs(temp_algns[cur_id].x.upper-dst.upper) + abs(temp_algns[cur_id].y.lower-src.lower) + abs(temp_algns[cur_id].y.upper-src.upper);
					min_id = cur_id;
					*is_x = true;
					res = temp_algns[cur_id].l_id;
				}
			}
		}	
  }

	if( min_diff >  width(src) + width(dst)) res = -1;

  free(temp_algns);
  free(st);
  return(res);
}

struct slist find_mapping_algn(char ori, struct I src, struct I dst, struct DotList *alg_self, int num_algn_self, int size, struct ID_List *suspend, int *num_suspend_pairs, FILE *fp, struct DotList *init_algns, struct kdnode *tree, struct perm_pt *p_pts)
{
	int cur_id = -1;
	int *num_x, *num_y;
	int num_ins_regs;
	int pred_op;
	int *add_info;
	struct slist res_list;
	int sp_state;
	bool *is_x;
	int threshold = LOOSE;
	struct I tmp_reg_x, tmp_reg_y;

	num_x = (int *) ckalloc(sizeof(int));
	num_y = (int *) ckalloc(sizeof(int));
	add_info = (int *) ckalloc(sizeof(int));
	is_x = (bool *) ckalloc(sizeof(bool));

	num_ins_regs = *num_suspend_pairs;
	if( dst.lower > src.lower ) {
		tmp_reg_x = assign_I(src.lower, src.upper);
		tmp_reg_y = assign_I(dst.lower, dst.upper);
		cur_id = search_algns_for_op(ori, tmp_reg_x, tmp_reg_y, alg_self, num_algn_self, is_x);
	}
	else {
		tmp_reg_x = assign_I(dst.lower, dst.upper);
		tmp_reg_y = assign_I(src.lower+width(dst), src.upper+width(dst));
		cur_id = search_algns_for_op(ori, tmp_reg_x, tmp_reg_y, alg_self, num_algn_self, is_x);
		if( (*is_x) == true ) *is_x = false; // y region is removed
		else *is_x = true; // x region is removed
	}
	if( cur_id == -1 ) fatalf("corresponding alignment not found for %c: %d-%d, %d-%d\n", ori, src.lower, src.upper, dst.lower, dst.upper);

	sp_state = check_ins_dup_copy(cur_id, suspend, num_ins_regs, alg_self, ANCESTRAL, add_info);
	if( sp_state == INSERTED_IN_BOTH )
	{
		pred_op = is_left_to_right_count(num_x, num_y, cur_id, num_algn_self, alg_self, threshold, suspend, num_ins_regs, TANDEM_CHECK, tree, p_pts, size, fp, init_algns);
    if( (pred_op == OVERLAP_IN_X) || (pred_op == NO_OVERLAP) )        {
    	if( pred_op == NO_OVERLAP )
      {
        if( (*add_info) == SP_OVERLAP_IN_Y) sp_state = INS_IN_X_TAN;
        else sp_state = INS_IN_Y_TAN;
      }
      else sp_state = INS_IN_Y_TAN;
    }
    else if( pred_op == OVERLAP_IN_Y )
    {
      sp_state = INS_IN_X_TAN;
    }
	}

	if( cur_id != -1 ) {
		alg_self[cur_id].x = assign_I(tmp_reg_x.lower, tmp_reg_x.upper);
		alg_self[cur_id].y = assign_I(tmp_reg_y.lower, tmp_reg_y.upper);
		res_list.id = cur_id;
		res_list.is_x = *is_x;
		res_list.val = alg_self[cur_id].identity;
		res_list.sp_state = sp_state;
	}
	else res_list.id = -1;

	free(is_x);
	free(num_x);
	free(num_y); 
	free(add_info);

	return(res_list);	
}

int make_ops_algns(struct DotList *ops_algns, struct ops_list *ops, int b, int e)
{
	int i, j;
	int dup_len;
	int num_ops_algns;
	int p1, p2;

	j = 0;
	for( i = b; i <= e; i++ )
	{
		if((ops[i].sign == '+') || (ops[i].sign == '-')) {
			ops_algns[j].y = assign_I(ops[i].dst_b, ops[i].dst_e);
			dup_len = abs(ops[i].dst_e - ops[i].dst_b);
			if( ops[i].src_b >= ops[i].dst_b ) {
				ops_algns[j].x = assign_I(ops[i].src_b + dup_len, ops[i].src_e + dup_len);
			}
			else ops_algns[j].x = assign_I(ops[i].src_b, ops[i].src_e);

			if(ops[i].sign == '+') ops_algns[j].sign = 0;
			else ops_algns[j].sign = 1;
			ops_algns[j].l_id = i; // jth alignment corresponds ith ops
			j++;
		}
	}
	num_ops_algns = j;

	if( num_ops_algns > 0 )
	{
		for( j = 0; j < num_ops_algns; j++ ) {
			p1 = cal_init_position(ops, ops_algns[j].l_id, b, ops_algns[j].x.lower);
			p2 = cal_init_position(ops, ops_algns[j].l_id, b, ops_algns[j].x.upper);
			ops_algns[j].x = assign_I(p1, p2);
			p1 = cal_init_position(ops, ops_algns[j].l_id, b, ops_algns[j].y.lower);
			p2 = cal_init_position(ops, ops_algns[j].l_id, b, ops_algns[j].y.upper);
			ops_algns[j].y = assign_I(p1, p2);
		}
	}

	return(num_ops_algns);
}
