#include "main.h"
#include "sec_round.h"
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
#include "apply_ops.h"
#include "rollback.h"
#include "pred_dels.h"

extern int debug_mode;

void reorder_dups(struct slist *sorted, int *num_alg_self, struct DotList *alg_self, int size, int threshold, struct DotList *init_algns, int num_init_algns, FILE *fp, float avg_pid, struct cv_list *cv, int *num_cv, struct ops_list *ops, int *num_ops, struct exons_list *exons, int *num_exons, struct exons_list *genes, int *num_genes, int *old_dups, int num_old_dups, int run_mode, int size1)
{
	struct slist rm_list;
	int num_ins_regs = 0;
	int *num_suspend_pairs;
	struct ID_List *suspend_list;
	struct DotList *candi_algns;
	struct DotList *pair_algns;
	int *num_pair;
	int num_candi = 0;
	int opt_id = -1, i = 0;
	int temp_count = 0;
	int org_id = 0;

	struct kdnode *tree;
	struct perm_pt *p_pts;
	struct kdnode *m_tree;
	struct perm_pt *m_p_pts;

	int mode = BEFORE_SP; // before inferring speciation
	int cur_num_cv = 0;
	int cur_num_exons = 0;
	int cur_num_genes = 0;
	int num_del = 0;
	struct gap_list *del_list;
	int h_pid = 0; // current highest identity
	struct I del_reg;

	del_reg = assign_I(0, 1);
	cur_num_cv = *num_cv;
	cur_num_exons = *num_exons;
	cur_num_genes = *num_genes;

	i = 0;
	h_pid = alg_self[sorted[i].id].identity;
	while( (i < (*num_alg_self)) && (width(alg_self[sorted[i].id].x) < EFFEC_DEL_LEN_TH) ) i++;

	if( (i >= (*num_alg_self)) || ((i < (*num_alg_self)) && ((h_pid-PID_DIFF) > alg_self[sorted[i].id].identity)) ) {
		h_pid = h_pid - PID_DIFF;
	}
	else h_pid = alg_self[sorted[i].id].identity-2;

	pair_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * (*num_alg_self));
	candi_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * (*num_alg_self));
	num_suspend_pairs = (int *) ckalloc(sizeof(int));
	suspend_list = (struct ID_List *) ckalloc(sizeof(struct ID_List) * (*num_alg_self));
	del_list = (struct gap_list *) ckalloc(sizeof(struct gap_list) * (*num_alg_self));
	num_pair = (int *) ckalloc(sizeof(int));

  initialize_algns(pair_algns, 0, *num_alg_self);
  initialize_algns(candi_algns, 0, *num_alg_self);

  for( i = 0; i < (*num_alg_self); i++ ) { // initialization for ops
    suspend_list[i].is_x = true;
    suspend_list[i].m_id = -1;
    suspend_list[i].left_id = -1;
    suspend_list[i].right_id = -1;
    suspend_list[i].f_is_x = true;
    suspend_list[i].is_t_ins = true;
		del_list[i].gid = -1;
		del_list[i].type = -1;
		del_list[i].id1 = -1;
		del_list[i].id2 = -1;
		del_list[i].x1 = 0;
		del_list[i].x2 = 1;
		del_list[i].y1 = 0;
		del_list[i].y2 = 1;
		del_list[i].offset = 0;
  }

  *num_ops = 0;
  *num_suspend_pairs = 0;

	if( debug_mode == TRUE ) printf("avg pid: %f\n", avg_pid);
	opt_id = 0;
	while ( (opt_id != NO_EXIST) && (mode == BEFORE_SP) && ((*num_alg_self) > 0) ) {
		temp_count++;

//		tree = (struct kdnode *) ckalloc(sizeof(struct kdnode));
		p_pts = (struct perm_pt *) ckalloc((*num_alg_self) * sizeof(struct perm_pt));

		assign_perm(p_pts, (*num_alg_self), alg_self, LEFT);
		tree = build_kd(p_pts, 0, (*num_alg_self)-1);

		num_del = 0;
	  num_del = pred_dels(alg_self, num_alg_self, tree, p_pts, size, del_list, h_pid);
	  if( num_del > 0 ) num_del = remove_false_dels(alg_self, del_list, num_del, suspend_list, num_ins_regs);
		
		if( num_del > 0 ) {
			if( debug_mode == TRUE ) {
				for( i = 0; i < num_del; i++ ) {
		      printf("Del candi: %d-%d, %d-%d and %d-%d, %d-%d\n", alg_self[del_list[i].id1].x.lower, alg_self[del_list[i].id1].x.upper, alg_self[del_list[i].id1].y.lower, alg_self[del_list[i].id1].y.upper, alg_self[del_list[i].id2].x.lower, alg_self[del_list[i].id2].x.upper, alg_self[del_list[i].id2].y.lower, alg_self[del_list[i].id2].y.upper);
				}
			}
    	rm_list = rollback_step_del(num_del, del_list, num_alg_self, alg_self, num_ops, ops, init_algns, fp, -1, -1, size1);
			del_reg = assign_I(rm_list.val_red, rm_list.val_red + rm_list.val);
			cur_num_cv = update_conv_del(del_reg, cv, cur_num_cv);
			cur_num_exons = update_exons_del(del_reg, exons, cur_num_exons);
			cur_num_genes = update_exons_del(del_reg, genes, cur_num_genes);
			rm_list.id = DEL_EXIST;
		}
		else {
			handle_tandem_dup(alg_self, num_alg_self, init_algns);
//		m_tree = (struct kdnode *) ckalloc(sizeof(struct kdnode));
			m_p_pts = (struct perm_pt *) ckalloc((*num_alg_self) * sizeof(struct perm_pt));
			assign_perm(m_p_pts, (*num_alg_self), alg_self, LEFT);  
			m_tree = build_kd(m_p_pts, 0, (*num_alg_self)-1);
			num_ins_regs = find_dup_copy(alg_self, num_alg_self, m_tree, m_p_pts, size, suspend_list, fp, init_algns);
			*num_suspend_pairs = num_ins_regs;

			for( i = 0; i < (*num_alg_self); i++ ) sorted[i].id = i;
			sort_by_pid(sorted, alg_self, *num_alg_self);
		
			num_candi = choose_candi_algns_cmp_pid(candi_algns, sorted, alg_self, *num_alg_self, pair_algns, num_pair, old_dups, num_old_dups);

			rm_list = latest_dup(candi_algns, num_candi, size, threshold, num_suspend_pairs, suspend_list, mode, tree, p_pts, alg_self, *num_alg_self, init_algns, num_init_algns, fp, pair_algns, *num_pair, exons, cur_num_exons, genes, cur_num_genes, avg_pid, ops, *num_ops);
			opt_id = rm_list.id;

		  if( (rm_list.id != NO_EXIST) &&( rm_list.id != SP_EVENT) )
 		 	{
 		   	if( num_ins_regs > 0 )
 		   	{
 		     	extend_slist(alg_self, num_alg_self, m_tree, m_p_pts, size, suspend_list, num_ins_regs, rm_list.id, rm_list.is_x, fp, init_algns);
 		   	}
 		 	}

		  if( (rm_list.id != SP_EVENT) && ((rm_list.id != NO_EXIST) && (rm_list.id != DEL_EXIST)) && ((rm_list.sp_state == INSERTED_IN_X) || (rm_list.sp_state == INSERTED_IN_Y) || (rm_list.sp_state == INS_IN_X_TAN) || (rm_list.sp_state == INS_IN_Y_TAN)) )
 		 	{ 
 		  	num_ins_regs = *num_suspend_pairs;
 		 	  rollback_ins_dup(rm_list.sp_state, rm_list.id, alg_self, suspend_list, num_ins_regs, init_algns, fp); // when the inferred duplication occurs in an inserted region of the suspend list
			}
    
 		  if( rm_list.id == NO_EXIST ) opt_id = NO_EXIST;
 		  else if( rm_list.id == DEL_EXIST )
	  	{
	      size = size + rm_list.val;
	    }
 	  	else
    	{ 
				if( rm_list.id != sorted[0].id ) h_pid = alg_self[sorted[0].id].identity - 2;
				else h_pid = alg_self[sorted[1].id].identity - 2;

				if( run_mode == INF_DUP ) {
					print_conv_on_dup((int)(avg_pid+0.5), alg_self[rm_list.id], cv, cur_num_cv, num_ops, ops, run_mode); // cv is sorted by similarity levels
				}
				cur_num_cv = update_conv(alg_self, rm_list.id, rm_list.is_x, cv, cur_num_cv);
				cur_num_exons = update_exons(alg_self, rm_list.id, rm_list.is_x, exons, cur_num_exons);
				cur_num_genes = update_exons(alg_self, rm_list.id, rm_list.is_x, genes, cur_num_genes);
   		  rollback_init_dots(alg_self, rm_list.id, rm_list.is_x, *num_alg_self, init_algns, num_init_algns, fp, *num_ops, ops, SECOND_RUN, size1);
				org_id = alg_self[rm_list.id].index;
				update_algns_sign(alg_self, rm_list.id, *num_alg_self, init_algns);
   		  predict_op(rm_list.is_x, rm_list.id, num_alg_self, alg_self, rm_list.sp_state, *num_ops, ops);
				ops[*num_ops].id = org_id; // ops[].id changed in #173, so it is assigned again
				if( ops[*num_ops].srcStart <= ops[*num_ops].dstStart ) {
					ops[*num_ops].ctg_id1 = init_algns[org_id].ctg_id1; // ops[].id changed in #173, so it is assigned again
					ops[*num_ops].ctg_id2 = init_algns[org_id].ctg_id2; // ops[].id changed in #173, so it is assigned again
				}
				else {
					ops[*num_ops].ctg_id1 = init_algns[org_id].ctg_id2; // ops[].id changed in #173, so it is assigned again
					ops[*num_ops].ctg_id2 = init_algns[org_id].ctg_id1; // ops[].id changed in #173, so it is assigned again
				}
   		  *num_ops = (*num_ops) + 1;
    	} 

			if( (*num_alg_self) == 0 ) {
				opt_id = NO_EXIST;
			}
			free(m_p_pts);
			free_kd(m_tree);
		}
		free(p_pts);
		free_kd(tree);
	}	

	*num_cv = cur_num_cv;
	*num_exons = cur_num_exons;
	*num_genes = cur_num_genes;
	free(del_list);
	free(num_pair);
	free(pair_algns);
	free(candi_algns);
	free(suspend_list);
	free(num_suspend_pairs);
}

int choose_candi_algns_pid(struct DotList *candi_algns, struct slist *sorted, struct DotList *algns, int num_algns, struct DotList *pair_algns, int *num_pair)
{
	int i, j;
	int highest_pid = 0;
	int num_candi = 0;
	int count = 0;

	j = 0;
	i = 0;
	while( (i < num_algns) && ((algns[sorted[i].id].sp_id == PAIR) || (algns[sorted[i].id].sign == DELETED)) ) {
		if( (algns[sorted[i].id].sp_id == PAIR) && (algns[sorted[i].id].sign != DELETED) ) {
			assign_algn(pair_algns, count, algns[sorted[i].id]);
			count++;
		}
		i++;
	}

	if( i >= num_algns ) num_candi = 0;
	else {
		highest_pid = algns[sorted[i].id].identity;
		while( (i < num_algns) && ((algns[sorted[i].id].identity) >= (highest_pid - PID_DIFF)) ) {
			if( (algns[sorted[i].id].sign != DELETED) && (algns[sorted[i].id].sp_id != PAIR) ) {
				assign_algn(candi_algns, j, algns[sorted[i].id]);
				candi_algns[j].l_id = sorted[i].id; // l_id is used for a different purpose at this time
				j++;
			}
			else if( (algns[sorted[i].id].sign != DELETED) && (algns[sorted[i].id].sp_id == PAIR) ) {
				assign_algn(pair_algns, count, algns[sorted[i].id]);
				count++;
			}
			i++;
		}
		num_candi = j;
	}

	*num_pair = count;
	return(num_candi);
}

int choose_candi_algns_cmp_pid(struct DotList *candi_algns, struct slist *sorted, struct DotList *algns, int num_algns, struct DotList *pair_algns, int *num_pair, int *old_dups, int num_old_dups)
{
	int i, j;
	int num_candi = 0;
	int count = 0;
	bool *is_exist;
	bool is_in_old_dup = false;

	is_exist = (bool *) ckalloc(sizeof(bool));
	*is_exist = false;
	j = 0;
	i = 0;
	while( i < num_algns ) {
		if( (algns[sorted[i].id].sp_id == PAIR) && (algns[sorted[i].id].sign != DELETED) ) {
			assign_algn(pair_algns, count, algns[sorted[i].id]);
			count++;
		}
		i++;
	}

	i = 0;
	while( i < num_algns ) {
/*
		if( (algns[sorted[i].id].sp_id != PAIR) && (algns[sorted[i].id].sign != DELETED) && (is_in(algns[sorted[i].id].index, old_dups, num_old_dups) == true) ) {
			if( debug_mode == TRUE ) {
				printf("%d-%d, %d-%d\n", algns[sorted[i].id].x.lower, algns[sorted[i].id].x.upper, algns[sorted[i].id].y.lower, algns[sorted[i].id].y.upper);
			}
		}
		else if( (algns[sorted[i].id].sp_id != PAIR) && (algns[sorted[i].id].sign != DELETED ) && ((is_overlap_ortho_algns(algns[sorted[i].id], true, pair_algns, count, is_exist) == false) || (is_overlap_ortho_algns(algns[sorted[i].id], false, pair_algns, count, is_exist) == false))) {
*/
		if( (algns[sorted[i].id].sp_id != PAIR) && (is_tandem(algns[sorted[i].id]) == true) ) {
		}
		else if( (algns[sorted[i].id].sp_id != PAIR) && (algns[sorted[i].id].sign != DELETED) && (is_in(algns[sorted[i].id].index, old_dups, num_old_dups) == true) ) {
			is_in_old_dup = true;
		}
		
		if( (algns[sorted[i].id].sp_id != PAIR) && (algns[sorted[i].id].sign != DELETED ) && ((is_overlap_ortho_algns(algns[sorted[i].id], true, pair_algns, count, is_exist, is_in_old_dup) == false) || (is_overlap_ortho_algns(algns[sorted[i].id], false, pair_algns, count, is_exist, is_in_old_dup) == false))) {
			assign_algn(candi_algns, j, algns[sorted[i].id]);
			candi_algns[j].l_id = sorted[i].id; // l_id is used for a different purpose at this time
			j++;
		}
		i++;
	}

	num_candi = j;
	*num_pair = count;
	free(is_exist);
	return(num_candi);
}

bool is_overlap_ortho_algns(struct DotList cur_algn, bool is_x, struct DotList *pair_algns, int num_pair, bool *is_exist, bool is_in_old_dup)
{
	bool res = false;
	int i;
	struct I reg, tmp;

	*is_exist = false;
	if( is_x == true ) tmp = assign_I(cur_algn.x.lower, cur_algn.x.upper);
	else tmp = assign_I(cur_algn.y.lower, cur_algn.y.upper);

	i = 0;
	while(i < num_pair) {
		
		if( cur_algn.sp_id == SELF1 ) reg = assign_I(pair_algns[i].x.lower, pair_algns[i].x.upper);
		else if( cur_algn.sp_id == SELF2 ) reg = assign_I(pair_algns[i].y.lower, pair_algns[i].y.upper);
		else fatalf("unexpected candidate alignment: %d-%d, %d-%d:%d\n", cur_algn.x.lower, cur_algn.x.upper, cur_algn.y.lower, cur_algn.y.upper, cur_algn.sp_id);

		if( ( debug_mode == true ) && ( proper_overlap(tmp, reg) == true ) ) {
			debug_mode = true; // nothing happens
		}

/*
		if( (strict_almost_equal(tmp, reg) == true) && ( cur_algn.identity <= pair_algns[i].identity) ) {
			res = true;
		}
		else if( (is_tandem(cur_algn) == false) && (strict_subset(tmp, reg) == true) ) {
			res = true;
		}
		else if( strict_almost_equal(tmp, reg) == true) {
			res = false;
		}
		else if( (f_loose_subset(reg, tmp, STRICT) == true) && ( cur_algn.identity <= pair_algns[i].identity )) {
			res = true;
		}
		else if( strict_subset(reg, tmp) == true ) {
			res = false;
		}
		else if( (is_tandem(cur_algn) == false) && (f_loose_subset(tmp, reg, STRICT) == true)) {
			res = true;
		}
		else if( f_loose_subset(reg, tmp, STRICT) == true ) {
			res = false;
		}
		else if( (f_loose_overlap(tmp, reg, STRICT) == true ) || ( (proper_overlap(tmp, reg) == true) && (width(intersect(tmp, reg)) >= (width(reg)/2)) && (width(intersect(tmp, reg)) >= (width(tmp)/2))) ) 
		{
			if( cur_algn.identity <= pair_algns[i].identity ) res = true;
			else {
				res = false;
			}
		}
		else {
			res = false;
		}
*/

		if( (strict_almost_equal(tmp, reg) == true) || (strict_almost_equal(reg, tmp) == true) ) {
			if( cur_algn.identity < pair_algns[i].identity ) {
				res = true;
				*is_exist = true;
			}
		}
		else if( f_loose_subset(reg, tmp, STRICT) == true ) {
			if( cur_algn.identity < pair_algns[i].identity ) {
				res = true;
				*is_exist = true;
			}
		}
		else if( (f_loose_overlap(tmp, reg, STRICT) == true ) || ( (proper_overlap(tmp, reg) == true) && (width(intersect(tmp, reg)) >= (width(reg)/2)) && (width(intersect(tmp, reg)) >= (width(tmp)/2))) ) 
		{
			if( is_in_old_dup == true ) {
				res = true;
				*is_exist = true;
			}
			else if( is_tandem(cur_algn) == false ) {
/*
					res = true;
					if( cur_algn.identity < pair_algns[i].identity ) {
						*is_exist = true;
					}
*/
				if( both_ends_overlap(tmp, reg, STRICT) == true ) {
					res = true;
					if( cur_algn.identity < pair_algns[i].identity ) {
						*is_exist = true;
					}
				}
				else {
					if( cur_algn.identity < pair_algns[i].identity ) {
						res = true;
						*is_exist = true;
					}
				}
			}
			else {
				if( cur_algn.identity < pair_algns[i].identity ) {
					res = true;
					*is_exist = true;
				}
			}
		}
		i++;
	}

	return(res);
}

struct slist latest_dup(struct DotList *candi_algns, int num_candi, int size, int threshold, int *num_suspend_pairs, struct ID_List *suspend, int mode, struct kdnode *tree, struct perm_pt *p_pts, struct DotList *alg_self, int num_algn_self, struct DotList *init_algns, int num_init_algns, FILE *fp, struct DotList *pair_algns, int num_pair, struct exons_list *exons, int num_exons, struct exons_list *genes, int num_genes, float avg_pid, struct ops_list *ops, int num_ops)
{
	struct slist *op_info;
	int i = 0, j = 0;
	int cur_id = 0;
	int *num_x = NULL, *num_y = NULL;
	int num_ins_regs = 0;
	int pred_op = NO_OVERLAP;
	int *add_info = NULL;
	int opt_id = -1;
	struct slist res_list = {0, 0, 0, NO_OVERLAP, NO_OVERLAP, true};
	int sp_state = NO_OVERLAP;
	int h_pid = 0;
	int side = TIE;
	int run_no = FIRST_RUN;
	int num_overlap = 0, pid_th = 0;
	bool *is_exist = NULL;
	int index = -1;
	int critical_side = TIE;

	if( num_candi > 0 ) {
		op_info = (struct slist *) ckalloc(sizeof(struct slist) * num_candi);
		initialize_slist(op_info, 0, num_candi);
		for( i = 0; i < num_candi; i++ ) {
			op_info[i].id = -1;
			op_info[i].val = -1;
			op_info[i].val_red = -1;
			op_info[i].sp_state = -1;
			op_info[i].add_sp_state = -1;
			op_info[i].is_x = true;
		}
	}
	num_x = (int *) ckalloc(sizeof(int));
	num_y = (int *) ckalloc(sizeof(int));
	add_info = (int *) ckalloc(sizeof(int));
	is_exist = (bool *) ckalloc(sizeof(bool));

	res_list.id = -1;
	res_list.val = 0;
	res_list.val_red = 0;
	res_list.sp_state = NO_OVERLAP;
	res_list.add_sp_state = NO_OVERLAP;
	res_list.is_x = true;

	*num_x = 0;
	*num_y = 0;
	*add_info = 0;

	num_ins_regs = *num_suspend_pairs;
//	h_pid = candi_algns[0].identity;
	h_pid = (int)avg_pid;
	i = 0;
	j = 0;
//	while( (i < num_candi) && (opt_id == -1) ) {
	while( i < num_candi ) {
		*is_exist = false;
		cur_id = candi_algns[i].l_id;
		alg_self[cur_id].identity = init_algns[alg_self[cur_id].index].identity;
		sp_state = check_ins_dup_copy(cur_id, suspend, num_ins_regs, alg_self, mode, add_info);
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
		else if( is_tandem(alg_self[cur_id]) ) {
			pred_op = is_left_to_right_count(num_x, num_y, cur_id, num_algn_self, alg_self, threshold, suspend, num_ins_regs, TANDEM_CHECK, tree, p_pts, size, fp, init_algns);
		}
		else pred_op = is_left_to_right_count(num_x, num_y, cur_id, num_algn_self, alg_self, threshold, suspend, num_ins_regs, GENERAL_CHECK, tree, p_pts, size, fp, init_algns); // a flag is temporary

  	if( (pred_op == OVERLAP_IN_X) && ( (sp_state == SP_OVERLAP_IN_X) || (sp_state == UNSUSPENDED) || (sp_state == INSERTED_IN_Y) || (sp_state == COVER_INS_SP_IN_Y) || (sp_state == INS_IN_Y_TAN) || (sp_state == INSERTED_IN_BOTH) || (sp_state == COVER_INS_SP_BOTH)) ) op_info[i].val_red = 0;
		else op_info[i].val_red = *num_y;

		if( is_overlap_ortho_algns(alg_self[cur_id], false, pair_algns, num_pair, is_exist, false) == true ) op_info[i].val_red = -1;

		if( (pred_op == OVERLAP_IN_Y) && ( (sp_state == SP_OVERLAP_IN_Y) || (sp_state == UNSUSPENDED) || (sp_state == INSERTED_IN_X) || (sp_state == COVER_INS_SP_IN_X) || (sp_state == INS_IN_X_TAN) || (sp_state == INSERTED_IN_BOTH) || (sp_state == COVER_INS_SP_BOTH)) ) op_info[i].val = 0;
		else op_info[i].val = *num_x;

		if( is_overlap_ortho_algns(alg_self[cur_id], true, pair_algns, num_pair, is_exist, false ) == true ) op_info[i].val = -1;

		if( (*is_exist) == false ) {
			if( (is_tandem(alg_self[cur_id])) && (op_info[i].val > 0) && (op_info[i].val_red > 0) ) {
				if( op_info[i].val >= op_info[i].val_red ) {
					op_info[i].val = NUM_OVERLAP;
					op_info[i].val_red = 0;
				}
				else {
					op_info[i].val_red = NUM_OVERLAP;
					op_info[i].val = 0;
				}
			}

			op_info[j].id = cur_id;
			op_info[j].val = op_info[i].val;
			op_info[j].val_red = op_info[i].val_red;
			op_info[j].sp_state = sp_state;
			op_info[j].add_sp_state = pred_op;
			j++;
		}

		if( (opt_id == -1) && (sp_state != SP_OVERLAP_IN_X)  && (op_info[i].val == 0)	&& (sp_state != SP_OVERLAP_IN_Y) && (op_info[i].val_red == 0) && (alg_self[cur_id].identity >= (h_pid-PID_DIFF_TH_STRICT))) {
			side = TIE;
			side = cmp_ortho_mappings(alg_self[cur_id], pair_algns, num_pair);
			if( side == TIE ) {
				side =  check_exons_bound(alg_self[cur_id], exons, num_exons, genes, num_genes);
				if( side == TIE ) {
					side = check_distance_prev_dup(alg_self[cur_id].sign, alg_self[cur_id].index, init_algns, ops, num_ops);
				}
			}

			if( side == TIE ) {
				res_list.is_x = false;
				ops[num_ops].dir = 0; // both regions can be a source or a target
			}
			else if( side == LEFT_SIDE ) {
				res_list.is_x = true; // LEFT means x region is removed as a duplicated region
				ops[num_ops].dir = 2; // y -> x
			}
			else if( side == RIGHT_SIDE ) {
				res_list.is_x = false;
				ops[num_ops].dir = 1; // x -> y
			}
			
			res_list.id = cur_id;
			opt_id = cur_id;
			index = i;
			res_list.sp_state = sp_state;
		}
		else if((opt_id == -1) && (sp_state != SP_OVERLAP_IN_X) && (op_info[i].val == 0) && (alg_self[cur_id].identity >= (h_pid-PID_DIFF_TH_STRICT))) 	
//		else if((sp_state != SP_OVERLAP_IN_X) && (op_info[i].val == 0))
		{
			opt_id = cur_id;
			index = i;
			res_list.id = cur_id;
			res_list.is_x = true;
			res_list.sp_state = sp_state;
			ops[num_ops].dir = 2; // y -> x
			critical_side = LEFT_SIDE;
		}
		else if( (opt_id == -1) && (sp_state != SP_OVERLAP_IN_Y) && (op_info[i].val_red == 0) && (alg_self[cur_id].identity >= (h_pid-PID_DIFF_TH_STRICT))) 
//		else if( (sp_state != SP_OVERLAP_IN_Y) && (op_info[i].val_red == 0) )
		{
			opt_id = cur_id;
			index = i;
			res_list.id = cur_id;
			res_list.is_x = false;
			res_list.sp_state = sp_state;
			ops[num_ops].dir = 1; // x -> y
			critical_side = RIGHT_SIDE;
		}
		i++;
	}

//// improve this part using machine learning algorithm or statistical method
	num_candi = j;
	num_overlap = NUM_OVERLAP;
	pid_th = h_pid - PID_DIFF_TH_STRICT;
	if( opt_id == -1 ) {
		i = 0;
		run_no = FIRST_RUN;
		while( (i < num_candi) && (opt_id == -1)  ) {
			if( (op_info[i].val != -1) && (op_info[i].val <= num_overlap) && (op_info[i].val_red != -1) && (op_info[i].val_red <= num_overlap) && (alg_self[op_info[i].id].identity >= pid_th)) 
			{
				opt_id = op_info[i].id;
				res_list.id = op_info[i].id;
				res_list.sp_state = op_info[i].sp_state;
				side = TIE;
				side = cmp_ortho_mappings(alg_self[opt_id], pair_algns, num_pair);
				if( side == TIE ) {
					side =  check_exons_bound(alg_self[cur_id], exons, num_exons, genes, num_genes);
					if( side == TIE ) side = check_distance_prev_dup(alg_self[cur_id].sign, alg_self[cur_id].index, init_algns, ops, num_ops);

					if( side == LEFT_SIDE ) {
						res_list.is_x = true;
						ops[num_ops].dir = 2;
					}
					else {
						if( side == TIE ) ops[num_ops].dir = 0;
						else ops[num_ops].dir = 1;
						res_list.is_x = false;
					}
				}
				else {
					if( side == LEFT_SIDE ) {
						res_list.is_x = true; // LEFT means x region is removed as a duplicated region
						ops[num_ops].dir = 2;
					}
					else if( side == RIGHT_SIDE ) {
						res_list.is_x = false;
						ops[num_ops].dir = 1;
					}
				}
			}
			else if( (op_info[i].val != -1) && (op_info[i].val <= num_overlap) && (alg_self[op_info[i].id].identity >= pid_th)) {
				opt_id = op_info[i].id;
				res_list.id = op_info[i].id;
				res_list.is_x = true;
				res_list.sp_state = op_info[i].sp_state;
				ops[num_ops].dir = 2;
			}
			else if( (op_info[i].val_red != -1) && ( op_info[i].val_red <= num_overlap) && (alg_self[op_info[i].id].identity >= pid_th)) {
				opt_id = op_info[i].id;
				res_list.id = op_info[i].id;
				res_list.is_x = false;
				res_list.sp_state = op_info[i].sp_state;
				ops[num_ops].dir = 1;
			}
			i++;
			if( (i >= num_candi) && (opt_id == -1) && (run_no == FIRST_RUN) ) {
				run_no = SECOND_RUN;
				i = 0;
				pid_th++;
				num_overlap++;
			}
		}
	}

	if( opt_id != -1 ) {
		res_list = replace_by_better(alg_self, res_list, op_info, num_candi, pid_th, num_overlap);

		if( critical_side == TIE ) 
		{
			if( (is_tandem(alg_self[res_list.id])) && (num_ops >= 1) ) {
				side = TIE;

				if( res_list.id != -1 ) {
					side = check_prev_tandem_dup(init_algns, num_init_algns, alg_self[res_list.id].index, ops, num_ops);
					if( side == LEFT_SIDE ) {
						res_list.is_x = true;
						ops[num_ops].dir = 2;
					}
					else if( side == RIGHT_SIDE ) {
						res_list.is_x = false;
						ops[num_ops].dir = 1;
					}
				}
			}
		}
	}

	if( num_candi > 0 ) {
		free(op_info);
	}
	free(num_x);
	free(num_y); 
	free(add_info);
	free(is_exist);

	if( opt_id == -1 ) res_list.id = -1;

	return(res_list);	
}

/// need to improve
int cmp_ortho_mappings(struct DotList cur_algn, struct DotList *ortho, int num_ortho)
{
	struct slist *st;
	int i;
	int mode;
	int s_loc, e_loc, s_d, e_d;
	float rate1, rate2;
	int *pid1, *pid2;
	int res = TIE;

	if( cur_algn.sp_id == SELF1 ) mode = INIT_SELF1;
	else if( cur_algn.sp_id == SELF2 ) mode = INIT_SELF2;
	else {
		fatalf("unexpected sp_id: %d in the second round\n", cur_algn.sp_id);
	}

	pid1 = (int *) ckalloc(sizeof(int));
	pid2 = (int *) ckalloc(sizeof(int));
	*pid1 = 0;
	*pid2 = 0;

	if( num_ortho > 0 ) {
		st = (struct slist *) ckalloc(sizeof(struct slist) * num_ortho);	

		initialize_slist(st, 0, num_ortho);
		for( i = 0; i < num_ortho; i++ ) {
			st[i].id = i;
			st[i].val = 0;
			st[i].val_red = 0;
			st[i].is_x = true;
		}
		sort_init_algns(st, ortho, num_ortho, mode);
	}

	if( num_ortho > 0 ) {
		s_loc = search_range_b(st, ortho, num_ortho, cur_algn.x.lower, cur_algn.sp_id);
		e_loc = search_range_e(st, ortho, num_ortho, cur_algn.x.upper, cur_algn.sp_id);
		rate1 = cal_cover_rate(cur_algn.x, st, s_loc, e_loc, ortho,  pid1, cur_algn.sp_id);
	}
	else {
		rate1 = (float) 0;
		*pid1 = 0;
	}

	if( num_ortho > 0 ) {
		s_d = search_range_b(st, ortho, num_ortho, cur_algn.y.lower, cur_algn.sp_id);
		e_d = search_range_e(st, ortho, num_ortho, cur_algn.y.upper, cur_algn.sp_id);
		rate2 = cal_cover_rate(cur_algn.y, st, s_d, e_d, ortho, pid2, cur_algn.sp_id);
	}
	else {
		rate2 = (float) 0;
		*pid2 = 0;
	}

	if( ((int)rate1+(*pid1)) > (((int)rate2+(*pid2))+1) ) res = RIGHT_SIDE;
	else if( (((int)rate1+(*pid1))+1) < ((int)rate2+(*pid2)) ) res = LEFT_SIDE;
	else res = TIE;

	free(pid1);
	free(pid2);

	if( num_ortho > 0 ) {
		free(st);
	}
	return(res);
}

// need to use dynamic programming
float cal_cover_rate(struct I reg, struct slist *st, int b, int e, struct DotList *ortho, int *pid, int sp_id)
{
	int i;
	int cur_id, max_id = -1;
	int len = 0;
	struct I tmp;
	float rate;

	for( i = b; i <= e; i++ ) {
		cur_id = st[i].id;
		if( sp_id == SELF1 ) tmp = ortho[cur_id].x;
		else if( sp_id == SELF2 ) tmp = ortho[cur_id].y;
		if( proper_overlap(tmp, reg) == true ) {
			if( width(intersect(tmp, reg)) > len ) {
				max_id = i;
				len = width(intersect(tmp, reg));
			}
		}
	}

	if( max_id != -1 ) {
		*pid = ortho[st[max_id].id].identity;
		rate = ((float)len / ((float)width(reg))) * ((float)100);
	}
	else {
		*pid = 0;
		rate = (float)0;
	}

	return(rate);
}

bool is_in(int id, int *list, int num_list) 
{
	int i = 0;
	bool res = false;

	while( (i < num_list) && (res == false) ) 
	{
		if( id == list[i] ) res = true;
		i++;
	}

	return(res);
}

int check_distance_prev_dup(int sign, int id, struct DotList *init_algns, struct ops_list *ops, int num_ops)
{
	int cid = -1;
	int old_id = -1;
	int b1 = 0, e1 = 0, b2 = 0, e2 = 0;
	int i = 0;
	int count = 0;
	int min_x1 = -1, min_x2 = -1, min_y1 = -1, min_y2 = -1;
	int min_x = -1, min_y = -1;
	int res = TIE;

  old_id = init_algns[id].index;
  cid = init_algns[id].c_id;
  while( cid != -1 ) {
    old_id = cid;
    cid = init_algns[cid].c_id;
  }
  if( (sign != DELETED) && (init_algns[id].sign == DELETED) ) {
    init_algns[id].sign = sign; // after rolling back a tandem duplication, a converted alignment in the step of handling tandem duplications could be deleted during the rollback process  
  }

  b1 = init_algns[id].x.lower + init_algns[id].xl_diff;
  e1 = init_algns[old_id].x.upper - init_algns[old_id].xr_diff;
  if( init_algns[id].sign == 0 ) {
    b2 = init_algns[id].y.lower + init_algns[id].yl_diff;
    e2 = init_algns[old_id].y.upper - init_algns[old_id].yr_diff;
  }
  else if( init_algns[id].sign == 1 ) {
    b2 = init_algns[old_id].y.lower + init_algns[old_id].yl_diff;
    e2 = init_algns[id].y.upper - init_algns[id].yr_diff;
  }

	if( num_ops == 0 ) {
		res = TIE;
	}
	else {
		for( i = 0; i < num_ops; i++ ) {
			if( init_algns[id].sp_id == ops[i].sp_id ) {
				if( ops[i].dstEnd <= (b1+(2*ERR_TH)) ) {
					if( (min_x1 == -1) || ( min_x1 > abs(b1 - ops[i].dstEnd) ) ) { 
						min_x1 = abs(b1 - ops[i].dstEnd);
					}	
				}

				if( ops[i].dstStart >= (e1-(2*ERR_TH)) ) {
					if( (min_x2 == -1) || ( min_x2 > abs(e1 - ops[i].dstStart) ) ) { 
						min_x2 = abs(e1 - ops[i].dstStart);
					}	
				}

				if( ops[i].dstEnd <= (b2+(2*ERR_TH)) ) {
					if( (min_y1 == -1) || ( min_y1 > abs(b2 - ops[i].dstEnd) ) ) { 
						min_y1 = abs(b2 - ops[i].dstEnd);
					}	
				}

				if( ops[i].dstStart >= (e2-(2*ERR_TH)) ) {
					if( (min_y2 == -1) || ( min_y2 > abs(e2 - ops[i].dstStart) ) ) { 
						min_y2 = abs(e2 - ops[i].dstStart);
					}	
				}
				count++;
			}
		}
	}

	if( count == 0 ) {
		res = TIE;
	}
	else {
		if( (min_x1 == -1) && (min_x2 == -1) ) {
			min_x = -1;
		}
		else if( min_x1 == -1 ) {
			min_x = min_x2;
		}
		else if( min_x2 == -1 ) {
			min_x = min_x1;
		}
		else {
			if( min_x1 < min_x2 ) min_x = min_x1;
			else min_x = min_x2;
		}

		if( (min_y1 == -1) && (min_y2 == -1) ) {
			min_y = -1;
		}
		else if( min_y1 == -1 ) {
			min_y = min_y2;
		}
		else if( min_y2 == -1 ) {
			min_y = min_y1;
		}
		else {
			if( min_y1 < min_y2 ) min_y = min_y1;
			else min_y = min_y2;
		}

		if( (min_x == -1) && (min_y == -1) ) {
			res = TIE;
		}
		else if( abs(min_x - min_y) <= (10*ERR_TH) ) {
			if( min_x <= min_y ) res = TIE;
			else res = TIE;
		}
		else if( min_x < min_y ) {
			res = RIGHT_SIDE;
		}
		else {
			res = LEFT_SIDE;
		}
	}

	return(res);
}

void update_algns_sign(struct DotList *alg_self, int id, int num_alg_self, struct DotList *init_algns)
{
	int i = 0;
	int index = 0;

	for( i = 0; i < num_alg_self; i++ ) {
		if( i != id ) {	
			index = alg_self[i].index;
			if( init_algns[index].sign == DELETED ) alg_self[i].sign = DELETED;
		}
	}
}

struct slist replace_by_better(struct DotList *alg_self, struct slist res_list, struct slist *op_info, int num_op, int pid_th, int num_overlap)
{
	struct I dup_reg;
	int opt_id = 0;
	int i = 0, cur_id = -1;
	struct slist res;

	res = assign_slist(res_list);
	opt_id = res_list.id;

	if( opt_id == -1 ) {
		fatalf("alignment id for duplication has not been determined %d\n", opt_id);
	}
	else if( res_list.is_x == true ) // x is a dup region
	{
		dup_reg = assign_I(alg_self[opt_id].x.lower, alg_self[opt_id].x.upper);
	}
	else {
		dup_reg = assign_I(alg_self[opt_id].y.lower, alg_self[opt_id].y.upper);
	}

	for( i = 0; i < num_op; i++ ) {
		cur_id = op_info[i].id;
		if((strict_almost_equal(alg_self[cur_id].x, dup_reg) == true) || (strict_almost_equal(alg_self[cur_id].y, dup_reg) == true) ) {}
		else if((strict_almost_equal(dup_reg, alg_self[cur_id].x) == true) || (strict_almost_equal(dup_reg, alg_self[cur_id].y) == true) ) {}
		else if( (f_loose_subset(dup_reg, alg_self[cur_id].x, STRICT) == true) && (op_info[i].val != -1) && (op_info[i].val <= num_overlap) && (alg_self[cur_id].identity >= pid_th)) {
			res = assign_slist(op_info[i]);
			res.is_x = true;
			dup_reg = assign_I(alg_self[cur_id].x.lower, alg_self[cur_id].x.upper);
		}
		else if( (f_loose_subset(dup_reg, alg_self[cur_id].y, STRICT) == true) && (op_info[i].val_red != -1) && (op_info[i].val_red <= num_overlap) && (alg_self[cur_id].identity >= pid_th)) {
			res = assign_slist(op_info[i]);
			res.is_x = false;
			dup_reg = assign_I(alg_self[cur_id].y.lower, alg_self[cur_id].y.upper);
		}
	}

	return(res);
}

int check_prev_tandem_dup(struct DotList *init_algns, int num_init_algns, int id, struct ops_list *ops, int num_ops)
{
	int res = TIE;
	int i = 0;
	struct I src, dst;
	int cmp_id = -1;
	int sp_id = PAIR;
	
	src = assign_I(0, 1);
	dst = assign_I(0, 1);

	sp_id = init_algns[id].sp_id;
	i = 0; 
	while( (i < num_ops) && (res == TIE) ) {
		cmp_id = ops[i].id;
		if( (cmp_id < 0) || (cmp_id >= num_init_algns) ) {}
		else {
			if( sp_id == init_algns[cmp_id].sp_id ) {
				src = assign_I(ops[i].src_b, ops[i].src_e);	
				dst = assign_I(ops[i].dst_b, ops[i].dst_e);	
				if( overlap(src, dst) == true ) { // tandem dup
					if( (overlap( init_algns[id].x, init_algns[cmp_id].x ) && overlap( init_algns[id].x, init_algns[cmp_id].x ) ) || ( overlap( init_algns[id].x, init_algns[cmp_id].y ) && overlap( init_algns[id].y, init_algns[cmp_id].x ) ) )
					{
						if( (ops[i].srcStart < ops[i].dstStart) && (ops[i].srcEnd < ops[i].dstEnd) ) 
						{
							res = RIGHT_SIDE;
						} 
						else {
							res = LEFT_SIDE;
						}
					} 
				}
			}
		}
		i++;
	}

	return(res);
}
