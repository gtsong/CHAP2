#include "main.h"
#include "id_inpar.h"
#include "pred_regions.h"
#include "util_gen.h"
#include "regions.h"
#include "kd_tree.h"
#include "find_dup_copy.h"
#include "ins_dup_copy.h"
#include "find_own_dup_copy.h"
#include "extend_slist.h"
#include "handle_tandem_dup.h"
#include "cmp_two_reg.h"
#include "cmp_pid.h"
#include "util.h"
#include "find_merging.h"

extern int debug_mode;

struct slist iden_inpar(int *num_alg_self, struct DotList *alg_self, int size, int threshold, int *num_suspend_pairs, struct ID_List *suspend, int mode, struct DotList *init_algns, FILE *fp)
{
	struct slist res_list;
	int num_del;
	int num_ins_regs;

	struct kdnode *tree;
	struct perm_pt *p_pts;
	struct kdnode *m_tree;
	struct perm_pt *m_p_pts;

	p_pts = (struct perm_pt *) ckalloc((*num_alg_self) * sizeof(struct perm_pt));
	m_p_pts = (struct perm_pt *) ckalloc((*num_alg_self) * sizeof(struct perm_pt));

	res_list.id = NO_EXIST;
	res_list.is_x = true;
	res_list.sp_state = UNSUSPENDED;
  res_list.val = 0;
  res_list.val_red = 0;
  res_list.sp_state = 0;
  res_list.add_sp_state = 0;

	assign_perm(p_pts, (*num_alg_self), alg_self, LEFT);

	tree = build_kd(p_pts, 0, (*num_alg_self)-1);

	handle_tandem_dup(alg_self, num_alg_self, init_algns);

	find_own_dup_copy(alg_self, num_alg_self, tree, p_pts, size); // for intra-posed duplication

	assign_perm(m_p_pts, (*num_alg_self), alg_self, LEFT);
	m_tree = build_kd(m_p_pts, 0, (*num_alg_self)-1);

	num_ins_regs = find_dup_copy(alg_self, num_alg_self, m_tree, m_p_pts, size, suspend, fp, init_algns);
	*num_suspend_pairs = num_ins_regs;
	num_del = 0;

	res_list = find_opt_inpar(*num_alg_self, alg_self, size, threshold, num_suspend_pairs, suspend, mode, m_tree, m_p_pts, init_algns, fp); 

	if( (res_list.id != NO_EXIST) &&( res_list.id != SP_EVENT) )
	{
		if( num_ins_regs > 0 )
		{
			extend_slist(alg_self, num_alg_self, m_tree, m_p_pts, size, suspend, num_ins_regs, res_list.id, res_list.is_x, fp, init_algns);
		}
	}

	free(m_p_pts);
	free_kd(m_tree);
	free(p_pts);
	free_kd(tree);

	return(res_list);
}

struct slist find_opt_inpar(int num_alg_self, struct DotList *alg_self, int size, int threshold, int *num_suspend_pairs, struct ID_List *suspend, int mode, struct kdnode *tree, struct perm_pt *p_pts, struct DotList *init_algns, FILE *fp)
{
	struct slist *alg_id; // the list of candidate inparalogous regions
	struct slist *ins_alg_id; // the list of candidate inparalogous regions
	int num_ins_alg = 0;
	int num_id = 0; // the number of candidate inparalogous regions
	int i, j;
	struct slist res_list;
	int pred_op;
	int sp_state; // suspend state
	int num_ins_regs;
	int cmp_res;
	struct DotList *s_alg_self;
	struct DotList *t_alg_self;
	int num_list;
	int x_val = 0, y_val = 0;
	int *num_x;
	int *num_y;
	int *add_info;
	int f_id;
	int rp_id;
	int final_id, alt_final_id;
	bool is_loose;

	res_list.id = NO_EXIST;
	res_list.is_x = true;
	res_list.sp_state = UNSUSPENDED;
  res_list.val = 0;
  res_list.val_red = 0;
  res_list.sp_state = 0;
  res_list.add_sp_state = 0;

	if( threshold > LOOSE_RUN ) is_loose = true;
	else is_loose = false;

	num_ins_regs = *num_suspend_pairs;
	alg_id = (struct slist *) ckalloc(sizeof(struct slist) * (num_alg_self));
	ins_alg_id = (struct slist *) ckalloc(sizeof(struct slist) * (num_alg_self));
	num_x = (int *) ckalloc(sizeof(int));
	num_y = (int *) ckalloc(sizeof(int));
	add_info = (int *) ckalloc(sizeof(int));

	for( i = 0; i < num_alg_self; i++ ) {
		alg_id[i].id = 0;
		alg_id[i].val = 0;
		alg_id[i].val_red = 0;
		alg_id[i].is_x = true;
		ins_alg_id[i].id = 0;
		ins_alg_id[i].val = 0;
		ins_alg_id[i].val_red = 0;
		ins_alg_id[i].is_x = true;
	}

	s_alg_self = (struct DotList *) ckalloc(sizeof(struct DotList) * (num_alg_self));
	t_alg_self = (struct DotList *) ckalloc(sizeof(struct DotList) * (2 * num_alg_self));
	add_symmetric_points(s_alg_self, alg_self, num_alg_self);
	make_into_one(num_alg_self, alg_self, s_alg_self, t_alg_self);
	num_list = num_alg_self;
	num_list = 2*num_list;

	for( i = 0; i < num_alg_self; i++ )
	{
		if( alg_self[i].pair_self == PAIR ) {
			if(debug_mode == TRUE) printf("PAIR: %d-%d, %d-%d: %d\n", alg_self[i].x.lower, alg_self[i].x.upper, alg_self[i].y.lower, alg_self[i].y.upper, alg_self[i].identity);
		}
		else if( alg_self[i].pair_self == SELF )
		{
			sp_state = check_ins_dup_copy(i, suspend, num_ins_regs, alg_self, mode, add_info);
			if( sp_state == INSERTED_IN_BOTH )
			{
				pred_op = is_left_to_right_count(num_x, num_y, i, num_alg_self, alg_self, threshold, suspend, num_ins_regs, TANDEM_CHECK, tree, p_pts, size, fp, init_algns);
				if( (pred_op == OVERLAP_IN_X) || (pred_op == NO_OVERLAP) ) 
				{
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
			else
			{
				if( is_tandem(alg_self[i]) ) pred_op = is_left_to_right_count(num_x, num_y, i, num_alg_self, alg_self, threshold, suspend, num_ins_regs, TANDEM_CHECK, tree, p_pts, size, fp, init_algns);
				else pred_op = is_left_to_right_count(num_x, num_y, i, num_alg_self, alg_self, threshold, suspend, num_ins_regs, GENERAL_CHECK, tree, p_pts, size, fp, init_algns); // a flag is temporary
			}

			if(debug_mode == TRUE) printf("STATE: %d-%d, %d-%d: %d, %d, %d\n", alg_self[i].x.lower, alg_self[i].x.upper, alg_self[i].y.lower, alg_self[i].y.upper, sp_state, pred_op, *add_info);

			if(debug_mode == TRUE) printf("       %d-%d, %d-%d: %d\n", init_algns[alg_self[i].index].x.lower, init_algns[alg_self[i].index].x.upper, init_algns[alg_self[i].index].y.lower, init_algns[alg_self[i].index].y.upper, init_algns[alg_self[i].index].identity);

			if( (pred_op == OVERLAP_IN_X) && ( (sp_state == SP_OVERLAP_IN_X) || (sp_state == UNSUSPENDED) || (sp_state == INSERTED_IN_Y) || (sp_state == COVER_INS_SP_IN_Y) || (sp_state == INS_IN_Y_TAN) || (sp_state == INSERTED_IN_BOTH) || (sp_state == COVER_INS_SP_BOTH)) )
			{
				y_val = 0;
				for( j = 0; j < num_list; j++ )
				{
					y_val = y_val + increase_count(alg_self, t_alg_self, i, j, false);
				}

				alg_id[num_id].id = i;
				alg_id[num_id].is_x = false;
				alg_id[num_id].sp_state = sp_state; // the flag to indicate the inserted copy into other alignments
				alg_id[num_id].val_red = y_val;
				alg_id[num_id].add_sp_state = *add_info;
				num_id++;
			}
			else if( (pred_op == OVERLAP_IN_Y) && ( (sp_state == SP_OVERLAP_IN_Y) || (sp_state == UNSUSPENDED) || (sp_state == INSERTED_IN_X) || (sp_state == COVER_INS_SP_IN_X) || (sp_state == INS_IN_X_TAN) || (sp_state == INSERTED_IN_BOTH) || (sp_state == COVER_INS_SP_BOTH)) )
			{
				x_val = 0;
				for( j = 0; j < num_list; j++ )
				{
					x_val = x_val + increase_count(alg_self, t_alg_self, i, j, true);
				}

				alg_id[num_id].id = i;
				alg_id[num_id].is_x = true;
				alg_id[num_id].sp_state = sp_state; // the flag to indicate the inserted copy into other alignments
				alg_id[num_id].val_red = x_val;
				alg_id[num_id].add_sp_state = *add_info;
				num_id++;
			}
			else if( sp_state == SUSPENDED )
			{
			}
			else if( pred_op == NO_OVERLAP )
			{
				if( (sp_state == UNSUSPENDED) || (sp_state == INSERTED_IN_BOTH) )
				{
					cmp_res = det_dup_reg_in_self(num_alg_self, alg_self, alg_self[i].x, alg_self[i].y, alg_self[i].sign);
				}
				else if( (sp_state == INSERTED_IN_X) || (sp_state == SP_OVERLAP_IN_Y) || (sp_state == COVER_INS_SP_IN_X) || (sp_state == INS_IN_X_TAN) )
				{
					cmp_res = RIGHT_SIDE;
				}
				else
				{
					cmp_res = LEFT_SIDE;
				}

				if( cmp_res == LEFT_SIDE ) // y region is removed
				{
					if( sp_state == INSERTED_IN_BOTH ) sp_state = INS_IN_Y_TAN;

					y_val = 0;
					for( j = 0; j < num_list; j++ )
					{
						y_val = y_val + increase_count(alg_self, t_alg_self, i, j, false);
					}

					alg_id[num_id].id = i;
					alg_id[num_id].is_x = false;
					alg_id[num_id].sp_state = sp_state; // the flag to indicate the inserted copy into other alignments
					alg_id[num_id].val_red = y_val;
					alg_id[num_id].add_sp_state = *add_info;
					num_id++;
				}
				else if( cmp_res == RIGHT_SIDE )
				{
					if( sp_state == INSERTED_IN_BOTH ) sp_state = INS_IN_X_TAN;

					x_val = 0;
					for( j = 0; j < num_list; j++ )
					{
						x_val = x_val + increase_count(alg_self, t_alg_self, i, j, true);
					}
					alg_id[num_id].id = i;
					alg_id[num_id].is_x = true;
					alg_id[num_id].sp_state = sp_state; // the flag to indicate the inserted copy into other alignments
					alg_id[num_id].val_red = x_val;
					alg_id[num_id].add_sp_state = *add_info;
					num_id++;
				}
				else
				{
					x_val = 0;
					for( j = 0; j < num_list; j++ )
					{
						x_val = x_val + increase_count(alg_self, t_alg_self, i, j, true);
					}

					y_val = 0;
					for( j = 0; j < num_list; j++ )
					{
						y_val = y_val + increase_count(alg_self, t_alg_self, i, j, false);
					}

					if( x_val > y_val )
					{
						alg_id[num_id].id = i;
						alg_id[num_id].is_x = true;
						alg_id[num_id].sp_state = sp_state; // the flag to indicate the inserted copy into other alignments
						alg_id[num_id].val_red = x_val;
						alg_id[num_id].add_sp_state = *add_info;
						num_id++;
					}
					else
					{
						alg_id[num_id].id = i;
						alg_id[num_id].is_x = false;
						alg_id[num_id].sp_state = sp_state; // the flag to indicate the inserted copy into other alignments
						alg_id[num_id].val_red = y_val;
						alg_id[num_id].add_sp_state = *add_info;
						num_id++;
					}	
				}
			}
		}
	}
	
	if( num_id == 0 ) 
	{
		res_list.id = NO_EXIST;
	}
	else
	{
		if( mode == BEFORE_SP ) {
		}

		j = 0;
		for( i = 0; i < num_id; i++ ) {
			if( alg_id[i].id != -1 ) {
				alg_id[j].id = alg_id[i].id;
				alg_id[j].is_x = alg_id[i].is_x;
				alg_id[j].sp_state = alg_id[i].sp_state;
				alg_id[j].val_red = alg_id[i].val_red;
				alg_id[j].add_sp_state = alg_id[i].add_sp_state;
				j++;
			}
		}
		num_id = j;

		num_ins_alg = 0;

		for( i = 0; i < num_id; i++ )
		{
			if( (alg_id[i].sp_state == INSERTED_IN_X) || (alg_id[i].sp_state == INSERTED_IN_Y) || (alg_id[i].sp_state == COVER_INS_SP_IN_X) || (alg_id[i].sp_state == COVER_INS_SP_IN_Y) || (alg_id[i].sp_state == INS_IN_X_TAN) || (alg_id[i].sp_state == INS_IN_Y_TAN))
			{
				if( (alg_id[i].add_sp_state == SUSPENDED) || ((alg_id[i].add_sp_state == SP_OVERLAP_IN_X) && ((alg_id[i].sp_state == INSERTED_IN_X) || (alg_id[i].sp_state == INS_IN_X_TAN))) || ((alg_id[i].add_sp_state == SP_OVERLAP_IN_Y) && ((alg_id[i].sp_state == INSERTED_IN_Y) || (alg_id[i].sp_state == INS_IN_Y_TAN))) ) {} 
				else if( alg_self[alg_id[i].id].pair_self != PAIR )
				{
					ins_alg_id[num_ins_alg].id = alg_id[i].id;
					ins_alg_id[num_ins_alg].is_x = alg_id[i].is_x;
					ins_alg_id[num_ins_alg].sp_state = alg_id[i].sp_state;
					ins_alg_id[num_ins_alg].val_red = alg_id[i].val_red;
					num_ins_alg++;
				}
			}
		}

		if( num_ins_alg == 0 )
		{
			for( i = 0; i < num_id; i++ )
			{
				if( (alg_id[i].sp_state == INSERTED_IN_X) || (alg_id[i].sp_state == INSERTED_IN_Y) || (alg_id[i].sp_state == COVER_INS_SP_IN_X) || (alg_id[i].sp_state == COVER_INS_SP_IN_Y) || (alg_id[i].sp_state == INS_IN_X_TAN) || (alg_id[i].sp_state == INS_IN_Y_TAN))
				{
					if( alg_self[alg_id[i].id].pair_self != PAIR )
					{
						ins_alg_id[num_ins_alg].id = alg_id[i].id;
						ins_alg_id[num_ins_alg].is_x = alg_id[i].is_x;
						ins_alg_id[num_ins_alg].sp_state = alg_id[i].sp_state;
						ins_alg_id[num_ins_alg].val_red = alg_id[i].val_red;
						num_ins_alg++;
					}
				}
			}
		}

		if( num_ins_alg > 0 )
		{
			sort_by_pid(ins_alg_id, alg_self, num_ins_alg);

			final_id = 0;
			while( (final_id < num_ins_alg) && ((alg_self[ins_alg_id[final_id].id].lock == BOTH_LOCK) || ((ins_alg_id[final_id].is_x == true) && (alg_self[ins_alg_id[final_id].id].lock == LEFT_LOCK)) || ((ins_alg_id[final_id].is_x == false) && (alg_self[ins_alg_id[final_id].id].lock == RIGHT_LOCK)) )) final_id++;

			if( final_id >= num_ins_alg ) final_id = 0;
			res_list.id = ins_alg_id[final_id].id;
			res_list.is_x = ins_alg_id[final_id].is_x;
			res_list.val = alg_self[ins_alg_id[final_id].id].identity;
			res_list.val_red = ins_alg_id[final_id].val_red;
			res_list.sp_state = ins_alg_id[final_id].sp_state;
			if( (res_list.sp_state == COVER_INS_SP_IN_X) || (res_list.sp_state == COVER_INS_SP_IN_Y) )
			{
				res_list.sp_state = UNSUSPENDED;
			}
			else 
			{
				res_list.sp_state = ins_alg_id[0].sp_state;
			}

			if( is_loose ) {
				sort_by_pid(alg_id, alg_self, num_id);
				alt_final_id = 0;
				while( (alt_final_id < num_id) && ((alg_self[alg_id[alt_final_id].id].lock == BOTH_LOCK) || ((alg_id[alt_final_id].is_x == true) && (alg_self[alg_id[alt_final_id].id].lock == LEFT_LOCK)) || ((alg_id[alt_final_id].is_x == false) && (alg_self[alg_id[alt_final_id].id].lock == RIGHT_LOCK)) )) alt_final_id++;

				if( alt_final_id >= num_id ) alt_final_id = 0;

				if( ((alg_self[ins_alg_id[final_id].id].identity - alg_self[alg_id[alt_final_id].id].identity) < 0) || ( ((alg_self[ins_alg_id[final_id].id].identity - alg_self[alg_id[alt_final_id].id].identity) == 0) && ( alg_id[alt_final_id].val_red > ins_alg_id[final_id].val_red ) ) ) {
					res_list.id = alg_id[alt_final_id].id;
					res_list.is_x = alg_id[alt_final_id].is_x;
					res_list.val = alg_self[alg_id[alt_final_id].id].identity;
					res_list.val_red = alg_id[alt_final_id].val_red;
					res_list.sp_state = alg_id[alt_final_id].sp_state;
					if( (res_list.sp_state == COVER_INS_SP_IN_X) || (res_list.sp_state == COVER_INS_SP_IN_Y) ) res_list.sp_state = UNSUSPENDED;
					else res_list.sp_state = alg_id[0].sp_state;
				}
			}

		}
		else if( threshold == STRICT )
		{
			if( mode == BEFORE_SP )
			{
				res_list.id = -1;
			}
			else
			{
				sort_by_pid(alg_id, alg_self, num_id);

				res_list.id = alg_id[0].id;
				res_list.is_x = alg_id[0].is_x;
				res_list.val = alg_self[alg_id[0].id].identity;
				res_list.sp_state = alg_id[0].sp_state;
				res_list.val_red = alg_id[0].val_red;
			}
		}
		else
		{
			sort_by_pid(alg_id, alg_self, num_id);

			rp_id = -1;
			if( alg_self[alg_id[0].id].pair_self == PAIR )
			{
				f_id = 1;
				while( (f_id < num_id ) && (alg_self[alg_id[f_id].id].pair_self == PAIR) ) f_id++;

				if( f_id == num_id )
				{
					rp_id = SP_EVENT;
				}
			}
			else
			{
				f_id = 0;
			}

			if( rp_id == -1 )
			{
				res_list.id = alg_id[f_id].id;
				res_list.is_x = alg_id[f_id].is_x;
				res_list.val = alg_self[alg_id[f_id].id].identity;
				res_list.sp_state = alg_id[f_id].sp_state;
				res_list.val_red = alg_id[f_id].val_red;
			}
			else if( rp_id == SP_EVENT )
			{
				res_list.id = SP_EVENT;
			}
		}
	}

	if(debug_mode == TRUE) {
		if( (res_list.id != NO_EXIST) && (res_list.id != SP_EVENT) ) printf("duplication in %d-%d %d-%d\n", alg_self[res_list.id].x.lower, alg_self[res_list.id].x.upper, alg_self[res_list.id].y.lower, alg_self[res_list.id].y.upper);
		else if( res_list.id == SP_EVENT ) printf("speciation event\n");
		else printf("no dup alignments exist\n");
	}

	free(alg_id);
	free(ins_alg_id);
	free(num_y);
	free(num_x);
	free(add_info);
	free(s_alg_self);
	free(t_alg_self);

	return(res_list);
}
