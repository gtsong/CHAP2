#include "main.h"
#include "regions.h"
#include "ins_dup_copy.h"
#include "adjust_algn.h"
#include "find_merging.h"
#include "update_init_algns.h"
#include "util.h"
#include "util_i.h"

extern int debug_mode;

int check_ins_dup_copy(int id, struct ID_List *dlist, int num_dup_copy, struct DotList *dots, int mode, int *add_info)
{
	int i = 0, j = 0;
	int res = UNSUSPENDED;
	int ins_res = SUSPENDED;
	int sp_res = UNSUSPENDED;
	int cover_res = UNSUSPENDED;
	int temp_res = UNSUSPENDED;
	struct I ins_reg, left_flank, right_flank, split_left, split_right;
	struct ID_List c_list;
	int temp_index = -1;

	while( (i < num_dup_copy) && (ins_res == SUSPENDED))
	{
		c_list = dlist[i];

		if(c_list.is_x == true )
		{
			ins_reg = assign_I(dots[c_list.m_id].x.lower, dots[c_list.m_id].x.upper);
		}
		else
		{
			ins_reg = assign_I(dots[c_list.m_id].y.lower, dots[c_list.m_id].y.upper);
		}

		if( id == dlist[i].m_id )
		{
			if( dlist[i].is_x == true ) ins_res = INSERTED_IN_X;
			else ins_res = INSERTED_IN_Y;

			temp_index = i;
		}
		i++;
	}

	for( i = 0; i < num_dup_copy; i++ )
	{
		temp_res = UNSUSPENDED;
		c_list = dlist[i];

		if(c_list.is_x == true )
		{
			ins_reg = assign_I(dots[c_list.m_id].x.lower, dots[c_list.m_id].x.upper);
		}
		else
		{
			ins_reg = assign_I(dots[c_list.m_id].y.lower, dots[c_list.m_id].y.upper);
		}

		if( c_list.f_is_x == true )
		{
			left_flank = assign_I(dots[c_list.left_id].x.lower, dots[c_list.left_id].x.upper);
			right_flank = assign_I(dots[c_list.right_id].x.lower, dots[c_list.right_id].x.upper);
			split_left = assign_I(dots[c_list.left_id].y.lower, dots[c_list.left_id].y.upper);
			split_right = assign_I(dots[c_list.right_id].y.lower, dots[c_list.right_id].y.upper);
		}
		else
		{
			left_flank = assign_I(dots[c_list.left_id].y.lower, dots[c_list.left_id].y.upper);
			right_flank = assign_I(dots[c_list.right_id].y.lower, dots[c_list.right_id].y.upper);
			split_left = assign_I(dots[c_list.left_id].x.lower, dots[c_list.left_id].x.upper);
			split_right = assign_I(dots[c_list.right_id].x.lower, dots[c_list.right_id].x.upper);
		}

		if( id == dlist[i].m_id ) {}
		else if( dots[dlist[i].m_id].pair_self == PAIR ) {}
		else if( check_later_dup(dlist[i], id, dots) == true )
		{
		}	
		else if( (ins_res != INSERTED_IN_X) && (ins_res != INSERTED_IN_Y) && ( strict_almost_equal(dots[id].x, ins_reg) == false ) && ( strict_almost_equal(dots[id].y, ins_reg) == false )  && ((temp_res = check_cover_ins_sp(dlist[i], id, dots)) != UNSUSPENDED ) )
		{
			if( cover_res == UNSUSPENDED )
			{
				cover_res = temp_res;
			}
			else if( cover_res == COVER_INS_SP_IN_X )
			{
				if( temp_res == COVER_INS_SP_IN_Y )
				{
					cover_res = COVER_INS_SP_BOTH;
				}
			}
			else if( cover_res == COVER_INS_SP_IN_Y )
			{
				if( temp_res == COVER_INS_SP_IN_X )
				{
					cover_res = COVER_INS_SP_BOTH;
				}
			}
		}
		else
		{
			if( (id == dlist[i].left_id) || (id == dlist[i].right_id) )
			{
				sp_res = SUSPENDED;
			}
			else if( (( strict_almost_equal(ins_reg, dots[id].x) == false) && (loose_overlap(ins_reg, dots[id].x) == true)) || (loose_overlap(left_flank, dots[id].x) == true) || (loose_overlap(right_flank, dots[id].x) == true) )
			{
				if( sp_res == SP_OVERLAP_IN_Y) sp_res = SUSPENDED;
				else if( sp_res == UNSUSPENDED ) sp_res = SP_OVERLAP_IN_X;
			}
			else if( (loose_overlap(split_left, dots[id].x) == true) || (loose_overlap(split_right, dots[id].x) == true) )
			{
				if( sp_res == SP_OVERLAP_IN_Y) sp_res = SUSPENDED;
				else if( sp_res == UNSUSPENDED ) sp_res = SP_OVERLAP_IN_X;
			}

			if( (id == dlist[i].left_id) || (id == dlist[i].right_id) )
			{
				sp_res = SUSPENDED;
			}
			else if( (( strict_almost_equal(ins_reg, dots[id].y) == false) && (loose_overlap(ins_reg, dots[id].y) == true)) || (loose_overlap(left_flank, dots[id].y) == true) || (loose_overlap(right_flank, dots[id].y) == true) )
			{
				if( sp_res == SP_OVERLAP_IN_X) sp_res = SUSPENDED;
				else if( sp_res == UNSUSPENDED ) sp_res = SP_OVERLAP_IN_Y;
			}
			else if( (loose_overlap(split_left, dots[id].y) == true) || (loose_overlap(split_right, dots[id].y) == true) )
			{
				if( sp_res == SP_OVERLAP_IN_X) sp_res = SUSPENDED;
				else if( sp_res == UNSUSPENDED ) sp_res = SP_OVERLAP_IN_Y;
			}
		}
	}

// if res is OVERLAP_IN_X, x region should not be removed
// else if res is OVERLAP_IN_Y, y region should not be removed
// else if res is SUSPENDED, the alignment can not be removed
	if( (ins_res == INSERTED_IN_X) || (ins_res == INSERTED_IN_Y) )
	{
		if( is_tandem(dots[id]) == true )
		{
			for( j = 0; j < num_dup_copy; j++ )
			{
				if( j == temp_index ) {}
				else if( dlist[temp_index].m_id == dlist[j].m_id )
				{
					if( dlist[j].is_t_ins == true )
					{
						temp_index = j;
					}
				}
			}
		}

		if( (dlist[temp_index].is_t_ins == true) && (is_tandem(dots[id]) == true) )
		{
			ins_res = INSERTED_IN_BOTH;
		}
	}

	if( mode == BEFORE_SP )
	{
		if( ins_res != SUSPENDED )
		{
			res = ins_res;
		}
		else
		{
			if( sp_res == UNSUSPENDED )
			{
				res = cover_res;
			}
			else
			{
				res = sp_res;
			}
		}
	}
	else if( mode == AFTER_SP )
	{
		if( sp_res == UNSUSPENDED )
		{
			if( ins_res != SUSPENDED )
			{
				res = ins_res;
			}
			else res = cover_res;
		}
		else
		{
			if( ins_res == SUSPENDED )
			{
				if( cover_res != UNSUSPENDED ) res = cover_res;
				else res = sp_res;
			}
			else
			{
				if( sp_res == SUSPENDED )
				{
					res = sp_res;
				}
				else if( sp_res == SP_OVERLAP_IN_X )
				{
					if( (ins_res == INSERTED_IN_Y) || (ins_res == INS_IN_Y_TAN) || (ins_res == INSERTED_IN_BOTH) ) res = ins_res;
					else res = sp_res;
				}
				else if( sp_res == SP_OVERLAP_IN_Y )
				{
					if( (ins_res == INSERTED_IN_X) || (ins_res == INS_IN_X_TAN) || (ins_res == INSERTED_IN_BOTH) ) res = ins_res;
					else res = sp_res;
				}
				else res = sp_res;
			}
		}
	}

	*add_info = sp_res;
	return(res);
}

void merge_ins_dup(int dup_copy, int id, struct DotList *dots, struct ID_List *dlist, int num_dup_copy, int *ch_index, int sp1, int sp2)
{
	int code_matrix[NUM_SP][NUM_SP];
	int i = 0;
	struct I temp;
	struct I cur;
	int cur_len = 0;
	int l_id = 0, r_id = 0;
	int c_fid = 0, d_fid = 0;
	int l_val = 0, r_val = 0;

	temp = assign_I(0,1);
	cur = assign_I(0,1);

	if( (dup_copy == INSERTED_IN_X) || (dup_copy == INS_IN_X_TAN) )
	{
		cur = assign_I(dots[id].x.lower, dots[id].x.upper);
	}
	else if( (dup_copy == INSERTED_IN_Y) || (dup_copy == INS_IN_Y_TAN) )
	{
		cur = assign_I(dots[id].y.lower, dots[id].y.upper);
	}
	cur_len = width(cur);

	for( i = 0; i < num_dup_copy; i++ )
	{
		if( dlist[i].is_x == true )
		{
			temp = assign_I(dots[dlist[i].m_id].x.lower, dots[dlist[i].m_id].x.upper);
		}
		else
		{
			temp = assign_I(dots[dlist[i].m_id].y.lower, dots[dlist[i].m_id].y.upper);
		}

		if( strict_almost_equal(cur, temp) == true )
		{
			l_id = dlist[i].left_id;
			r_id = dlist[i].right_id;
			c_fid = (dots[l_id].index)/BASE_NUM;
			d_fid = (dots[r_id].index)/BASE_NUM;
			if( dots[l_id].sp_id == code_matrix[sp1][sp2] ) ch_index[c_fid] = d_fid;
			
			if( dots[l_id].x.lower <= dots[r_id].x.lower ) l_val = dots[l_id].x.lower;
			else l_val = dots[r_id].x.lower;

			if( dots[l_id].x.upper >= dots[r_id].x.upper ) r_val = dots[l_id].x.upper;
			else r_val = dots[r_id].x.upper;
			dots[l_id].x = assign_I(l_val, r_val);

			if( dots[l_id].y.lower <= dots[r_id].y.lower ) l_val = dots[l_id].y.lower;
			else l_val = dots[r_id].y.lower;

			if( dots[l_id].y.upper >= dots[r_id].y.upper ) r_val = dots[l_id].y.upper;
			else r_val = dots[r_id].y.upper;
			dots[l_id].y = assign_I(l_val, r_val);
			if( debug_mode == TRUE) {
				printf("merged %d-%d, %d-%d into %d-%d, %d-%d\n", dots[r_id].x.lower, dots[r_id].x.upper, dots[r_id].y.lower, dots[r_id].y.upper, dots[l_id].x.lower, dots[l_id].x.upper, dots[l_id].y.lower, dots[l_id].y.upper);
			}
			dots[r_id].sign = DELETED;
		}
	}
}

void rollback_ins_dup(int dup_copy_flag, int id, struct DotList *dots, struct ID_List *dlist, int num_dup_copy, struct DotList *init_dots, FILE *fp)
{
	int i;
	struct I temp;
	struct I cur;
	struct I dup_reg;
	int cur_len = 0;
	bool is_x = true;
	struct I left, right;

	int old_id = 0;
	int cid = 0, org_id = 0;
	int b = 0, e = 1;

	temp = assign_I(0,1);
	dup_reg = assign_I(0,1);
	cur = assign_I(0,1);
	left = assign_I(0,1);
	right = assign_I(0,1);

	if( (dup_copy_flag == INSERTED_IN_X) || (dup_copy_flag == INS_IN_X_TAN) )
	{
		cur = assign_I(dots[id].x.lower, dots[id].x.upper);
		is_x = true;
	}
	else if( (dup_copy_flag == INSERTED_IN_Y) || (dup_copy_flag == INS_IN_Y_TAN) )
	{
		cur = assign_I(dots[id].y.lower, dots[id].y.upper);
		is_x = false;
	}
	cur_len = width(cur);

	org_id = dots[id].index;
  old_id = init_dots[org_id].index;
  cid = init_dots[org_id].c_id;
  while( cid != -1 ) {
    old_id = cid;
    cid = init_dots[cid].c_id;
  }

  if( (dots[id].sign != DELETED) && (init_dots[org_id].sign == DELETED) ) {
    init_dots[org_id].sign = dots[id].sign; // after rolling back a tandem duplication, a converted alignment in the step of handling tandem duplications could be deleted during the rollback process  
  }

  if( is_x == true ) {
    b = init_dots[org_id].x.lower + init_dots[org_id].xl_diff;
    e = init_dots[old_id].x.upper - init_dots[old_id].xr_diff;
  }
  else {
    if( init_dots[org_id].sign == 0 ) {
      b = init_dots[org_id].y.lower + init_dots[org_id].yl_diff;
      e = init_dots[old_id].y.upper - init_dots[old_id].yr_diff;
    }
    else if( init_dots[org_id].sign == 1 ) {
      b = init_dots[old_id].y.lower + init_dots[old_id].yl_diff;
      e = init_dots[org_id].y.upper - init_dots[org_id].yr_diff;
    }
    else {
      fatalf("update_init_algns: invalid sign in %d and %d\n", org_id, old_id);
    }
  }

	if( b > e ) {
    fatalf("update_init_algns: invalid interval [%d, %d]\n", b, e);
	}
	else {
		dup_reg = assign_I(b, e);
	}

	if( debug_mode == TRUE ) printf("Interval in the original dot plot: [%d, %d]\n", b, e);

	for( i = 0; i < num_dup_copy; i++ )
	{
		if( dlist[i].is_x == true )
		{
			temp = assign_I(dots[dlist[i].m_id].x.lower, dots[dlist[i].m_id].x.upper);
		}
		else
		{
			temp = assign_I(dots[dlist[i].m_id].y.lower, dots[dlist[i].m_id].y.upper);
		}

		if( dlist[i].f_is_x == true ) {
			left = assign_I(dots[dlist[i].left_id].x.lower, dots[dlist[i].left_id].x.upper);
			right = assign_I(dots[dlist[i].right_id].x.lower, dots[dlist[i].right_id].x.upper);
		}
		else {
			left = assign_I(dots[dlist[i].left_id].y.lower, dots[dlist[i].left_id].y.upper);
			right = assign_I(dots[dlist[i].right_id].y.lower, dots[dlist[i].right_id].y.upper);
		}

		if( strict_almost_equal(cur, temp) == true )
		{
			if( (dots[dlist[i].left_id].sign == DELETED ) || (dots[dlist[i].right_id].sign == DELETED ) ) {}
			else if( (dots[dlist[i].left_id].sign == 20 ) || (dots[dlist[i].right_id].sign == 20 ) ) {}
			else if( (dots[dlist[i].left_id].sign == 21 ) || (dots[dlist[i].right_id].sign == 21 ) ) {}
			else if( strict_subset(left, cur) || strict_subset(right, cur) ) {}
			else {
				update_init_algn(dup_reg, dots, id, is_x, dlist[i].left_id, init_dots, fp, INS);
				update_init_algn(dup_reg, dots, id, is_x, dlist[i].right_id, init_dots, fp, INS);
				mark_chain(dots, dlist[i].left_id, dlist[i].right_id, init_dots);
				merge_two_algns(dots, dlist[i].left_id, dlist[i].right_id, dlist[i].f_is_x, cur, dup_copy_flag);

				if( (i < (num_dup_copy-1)) && (dlist[i+1].m_id == dlist[i].m_id) && (dlist[i+1].left_id == dlist[i].left_id) && (dlist[i+1].right_id == dlist[i].right_id) ) i++; 
			}
		}
	}
}

void merge_two_algns(struct DotList *dots, int left_id, int right_id, bool is_x, struct I cur, int dup_copy)
{
	int len = 0;
	int id1 = 0, id2 = 0;
	int b = 0, e = 1;

	if( (dup_copy == INS_IN_X_TAN) || (dup_copy == INS_IN_Y_TAN) )
	{
		if( is_x == true )
		{
			if( proper_overlap( dots[left_id].x, cur ) == true )
			{
				adjust_alignment(dots, left_id, cur);
			}

			if( proper_overlap( dots[right_id].x, cur ) == true )
			{
				adjust_alignment(dots, right_id, cur);
			}
		}
		else
		{
			if( proper_overlap( dots[left_id].y, cur ) == true )
			{
				adjust_alignment(dots, left_id, cur);
			}

			if( proper_overlap( dots[right_id].y, cur ) == true )
			{
				adjust_alignment(dots, right_id, cur);
			}
		}
	}

	if( is_x == true ) {
		if( dots[left_id].x.lower <= dots[right_id].x.lower )
		{
			id1 = left_id;
			id2 = right_id;
		}
		else 
		{
			id1 = right_id;
			id2 = left_id;
		}
	}
	else
	{
		if( dots[left_id].y.lower <= dots[right_id].y.lower )
		{
			if( dots[left_id].init_sign == 0 ) {
				id1 = left_id;
				id2 = right_id;
			}
			else {
				id1 = right_id;
				id2 = left_id;
			}
		}
		else 
		{
			if( dots[left_id].init_sign == 0 ) {
				id1 = right_id;
				id2 = left_id;
			}
			else {
				id1 = left_id;
				id2 = right_id;
			}
		}
	}

	len = width(cur);
	if( (width(dots[id1].x) >= MIN_LEN_FOR_RECAL_PID) && (width(dots[id2].x) >= MIN_LEN_FOR_RECAL_PID) )
	{
		dots[id1].identity = ((dots[id1].identity * width(dots[id1].x)) + (dots[id2].identity * width(dots[id2].x))) / (width(dots[id1].x) + width(dots[id2].x));
	}
	else if( (width(dots[id1].x) <= MIN_LEN_FOR_RECAL_PID) && (width(dots[id2].x) <= MIN_LEN_FOR_RECAL_PID) )
	{
		dots[id1].identity = ((dots[id1].identity * width(dots[id1].x)) + (dots[id2].identity * width(dots[id2].x))) / (width(dots[id1].x) + width(dots[id2].x));
	}
	else if( width(dots[id2].x) >= MIN_LEN_FOR_RECAL_PID )
	{
		dots[id1].identity = dots[id2].identity;
	}

	if( dots[id1].init_sign == 0 )
	{
		if( is_x == true ) {
			b = dots[id1].x.lower;
			e = dots[id2].x.upper - len;
			if( b >= e ) {
				e = dots[id1].x.upper;
			}
		}
		else {
			b = dots[id1].x.lower;
			e = dots[id2].x.upper;
		}

		if( b >= e ) {
			fatalf("merging error: [%d-%d], [%d-%d]\n", dots[id1].x.lower, dots[id1].x.upper, dots[id2].x.lower, dots[id2].x.upper);
		}
		else {
			dots[id1].x = assign_I(b, e);
		}

		if( is_x == true ) {
			b = dots[id1].y.lower;
			e = dots[id2].y.upper;
		}
		else {
			b = dots[id1].y.lower;
			e = dots[id2].y.upper - len;
			if( b >= e ) {
				e = dots[id1].y.upper;
			}
		}

		if( b >= e ) {
			fatalf("merging error: [%d-%d], [%d-%d]\n", dots[id1].y.lower, dots[id1].y.upper, dots[id2].y.lower, dots[id2].y.upper);
		}
		else {
			dots[id1].y = assign_I(b, e);
		}

		dots[id1].sign = 20; // mark to be merged
		dots[id2].sign = 2; // the right alignment is deleted since it was merged into the left alignment
	}
	else if( dots[id1].init_sign == 1 )
	{
		if( is_x == true ) {
			b = dots[id1].x.lower;
			e = dots[id2].x.upper - len;
			if( b >= e ) {
				e = dots[id1].x.upper;
			}
		}
		else {
			b = dots[id1].x.lower;
			e = dots[id2].x.upper;
		}

		if( b >= e ) {
			fatalf("merging error: [%d-%d], [%d-%d]\n", dots[id1].x.lower, dots[id1].x.upper, dots[id2].x.lower, dots[id2].x.upper);
		}
		else {
			dots[id1].x = assign_I(b, e);
		}

		if( is_x == true ) {
			b = dots[id2].y.lower;
			e = dots[id1].y.upper;
		}
		else {
			b = dots[id2].y.lower;
			e = dots[id1].y.upper - len;
			if( b >= e ) {
				e = dots[id2].y.upper;
			}
		}

		if( b >= e ) {
			fatalf("merging error: [%d-%d], [%d-%d]\n", dots[id1].y.lower, dots[id1].y.upper, dots[id2].y.lower, dots[id2].y.upper);
		}
		else {
			dots[id1].y = assign_I(b, e);
		}
		dots[id1].sign = 21; // mark to be merged
		dots[id2].sign = 2; // the right alignment is deleted since it was merged into the left alignment
	}	
}

bool check_later_dup(struct ID_List c_list, int id, struct DotList *dots)
{
	bool res = false;
	struct I ins_reg, left_flank, right_flank;
	struct I split_left, split_right;

	if(c_list.is_x == true )
	{
		ins_reg = assign_I(dots[c_list.m_id].x.lower, dots[c_list.m_id].x.upper);
	}
	else
	{
		ins_reg = assign_I(dots[c_list.m_id].y.lower, dots[c_list.m_id].y.upper);
	}

	if( c_list.f_is_x == true )
	{
		left_flank = assign_I(dots[c_list.left_id].x.lower, dots[c_list.left_id].x.upper);
		right_flank = assign_I(dots[c_list.right_id].x.lower, dots[c_list.right_id].x.upper);
		split_left = assign_I(dots[c_list.left_id].y.lower, dots[c_list.left_id].y.upper);
		split_right = assign_I(dots[c_list.right_id].y.lower, dots[c_list.right_id].y.upper);
	}
	else
	{
		left_flank = assign_I(dots[c_list.left_id].y.lower, dots[c_list.left_id].y.upper);
		right_flank = assign_I(dots[c_list.right_id].y.lower, dots[c_list.right_id].y.upper);
		split_left = assign_I(dots[c_list.left_id].x.lower, dots[c_list.left_id].x.upper);
		split_right = assign_I(dots[c_list.right_id].x.lower, dots[c_list.right_id].x.upper);
	}

	if( (loose_subset(ins_reg, dots[id].x) == true) && (loose_subset(left_flank, dots[id].x) == true) && (loose_subset(right_flank, dots[id].x) == true) )
	{
		res = true;
	}

	if( (loose_subset(ins_reg, dots[id].y) == true) && (loose_subset(left_flank, dots[id].y) == true) && (loose_subset(right_flank, dots[id].y) == true) ) 
	{
		res = true;
	}
	
	if( (loose_subset(split_left, dots[id].x) == true) && (loose_subset(split_right, dots[id].x) == true) ) 
	{
		res = true;
	}

	if( (loose_subset(split_left, dots[id].y) == true) && (loose_subset(split_right, dots[id].y) == true) ) 
	{
		res = true;
	}

	return(res);
}

int check_cover_ins_sp(struct ID_List c_list, int id, struct DotList *dots)
{
	int temp_res = UNSUSPENDED;
	struct I ins_reg, left_flank, right_flank;
	struct I split_left, split_right;
	struct I merged_left, merged_right;
	int left_id, right_id;

	if(c_list.is_x == true )
	{
		ins_reg = assign_I(dots[c_list.m_id].x.lower, dots[c_list.m_id].x.upper);
	}
	else
	{
		ins_reg = assign_I(dots[c_list.m_id].y.lower, dots[c_list.m_id].y.upper);
	}

	if( dots[c_list.left_id].x.lower < dots[c_list.right_id].x.lower )
	{
		left_id = c_list.left_id;
		right_id = c_list.right_id;
	}
	else if( dots[c_list.left_id].x.lower > dots[c_list.right_id].x.lower )
	{
		right_id = c_list.left_id;
		left_id = c_list.right_id;
	}
	else
	{
		if( dots[c_list.left_id].x.upper <= dots[c_list.right_id].x.upper )
		{
			left_id = c_list.left_id;
			right_id = c_list.right_id;
		}
		else 
		{
			right_id = c_list.left_id;
			left_id = c_list.right_id;
		}
	}

	if( (id == c_list.left_id) || (id == c_list.right_id) )
	{
		temp_res = UNSUSPENDED;
		return(temp_res);	
	}
	else if( c_list.f_is_x == true )
	{
		left_flank = assign_I(dots[left_id].x.lower, dots[left_id].x.upper);
		right_flank = assign_I(dots[right_id].x.lower, dots[right_id].x.upper);
		
		if( dots[left_id].sign == 0 )
		{
			split_left = assign_I(dots[left_id].y.lower, dots[left_id].y.upper);
			split_right = assign_I(dots[right_id].y.lower, dots[right_id].y.upper);
		}
		else if( dots[left_id].sign == 1 )
		{
			split_right = assign_I(dots[left_id].y.lower, dots[left_id].y.upper);
			split_left = assign_I(dots[right_id].y.lower, dots[right_id].y.upper);
		}
	}
	else
	{
		if( dots[left_id].sign == 0 )
		{
			left_flank = assign_I(dots[left_id].y.lower, dots[left_id].y.upper);
			right_flank = assign_I(dots[right_id].y.lower, dots[right_id].y.upper);
		}
		else if( dots[left_id].sign == 1 )
		{
			right_flank = assign_I(dots[left_id].y.lower, dots[left_id].y.upper);
			left_flank = assign_I(dots[right_id].y.lower, dots[right_id].y.upper);
		}
		split_left = assign_I(dots[left_id].x.lower, dots[left_id].x.upper);
		split_right = assign_I(dots[right_id].x.lower, dots[right_id].x.upper);
	}

	if( (left_flank.lower < ins_reg.upper) && (ins_reg.lower < right_flank.upper) )
	{
		merged_left = assign_I(left_flank.lower, ins_reg.upper);
		merged_right = assign_I(ins_reg.lower, right_flank.upper);
	}
	else if( (right_flank.lower < ins_reg.upper) && (ins_reg.lower < left_flank.upper) )
	{
		merged_left = assign_I(right_flank.lower, ins_reg.upper);
		merged_right = assign_I(ins_reg.lower, left_flank.upper);
	}
	else 
	{
		temp_res = UNSUSPENDED;
		return temp_res;
	}

	if( (strict_almost_equal(dots[id].x, merged_left) == true ) || (strict_almost_equal(dots[id].x, merged_right) == true ) )
	{
		temp_res = COVER_INS_SP_IN_X;
	}

	if( (strict_almost_equal(dots[id].y, merged_left) == true ) || (strict_almost_equal(dots[id].y, merged_right) == true ) )
	{
		temp_res = COVER_INS_SP_IN_Y;
	}

	return(temp_res);
}
