#include "main.h"
#include "regions.h"
#include "update_init_algns.h"
#include "read_algn.h"
#include "util_gen.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "find_merging.h"
#include "id_ortho.h"
#include "apply_ops.h"
#include "util_ops.h"
#include "util_algns.h"

extern int debug_mode;

/* the first item in each chain should be the leftmost local alignment */
void mark_chain(struct DotList *dots, int id1, int id2, struct DotList *init_dots)
{
	int left_id, right_id;
	int old_id, cid;
	int l_pos1, l_pos2;

	if( dots[id1].x.lower < dots[id2].x.lower )
  {
    left_id = dots[id1].index;
    right_id = dots[id2].index;
  }
  else if( dots[id1].x.lower > dots[id2].x.lower )
  {
		left_id = dots[id2].index;
		right_id = dots[id1].index;
  }
  else
  {
    left_id = dots[id1].index;
    right_id = dots[id2].index; 
	}

	l_pos1 = init_dots[left_id].x.lower + init_dots[left_id].xl_diff;
	l_pos2 = init_dots[right_id].x.lower + init_dots[right_id].xl_diff;

	if( debug_mode == TRUE ) {
		if( l_pos1 > l_pos2 ) {
				printf("mismatched coordinates in %d-%d and %d-%d\n", dots[id1].x.lower, dots[id2].x.lower, l_pos1, l_pos2);
		}

		if( init_dots[left_id].m_id != -1 ) {
			printf("%d:%d-%d, %d-%d is not the beginning of the chain\n", init_dots[left_id].index, init_dots[left_id].x.lower, init_dots[left_id].x.upper, init_dots[left_id].y.lower, init_dots[left_id].y.upper);
		}

		if( init_dots[right_id].m_id != -1 ) {
			printf("%d:%d-%d, %d-%d is not the beginning of the chain\n", init_dots[right_id].index, init_dots[right_id].x.lower, init_dots[right_id].x.upper, init_dots[right_id].y.lower, init_dots[right_id].y.upper);
		}
	}

	old_id = init_dots[left_id].index;
	cid = init_dots[left_id].c_id;	
	while( cid != -1 ) {
		old_id = cid;
		cid = init_dots[cid].c_id;
	}

	if( debug_mode == TRUE ) {
		l_pos1 = init_dots[old_id].x.upper - init_dots[old_id].xr_diff;
		l_pos2 = init_dots[right_id].x.lower + init_dots[right_id].xl_diff;
		if( l_pos1 > l_pos2 ) {
			printf("%d-%d,%d-%d and %d-%d,%d-%d are not supposed to overlap: %d, %d\n", init_dots[old_id].x.lower, init_dots[old_id].x.upper, init_dots[old_id].y.lower, init_dots[old_id].y.upper, init_dots[right_id].x.lower, init_dots[right_id].x.upper, init_dots[right_id].y.lower, init_dots[right_id].y.upper, l_pos1, l_pos2);
		}

		l_pos1 = init_dots[old_id].y.upper - init_dots[old_id].yr_diff;
		l_pos2 = init_dots[right_id].y.lower + init_dots[right_id].yl_diff;
		if( l_pos1 > l_pos2 ) {
			printf("%d-%d,%d-%d and %d-%d,%d-%d are not supposed to overlap: %d, %d\n", init_dots[old_id].x.lower, init_dots[old_id].x.upper, init_dots[old_id].y.lower, init_dots[old_id].y.upper, init_dots[right_id].x.lower, init_dots[right_id].x.upper, init_dots[right_id].y.lower, init_dots[right_id].y.upper, l_pos1, l_pos2);
		}
	}

	if( old_id != init_dots[right_id].index ) {
		init_dots[old_id].c_id = init_dots[right_id].index;
		init_dots[right_id].m_id = init_dots[old_id].index;
	}
}

// ALL mode: when rolling-back for all local alignments
// INS mode: when rolling-back for an inserted copy
void update_init_algn(struct I reg, struct DotList *dots, int id, bool is_x, int cmp_id, struct DotList *init_dots, FILE *fp, int mode)
{
	int count = 1; // the number of the alignments linked in a chain
	int *cid_list;
	int b = -1, e = -1;
	int temp_b = -1, temp_e = -1;
	int org_id = dots[id].index;
	int cid = -1, old_id = -1;
	int cmp_org_id = dots[cmp_id].index;
	int cmp_cur_id = -1;
	struct I dup_reg;
	struct DotList cmp_algn;
	int case_flag; // CASE_1 - cmp_id is SELF, CASE_2 - cmp_id is PAIR and org_id is SELF1, and CASE_3 - cmp_id is PAIR and org_id is SELF2
	int temp_id = -1;
	int i = 0;

// if dots[id] is a chained alignment, init_dots[org_id] should be the beginning of the chained alignment.

	if( debug_mode == TRUE ) {
		if( init_dots[org_id].m_id != -1 ) {
			printf("id %d:%d-%d,%d-%d should be the beginning or the chained alignment\n", org_id, init_dots[org_id].x.lower, init_dots[org_id].x.upper, init_dots[org_id].y.lower, init_dots[org_id].y.upper);
		}

		if( init_dots[cmp_org_id].m_id != -1 ) {
			printf("cmp id %d:%d-%d,%d-%d should be the beginning or the chained alignment\n", cmp_org_id, init_dots[cmp_org_id].x.lower, init_dots[cmp_org_id].x.upper, init_dots[cmp_org_id].y.lower, init_dots[cmp_org_id].y.upper);
		}
	}

  count = 0;
  temp_id = init_dots[dots[cmp_id].index].c_id;
  while( temp_id != -1 )    
  {
    count++;
    temp_id = init_dots[temp_id].c_id;
  } 
  count++;

	cid_list = (int *) ckalloc(count * sizeof(int));

	for( i = 0; i < count; i++ ) cid_list[i] = 0;

	count = 0;
	cid_list[count] = dots[cmp_id].index;
	temp_id = init_dots[dots[cmp_id].index].c_id;
	while( temp_id != -1 )		
	{
		count++;
		cid_list[count] = temp_id;
		temp_id = init_dots[temp_id].c_id;
	}
	count++;

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

	if( (e - b) <= 0 ) {
		fatalf("update_init_algns: empty interval in (%d, %d)\n", b, e);
	}
	else dup_reg = assign_I(b, e);

	if( (mode == ALL) && (is_tandem(dots[id]) == true) ) {
		dup_reg = assign_I(reg.lower, reg.upper);
	}

	if( dots[cmp_id].sign == DELETED ) {
		if( mode == ALL ) {
			i = 0;
			cmp_org_id = dots[cmp_id].index;
			while( (i < count) && (cmp_org_id != -1) ) {
				temp_id = init_dots[cmp_org_id].c_id;
				if( check_in_cid_list(cmp_org_id, cid_list, count) == false ) {
					init_dots[cmp_org_id].sign = DELETED;
					if( (init_dots[cmp_org_id].c_id != -1) && (init_dots[cmp_org_id].m_id != -1) ) 
					{
						init_dots[init_dots[cmp_org_id].m_id].c_id = init_dots[cmp_org_id].c_id;
						init_dots[init_dots[cmp_org_id].c_id].m_id = init_dots[cmp_org_id].m_id;
						init_dots[cmp_org_id].m_id = -1;
						init_dots[cmp_org_id].c_id = -1;
					}
				}
				cmp_org_id = temp_id;
				i++;
			}
		}
	}
	else if( (init_dots[cmp_org_id].sp_id != PAIR) && (init_dots[cmp_org_id].sp_id != init_dots[org_id].sp_id) ) {
		fatalf("update_init_algn: incorrect suspend_list (%d, %d) and (%d, %d)\n", init_dots[cmp_org_id].x.lower, init_dots[cmp_org_id].x.upper, init_dots[cmp_org_id].y.lower, init_dots[cmp_org_id].y.upper);
	}
	else cmp_algn.sign = init_dots[cmp_org_id].sign;

	cmp_org_id = dots[cmp_id].index;
	temp_b = init_dots[cmp_org_id].x.lower + init_dots[cmp_org_id].xl_diff;
	temp_e = init_dots[cmp_org_id].x.upper - init_dots[cmp_org_id].xr_diff;
	if( (temp_e - temp_b) <= 0 ) {
		dots[cmp_id].sign = DELETED;
		init_dots[cmp_org_id].sign = DELETED;
		if( (init_dots[cmp_org_id].c_id != -1) && (init_dots[cmp_org_id].m_id != -1) ) {
			init_dots[init_dots[cmp_org_id].m_id].c_id = init_dots[cmp_org_id].c_id;
			init_dots[init_dots[cmp_org_id].c_id].m_id = init_dots[cmp_org_id].m_id;
			init_dots[cmp_org_id].m_id = -1;
			init_dots[cmp_org_id].c_id = -1;
		}
	}
	else cmp_algn.x = assign_I(temp_b, temp_e);
	
	temp_b = init_dots[cmp_org_id].y.lower + init_dots[cmp_org_id].yl_diff;
	temp_e = init_dots[cmp_org_id].y.upper - init_dots[cmp_org_id].yr_diff;

	if( (temp_e - temp_b) <= 0 ) {
		dots[cmp_id].sign = DELETED;
		init_dots[cmp_org_id].sign = DELETED;
		if( (init_dots[cmp_org_id].c_id != -1) && (init_dots[cmp_org_id].m_id != -1) ) {
			init_dots[init_dots[cmp_org_id].m_id].c_id = init_dots[cmp_org_id].c_id;
			init_dots[init_dots[cmp_org_id].c_id].m_id = init_dots[cmp_org_id].m_id;
			init_dots[cmp_org_id].m_id = -1;
			init_dots[cmp_org_id].c_id = -1;
		}
	}
	else cmp_algn.y = assign_I(temp_b, temp_e);

	i = 0;
//	if( mode == ALL ) count = 0;
	
	cmp_cur_id = dots[cmp_id].index;
	while ( i < count ) {
		if( dots[cmp_id].sign == DELETED ) {
			case_flag = -1;
		}
		else if( init_dots[cmp_cur_id].pair_self == SELF) case_flag = CASE_1;
		else if( (init_dots[cmp_cur_id].pair_self == PAIR) && (init_dots[org_id].sp_id == SELF1) ) {
			case_flag = CASE_2;
		}
		else if( (init_dots[cmp_cur_id].pair_self == PAIR) && (init_dots[org_id].sp_id == SELF2) ) {
			case_flag = CASE_3;
		}

		if( case_flag == CASE_1) {
			if(f_loose_subset(cmp_algn.x, dup_reg, LOOSE) || f_loose_subset(cmp_algn.y, dup_reg, LOOSE)) {
				init_dots[cmp_cur_id].sign = DELETED;
				if( (init_dots[cmp_cur_id].c_id != -1) && (init_dots[cmp_cur_id].m_id != -1) ) {
					init_dots[init_dots[cmp_cur_id].m_id].c_id = init_dots[cmp_cur_id].c_id;
					init_dots[init_dots[cmp_cur_id].c_id].m_id = init_dots[cmp_cur_id].m_id;
					init_dots[cmp_cur_id].m_id = -1;
					init_dots[cmp_cur_id].c_id = -1;
				}
			}
			else {
				if( proper_overlap(cmp_algn.x, dup_reg) || proper_overlap(cmp_algn.y, dup_reg) ) 
				{
					if( proper_overlap(cmp_algn.x, dup_reg) ) {
						adjust_init_algn(init_dots, cmp_cur_id, dup_reg, true, fp, DUP, is_x);
					}
					
					temp_b = init_dots[cmp_cur_id].y.lower + init_dots[cmp_cur_id].yl_diff;
					temp_e = init_dots[cmp_cur_id].y.upper - init_dots[cmp_cur_id].yr_diff;
					if( (temp_e - temp_b) <= 0 ) {
						init_dots[cmp_cur_id].sign = DELETED;
						if( (init_dots[cmp_cur_id].c_id != -1) && (init_dots[cmp_cur_id].m_id != -1) ) {
							init_dots[init_dots[cmp_cur_id].m_id].c_id = init_dots[cmp_cur_id].c_id;
							init_dots[init_dots[cmp_cur_id].c_id].m_id = init_dots[cmp_cur_id].m_id;
							init_dots[cmp_cur_id].m_id = -1;
							init_dots[cmp_cur_id].c_id = -1;
						}
					}
					else {
						cmp_algn.y = assign_I(temp_b, temp_e);
						if( proper_overlap(cmp_algn.y, dup_reg) ) {
							adjust_init_algn(init_dots, cmp_cur_id, dup_reg, false, fp, DUP, is_x);
						}
					}
				}
//				else if( f_loose_subset(dup_reg, cmp_algn.x, LOOSE) || f_loose_subset(dup_reg, cmp_algn.y, LOOSE) ) {
				else if( subset(dup_reg, cmp_algn.x) || subset(dup_reg, cmp_algn.y) ) {
					fatalf("update_init_algns: incorrect inferrence 1 in (%d, %d)\n", dup_reg.lower, dup_reg.upper);
				}
			}
		}
		else if( case_flag == CASE_2) {
			if(f_loose_subset(cmp_algn.x, dup_reg, LOOSE)) {
				init_dots[cmp_cur_id].sign = DELETED;
				if( (init_dots[cmp_cur_id].c_id != -1) && (init_dots[cmp_cur_id].m_id != -1))  {
					init_dots[init_dots[cmp_cur_id].m_id].c_id = init_dots[cmp_cur_id].c_id;
					init_dots[init_dots[cmp_cur_id].c_id].m_id = init_dots[cmp_cur_id].m_id;
					init_dots[cmp_cur_id].m_id = -1;
					init_dots[cmp_cur_id].c_id = -1;
				}
			}
			else if( proper_overlap(cmp_algn.x, dup_reg) ) {
				adjust_init_algn(init_dots, cmp_cur_id, dup_reg, true, fp, DUP, is_x);
			}
		}
		else if( case_flag == CASE_3) {
			if(f_loose_subset(cmp_algn.y, dup_reg, LOOSE)) {
				init_dots[cmp_cur_id].sign = DELETED;
				if( (init_dots[cmp_cur_id].c_id != -1) && (init_dots[cmp_cur_id].m_id != -1) ) {
					init_dots[init_dots[cmp_cur_id].m_id].c_id = init_dots[cmp_cur_id].c_id;
					init_dots[init_dots[cmp_cur_id].c_id].m_id = init_dots[cmp_cur_id].m_id;
					init_dots[cmp_cur_id].m_id = -1;
					init_dots[cmp_cur_id].c_id = -1;
				}
			}
			else if( proper_overlap(cmp_algn.y, dup_reg) ) {
				adjust_init_algn(init_dots, cmp_cur_id, dup_reg, false, fp, DUP, is_x);
			}
		}
		i++;

		if((count > 1) && (i < count)) {
			cmp_cur_id = cid_list[i];
			if( init_dots[cmp_cur_id].sign == DELETED ) {}
			else if( (init_dots[cmp_cur_id].sp_id != PAIR) && (init_dots[cmp_cur_id].sp_id != init_dots[org_id].sp_id) ) {
				fatalf("update_init_algn: incorrect suspend_list (%d, %d) and (%d, %d)\n", init_dots[cmp_cur_id].x.lower, init_dots[cmp_cur_id].x.upper, init_dots[cmp_cur_id].y.lower, init_dots[cmp_cur_id].y.upper);
			}
			else cmp_algn.sign = init_dots[cmp_cur_id].sign;

			temp_b = init_dots[cmp_cur_id].x.lower + init_dots[cmp_cur_id].xl_diff;
			temp_e = init_dots[cmp_cur_id].x.upper - init_dots[cmp_cur_id].xr_diff;
			if( (temp_e - temp_b) <= 0 ) {
				init_dots[cmp_cur_id].sign = DELETED;
				fatalf("update_init_algns: empty interval 1 in (%d, %d)\n", temp_b, temp_e);
			}
			else cmp_algn.x = assign_I(temp_b, temp_e);
	
			temp_b = init_dots[cmp_cur_id].y.lower + init_dots[cmp_cur_id].yl_diff;
			temp_e = init_dots[cmp_cur_id].y.upper - init_dots[cmp_cur_id].yr_diff;

			if( (temp_e - temp_b) <= 0 ) {
				init_dots[cmp_cur_id].sign = DELETED;
				fatalf("update_init_algns: empty interval 2 in (%d, %d)\n", temp_b, temp_e);
			}
			else cmp_algn.y = assign_I(temp_b, temp_e);
		}
	}

	free(cid_list);
}

void adjust_init_algn(struct DotList *init_dots, int id, struct I dup_reg, bool is_x, FILE *fp, int mode, bool dup_is_x) 
{
	int cur = 0, old = 0;
	int loc = 0;
	int old_diff = 0, new_diff = 0;
	struct I cur_reg = {0, 1};
	int b = 0, e = 1;
	int threshold = 0;
	int cut_point = 0;
	
	if( is_x == true ) {
		b = init_dots[id].x.lower + init_dots[id].xl_diff; 
		e = init_dots[id].x.upper - init_dots[id].xr_diff;
		if( b >= e ) {
			if( init_dots[id].sign != DELETED ) init_dots[id].sign = DELETED;
		}
		else {
			cur_reg = assign_I(b,e);

			if( f_loose_subset(cur_reg, dup_reg, LOOSE) == true ) {
				if( mode == DEL ) {
					fatalf("incorrect inferrence of deletion in %d-%d\n", dup_reg.lower, dup_reg.upper);	
				}
				else init_dots[id].sign = DELETED;
			}
			else if( ((mode == DUP) && (f_loose_subset(dup_reg, cur_reg, STRICT) == true)) ) // for tandem duplication
			{
				threshold = 5 * LOOSE * UNIT_THRESHOLD; 	

				if((abs(dup_reg.lower - cur_reg.lower) <= threshold) && (abs(dup_reg.upper - cur_reg.upper) <= threshold)) {
					init_dots[id].sign = DELETED;
				}
				else if(abs(dup_reg.lower - cur_reg.lower) <= threshold) {
					cur_reg = assign_I(dup_reg.lower+1, e);	
				}
				else if( abs(dup_reg.upper - cur_reg.upper) <= threshold) {
					cur_reg = assign_I(b, dup_reg.upper-1);
				}
				else if( (dup_is_x == false) && (abs(dup_reg.lower - cur_reg.lower) <= width(dup_reg)) ) 
				{
					cur_reg = assign_I(dup_reg.lower+1, e);	
				}
				else if( (dup_is_x == true) && (abs(dup_reg.upper - cur_reg.upper) <= width(dup_reg))) {
					cur_reg = assign_I(b, dup_reg.upper-1);
				}
				else {
					init_dots[id].sign = DELETED;
//      		fatalf("update_init_algns: %d-%d contains a dup region %d-%d\n", cur_reg.lower, cur_reg.upper, dup_reg.lower, dup_reg.upper);
				}
			}
			else if( (mode == DEL) && (subset(dup_reg, cur_reg) == true) )  {
				if( (abs(cur_reg.upper-dup_reg.upper) <= DEL_TH) ) {
					cur_reg = assign_I(cur_reg.lower, dup_reg.upper-1);
				}
				else if(abs(dup_reg.lower-cur_reg.lower) <= DEL_TH) {
					cur_reg = assign_I(dup_reg.lower+1, cur_reg.upper);
				}
				else {
					fatalf("incorrect inferrence of deletion in %d-%d\n", dup_reg.lower, dup_reg.upper);
				}
			}

			if( mode == DEL ) cut_point = (dup_reg.lower + dup_reg.upper)/2;

			if( ((mode == DUP) && (cur_reg.lower > dup_reg.lower)) || ((mode == DEL) && (in(cur_reg.lower, dup_reg) == true) && (cur_reg.lower < cut_point) == true)  ) {
				old_diff = init_dots[id].xl_diff;
				if( mode == DUP) new_diff = dup_reg.upper - init_dots[id].x.lower;
				else if( mode == DEL ) new_diff = cut_point - init_dots[id].x.lower;
				else new_diff = 0;

				if( new_diff > old_diff ) { 
					init_dots[id].xl_diff = new_diff;
					if( init_dots[id].sign == 0 ) {
						if( init_dots[id].xl_offset == 0 ) {
							old = init_dots[id].y.lower; // y.lower - yl_offset is the original coordinate
							cur = find_yloc_one(init_dots[id], fp, init_dots[id].xl_diff, NO_GAP_INC);
						}
						else {
							old = find_yloc_one(init_dots[id], fp, init_dots[id].xl_offset, NO_GAP_INC);
							cur = find_yloc_one(init_dots[id], fp, init_dots[id].xl_diff + init_dots[id].xl_offset, NO_GAP_INC);
						}

						if( cur > old ) {
							init_dots[id].yl_diff = cur - old;
						}
						else {
							if( debug_mode == TRUE ) printf("1:a new offset is negative: %d\n", abs(old - cur));
						}
					}
					else if( init_dots[id].sign == 1 ) {
						if( init_dots[id].xl_offset == 0 ) {
							old = init_dots[id].y.upper; // y.upper - yr_offset is original
							cur = find_yloc_one(init_dots[id], fp, init_dots[id].xl_diff, NO_GAP_INC);
						}
						else {
							old = find_yloc_one(init_dots[id], fp, init_dots[id].xl_offset, NO_GAP_INC);
							cur = find_yloc_one(init_dots[id], fp, init_dots[id].xl_diff+init_dots[id].xl_offset, NO_GAP_INC);
						}

						if( cur < old ) init_dots[id].yr_diff = old - cur;
						else {
							if( debug_mode == TRUE ) printf("2:a new offset is negative: %d\n", abs(old - cur));
						}
					}
				}
			}
			else if( ((mode == DUP) && (cur_reg.upper < dup_reg.upper)) || ((mode == DEL) && (in(cur_reg.upper, dup_reg) == true) && (cur_reg.upper > cut_point)) ) {
				old_diff = init_dots[id].xr_diff;
				if( mode == DUP ) {
					new_diff = init_dots[id].x.upper - dup_reg.lower;
					cut_point = dup_reg.lower;
				}	
				else if( mode == DEL ) new_diff = cur_reg.upper - cut_point;

				if( new_diff > old_diff ) {
					init_dots[id].xr_diff = new_diff;
					if( init_dots[id].sign == 0 ) {
						if( (init_dots[id].xl_offset == 0) && (init_dots[id].xr_offset == 0) ) {
							old = init_dots[id].y.upper;
							cur = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].x.lower, NO_GAP_INC);
						}
						else {
							old = find_yloc_one(init_dots[id], fp, init_dots[id].x.upper - init_dots[id].x.lower + init_dots[id].xl_offset, NO_GAP_INC);
							cur = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].x.lower + init_dots[id].xl_offset, NO_GAP_INC);
						}

						if( cur < old ) init_dots[id].yr_diff = old - cur;
						else {
							if( debug_mode == TRUE ) printf("3:a new offset is negative: %d\n", abs(old - cur));
						}
					}
					else if( init_dots[id].sign == 1 ) {
						if( (init_dots[id].xl_offset == 0) && (init_dots[id].xr_offset == 0) ) {
							old = init_dots[id].y.lower;
							cur = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].x.lower, NO_GAP_INC);
						}
						else {
							old = find_yloc_one(init_dots[id], fp, init_dots[id].x.upper - init_dots[id].x.lower + init_dots[id].xl_offset, NO_GAP_INC);
							cur = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].x.lower + init_dots[id].xl_offset, NO_GAP_INC);
						}

						if( cur > old) {
							init_dots[id].yl_diff = cur - old;
						}
						else {
							if( debug_mode == TRUE ) printf("4:a new offset is negative: %d\n", abs(old - cur));
						}
					}
				}
			}
		}
	}
	else {
		b = init_dots[id].y.lower + init_dots[id].yl_diff; 
		e = init_dots[id].y.upper - init_dots[id].yr_diff;
		if( b >= e ) {
			if( init_dots[id].sign != DELETED ) init_dots[id].sign = DELETED;
		}
		else {
			cur_reg = assign_I(b,e);

			if( f_loose_subset(cur_reg, dup_reg, LOOSE) == true ) {
				if( mode == DEL ) {
					fatalf("incorrect inferrence of deletion in %d-%d\n", dup_reg.lower, dup_reg.upper);	
				}
				else init_dots[id].sign = DELETED;
			}
			else if( (mode == DUP) && (f_loose_subset(dup_reg, cur_reg, STRICT) == true) ) {
				threshold = 5 * LOOSE * UNIT_THRESHOLD; 	

				if((abs(dup_reg.lower - cur_reg.lower) <= threshold) && (abs(dup_reg.upper - cur_reg.upper) <= threshold)) {
					init_dots[id].sign = DELETED;
				}
				else if(abs(dup_reg.lower - cur_reg.lower) <= threshold) {
					cur_reg = assign_I(dup_reg.lower+1, e);	
				}
				else if(abs(dup_reg.upper - cur_reg.upper) <= threshold) {
					cur_reg = assign_I(b, dup_reg.upper-1);
				}
				else if( (dup_is_x == false) && (abs(dup_reg.lower - cur_reg.lower) <= width(dup_reg)) ) 
				{
					cur_reg = assign_I(dup_reg.lower+1, e);	
				}
				else if( (dup_is_x == true) && (abs(dup_reg.upper - cur_reg.upper) <= width(dup_reg))) {
					cur_reg = assign_I(b, dup_reg.upper-1);
				}
				else {
					init_dots[id].sign = DELETED;
//      		fatalf("update_init_algns: %d-%d contains a dup region %d-%d\n", cur_reg.lower, cur_reg.upper, dup_reg.lower, dup_reg.upper);
				}
			}
			else if( (mode == DEL) && (subset(dup_reg, cur_reg) == true) )  {
				fatalf("incorrect inferrence of deletion in %d-%d\n", dup_reg.lower, dup_reg.upper);
			}

			if( mode == DEL ) cut_point = (dup_reg.lower + dup_reg.upper)/2;

			if( ((mode == DUP) && (cur_reg.lower < dup_reg.lower)) || ((mode == DEL) && (in(cur_reg.upper, dup_reg) == true) && (cur_reg.upper > cut_point)) ) {
				old_diff = init_dots[id].yr_diff;
				if( mode == DUP ) {	
					new_diff = init_dots[id].y.upper - dup_reg.lower;
					cut_point = dup_reg.lower;
				}
				else if( mode == DEL ) new_diff = init_dots[id].y.upper - cut_point;
				
				if( new_diff > old_diff ) {
					init_dots[id].yr_diff = new_diff;
					if( init_dots[id].sign == 0 ) {
						if( (init_dots[id].yl_offset == 0 ) && (init_dots[id].yr_offset == 0) ) {
							old = init_dots[id].x.upper;
							loc = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].y.lower, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}
						else {
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].y.upper - init_dots[id].y.lower + init_dots[id].yl_offset, GAP_INC_IN_Y);
							old = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
							loc = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].y.lower + init_dots[id].yl_offset, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}

						if( old > cur ) init_dots[id].xr_diff = old - cur;
						else {
							if( debug_mode == TRUE ) printf("5:a new offset is negative: %d\n", abs(old - cur));
						}
					}
					else if( init_dots[id].sign == 1 ) {
						if( (init_dots[id].yl_offset == 0 ) && (init_dots[id].yr_offset == 0) ) {
							old = init_dots[id].x.lower;
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].y.upper - cut_point, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}
						else {
							loc = find_yloc_one(init_dots[id], fp, abs(init_dots[id].yr_offset), GAP_INC_IN_Y);
							old = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].y.upper - cut_point + abs(init_dots[id].yr_offset), GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}

						if( cur > old ) init_dots[id].xl_diff = cur - old;
						else {
							if( debug_mode == TRUE ) printf("6:a new offset is negative: %d\n", abs(old - cur));
						}
					}
				}
			}
			else if( ((mode == DUP) && (cur_reg.upper > dup_reg.upper)) || ((mode == DEL) && (in(cur_reg.lower, dup_reg) == true ) && (cur_reg.lower < cut_point)) ) {
				old_diff = init_dots[id].yl_diff;
				if( mode == DUP ) {	
					new_diff = dup_reg.upper - init_dots[id].y.lower;
					cut_point = dup_reg.upper;
				}
				else if( mode == DEL ) new_diff = cut_point - init_dots[id].y.lower;

				if( new_diff > old_diff ) {
					init_dots[id].yl_diff = new_diff;
					if( init_dots[id].sign == 0 )  {
						if( (init_dots[id].yl_offset == 0 ) && (init_dots[id].yr_offset == 0) ) {
							old = init_dots[id].x.lower;
							loc = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].y.lower, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}
						else {
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].yl_offset, GAP_INC_IN_Y);
							old = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
							loc = find_yloc_one(init_dots[id], fp, cut_point - init_dots[id].y.lower + init_dots[id].yl_offset, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}

						if( cur > old ) init_dots[id].xl_diff = cur - old;
						else {
							if( debug_mode == TRUE ) printf("7:a new offset is negative: %d\n", abs(old - cur));
						}
					}
					else if( init_dots[id].sign == 1 ) {
						if( (init_dots[id].yl_offset == 0 ) && (init_dots[id].yr_offset == 0) ) {
							old = init_dots[id].x.upper;
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].y.upper - cut_point, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}
						else {
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].y.upper - init_dots[id].yr_offset - init_dots[id].y.lower, GAP_INC_IN_Y);
							old = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
							loc = find_yloc_one(init_dots[id], fp, init_dots[id].y.upper - init_dots[id].yr_offset - cut_point, GAP_INC_IN_Y);
							cur = find_xloc_one(init_dots[id], fp, loc, GAP_INC);
						}

						if( old > cur ) init_dots[id].xr_diff = old - cur;
						else {
							if( debug_mode == TRUE ) printf("8:a new offset is negative: %d\n", abs(old - cur));
						}
					}
				}
			}
		}
	}	
}

void rollback_init_dots(struct DotList *algns, int id, bool is_x, int num_algns, struct DotList *init_algns, int num_init_algns, FILE *fp, int num_ops, struct ops_list *ops, int mode, int size)
{  
	int i;  
	bool is_skip = false;
	struct I to, from, dup_reg;
	int cid = -1;
	int cur_id = -1, tmp_id = -1;
	int b = -1, e = -1;
	int src_b = -1, src_e = -1;
	int old_id = -1, org_id = -1;
	int cur_num = 0; 
	char sign;
	int pid = 0;
	int num_temp_ops = 1;
	struct ops_list *temp_ops;
	struct ops_list *new_ops;

	cur_num = num_ops;
	temp_ops = (struct ops_list *) ckalloc(sizeof(struct ops_list));
	
	if( ( id < 0 ) || ( id >= num_algns ) ) fatalf("rollback failure for %dth alignment", id);

	if( is_x ) {
		to = assign_I(algns[id].x.lower, algns[id].x.upper);
		from = assign_I(algns[id].y.lower, algns[id].y.upper);
	}
	else {
		from = assign_I(algns[id].x.lower, algns[id].x.upper);
		to = assign_I(algns[id].y.lower, algns[id].y.upper);
	}
	
  org_id = algns[id].index;  
	old_id = init_algns[org_id].index;  
	cid = init_algns[org_id].c_id;  
	while( cid != -1 ) {    
		old_id = cid;    
		cid = init_algns[cid].c_id;  
	}  
	if( (algns[id].sign != DELETED) && (init_algns[org_id].sign == DELETED) ) {    
		init_algns[org_id].sign = algns[id].sign; // after rolling back a tandem duplication, a converted alignment in the step of handling tandem duplications could be deleted during the rollback process  
  }

	if( init_algns[org_id].sign == 0 ) sign = '+';
	else if( init_algns[org_id].sign == 1 ) sign = '-';
	else {
		fatalf("unsupported sign %c\n", init_algns[org_id].sign);
	}
	pid = init_algns[org_id].identity;

  if( is_x == true ) {
    b = init_algns[org_id].x.lower + init_algns[org_id].xl_diff;
    e = init_algns[old_id].x.upper - init_algns[old_id].xr_diff;
    if( init_algns[org_id].sign == 0 ) {
      src_b = init_algns[org_id].y.lower + init_algns[org_id].yl_diff;
      src_e = init_algns[old_id].y.upper - init_algns[old_id].yr_diff;
    }
    else if( init_algns[org_id].sign == 1 ) {
      src_b = init_algns[old_id].y.lower + init_algns[old_id].yl_diff;
      src_e = init_algns[org_id].y.upper - init_algns[org_id].yr_diff;
    }
  }
  else {
    if( init_algns[org_id].sign == 0 ) {
      b = init_algns[org_id].y.lower + init_algns[org_id].yl_diff;
      e = init_algns[old_id].y.upper - init_algns[old_id].yr_diff;
    	src_b = init_algns[org_id].x.lower + init_algns[org_id].xl_diff;
    	src_e = init_algns[old_id].x.upper - init_algns[old_id].xr_diff;
    }
    else if( init_algns[org_id].sign == 1 ) {
      b = init_algns[old_id].y.lower + init_algns[old_id].yl_diff;
      e = init_algns[org_id].y.upper - init_algns[org_id].yr_diff;
    	src_b = init_algns[org_id].x.lower + init_algns[org_id].xl_diff;
    	src_e = init_algns[old_id].x.upper - init_algns[old_id].xr_diff;
    }
    else {
      fatalf("update_init_algns: invalid sign in %d and %d\n", org_id, old_id);
    }
  }

  if( debug_mode == TRUE ) printf("Interval in the original dot plot: [%d, %d]\n", b, e);

	if( (mode == SECOND_RUN) && (is_tandem(algns[id]) == false)) {
		if( (b != -1) && (e != -1) && (src_b != -1) && (src_e != -1) ) {
			temp_ops[0].srcStart = b;
			temp_ops[0].srcEnd = e;
			temp_ops[0].dstStart = src_b;
			temp_ops[0].dstEnd = src_e;
			temp_ops[0].src_b = b;
			temp_ops[0].src_e = e;
			temp_ops[0].dst_b = src_b;
			temp_ops[0].dst_e = src_e;
			temp_ops[0].id = org_id;
			temp_ops[0].sign = sign;
			temp_ops[0].pid = pid;
			redo_dups_for_mtom(num_temp_ops, temp_ops, num_init_algns, init_algns, fp, init_algns[org_id].sp_id);
		}
	}

	ops[cur_num].id = org_id;
	ops[cur_num].srcStart = src_b;
	ops[cur_num].srcEnd = src_e;
	ops[cur_num].dstStart = b;
	ops[cur_num].dstEnd = e;
	ops[cur_num].src_b = from.lower;
	ops[cur_num].src_e = from.upper;
	ops[cur_num].dst_b = to.lower;
	ops[cur_num].dst_e = to.upper;
	ops[cur_num].sign = sign;
	ops[cur_num].pid = pid;
	ops[cur_num].sp_id = algns[id].sp_id;

	new_ops = (struct ops_list *) ckalloc((cur_num+1) * sizeof(struct ops_list));
	init_ops(new_ops, 0, cur_num+1);
	if( algns[id].sp_id == SELF1 ) {
		cur_num = cal_cur_pos_ops(cur_num+1, ops, new_ops, SELF1, 0);
	}
	else if( algns[id].sp_id == SELF2 ) {
		cur_num = cal_cur_pos_ops(cur_num+1, ops, new_ops, SELF2, size);
	}
	else {
		fatalf("unexpected species id: %d\n", algns[id].sp_id);
	}

///// if the alignment is chained one, the endpoints of a region need to be traced 
	cur_id = algns[id].index;
	dup_reg = assign_I(new_ops[cur_num-1].dstStart, new_ops[cur_num-1].dstEnd);
	for( i = 0 ; i < num_algns; i++ ) { // Excluding the test for split region newly  
		tmp_id = algns[i].index;
		is_skip = false;
		if((i != id) && is_tandem(algns[id]) && is_tandem(algns[i])) {
//			if(f_loose_subset(algns[i].x, to, LOOSE) || f_loose_subset(algns[i].y, to, LOOSE)) { /* Apr 23, 2015
			if(((is_x == true ) && f_loose_subset(algns[i].x, to, LOOSE)) || ((is_x == false) && f_loose_subset(algns[i].y, to, LOOSE))) {
				is_skip = true;
			}
		}

		if( (i != id) && (!is_skip) && ((algns[i].sp_id == PAIR) || (algns[i].sp_id == algns[id].sp_id))) {
			if( algns[i].sp_id == PAIR )
			{
				if( debug_mode == TRUE ) {
					printf("this is an inter-species alignment\n");
				}
			}
			update_init_algn(dup_reg, algns, id, is_x, i, init_algns, fp, ALL);
		}

 		b = init_algns[tmp_id].x.lower + init_algns[tmp_id].xl_diff;
    e = init_algns[tmp_id].x.upper - init_algns[tmp_id].xr_diff;
		if( (e - b) <= 0 ) {
			init_algns[tmp_id].sign = DELETED;
			if( (init_algns[tmp_id].c_id != -1) && (init_algns[tmp_id].m_id != -1)) {
				init_algns[init_algns[tmp_id].m_id].c_id = init_algns[tmp_id].c_id;
				init_algns[init_algns[tmp_id].c_id].m_id = init_algns[tmp_id].m_id;
				init_algns[tmp_id].m_id = -1;
				init_algns[tmp_id].c_id = -1;
			}
		}	

 		b = init_algns[tmp_id].y.lower + init_algns[tmp_id].yl_diff;
    e = init_algns[tmp_id].y.upper - init_algns[tmp_id].yr_diff;
		if( (e - b) <= 0 ) {
			init_algns[tmp_id].sign = DELETED;
			if( (init_algns[tmp_id].c_id != -1) && (init_algns[tmp_id].m_id != -1)) {
				init_algns[init_algns[tmp_id].m_id].c_id = init_algns[tmp_id].c_id;
				init_algns[init_algns[tmp_id].c_id].m_id = init_algns[tmp_id].m_id;
				init_algns[tmp_id].m_id = -1;
				init_algns[tmp_id].c_id = -1;
			}
		}	
	}	

	cur_id = algns[id].index;
	init_algns[cur_id].sign = DELETED;
	cid = init_algns[cur_id].c_id; 
	while( cid != -1 ) {
		init_algns[cid].sign = DELETED;
		cur_id = cid;
		cid = init_algns[cid].c_id;
		init_algns[cur_id].c_id = -1;
		init_algns[cur_id].m_id = -1;
	}
	init_algns[cur_id].c_id = -1;
	init_algns[cur_id].m_id = -1;

	free(new_ops);
	free(temp_ops);
}

// flag is SELF1, SELF2 or PAIR
void write_init_maf(FILE *output_f, struct DotList *init_algns, int num_algns, struct n_pair *contigs1, struct n_pair *contigs2, int len1, int len2, FILE *fp, int flag, char *species, char *species2)
{
	char S2[BIG], T2[BIG];
	struct b_list *a_info;
	int i = 0, j = 0;
	int beg = 0, end = 0;
	int pid = 0;
	int y_b = 0, y_e = 0;
	int e1 = 0, e2 = 0;
	int gap1 = 0, gap2 = 0;
	int skip_nu = 0;
	int len = 0;
	int ctg_id1 = -1, ctg_id2 = -1;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;

  fprintf(output_f, "##maf version=1 scoring=lastz-pid %s %s\n", species, species2);

	for( i = 0; i < num_algns; i++ ) 
	{
		if( (init_algns[i].sp_id == flag) && (init_algns[i].sign != DELETED)) 
		{
			end = 0;
			if( init_algns[i].rp1_id != -1 ) {
				beg = init_algns[i].xl_diff + init_algns[i].xl_offset;
				skip_nu = find_xloc_one(init_algns[i], fp, init_algns[i].xl_diff + init_algns[i].xl_offset, NO_GAP_INC);
				end = find_end_xloc(init_algns[i], fp, init_algns[i].x.upper - init_algns[i].x.lower - init_algns[i].xr_diff - init_algns[i].xr_offset);
			}
			else {
				beg = init_algns[i].xl_diff;
				skip_nu = find_xloc_one(init_algns[i], fp, init_algns[i].xl_diff, NO_GAP_INC);
				end = find_end_xloc(init_algns[i], fp, init_algns[i].x.upper - init_algns[i].x.lower - init_algns[i].xr_diff);
			}
			(*a_info).b1 = init_algns[i].x.lower + init_algns[i].xl_diff - 1;
			y_b = init_algns[i].y.lower + init_algns[i].yl_diff;
			y_e = init_algns[i].y.upper - init_algns[i].yr_diff;
      ctg_id1 = init_algns[i].ctg_id1;
      ctg_id2 = init_algns[i].ctg_id2;

			strcpy(S2, "");
			strcpy(T2, "");
			get_nth_algn(S2, T2, init_algns[i].fid, beg, fp, a_info, REG);

      if( ctg_id1 == -1 ) (*a_info).len1 = len1;
      else (*a_info).len1 = contigs1[ctg_id1].len;

      if( (*a_info).len1 != init_algns[i].len1 ) {
        fatalf("Contig length not match: %s %d-%d\n", contigs1[ctg_id1].name2, len1, init_algns[i].len1);
      }

      if( ctg_id2 == -1 ) (*a_info).len2 = len2;
      else (*a_info).len2 = contigs2[ctg_id2].len;

      if( (*a_info).len2 != init_algns[i].len2 ) {
        fatalf("Contig length not match: %s %d-%d\n", contigs2[ctg_id2].name2, len2, init_algns[i].len2);
      }

      end = end - skip_nu;
      if( end <= 0 ) {
        end = 0;
      }
      S2[end] = '\0';
      T2[end] = '\0';

/*
			if( init_algns[i].sign == 0 ) {
				(*a_info).b2 = y_b-1;
				(*a_info).strand = '+';
			}
			else if( init_algns[i].sign == 1 ) {
				(*a_info).b2 = (a_info->len2) - y_e + 1;
				(*a_info).strand = '-';
			}
*/
      if( init_algns[i].sign == 1 ) {
        (*a_info).b2 = (a_info->len2) - (a_info->b2) + 2;
        (*a_info).strand = '-';
      }

		  e1 = 0;
		  e2 = 0;
		  gap1 = 0;
		  gap2 = 0;
		  for( j = 0; (j < end) && ((S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n')); j++ ) {
    		if( strchr("ACGTN", toupper(S2[j])) ) e1++;
				if( S2[j] == '-' ) gap1++;
    		if( strchr("ACGTN", toupper(T2[j])) ) e2++;
				if( T2[j] == '-' ) gap2++;
  		}

			if( (e1 >= MIN_GAP) && (e2 >= MIN_GAP) ) {
				(*a_info).e1 = e1;
				(*a_info).e2 = e2;
			
   	 	  pid = (int)(cal_pid_maf(S2, T2, end) + 0.5);

        if( (((a_info->b1)-1+(a_info->e1)) > (a_info->len1)) || ((a_info->b1-1) < 0 ) ) {
					fatalf("line 930: Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
				}

     		fprintf(output_f, "a srcblock=%d\n", init_algns[i].indiv_fid);
        fprintf(output_f, "s %s %d %d + %d %s", init_algns[i].name1, (a_info->b1)-1, a_info->e1, a_info->len1, S2);
        len = strlen(S2);
        if( S2[len-1] != '\n' ) fprintf(output_f, "\n");

        fprintf(output_f, "s %s %d %d %c %d %s", init_algns[i].name2, (a_info->b2)-1, a_info->e2, a_info->strand, a_info->len2, T2);
        len = strlen(T2);

      	if( T2[len-1] != '\n' ) fprintf(output_f, "\n");
      	fprintf(output_f, "\n");
			}
		}
	}

	free(a_info);
}

bool check_in_cid_list(int index, int *list, int num_list) 
{
	bool res = false;
	int i;

	for( i = 0; i < num_list; i++ ) {
		if( index == list[i] ) res = true;
	}

	return res;
}

void update_init_algn_del(struct DotList *dots, int cmp_id, struct I del_reg, bool is_x, FILE *fp, struct DotList *init_dots)
{
	int count = 1; // the number of the alignments linked in a chain
	int *cid_list;
	int temp_b, temp_e;
	int index = dots[cmp_id].index;
	int cmp_org_id = index;
	int cmp_cur_id;
	struct DotList cmp_algn;
	int temp_id;
	int i;

// if dots[id] is a chained alignment, init_dots[org_id] should be the beginning of the chained alignment.

	if( debug_mode == TRUE ) {
		if( init_dots[cmp_org_id].m_id != -1 ) {
			printf("cmp id %d:%d-%d,%d-%d should be the beginning or the chained alignment\n", cmp_org_id, init_dots[cmp_org_id].x.lower, init_dots[cmp_org_id].x.upper, init_dots[cmp_org_id].y.lower, init_dots[cmp_org_id].y.upper);
		}
	}

  count = 0;
  temp_id = init_dots[index].c_id;
  while( temp_id != -1 )    
  {
    count++;
    temp_id = init_dots[temp_id].c_id;
  } 
  count++;

	cid_list = (int *) ckalloc(count * sizeof(int));

	for( i = 0; i < count; i++ ) cid_list[i] = 0;

	count = 0;
	cid_list[count] = index;
	temp_id = init_dots[index].c_id;
	while( temp_id != -1 )		
	{
		count++;
		cid_list[count] = temp_id;
		temp_id = init_dots[temp_id].c_id;
	}
	count++;

	cmp_algn.sign = init_dots[cmp_org_id].sign;

	temp_b = init_dots[cmp_org_id].x.lower + init_dots[cmp_org_id].xl_diff;
	temp_e = init_dots[cmp_org_id].x.upper - init_dots[cmp_org_id].xr_diff;
	if( (temp_e - temp_b) <= 0 ) {
		dots[cmp_id].sign = DELETED;
		init_dots[cmp_org_id].sign = DELETED;
		if( (init_dots[cmp_org_id].c_id != -1) && (init_dots[cmp_org_id].m_id != -1) ) {
			init_dots[init_dots[cmp_org_id].m_id].c_id = init_dots[cmp_org_id].c_id;
			init_dots[init_dots[cmp_org_id].c_id].m_id = init_dots[cmp_org_id].m_id;
			init_dots[cmp_org_id].m_id = -1;
			init_dots[cmp_org_id].c_id = -1;
		}
	}
	else cmp_algn.x = assign_I(temp_b, temp_e);
	
	temp_b = init_dots[cmp_org_id].y.lower + init_dots[cmp_org_id].yl_diff;
	temp_e = init_dots[cmp_org_id].y.upper - init_dots[cmp_org_id].yr_diff;

	if( (temp_e - temp_b) <= 0 ) {
		dots[cmp_id].sign = DELETED;
		init_dots[cmp_org_id].sign = DELETED;
		if( (init_dots[cmp_org_id].c_id != -1) && (init_dots[cmp_org_id].m_id != -1) ) {
			init_dots[init_dots[cmp_org_id].m_id].c_id = init_dots[cmp_org_id].c_id;
			init_dots[init_dots[cmp_org_id].c_id].m_id = init_dots[cmp_org_id].m_id;
			init_dots[cmp_org_id].m_id = -1;
			init_dots[cmp_org_id].c_id = -1;
		}
	}
	else cmp_algn.y = assign_I(temp_b, temp_e);

	i = 0;
	if( dots[cmp_id].sign == DELETED ) count = 0;
	
	cmp_cur_id = cmp_org_id;
	while ( i <= count ) {
		if( is_x == true) {
			if(f_loose_subset(cmp_algn.x, del_reg, LOOSE)) {
				init_dots[cmp_cur_id].sign = DELETED;
				if( (init_dots[cmp_cur_id].c_id != -1) && (init_dots[cmp_cur_id].m_id != -1))  {
					init_dots[init_dots[cmp_cur_id].m_id].c_id = init_dots[cmp_cur_id].c_id;
					init_dots[init_dots[cmp_cur_id].c_id].m_id = init_dots[cmp_cur_id].m_id;
					init_dots[cmp_cur_id].m_id = -1;
					init_dots[cmp_cur_id].c_id = -1;
				}
			}
			else if( proper_overlap(cmp_algn.x, del_reg) ) { // actual cut point is (del_reg.lower+del_reg.upper/2)
				if( ( in(cmp_algn.x.lower, del_reg) == true ) && ( in(cmp_algn.x.upper, del_reg) == true ) ) {
					fatalf("1: %d-%d has a conflict at the deletion cutpoint %d\n", cmp_algn.x.lower, cmp_algn.x.upper, (del_reg.lower + del_reg.upper)/2);
				}
				else if( in(cmp_algn.x.lower, del_reg) == true ) 
				{
					adjust_init_algn(init_dots, cmp_cur_id, del_reg, true, fp, DEL, is_x);
				}
				else if( in(cmp_algn.x.upper, del_reg) == true ) 
				{
					adjust_init_algn(init_dots, cmp_cur_id, del_reg, true, fp, DEL, is_x);
				}
				else if( subset(del_reg, cmp_algn.x) == true ) {
					if( (abs(cmp_algn.x.upper-del_reg.upper) <= DEL_TH) || (abs(del_reg.lower-cmp_algn.x.lower) <= DEL_TH) ) {
						adjust_init_algn(init_dots, cmp_cur_id, del_reg, true, fp, DEL, is_x);
					}
					else {
						fatalf("2: %d-%d has a conflict at the deletion cutpoint %d\n", cmp_algn.x.lower, cmp_algn.x.upper, (del_reg.lower + del_reg.upper)/2);
					}
				}
			}
		}
		else {
			if(f_loose_subset(cmp_algn.y, del_reg, LOOSE)) {
				init_dots[cmp_cur_id].sign = DELETED;
				if( (init_dots[cmp_cur_id].c_id != -1) && (init_dots[cmp_cur_id].m_id != -1) ) {
					init_dots[init_dots[cmp_cur_id].m_id].c_id = init_dots[cmp_cur_id].c_id;
					init_dots[init_dots[cmp_cur_id].c_id].m_id = init_dots[cmp_cur_id].m_id;
					init_dots[cmp_cur_id].m_id = -1;
					init_dots[cmp_cur_id].c_id = -1;
				}
			}
			else if( proper_overlap(cmp_algn.y, del_reg) ) {
				if( ( in(cmp_algn.y.lower, del_reg) == true ) && ( in(cmp_algn.y.upper, del_reg) == true ) ) {
					fatalf("3: %d-%d has conflict at the deletion cutpoint %d\n", cmp_algn.y.lower, cmp_algn.y.upper, (del_reg.lower + del_reg.upper)/2);
				}
				else if( in(cmp_algn.y.lower, del_reg) == true ) 
				{
					adjust_init_algn(init_dots, cmp_cur_id, del_reg, false, fp, DEL, is_x);
				}
				else if( in(cmp_algn.y.upper, del_reg) == true ) 
				{
					adjust_init_algn(init_dots, cmp_cur_id, del_reg, false, fp, DEL, is_x);
				}
				else if( subset(del_reg, cmp_algn.y) == true ) {
					if( (abs(cmp_algn.y.upper-del_reg.upper) <= DEL_TH) || (abs(del_reg.lower-cmp_algn.y.lower) <= DEL_TH) ) {
						adjust_init_algn(init_dots, cmp_cur_id, del_reg, false, fp, DEL, is_x);
					}
					else {
						fatalf("4: %d-%d has conflict at the deletion cutpoint %d (%d-%d)\n", cmp_algn.y.lower, cmp_algn.y.upper, (del_reg.lower + del_reg.upper)/2, del_reg.lower, del_reg.upper);
					}
				}
			}
		}
		i++;

		if( (count > 1) && (i < count)) {
			cmp_cur_id = cid_list[i];
			if( init_dots[cmp_cur_id].sign == DELETED ) {}
			cmp_algn.sign = init_dots[cmp_cur_id].sign;

			temp_b = init_dots[cmp_cur_id].x.lower + init_dots[cmp_cur_id].xl_diff;
			temp_e = init_dots[cmp_cur_id].x.upper - init_dots[cmp_cur_id].xr_diff;
			if( (temp_e - temp_b) <= 0 ) {
				init_dots[cmp_cur_id].sign = DELETED;
				fatalf("update_init_algns: empty interval 1 in (%d, %d)\n", temp_b, temp_e);
			}
			else cmp_algn.x = assign_I(temp_b, temp_e);
	
			temp_b = init_dots[cmp_cur_id].y.lower + init_dots[cmp_cur_id].yl_diff;
			temp_e = init_dots[cmp_cur_id].y.upper - init_dots[cmp_cur_id].yr_diff;

			if( (temp_e - temp_b) <= 0 ) {
				init_dots[cmp_cur_id].sign = DELETED;
				fatalf("update_init_algns: empty interval 2 in (%d, %d)\n", temp_b, temp_e);
			}
			else cmp_algn.y = assign_I(temp_b, temp_e);
		}
	}

	free(cid_list);
}

void write_init_maf_stdout(struct DotList *init_algns, int num_algns, struct n_pair *contigs1, struct n_pair *contigs2, int len1, int len2, FILE *fp, int flag, char *species, char *species2)
{
	char S2[BIG], T2[BIG];
	struct b_list *a_info;
	int i = 0, j = 0;
	int beg = 0, end = 0;
	int pid = 0;
	int y_b = 0, y_e = 0;
	int e1 = 0, e2 = 0;
	int gap1 = 0, gap2 = 0;
	int skip_nu = 0;
	int len = 0;
	int ctg_id1 = -1, ctg_id2 = -1;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;

  printf("##maf version=1 scoring=lastz-pid %s %s\n", species, species2);

	for( i = 0; i < num_algns; i++ ) 
	{
		if( (init_algns[i].sp_id == flag) && (init_algns[i].sign != DELETED)) 
		{
			if( init_algns[i].rp1_id != -1 ) {
				beg = init_algns[i].xl_diff + init_algns[i].xl_offset;
				skip_nu = find_xloc_one(init_algns[i], fp, init_algns[i].xl_diff + init_algns[i].xl_offset, NO_GAP_INC);
				end = find_end_xloc(init_algns[i], fp, init_algns[i].x.upper - init_algns[i].x.lower - init_algns[i].xr_diff - init_algns[i].xr_offset);
			}
			else {
				beg = init_algns[i].xl_diff;
				skip_nu = find_xloc_one(init_algns[i], fp, init_algns[i].xl_diff, NO_GAP_INC);
				end = find_end_xloc(init_algns[i], fp, init_algns[i].x.upper - init_algns[i].x.lower - init_algns[i].xr_diff);
			}
			(*a_info).b1 = init_algns[i].x.lower + init_algns[i].xl_diff - 1;
			y_b = init_algns[i].y.lower + init_algns[i].yl_diff;
			y_e = init_algns[i].y.upper - init_algns[i].yr_diff;
			ctg_id1 = init_algns[i].ctg_id1;
			ctg_id2 = init_algns[i].ctg_id2;

			get_nth_algn(S2, T2, init_algns[i].fid, beg, fp, a_info, REG);

			if( ctg_id1 == -1 ) (*a_info).len1 = len1;	
			else (*a_info).len1 = contigs1[ctg_id1].len;

      if( (*a_info).len1 != init_algns[i].len1 ) {
	      fatalf("Contig length not match: %s %d-%d\n", contigs1[ctg_id1].name2, len1, init_algns[i].len1);
      }

			if( ctg_id2 == -1 ) (*a_info).len2 = len2;	
			else (*a_info).len2 = contigs2[ctg_id2].len;

      if( (*a_info).len2 != init_algns[i].len2 ) {
	      fatalf("Contig length not match: %s %d-%d\n", contigs2[ctg_id2].name2, len2, init_algns[i].len2);
      }

			end = end - skip_nu;
			if( end <= 0 ) {
				end = 0;
			}
			S2[end] = '\0';
			T2[end] = '\0';

/*
			if( init_algns[i].sign == 0 ) {
				(*a_info).b2 = y_b-1;
				(*a_info).strand = '+';
			}
			else if( init_algns[i].sign == 1 ) {
				(*a_info).b2 = (a_info->len2) - y_e + 1;
				(*a_info).strand = '-';
			}
*/
			if( init_algns[i].sign == 1 ) {
				(*a_info).b2 = (a_info->len2) - (a_info->b2) + 2;
        (*a_info).strand = '-';
			}

		  e1 = 0;
		  e2 = 0;
		  gap1 = 0;
		  gap2 = 0;
		  for( j = 0; (j < end) && ((S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n')); j++ ) {
    		if( strchr("ACGTN", toupper(S2[j])) ) e1++;
				if( S2[j] == '-' ) gap1++;
    		if( strchr("ACGTN", toupper(T2[j])) ) e2++;
				if( T2[j] == '-' ) gap2++;
  		}
			S2[j] = '\0';
      T2[j] = '\0';

			if( (e1 >= MIN_GAP) && (e2 >= MIN_GAP) ) {
				(*a_info).e1 = e1;
				(*a_info).e2 = e2;
			
   	 	  pid = (int)(cal_pid_maf(S2, T2, end) + 0.5);

				if( (((a_info->b1)-1+(a_info->e1)) > (a_info->len1)) || ((a_info->b1-1) < 0 ) ) {
					fatalf("line 1252: Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
				}

     		printf("a srcblock=%d\n", init_algns[i].indiv_fid);
      	printf("s %s %d %d + %d %s", init_algns[i].name1, (a_info->b1)-1, a_info->e1, a_info->len1, S2);
      	len = strlen(S2);
      	if( S2[len-1] != '\n' ) printf("\n");

      	printf("s %s %d %d %c %d %s", init_algns[i].name2, (a_info->b2)-1, a_info->e2, a_info->strand, a_info->len2, T2);
      	len = strlen(T2);
      	if( T2[len-1] != '\n' ) printf("\n");
      	printf("\n");
			}
		}
	}

	free(a_info);
}

void update_pid_init_algns(int num_init_algns, struct DotList *init_algns, FILE *f, float avg_pid)
{
	int i = 0;
	int pid = 0;

	pid = (int)(avg_pid-0.5)-PID_DIFF;

	for( i = 0; i < num_init_algns; i++ ) {
		if( (init_algns[i].sign != DELETED) && (init_algns[i].sp_id == PAIR) ) {
 	  	if( ((init_algns[i].x.upper-init_algns[i].xr_diff) - (init_algns[i].x.lower+init_algns[i].xl_diff)) < ERR_SM_TH ) {
     		init_algns[i].sign = DELETED;
     	}
 	  	else if( ((init_algns[i].y.upper-init_algns[i].yr_diff) - (init_algns[i].y.lower+init_algns[i].yl_diff)) < ERR_SM_TH ) {
     		init_algns[i].sign = DELETED;
			}
   	  else {
   	    init_algns[i].identity = cal_pid_part_algn(init_algns, i, init_algns[i].xl_diff, init_algns[i].xr_diff, f, SELF1);
				if( init_algns[i].identity < pid ) init_algns[i].sign = DELETED;
				else {
					init_algns[i].m_x = assign_I(init_algns[i].x.lower+init_algns[i].xl_diff, init_algns[i].x.upper-init_algns[i].xr_diff);
					init_algns[i].m_y = assign_I(init_algns[i].y.lower+init_algns[i].yl_diff, init_algns[i].y.upper-init_algns[i].yr_diff);
				}
      }
    }
  }
}
