#include "main.h"
#include "regions.h"
#include "cmp_pid.h"
#include "get_algn.h"
#include "read_maf.h"
#include "util_i.h"

extern int debug_mode;

void update_inpar_list(int id, struct slist *alg_id, int num_algns, struct DotList *algns, struct DotList *init_algns, FILE *fp)
{
	int diff;
	int i = 0;
	int cur_id;

	diff = 0;
	cur_id = alg_id[id].id;
	while( (i < num_algns) && (diff >= 0) ) { 
		if( (i != cur_id) && (algns[i].sign != DELETED ) && (algns[i].pair_self == PAIR) ) 
		{
			if( is_cand_cmp(algns[cur_id].x, algns[i].x) == true) diff = cmp_pid(cur_id, true, i, true, algns, init_algns, fp);
			else if( is_cand_cmp(algns[cur_id].x, algns[i].y) == true ) diff = cmp_pid(cur_id, true, i, false, algns, init_algns, fp);
			else if( is_cand_cmp(algns[cur_id].y, algns[i].x) == true ) diff = cmp_pid(cur_id, false, i, true, algns, init_algns, fp);
			else if( is_cand_cmp(algns[cur_id].y, algns[i].y) == true ) diff = cmp_pid(cur_id, false, i, false, algns, init_algns, fp);
			else diff = 0;
		}
		else diff = 0;
		i++;
	}

	if( diff < 0 ) alg_id[cur_id].id = -1;
}

int cmp_pid(int cur_id, bool is_x_cur, int cmp_id, bool is_x_cmp, struct DotList *algns, struct DotList *init_algns, FILE *fp)
{
	char S1[BIG], T1[BIG];
	char S2[BIG], T2[BIG];
	int b1, b2, e1, e2;
	int b1_pos, b2_pos;
	int e1_pos, e2_pos;
	int res = 0;
	int b_diff = 0, e_diff = 0;
	int id1, id2;
	int i, count;
	int col_len1, col_len2;
	int num_nu1_S, num_nu1_T, num_nu2_S, num_nu2_T;
	int pid1, pid2;

	id1 = algns[cur_id].index;
	id2 = algns[cmp_id].index;

	get_nth_algn_ch(S1, T1, algns[cur_id].index, init_algns, fp);
	get_nth_algn_ch(S2, T2, algns[cmp_id].index, init_algns, fp);

	if( S1[strlen(S1)-1] == '\n' ) col_len1 = strlen(S1)-1;  
	else col_len1 = strlen(S1);  

	if( S2[strlen(S2)-1] == '\n' ) col_len2 = strlen(S2) - 1;  
	else col_len2 = strlen(S2);

	num_nu1_S = 0;
	num_nu1_T = 0;
	for( i = 0; i < col_len1; i++ ) {
		if(strchr("ACGTN", toupper(S1[i]))) num_nu1_S++;
		if(strchr("ACGTN", toupper(T1[i]))) num_nu1_T++;
	}

	num_nu2_S = 0;
	num_nu2_T = 0;
	for( i = 0; i < col_len2; i++ ) {
		if(strchr("ACGTN", toupper(S2[i]))) num_nu2_S++;
		if(strchr("ACGTN", toupper(T2[i]))) num_nu2_T++;
	}

	if( is_x_cur == true ) {
		b1 = init_algns[id1].x.lower + init_algns[id1].xl_diff;
		e1 = b1 + num_nu1_S;
	}
	else {
		if( init_algns[id1].sign == 0 ) {
			b1 = init_algns[id1].y.lower + init_algns[id1].yl_diff;
			e1 = b1 + num_nu1_T;
		}
		else {
			e1 = init_algns[id1].y.upper + init_algns[id1].yr_diff;
			b1 = e1 - num_nu1_T;
		}
	}

	if( is_x_cmp == true ) {
		b2 = init_algns[id2].x.lower + init_algns[id2].xl_diff;
		e2 = b2 + num_nu2_S;
	}
	else {
		if( init_algns[id2].sign == 0 ) {
			b2 = init_algns[id2].y.lower + init_algns[id2].yl_diff;
			e2 = b2 + num_nu2_T;
		}
		else {
			e2 = init_algns[id2].y.upper + init_algns[id2].yr_diff;
			b2 = e2 - num_nu2_T;
		}
	}

	b_diff = b1 - b2;
	e_diff = e1 - e2;

	b1_pos = 0;
	b2_pos = 0;
	e1_pos = col_len1;
	e2_pos = col_len2;

	count = 0;
	if( b_diff > 0 ) { // b1 is greater than b2, so b2 part should be cut off
		if( is_x_cur == true ) b1_pos = 0;	
		else {
			if( init_algns[id1].sign == 0 ) b1_pos = 0;
			else e1_pos = col_len1;
		}

		if( is_x_cmp == true ) {
			i = 0;
			while( (i < col_len2) && (count < abs(b_diff)) ) {
				if(strchr("ACGTN", toupper(S2[i]))) count++;
				i++;
			}
			b2_pos = i;
		}
		else {
			if( init_algns[id2].sign == 0 ) {
				i = 0;
				while( (i < col_len2) && (count < abs(b_diff)) ) {
					if(strchr("ACGTN", toupper(T2[i]))) count++;
					i++;
				}
				b2_pos = i;
			}
			else if( init_algns[id2].sign == 1 ) {
				i = col_len2;
				while( (i > 0 ) && (count < abs(b_diff)) ) {
					if(strchr("ACGTN", toupper(T2[i-1]))) count++;
					i--;
				}
				e2_pos = i;
			}
		}
	}
	else if( b_diff < 0 ) {
		if( is_x_cmp == true ) b2_pos = 0;	
		else {
			if( init_algns[id2].sign == 0 ) b2_pos = 0;
			else e2_pos = col_len2;
		}

		if( is_x_cur == true ) {
			i = 0;
			while( (i < col_len1) && (count < abs(b_diff)) ) {
				if(strchr("ACGTN", toupper(S1[i]))) count++;
				i++;
			}
			b1_pos = i;
		}
		else {
			if( init_algns[id1].sign == 0 ) {
				i = 0;
				while( (i < col_len1) && (count < abs(b_diff)) ) {
					if(strchr("ACGTN", toupper(T1[i]))) count++;
					i++;
				}
				b1_pos = i;
			}
			else if( init_algns[id1].sign == 1 ) {
				i = col_len1;
				while( (i > 0) && (count < abs(b_diff)) ) {
					if(strchr("ACGTN", toupper(T1[i-1]))) count++;
					i--;
				}
				e1_pos = i;
			}
		}
	}
	else {
		if( is_x_cmp == true ) b2_pos = 0;	
		else {
			if( init_algns[id2].sign == 0 ) b2_pos = 0;
			else e2_pos = col_len2;
		}

		if( is_x_cur == true ) b1_pos = 0;	
		else {
			if( init_algns[id1].sign == 0 ) b1_pos = 0;
			else e1_pos = col_len1;
		}
	}

	count = 0;
	if( e_diff > 0 ) { // e1 is greater than e2, so e1 is cut off
		if( is_x_cmp == true ) e2_pos = col_len2;	
		else {
			if( init_algns[id2].sign == 0 ) e2_pos = col_len2;
			else b2_pos = 0;
		}

		if( is_x_cur == true ) {
			i = col_len1;
			while( ( i > 0 ) && (count < abs(e_diff)) ) {
				if(strchr("ACGTN", toupper(S1[i-1]))) count++;
				i--;
			}
			e1_pos = i;
		}
		else {
			if( init_algns[id1].sign == 0 ) {
				i = col_len1;
				while( (i > 0 ) && (count < abs(e_diff)) ) {
					if(strchr("ACGTN", toupper(T1[i-1]))) count++;
					i--;
				}
				e1_pos = i;
			}
			else if( init_algns[id1].sign == 1 ) {
				i = 0;
				while( (i < col_len1) && (count < abs(e_diff)) ) {
					if(strchr("ACGTN", toupper(T1[i]))) count++;
					i++;
				}
				b1_pos = i;
			}
		}
	}
	else if( e_diff < 0 ) {
		if( is_x_cur == true ) e1_pos = col_len1;	
		else {
			if( init_algns[id1].sign == 0 ) e1_pos = col_len1;
			else b1_pos = 0;
		}

		if( is_x_cmp == true ) {
			i = col_len2;
			while( ( i > 0 ) && (count < abs(e_diff)) ) {
				if(strchr("ACGTN", toupper(S2[i-1]))) count++;
				i--;
			}
			e2_pos = i;
		}
		else {
			if( init_algns[id2].sign == 0 ) {
				i = col_len2;
				while( (i > 0) && (count < abs(e_diff)) ) {
					if(strchr("ACGTN", toupper(S2[i-1]))) count++;
					i--;
				}
				e2_pos = i;
			}
			else if( init_algns[id2].sign == 1 ) {
				i = 0;
				while( ( i < col_len2) && (count < abs(b_diff)) ) {
					if(strchr("ACGTN", toupper(T2[i]))) count++;
					i++;
				}
				b2_pos = i;
			}
		}
	}
	else {
		if( is_x_cmp == true ) e2_pos = col_len2;	
		else {
			if( init_algns[id2].sign == 0 ) e2_pos = col_len2;
			else b2_pos = 0;
		}

		if( is_x_cur == true ) e1_pos = col_len1;	
		else {
			if( init_algns[id1].sign == 0 ) e1_pos = col_len1;
			else b1_pos = 0;
		}
	}

	if( (abs(e1_pos - b1_pos) < ALT_EFFEC_VALUE) || (abs(e2_pos - b2_pos) < ALT_EFFEC_VALUE) ) res = 0;
	else {
		pid1 = cal_pid_maf_beg(S1, T1, b1_pos, e1_pos);
		pid2 = cal_pid_maf_beg(S2, T2, b2_pos, e2_pos);
		res = pid1-pid2;
	}

	return(res);
}

bool is_cand_cmp(struct I cur_reg, struct I cmp_reg) // check if two intervals are overlapped over OP_RATIO_TH (50%)
{
	int len;

	len = width(cur_reg);
	if( (proper_overlap(cur_reg, cmp_reg) == true) && (width(intersect(cur_reg, cmp_reg)) >= (int)(OP_RATIO_TH * ((float)len)) ) ) return true;
	else return false;
}
