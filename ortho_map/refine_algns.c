#include "main.h"
#include "refine_algns.h"
#include "util_gen.h"
#include "util_i.h"
#include "util_input.h"
#include "util.h"
#include "id_ortho_conv.h"
#include "regions.h"
#include "update_init_algns.h"
#include "read_algn.h"
#include "read_maf.h"
#include "util_algns.h"

int adjust_all_algns_diff(struct DotList *ortho_algns, int num_ortho_algns, struct DotList *algns, int num_algns)
{
  int i = 0, j = 0;
  bool is_in = false;
  int cur_id = -1;
	int count = num_ortho_algns;

  for( i = 0; i < num_algns; i++ ) {
    j = 0;
    cur_id = -1;
    is_in = false;
    while( j < num_algns ) {      
			if( algns[i].indiv_fid == ortho_algns[j].indiv_fid ) {        
				is_in = true;
				ortho_algns[j].index = j;
        ortho_algns[j].fid = algns[i].fid;        
				ortho_algns[j].xl_diff = abs(ortho_algns[j].x.lower-algns[i].x.lower);        
				ortho_algns[j].xr_diff = abs(ortho_algns[j].x.upper-algns[i].x.upper);
        ortho_algns[j].yl_diff = abs(ortho_algns[j].y.lower-algns[i].y.lower);
        ortho_algns[j].yr_diff = abs(ortho_algns[j].y.upper-algns[i].y.upper);
        ortho_algns[j].x = assign_I(algns[i].x.lower, algns[i].x.upper); 
				ortho_algns[j].y = assign_I(algns[i].y.lower, algns[i].y.upper);
				ortho_algns[j].sp_id = PAIR;
      }
      j++;    
		}

    if( is_in == false ) {
			assign_algn(ortho_algns, count, algns[i]);
			ortho_algns[count].index = count;
			ortho_algns[count].sign = DELETED;
			ortho_algns[count].sp_id = PAIR;
			count++;
    }    
  }

	return(count);
}

void undo_events(struct I reg, struct DotList *algns, int num_algns, int mode, FILE *f)
{
	int i = 0;
	struct I cmp_reg;
	bool is_in = false;

	cmp_reg = assign_I(0, 1);

	for( i = 0; i < num_algns; i++ ) {
		if( mode == SELF1 ) {
			cmp_reg = assign_I(algns[i].x.lower+algns[i].xl_diff, algns[i].x.upper-algns[i].xr_diff);
		}
		else if( mode == SELF2 ) {
			cmp_reg = assign_I(algns[i].y.lower+algns[i].yl_diff, algns[i].y.upper-algns[i].yr_diff);
		}
		else {
			fatalf("mode has to be either SELF1 or SELF2: %d\n", mode);
		}

		is_in = false;
		if( algns[i].sign == DELETED ) {}
		else if( (strict_almost_equal(cmp_reg, reg) == true) || (strict_almost_equal(reg, cmp_reg) == true) || (f_loose_subset(cmp_reg, reg, LOOSE) == true)) {
			is_in = true;
		}
		else if( f_loose_subset(reg, cmp_reg, LOOSE) == true) {}
		else if( proper_overlap(cmp_reg, reg) == true ) {
			is_in = true;
		}

		if( is_in == true ) {
			if( mode == SELF1 ) {
				adjust_init_algn(algns, i, reg, true, f, DUP, true);
			}
			else {
				adjust_init_algn(algns, i, reg, false, f, DUP, false);
			}
		}
	}
}

void remove_false_mtom(struct DotList *algns, int num_algns, struct ops_list *ops, int num_ops, FILE *f)
{
	struct slist *sorted;
	int i = 0, j = 0, cur_id = 0, cmp_id = 0;
	int *b, *e;
	struct I src, dst;
	int loc = -1;
	struct I reg1, reg2;

	b = (int *) ckalloc(sizeof(int));
	e = (int *) ckalloc(sizeof(int));
	*b = -1;
	*e = -1;
	src = assign_I(0, 1);
	dst = assign_I(0, 1);
	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);
	sorted = (struct slist *) ckalloc(num_algns * sizeof(struct slist));
	initialize_slist(sorted, 0, num_algns);
	sort_init_algns(sorted, algns, num_algns, INIT_SELF2);
	for( i = 0; i < (num_algns-1); i++ ) {
		cur_id = sorted[i].id;
		if( algns[cur_id].sign == DELETED ) {}
		else {
			reg1 = assign_I(algns[cur_id].y.lower+algns[cur_id].yl_diff, algns[cur_id].y.upper-algns[cur_id].yr_diff);
			j = i+1;
			cmp_id = sorted[j].id;
			while( ( j < num_algns) && (algns[cmp_id].sign == DELETED) ) {
				j++;		
				if( j < num_algns ) cmp_id = sorted[j].id;
			}

			if( j < num_algns ) {
				reg2 = assign_I(algns[cmp_id].y.lower+algns[cmp_id].yl_diff, algns[cmp_id].y.upper-algns[cmp_id].yr_diff);
			}	
			while( ( j < num_algns ) && ( algns[cmp_id].sign != DELETED ) && (f_loose_overlap(reg1, reg2, STRICT) == true) ) {
				find_overlapping_ends(reg2, reg1, SELF2, algns, cur_id, f, b, e);
				src = assign_I(*b, *e);
				find_overlapping_ends(reg1, reg2, SELF2, algns, cmp_id, f, b, e);
				dst = assign_I(*b, *e);
				loc = -1;
				loc = where_in_ops(src, dst, ops, num_ops);

				if( loc == -1 ) {
					remove_algn_overlap(algns, num_algns, cur_id, cmp_id, src, dst, ops, num_ops, f);
					if( algns[cur_id].sign == DELETED ) {
						j = num_algns;
					}
					else if( ((algns[cur_id].y.upper-algns[cur_id].yr_diff)-(algns[cur_id].y.lower+algns[cur_id].yl_diff)) <= DEL_TH ) {
						j = num_algns;
						algns[cur_id].sign = DELETED;
					}
					else if( ((algns[cur_id].x.upper-algns[cur_id].xr_diff)-(algns[cur_id].x.lower+algns[cur_id].xl_diff)) <= DEL_TH ) {
						j = num_algns;
						algns[cur_id].sign = DELETED;
					}
					else {
						reg1 = assign_I(algns[cur_id].y.lower+algns[cur_id].yl_diff, algns[cur_id].y.upper-algns[cur_id].yr_diff);
					}
				}

				j++;
				if( j < num_algns ) cmp_id = sorted[j].id;
				while( ( j < num_algns) && (((algns[cur_id].x.upper-algns[cur_id].xr_diff)-(algns[cur_id].x.lower+algns[cur_id].xl_diff)) > DEL_TH) && (((algns[cur_id].y.upper-algns[cur_id].yr_diff)-(algns[cur_id].y.lower+algns[cur_id].yl_diff)) > DEL_TH) && (algns[cmp_id].sign == DELETED) ) {
					j++;		
					if( j < num_algns ) cmp_id = sorted[j].id;
				}
				if( j < num_algns ) {
					reg2 = assign_I(algns[cmp_id].y.lower+algns[cmp_id].yl_diff, algns[cmp_id].y.upper-algns[cmp_id].yr_diff);
				}
			}
		}
	}

	free(b);
	free(e);
	free(sorted);
}

int where_in_ops(struct I src, struct I dst, struct ops_list *ops, int num_ops)
{
	int i = 0;
	struct I reg1, reg2;
	int loc = -1;

	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);

	i = 0;
	while( (i < num_ops) && (loc == -1) ) {
		if( (ops[i].sign == 'c') || (ops[i].sign == 'v')) {}
		else {
			reg1 = assign_I(ops[i].srcStart, ops[i].srcEnd);
			reg2 = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( ((f_loose_subset(src, reg1, STRICT) == true) && (f_loose_subset(dst, reg2, STRICT) == true)) || (((f_loose_subset(src, reg2, STRICT) == true) && (f_loose_subset(dst, reg1, STRICT) == true) )) ) {
				loc = i;
			}
			else if( ((f_loose_subset(reg1, src, STRICT) == true) && (f_loose_subset(reg2, dst, STRICT) == true)) || ((f_loose_subset(reg2, src, STRICT) == true) && (f_loose_subset(reg1, dst, STRICT) == true) ) ) {
				loc = i;
			}
		}
		i++;
	}

	return(loc);
}

void remove_algn_overlap(struct DotList *algns, int num_algns, int cur_id, int cmp_id, struct I src, struct I dst, struct ops_list * ops, int num_ops, FILE *f)
{
	struct I ops_src;
	int i = 0;
	bool is_src = false;
	bool is_left = false; // Should the left alignment (cur_id) be remained?
	bool is_right = false;
	struct I ov; 
	struct I reg1, reg2;
	int loc = -1;
	bool *is_first;
	int ops_id1 = -1, ops_id2 = -1;

	is_first = (bool *) ckalloc(sizeof(bool));

	*is_first = true;
	ov = assign_I(0, 1);
	ops_src = assign_I(0, 1);
	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);
	reg1 = assign_I(algns[cur_id].y.lower+algns[cur_id].yl_diff, algns[cur_id].y.upper-algns[cur_id].yr_diff);
	reg2 = assign_I(algns[cmp_id].y.lower+algns[cmp_id].yl_diff, algns[cmp_id].y.upper-algns[cmp_id].yr_diff);
	i = 0;
	while( i < num_ops ) {
		ops_src = assign_I(ops[i].srcStart, ops[i].srcEnd);
		if( width(ops_src) < ERR_SM_TH ) {} 
		else if( (ops[i].sign == 'c') || (ops[i].sign == 'v') ) {}
		else if( (f_loose_subset(src, ops_src, STRICT) == true) || (f_loose_subset(ops_src, src, STRICT) == true)) {
			is_left = true;	
			is_src = true;
			ops_id1 = i;
		}
		else if( (f_loose_subset(dst, ops_src, STRICT) == true) || (f_loose_subset(ops_src, dst, STRICT) == true)) {
			is_right = true;
			is_src = true;
			ops_id2 = i;
		}
		i++;
	}

	if( (is_src == true) && (is_left == true) && (is_right == true ) ) {
		if( ops_id1 > ops_id2 ) {
			is_left = true;
			is_right = false;
		}
		else {
			is_left = false;
			is_right = true;
		}
	}
	else if( is_src == false ) {
		is_left = false;
		is_right = false;
		if(is_only_algn(algns, cur_id, src, num_algns) == true) {
			if( is_only_algn(algns, cmp_id, dst, num_algns) == true ) {
			}
			else {
				is_src = true;
				is_left = true;
			}
		}
		else if(is_only_algn(algns, cmp_id, dst, num_algns) == true ) {
			is_src = true;
			is_left = false;
		}
	}

	if( (f_loose_subset(reg2, reg1, STRICT) == true) || (f_loose_subset(reg1, reg2, STRICT) == true) ) {
		if( is_src == false ) {
			ov = intersect(reg1, reg2);
			if( (strict_almost_equal(reg1, reg2) == true ) || (strict_almost_equal(reg2, reg1) == true) ) {
				if( algns[cur_id].identity >= algns[cmp_id].identity ) {
					is_src = true;		
					is_left = true;
				}
				else {
					is_src = true;		
					is_left = false;
				}
			}
			else if( f_loose_subset(reg2, reg1, STRICT) == true) {
				is_src = true;
				is_left = true;
			}
			else {
				is_src = true;
				is_left = false;
			}
		}
	}
	else {
		ov = intersect(reg1, reg2);
		loc = -1;
		loc = find_loc_both_chained(algns, cur_id, cmp_id, ov, f, SELF2, is_first);
		if( loc != -1 ) {
			if(is_src == false) {
				is_src = true;
				is_left = *is_first;
			}
			else {
			}

			ov = assign_I(ov.lower, ov.lower+loc);
		}
		else {
			if( is_src == false ) {
				is_src = true;
				is_left = *is_first;
			}
		}
	}

	if( is_src == false ) {
		fatalf("Removing which alignment has not been determined yet, %d-%d, %d-%d\n, src.lower, src.upper, dst.lower, dst.upper");
	}
	else {
		if( is_left == true ) {
			if( f_loose_subset(reg2, reg1, STRICT) == true ) {
				algns[cmp_id].sign = DELETED;
			}
			else {
				if( loc != -1 ) {
					update_algn_diff(algns, cmp_id, ov.upper, f, YR_DIFF);
					update_algn_diff(algns, cur_id, ov.upper, f, YL_DIFF);
				}
				else if( algns[cmp_id].y.lower > algns[cur_id].y.lower ) {
					update_algn_diff(algns, cmp_id, ov.upper, f, YL_DIFF);
				}
				else {
					update_algn_diff(algns, cmp_id, ov.lower, f, YR_DIFF);
				}
			}
		}
		else {
			if( f_loose_subset(reg1, reg2, STRICT) == true ) {
				algns[cur_id].sign = DELETED;
			}
			else {
				if( loc != -1 ) {
					update_algn_diff(algns, cmp_id, ov.upper, f, YL_DIFF);
					update_algn_diff(algns, cur_id, ov.upper, f, YR_DIFF);
				}
				else if( algns[cur_id].y.lower > algns[cmp_id].y.lower ) {
					update_algn_diff(algns, cur_id, ov.upper, f, YL_DIFF);
				}
				else {
					update_algn_diff(algns, cur_id, ov.lower, f, YR_DIFF);
				}
			}
		}
	}

	if( ((algns[cur_id].x.upper - algns[cur_id].xr_diff) < (algns[cur_id].x.lower + algns[cur_id].xl_diff ))  || ( (algns[cur_id].y.upper - algns[cur_id].yr_diff) < (algns[cur_id].y.lower + algns[cur_id].yl_diff ) ))  {
		algns[cur_id].sign = DELETED;
	}

	if( ((algns[cmp_id].x.upper - algns[cmp_id].xr_diff) < (algns[cmp_id].x.lower + algns[cmp_id].xl_diff ))  || ( (algns[cmp_id].y.upper - algns[cmp_id].yr_diff) < (algns[cmp_id].y.lower + algns[cmp_id].yl_diff ) ))  {
		algns[cmp_id].sign = DELETED;
	}
	free(is_first);
}

bool is_only_algn(struct DotList *algns, int id, struct I reg, int num_algns)
{
	bool res = true;
	int i = 0;

	while( (i < num_algns) && (res == true) )
	{
		if( i == id ) {}
		else if( f_loose_subset(reg, algns[i].x, STRICT) == true ) {
			res = false;
		}
		i++;
	}

	return(res);
}

// if the first of the first alignment is taken, true
int find_loc_both_chained(struct DotList *algns, int id1, int id2, struct I ov, FILE *f, int sp_id, bool *is_first)
{ 
	int *beg1, *end1;
	int *beg2, *end2;
	struct b_list *a_info1, *a_info2;
	int res = -1;
	int count1 = 0, count2 = 0;
	char S3[BIG], T3[BIG];
	char S4[BIG], T4[BIG];

  a_info1 = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info2 = (struct b_list *) ckalloc(sizeof(struct b_list));
	beg1 = (int *) ckalloc(sizeof(int));
	beg2 = (int *) ckalloc(sizeof(int));
	end1 = (int *) ckalloc(sizeof(int));
	end2 = (int *) ckalloc(sizeof(int));

	cal_algn_beg_end(algns, id1, ov, sp_id, f, beg1, end1);
	cal_algn_beg_end(algns, id2, ov, sp_id, f, beg2, end2);

	strcpy(S3, "");
	strcpy(T3, "");
	strcpy(S4, "");
	strcpy(T4, "");

  get_nth_algn(S3, T3, algns[id1].fid, 0, f, a_info1, REG);
  get_nth_algn(S4, T4, algns[id2].fid, 0, f, a_info2, REG);

	*is_first = true;
	if( sp_id == SELF1 ) {
		res = find_loc_both_chained_sp1(S3, T3, S4, T4, *beg1, *beg2, *end1, *end2, is_first);
	}
	else if( sp_id == SELF2 ) {
		res = find_loc_both_chained_sp2(S3, T3, S4, T4, *beg1, *beg2, *end1, *end2, algns[id1].init_sign, algns[id2].init_sign, is_first);
	}

	if( res == -1 ) {
		count1 = cal_match_maf_beg(S3, T3, *beg1, *end1);						
		count2 = cal_match_maf_beg(S4, T4, *beg2, *end2);						
		if( count1 >= count2 )  {
			*is_first = true;
		}
		else {
			*is_first = false;
		}
	}

	free(a_info1);
	free(a_info2);
	free(beg1);
	free(end1);
	free(beg2);
	free(end2);
	return(res);
}

int find_loc_both_chained_sp1(char *seq1, char *t_seq1, char *seq2, char *t_seq2, int b1, int b2, int e1, int e2, bool *is_first)
{
	int i = 0, j = 0;
	int cur_N1 = 0, cur_N2 = 0;
	int nu = 0;
	int res = -1;
	int beg_count1 = 0, beg_count2 = 0;
	int end_count1 = 0, end_count2 = 0;
	int temp_count = 0;

	i = b1;
	j = e1;

	while( (res == -1) && (i < e1) && (j < e2) )
	{
		while( (!strchr("ACGTN", toupper(t_seq1[i]))) && (i < e1) ) {
			if( seq1[i] == 'N' ) cur_N1 = nu;	
			i++;
		}

		while( (!strchr("ACGTN", toupper(t_seq2[j]))) && (j < e2) ) {
			if( seq2[j] == 'N' ) cur_N2 = nu;	
			j++;
		}

		if( (i >= e1) || (j >= e2)) {
		}
		else if( strchr("ACGTN", toupper(t_seq1[i])) && (strchr("ACGTN", toupper(t_seq2[j])) ) ) 
		{
			if( t_seq1[i] == 'N' ) cur_N1 = nu;	
			if( t_seq2[j] == 'N' ) cur_N2 = nu;	
			nu++;
		}
		else {
			fatalf("no gaps expected %c and %c\n", t_seq1[i], t_seq2[j]);
		}

		temp_count = nu;
		if( (cur_N1 > ERR_TH) && (cur_N2 > ERR_TH) && (abs(cur_N1 - cur_N2) < DEL_TH) ) 
		{
			while((i < e1) || ( (seq1[i] == 'N') || (t_seq1[i] == 'N') ) ) {
				if(strchr("ACGTN", toupper(t_seq1[i]))) temp_count++;
				i++;
			}

			while((j < e2) || ( (seq2[j] == 'N') || (t_seq2[j] == 'N') ) ) {
				if(strchr("ACGTN", toupper(t_seq2[j]))) nu++;
				j++;
			}

			if( nu < temp_count ) {
				while((j < e2) && ( nu < temp_count) ) {
					if(strchr("ACGTN", toupper(t_seq2[j]))) nu++;
					j++;
				}
			}
			else if( temp_count < nu ) {
				while( (i < e1) && (temp_count < nu) ) {
					if(strchr("ACGTN", toupper(t_seq1[i]))) temp_count++;
					i++;
				}
			}

			beg_count1 = cal_match_maf_beg(seq1, t_seq1, b1, i);						
			end_count1 = cal_match_maf_beg(seq1, t_seq1, i, e1);						
			beg_count2 = cal_match_maf_beg(seq2, t_seq2, b2, j);						
			end_count2 = cal_match_maf_beg(seq2, t_seq2, j, e2);						

			if( (beg_count1 > beg_count2) && (end_count1 < end_count2) ) {
				if( (abs(beg_count1 - beg_count2) > ERR_SM_TH ) || (abs(end_count1-end_count2) > ERR_SM_TH ) ) {
					res = nu;
					*is_first = true;
				}
			}
			else if( (beg_count1 < beg_count2) && (end_count1 > end_count2) ) {
				if( (abs(beg_count1 - beg_count2) > ERR_SM_TH ) || (abs(end_count1-end_count2) > ERR_SM_TH ) ) {
					res = nu;
					*is_first = false;
				}
			}
		}

		i++;
		j++;
	}

	return(res);

}

int find_loc_both_chained_sp2(char *seq1, char *t_seq1, char *seq2, char *t_seq2, int b1, int b2, int e1, int e2, int sign1, int sign2, bool *is_first)
{
	int i = 0, j = 0;
	int cur_N1 = 0, cur_N2 = 0;
	int nu = 0;
	int res = -1;
	int beg_count1 = 0, beg_count2 = 0;
	int end_count1 = 0, end_count2 = 0;
	int temp_count = 0;

	if( sign1 == 0 ) {
		i = b1;
	}
	else {
		i = e1-1;
	}
	
	if( sign2 == 0 ) {
		j = b2;
	}
	else {
		j = e2-1;
	}

	while( (res == -1) && (((sign1 == 0) && (i < e1)) || ((sign1 == 1) && (i >= b1))) && ( ((sign2 == 0) && (j < e2)) || ((sign2 == 1) && (j >= b2))) )
	{
		while( (!strchr("ACGTN", toupper(t_seq1[i]))) && ((((sign1 == 0) && (i < e1)) || ((sign1 == 1) && (i >= b1)))) ) {
			if( seq1[i] == 'N' ) cur_N1 = nu;	
			if( sign1 == 0 ) i++;
			else i--;
		}

		while( (!strchr("ACGTN", toupper(t_seq2[j]))) && ((((sign2 == 0) && (j < e2)) || ((sign2 == 1) && (j >= b2)))) ) {
			if( seq2[j] == 'N' ) cur_N2 = nu;	
			if( sign2 == 0 ) j++;
			else j--;
		}

		if( ((sign1 == 0) && (i >= e1)) || ((sign1 == 1) && (i < b1)) || ((sign2 == 0) && (j >= e2)) || ((sign2 == 1) && (j < b2)) ) {
		}
		else if( strchr("ACGTN", toupper(t_seq1[i])) && strchr("ACGTN", toupper(t_seq2[j])) ) 
		{
			if( t_seq1[i] == 'N' ) cur_N1 = nu;	
			if( t_seq2[j] == 'N' ) cur_N2 = nu;	
			nu++;
		}
		else {
			fatalf("no gaps expected %c and %c\n", t_seq1[i], t_seq2[j]);
		}

		temp_count = nu;
		if( (cur_N1 > ERR_TH) && (cur_N2 > ERR_TH) && (abs(cur_N1 - cur_N2) < DEL_TH) ) 
		{
			while((((sign1 == 0) && (i < e1)) || ((sign1 == 1) && (i >= b1))) && ( (seq1[i] == 'N') || (t_seq1[i] == 'N') ) ) {
				if(strchr("ACGTN", toupper(t_seq1[i]))) temp_count++;
				if( sign1 == 0 ) i++;
				else i--;
			}

			while((((sign2 == 0) && (j < e2)) || ((sign2 == 1) && (j >= b2))) && ( (seq2[j] == 'N') || (t_seq2[j] == 'N') ) ) {
				if(strchr("ACGTN", toupper(t_seq2[j]))) nu++;
				if( sign2 == 0 ) j++;
				else j--;
			}

			if( nu < temp_count ) {
				while((((sign2 == 0) && (j < e2)) || ((sign2 == 1) && (j >= b2))) && ( nu < temp_count) ) {
					if(strchr("ACGTN", toupper(t_seq2[j]))) nu++;
					if( sign2 == 0 ) j++;
					else j--;
				}
			}
			else if( temp_count < nu ) {
				while((((sign1 == 0) && (i < e1)) || ((sign1 == 1) && (i >= b1))) && (temp_count < nu  ) ) {
					if(strchr("ACGTN", toupper(t_seq1[i]))) temp_count++;
					if( sign1 == 0 ) i++;
					else i--;
				}
			}

			if( sign1 == 0 ) {
				beg_count1 = cal_match_maf_beg(seq1, t_seq1, b1, i);						
				end_count1 = cal_match_maf_beg(seq1, t_seq1, i, e1);						
			}
			else {
				beg_count1 = cal_match_maf_beg(seq1, t_seq1, i, e1);						
				end_count1 = cal_match_maf_beg(seq1, t_seq1, b1, i);						
			}

			if( sign2 == 0 ) {
				beg_count2 = cal_match_maf_beg(seq2, t_seq2, b2, j);						
				end_count2 = cal_match_maf_beg(seq2, t_seq2, j, e2);						
			}
			else {
				beg_count2 = cal_match_maf_beg(seq2, t_seq2, j, e2);						
				end_count2 = cal_match_maf_beg(seq2, t_seq2, b2, j);						
			}

			if( (beg_count1 > beg_count2) && (end_count1 < end_count2) ) {
				if( (abs(beg_count1 - beg_count2) > ERR_SM_TH ) || (abs(end_count1-end_count2) > ERR_SM_TH ) ) {
					res = nu;
					*is_first = true;
				}
			}
			else if( (beg_count1 < beg_count2) && (end_count1 > end_count2) ) {
				if( (abs(beg_count1 - beg_count2) > ERR_SM_TH ) || (abs(end_count1-end_count2) > ERR_SM_TH ) ) {
					res = nu;
					*is_first = false;
				}
			}
		}

		if( sign1 == 0 ) i++;
		else i--;

		if( sign2 == 0 ) j++;
		else j--;
	}

	return(res);
}
