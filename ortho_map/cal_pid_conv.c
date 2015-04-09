#include "main.h"
#include "util_gen.h"
#include "util_algns.h"
#include "util_ops.h"
#include "util.h"
#include "util_i.h"
#include "regions.h"
#include "sec_round.h"
#include "find_merging.h"

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

void cal_pid_conv(struct DotList *algns, int num_algns, struct cv_list *cv, int num_cv, FILE *f)
{
	int i = 0, j = 0, k = 0;
	struct DotList *temp_algns;
	struct DotList *candi_algns;
	int num_candi = 0;
	struct I src, dst;
	int sign = DELETED;
	float pid = (float)(-1);

	temp_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);
	candi_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);

	for( i = 0 ; i < num_cv; i++ )
	{
//		src = assign_I(cv[i].s1, cv[i].s2);
//		dst = assign_I(cv[i].t1, cv[i].t2);
		src = assign_I(cv[i].a1, cv[i].a2);
		dst = assign_I(cv[i].b1, cv[i].b2);

		initialize_algns(temp_algns, 0, num_algns);
		initialize_algns(candi_algns, 0, num_algns);
		pid = (float)(-1);
		k = 0;
		for( j = 0; j < num_algns; j++ ) {
			if( cv[i].ori == '+' ) {
				if( (algns[j].sign == 0) && (algns[j].sp_id == cv[i].sp_id) ) {
					assign_algn(temp_algns, k, algns[j]);
					k++;
				}
			}
			else if( cv[i].ori == '-' ) {
				if( (algns[j].sign == 1) && (algns[j].sp_id == cv[i].sp_id)) {
					assign_algn(temp_algns, k, algns[j]);
					k++;
				}
			}	
			else {
				fatalf("unsupported orientation sign: %c\n", cv[i].ori);
			}
		}

		if( k > 0 ) {
			num_candi = search_candi_algns(temp_algns, k, src, dst, candi_algns);			
		}
		else {
			if( debug_mode == TRUE ) {
				printf("no candi alignments\n");
			}
		}

		if( num_candi > 0 ) {
			sign = candi_algns[0].sign;
			if( debug_mode == TRUE ) {
				for( j = 0; j < num_candi; j++ ) {
					printf("%d-%d %d-%d: %d\n", candi_algns[j].x.lower, candi_algns[j].x.upper, candi_algns[j].y.lower, candi_algns[j].y.upper, candi_algns[j].identity); 
				}
			}
			pid = pick_sim_level(src, dst, candi_algns, num_candi, f);	
			cv[i].conv_pid = pid;
			if( debug_mode == TRUE ) printf("pid for %d-%d, %d-%d: %f\n", src.lower, src.upper, dst.lower, dst.upper, pid);
		}
		else {
			if( debug_mode == TRUE ) {
				printf("alignment not found for %d-%d, %d-%d\n", src.lower, src.upper, dst.lower, dst.upper);
			}
		}
	}
	
	quick_sort_dec_conv(cv, 0, num_cv-1, PID_BASE);

	if( debug_mode == TRUE ) {
		for( i = 0; i < num_cv; i++ ) {
			if( cv[i].dir == 0 ) printf("con %c %d %d %d %d ?", cv[i].ori, cv[i].a1, cv[i].a2, cv[i].b1, cv[i].b2);
			else if( cv[i].dir == 1 ) printf("con %c %d %d %d %d", cv[i].ori, cv[i].b1, cv[i].b2, cv[i].a1, cv[i].a2);
			else if( cv[i].dir == 2 ) printf("con %c %d %d %d %d", cv[i].ori, cv[i].a1, cv[i].a2, cv[i].b1, cv[i].b2);
			else {
				fatalf("unsupported direction sign: %d\n", cv[i].dir);
			}
			if( debug_mode == TRUE ) {
				printf(" %f\n", cv[i].conv_pid);
			}
			else printf("\n");
		}
	}	

	free(candi_algns);
	free(temp_algns);
}

int make_old_dup_list(struct DotList *algns, int num_algns, struct cv_list *init_dup_cv, int num_dup_cv, int *old_dups)
{
	int i = 0, j = 0, k = 0;
	struct DotList *temp_algns;
	struct DotList *candi_algns;
	struct DotList *pair_algns;
	int count = 0;
	int num_candi = 0;
	struct I src, dst;
	int num_list = 0;
	int min_id = -1, min_val = -1, temp_val = 0;

	src = assign_I(0, 1);
	dst = assign_I(0, 1);
	temp_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);
	candi_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);
	pair_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * num_algns);

	for( i = 0; i < num_algns; i++ ) {
		if( (algns[i].sp_id == PAIR) && (algns[i].sign != DELETED) ) {
			assign_algn(pair_algns, count, algns[i]);
			count++;
		}
	}

	for( i = 0; i < num_dup_cv; i++ )
	{
		src = assign_I(init_dup_cv[i].a1, init_dup_cv[i].a2);
		dst = assign_I(init_dup_cv[i].b1, init_dup_cv[i].b2);

		k = 0;
		min_id = -1;
		min_val = -1;
		num_candi = 0;
		temp_val = 0;
		for( j = 0; j < num_algns; j++ ) {
			if( init_dup_cv[i].ori == '+' ) {
				if( (algns[j].sign == 0) && (algns[j].sp_id == init_dup_cv[i].sp_id) ) {
					assign_algn(temp_algns, k, algns[j]);
					k++;
				}
			}
			else if( init_dup_cv[i].ori == '-' ) {
				if( (algns[j].sign == 1) && (algns[j].sp_id == init_dup_cv[i].sp_id)) {
					assign_algn(temp_algns, k, algns[j]);
					k++;
				}
			}	
			else {
				fatalf("unsupported orientation sign: %c\n", init_dup_cv[i].ori);
			}
		}

		if( k > 0 ) {
			num_candi = search_candi_algns(temp_algns, k, src, dst, candi_algns);			
		}
		else {
			if( debug_mode == TRUE ) {
				printf("no candi alignments\n");
			}
		}

		if( num_candi > 0 ) {
			for( j = 0; j < num_candi; j++ ) {
				if( debug_mode == TRUE ) {
					printf("%d-%d %d-%d: %d\n", candi_algns[j].x.lower, candi_algns[j].x.upper, candi_algns[j].y.lower, candi_algns[j].y.upper, candi_algns[j].identity); 
				}

				if( (width(src) >= (int)((0.7) * (float)width(candi_algns[j].x)) ) && (strict_subset(src, candi_algns[j].x) == true) && (strict_subset(dst, candi_algns[j].y) == true) ) {
					if( candi_algns[j].sign == 0 ) {
						temp_val = abs( abs(candi_algns[i].x.lower - src.lower) - abs(candi_algns[i].y.lower - dst.lower) );
					}
					else if( candi_algns[j].sign == 1 ) {
						temp_val = abs( abs(candi_algns[i].x.lower - src.lower) - abs(candi_algns[i].y.upper - dst.upper) );
					}

					if( min_id == -1 ) {
						min_id = j;
						min_val = temp_val;
					}
					else if( min_val > temp_val ) {
						min_id = j;
						min_val = temp_val;
					}
				}
			}
			
			if( min_id != -1 ) {
/*
				if( ( is_overlap_ortho_algns(candi_algns[min_id], true, pair_algns, count) == false) || ( is_overlap_ortho_algns(candi_algns[min_id], false, pair_algns, count) == false) ) {
					if( debug_mode == TRUE ) {
						printf("old dup conversion but a duplication happened after the speciation potentially: [%d,%d], [%d,%d]\n", candi_algns[min_id].x.lower, candi_algns[min_id].x.upper, candi_algns[min_id].y.lower, candi_algns[min_id].y.upper);	
					}
				}
				else {
					old_dups[num_list] = candi_algns[min_id].index;
					num_list++;
				}
*/
				if( is_tandem(candi_algns[min_id]) == false ) {
					old_dups[num_list] = candi_algns[min_id].index;
					num_list++;
				}
			}
		}
	}
	
	free(pair_algns);
	free(candi_algns);
	free(temp_algns);

	return(num_list);
}
