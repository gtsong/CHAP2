#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_input.h"
#include "util_gen.h"
#include "util_algns.h"
#include "util_i.h"
#include "id_ortho_conv.h"
#include "regions.h"
#include "read_algn.h"
#include "update_init_algns.h"
#include "contigs_op.h"

#define ORTHO_CONTENT 1
#define ORTHO_POSITION 2
#define SIMPLE_CUT 3
#define SYM_PID_TH 3

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];
void make_symmetry(struct DotList **ortho_algns, int *num_ortho, struct DotList **symm_ortho_algns, int *num_symm_ortho, struct DotList *symm_algns, int *num_symm, int ortho_mode, FILE *f, FILE *g, float avg_pid);
int check_algn_in(struct DotList *ortho_algn, int id, struct slist *st, int s_loc, int e_loc, struct DotList *symm_algns, FILE *g, FILE *f, struct short_alist *list);
void adjust_algn_diff(struct DotList *ortho_algns, int num_ortho_algns, struct DotList *algns, int num_algns);
int toggle_list(struct DotList *ortho_algns, int id, struct short_alist *list, int num_list);
int cut_algns(struct DotList **ortho_algns, int id, int num1, struct short_alist *list, int num_list, FILE *f);
int add_algns(struct short_alist *add_list, int num_add, struct DotList *algns, int num_algns, FILE *f, struct DotList *added_algns, int num_added, int sign);
bool is_pid_higher(struct short_alist cur_list, struct DotList *ortho_algns, int num_ortho, FILE *f, float avg_pid, float cur_pid, int id); 
bool is_interval_bigger(struct short_alist cur_list, struct DotList *ortho_algns, int num_ortho, FILE *f, float avg_pid, float cur_pid, int id);
void chain_back(struct DotList *ortho_algns, int num_ortho);

int main(int argc, char **argv)
{
	struct DotList *ortho_algns;
	struct DotList *symm_ortho_algns;
	struct DotList *algns, *symm_algns; // alignments from chaining
	int *num_ortho, *num_symm_ortho, *num_algns, *num_symm;
	int *size1, *size2;
	int count = 0;
	int ortho_mode = ORTHO_CONTENT; 
	char species[100], species2[100];
	int i = 0;
	FILE *f, *g, *h, *fp;
	float avg_pid = (float)0; 
	struct n_pair *contigs1, *contigs2;
	int *num_contigs1, *num_contigs2;
	int *len_sum1, *len_sum2;
	int *num_alloc1, *num_alloc2; 
	char name1[LEN_NAME], name2[LEN_NAME];
	int len1 = 0, len2 = 0;
	int num_org1 = 0, num_org2 = 0;
	char buf[1000];

	debug_mode = FALSE;
	ortho_mode = ORTHO_CONTENT;
	if( (argc != 7) && (argc != 9) ) {
		fatal("args: ortho_maf_file symmetric_ortho_maf_file initial_maf_file symmetric_initial_maf_file output_maf_file symmetric_output_maf_file\n");
	}

	num_contigs1 = (int *) ckalloc(sizeof(int));
	num_contigs2 = (int *) ckalloc(sizeof(int));
	num_alloc1 = (int *) ckalloc(sizeof(int));
	num_alloc2 = (int *) ckalloc(sizeof(int));

	if( argc == 9 ) {
		if( (fp = ckopen(argv[7], "r")) == NULL ) {
			fatalf("file %s not exist\n", argv[7]);
		}	  
		else {
			count = count_lines(fp);
			contigs1 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * (count+1));
			len_sum1 = (int *) ckalloc(sizeof(int) * (count+1));
			init_n_pair(contigs1, 0, count);
			for( i = 0; i < (count+1); i++ ) len_sum1[i] = 0;
			i = 0;
			fseek(fp, 0, SEEK_SET);
			while( fgets(buf, 1000, fp) ) {
				if( i >= count ) {
					fatalf("counting error in %s\n", argv[7]);
				}
				sscanf(buf, "%s %s %d %d", name1, name2, &len1, &len2);	
				strcpy(contigs1[i].name1, name1);
				strcpy(contigs1[i].name2, name2);
				contigs1[i].len = len1;
				contigs1[i].id = i;
				len_sum1[i] = len2;
				i++;
			}
			*num_contigs1 = count;
			num_org1 = count;
			*num_alloc1 = count+1;
			fclose(fp);
		}

		if( (fp = ckopen(argv[8], "r")) == NULL ) {
			fatalf("file %s not exist\n", argv[8]);
		} 
		else {
			count = count_lines(fp);
			contigs2 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * (count+1));
			len_sum2 = (int *) ckalloc(sizeof(int) * (count+1));
			init_n_pair(contigs2, 0, count);
			for( i = 0; i < (count+1); i++ ) len_sum2[i] = 0;
			i = 0;
			fseek(fp, 0, SEEK_SET);
			while( fgets(buf, 1000, fp) ) {
				if( i >= count ) {
					fatalf("counting error in %s \n", argv[8]);
				}
				sscanf(buf, "%s %s %d %d", name1, name2, &len1, &len2);	
				strcpy(contigs2[i].name1, name1);
				strcpy(contigs2[i].name2, name2);
				contigs2[i].len = len1;
				contigs2[i].id = i;
				len_sum2[i] = len2;
				i++;
			}
			*num_contigs2 = count;
			*num_alloc2 = count+1;
			num_org2 = count;
			fclose(fp);
		}
	}
	else {
		count = ALLOC_UNIT;
		contigs1 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * count);
		contigs2 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * count);
		len_sum1 = (int *) ckalloc(sizeof(int) * count);
		len_sum2 = (int *) ckalloc(sizeof(int) * count);

		init_n_pair(contigs1, 0, count-1);
		init_n_pair(contigs2, 0, count-1);

		for( i = 0; i < count; i++ ) {
			len_sum1[i] = 0;
			len_sum2[i] = 0;
		}
		*num_contigs1 = 0;
		*num_contigs2 = 0;
		*num_alloc1 = count;
		*num_alloc2 = count;
		num_org1 = 0;
		num_org2 = 0;
	}

	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_ortho = (int *) ckalloc(sizeof(int));
	num_symm_ortho = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	num_symm = (int *) ckalloc(sizeof(int));

	*size1 = 0;
	*size2 = 0;
	*num_ortho = 0;
	*num_symm_ortho = 0;
	*num_algns = 0;
	*num_symm = 0;
	count = count_local_algns(argv[1], species, species2);
	if( count > 0 ) {
		ortho_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(ortho_algns, 0, count);

		if( ortho_mode == SIMPLE_CUT ) {
		  read_maf(argv[1], G_MODE, ortho_algns, num_ortho, size1, size2);
		}
		else if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
		  read_maf(argv[1], O_MODE, ortho_algns, num_ortho, size1, size2);
		}
	}

	if( count != (*num_ortho) ) {
		fatalf("counting error 1: %d vs. %d\n", count, *num_ortho);
	}

	count = count_local_algns(argv[2], species, species2);
	if( count > 0 ) {
		symm_ortho_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(symm_ortho_algns, 0, count);

		if( ortho_mode == SIMPLE_CUT) {
	 		read_maf(argv[2], G_MODE, symm_ortho_algns, num_symm_ortho, size1, size2);
		}
		else if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
	 		read_maf(argv[2], O_MODE, symm_ortho_algns, num_symm_ortho, size1, size2);
		}
	}

	if( count != (*num_symm_ortho) ) {
		fatalf("counting error 2: %d vs. %d\n", count, *num_symm_ortho);
	}

	count = count_local_algns(argv[3], species, species2);
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(algns, 0, count);
	  read_maf(argv[3], G_MODE, algns, num_algns, size1, size2);
		
/*
		if( (*num_algns) > 0 ) {
			contigs1 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * (*num_algns));
			contigs2 = (struct n_pair *) ckalloc(sizeof(struct n_pair) * (*num_algns));
			init_n_pair(contigs1, 0, (*num_algns)-1);
			init_n_pair(contigs2, 0, (*num_algns)-1);
		}
		*num_contigs1 = 0;
		*num_contigs2 = 0;
		*num_alloc1 = *num_algns;
		*num_alloc2 = *num_algns;
			
		adjust_multi_contig_pos(algns, *num_algns, size1, size2, contigs1, num_contigs1, contigs2, num_contigs2);

		if( (*num_contigs1) > 0 ) len_sum1 = (int *) ckalloc(sizeof(int) * (*num_contigs1));
		else len_sum1 = (int *) ckalloc(sizeof(int));

		if( (*num_contigs2) > 0 ) len_sum2 = (int *) ckalloc(sizeof(int) * (*num_contigs2));
		else len_sum2 = (int *) ckalloc(sizeof(int));

		cal_length_sum(len_sum1, contigs1, *num_contigs1);
		cal_length_sum(len_sum2, contigs2, *num_contigs2);
*/
		merge_sep_contigs(algns, *num_algns, &contigs1, num_contigs1, num_alloc1, &len_sum1, &contigs2, num_contigs2, num_alloc2, &len_sum2);

		if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
			adjust_algn_diff(ortho_algns, *num_ortho, algns, *num_algns);
		}
	}

	if( count != (*num_algns) ) {
		fatalf("counting error 3: %d vs. %d\n", count, *num_algns);
	}

	merge_sep_contigs(ortho_algns, *num_ortho, &contigs1, num_contigs1, num_alloc1, &len_sum1, &contigs2, num_contigs2, num_alloc2, &len_sum2);
	merge_sep_contigs(symm_ortho_algns, *num_symm_ortho, &contigs2, num_contigs2, num_alloc2, &len_sum2, &contigs1, num_contigs1, num_alloc1, &len_sum1);
	count = count_local_algns(argv[4], species, species2);
	if( count > 0 ) {
		symm_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(symm_algns, 0, count);
	  read_maf(argv[4], G_MODE, symm_algns, num_symm, size1, size2);
		merge_sep_contigs(symm_algns, *num_symm, &contigs2, num_contigs2, num_alloc2, &len_sum2, &contigs1, num_contigs1, num_alloc1, &len_sum1);

		if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION)) {
			adjust_algn_diff(symm_ortho_algns, *num_symm_ortho, symm_algns, *num_symm);
		}
	}

	if( count != (*num_symm) ) {
		fatalf("counting error 4: %d vs. %d\n", count, *num_symm);
	}

	if( (*num_ortho) > 0 ) {
		avg_pid = cal_avg_pid(ortho_algns, *num_ortho);
	}
	else if( (*num_symm_ortho) > 0 ) {
		avg_pid = cal_avg_pid(symm_ortho_algns, *num_symm_ortho);
	}

	if( ortho_mode == SIMPLE_CUT ) {
		f = ckopen(argv[1], "r");
		g = ckopen(argv[2], "r");
	}
	else if( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
		f = ckopen(argv[3], "r");
		g = ckopen(argv[4], "r");
	}

	make_symmetry(&ortho_algns, num_ortho, &symm_ortho_algns, num_symm_ortho, symm_algns, num_symm, ortho_mode, f, g, avg_pid); // if ortho_mode == ORTHO_CONTENT, this function adds some alignments to symm_ortho_algns and cuts other alignments from ortho_algns
	make_symmetry(&symm_ortho_algns, num_symm_ortho, &ortho_algns, num_ortho, algns, num_algns, ortho_mode, g, f, avg_pid);

	for( i = 0; i < (*num_ortho); i++ ) {
		ortho_algns[i].sp_id = PAIR;
	}

	for( i = 0; i < (*num_symm_ortho); i++ ) {
		symm_ortho_algns[i].sp_id = PAIR;
	}	

	chain_back(ortho_algns, *num_ortho);
	chain_back(symm_ortho_algns, *num_symm_ortho);

	if( argc == 9 ) {

		fp = ckopen(argv[7], "a+");
		print_contig_list(fp, contigs1, len_sum1, num_org1, *num_contigs1);
		fclose(fp);
		fp = ckopen(argv[8], "a+");
		print_contig_list(fp, contigs2, len_sum2, num_org2, *num_contigs2);
		fclose(fp);
	}

//	replace_contigs_len(contigs1, *num_contigs1);
//	replace_contigs_len(contigs2, *num_contigs2);

	h = ckopen(argv[5], "w");
	write_init_maf(h, ortho_algns, *num_ortho, contigs1, contigs2, *num_contigs1, *num_contigs2, f, PAIR);
	fclose(h);
	h = ckopen(argv[6], "w");
	write_init_maf(h, symm_ortho_algns, *num_symm_ortho, contigs2, contigs1, *num_contigs2, *num_contigs1, g, PAIR);
	fclose(h);

	fclose(g);
	fclose(f);

	free(size1);
	free(size2);

	if( (*num_alloc1) > 0 ) {
		free(contigs1);
		free(len_sum1);
	}

	if( (*num_alloc2) > 0 ) {
		free(contigs2);
		free(len_sum2);
	}

	if( (*num_algns) > 0 ) {
		free(algns);
	}

	if( (*num_symm_ortho) > 0 ) {
		free(symm_ortho_algns);
	}

	if( (*num_symm) > 0 ) {
		free(symm_algns);
	}

	if( (*num_ortho) > 0 ) {
		free(ortho_algns);
	}

	free(num_alloc1);
	free(num_alloc2);
	free(num_contigs1);
	free(num_contigs2);
	free(num_ortho);
	free(num_symm_ortho);
	free(num_symm);
	free(num_algns);
	return(EXIT_SUCCESS);
}

void make_symmetry(struct DotList **ortho_algns, int *num_ortho, struct DotList **symm_ortho_algns, int *num_symm_ortho, struct DotList *symm_algns, int *num_symm, int ortho_mode, FILE *f, FILE *g, float avg_pid)
{
	int num1 = 0, num2 = 0, num3 = 0;
	struct slist *st1, *st2;
	int i = 0, j = 0, k = 0;
	struct I query;
	int s_loc = 0, e_loc = 0;
	struct short_alist *list, *add_list, *cut_list;
	struct short_alist temp_list;
	int num_list = 0, num_add = 0, num_cut = 0;
	int org_num1 = 0;
	float pid = (float)0;
	struct DotList *added_algns;
	int num_added = 0;

	num1 = *num_ortho;
	num2 = *num_symm_ortho;
	num3 = *num_symm;
	org_num1 = num1;
	
	st1 = (struct slist *) ckalloc(num1 * sizeof(struct slist));
	st2 = (struct slist *) ckalloc(num2 * sizeof(struct slist));

	list = (struct short_alist *) ckalloc(num2 * sizeof(struct short_alist));
	add_list = (struct short_alist *) ckalloc(num2 * sizeof(struct short_alist));
	cut_list = (struct short_alist *) ckalloc(num2 * sizeof(struct short_alist));
	added_algns = (struct DotList *) ckalloc(num1 * sizeof(struct DotList));
	
	temp_list.id = -1;
	temp_list.x = assign_I(0,1);
	temp_list.y = assign_I(0,1);
	temp_list.val = 0;

	initialize_alist(list, 0, num2);
	initialize_alist(add_list, 0, num2);
	initialize_alist(cut_list, 0, num2);
	initialize_algns(added_algns, 0, num1);

	initialize_slist(st1, 0, num1);
	initialize_slist(st2, 0, num2);

	sort_init_algns(st1, *ortho_algns, num1, SELF1);
	sort_init_algns(st2, *symm_ortho_algns, num2, SELF2);

	for( k = 0; k < org_num1; k++) {
		i = st1[k].id;	
		if( (*ortho_algns)[i].sign == DELETED ) {

		}
		else {
			query = assign_I((*ortho_algns)[i].x.lower + (*ortho_algns)[i].xl_diff, (*ortho_algns)[i].x.upper - (*ortho_algns)[i].xr_diff);
		  s_loc = search_range_b(st2, *symm_ortho_algns, num2, query.lower, SELF2);
		  e_loc = search_range_e(st2, *symm_ortho_algns, num2, query.upper, SELF2);
			initialize_alist(list, 0, num2);

			num_list = 0;
			num_list = check_algn_in(*ortho_algns, i, st2, s_loc, e_loc, *symm_ortho_algns, g, f, list);

			if( ortho_mode == SIMPLE_CUT ) {
				if( num_list <= 0 ) {
					(*ortho_algns[i]).sign = DELETED;
				}
				else if( num_list == 1 ) {
					if( strict_subset((*ortho_algns)[i].x, list[0].x) == true ) {
						if( strict_subset((*ortho_algns)[i].y, list[0].y) == false ) {
							printf("warning: [%d,%d] [%d,%d] found incorrectly for [%d,%d] [%d,%d]\n", list[0].x.lower, list[0].x.upper, list[0].y.lower, list[0].y.upper, (*ortho_algns)[i].x.lower, (*ortho_algns)[i].x.upper, (*ortho_algns)[i].y.lower, (*ortho_algns)[i].y.upper);
						}
					}
					else {
						update_algn_diff(*ortho_algns, i, list[0].y.lower, f, YL_DIFF);
						update_algn_diff(*ortho_algns, i, list[0].y.upper, f, YR_DIFF);
					}
				}
				else { // num_list > 1
					num1 = cut_algns(ortho_algns, i, num1, list, num_list, f);
				}
			}
			else if ( (ortho_mode == ORTHO_CONTENT) || (ortho_mode == ORTHO_POSITION) ) {
				initialize_alist(add_list, 0, num2);
				initialize_alist(cut_list, 0, num2);
				num_add = 0;
				num_cut = 0;
				num_list = toggle_list(*ortho_algns, i, list, num_list);
				for( j = 0; j < num_list; j++ ) {
					pid = cal_pid_part_algn(*ortho_algns, i, abs(list[j].y.lower-(*ortho_algns)[i].y.lower), abs(list[j].y.upper-(*ortho_algns)[i].y.upper), f, SELF2);	
					temp_list = assign_alist(list[j]);
					temp_list.x = assign_I(list[j].y.lower, list[j].y.upper);
					temp_list.y = assign_I(list[j].x.lower, list[j].x.upper);

					if( ortho_mode == ORTHO_CONTENT ) {
						if( (is_pid_higher(list[j], *ortho_algns, num1, f, avg_pid, pid, i) == true) && (is_pid_higher(temp_list, *symm_ortho_algns, num2, g, avg_pid, pid, -1) == true) )
						{
							add_list[num_add] = assign_alist(list[j]);
							num_add++;
						}
						else {
							cut_list[num_cut] = assign_alist(list[j]);
							num_cut++;
						}
					}
					else if( ortho_mode == ORTHO_POSITION ) {
						if( (is_interval_bigger(list[j], *ortho_algns, num1, f, avg_pid, pid, i) == true) && (is_interval_bigger(temp_list, *symm_ortho_algns, num2, g, avg_pid, pid, -1) == true) )
						{
							add_list[num_add] = assign_alist(list[j]);
							num_add++;
						}
						else {
							cut_list[num_cut] = assign_alist(list[j]);
							num_cut++;
						}
					}
				}

				if( num_cut > 0 ) {
					num_cut = toggle_list(*ortho_algns, i, cut_list, num_cut);
					if( num_cut == 0 ) {
						(*ortho_algns)[i].sign = DELETED;
					}
					else {
						num1 = cut_algns(ortho_algns, i, num1, cut_list, num_cut, f);
					}
				}

				if( num_add > 0 ) {
					num_added = add_algns(add_list, num_add, symm_algns, num3, g, added_algns, num_added, (*ortho_algns)[i].sign);
				}
			}
		}
	}

	if( num_added > 0 ) {
		*symm_ortho_algns = (struct DotList *) ckrealloc(*symm_ortho_algns, (num2+num_added) * sizeof(struct DotList));
		initialize_algns(*symm_ortho_algns, num2, num2+num_added);
		for( i = 0; i < num_added; i++ ) {
			assign_algn(*symm_ortho_algns, num2+i, added_algns[i]);
		}
		num2 = num2 + num_added;
	}

	*num_ortho = num1;
	*num_symm_ortho = num2;
	*num_symm = num3;
	free(added_algns);
	free(add_list);
	free(cut_list);
	free(list);
	free(st1);
	free(st2);
}

int check_algn_in(struct DotList *ortho_algns, int id, struct slist *st, int s_loc, int e_loc, struct DotList *symm_algns, FILE *g, FILE *f, struct short_alist *list)
{
	int i = 0, j = 0;
	struct I query_x, query_y, cmp, reg, tmp;
	int *res_b, *res_e;
	bool is_done = false;
	int b = 0, e = 0;
	int num_list = 0;
	struct I temp_reg; 

	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));
	*res_b = -1;
	*res_e = -1;
	cmp = assign_I(0, 1);
	reg = assign_I(0, 1);
	tmp = assign_I(0, 1);
	temp_reg = assign_I(0, 1);
	query_x = assign_I(ortho_algns[id].x.lower+ortho_algns[id].xl_diff, ortho_algns[id].x.upper-ortho_algns[id].xr_diff);
	query_y = assign_I(ortho_algns[id].y.lower+ortho_algns[id].yl_diff, ortho_algns[id].y.upper-ortho_algns[id].yr_diff);
	i = s_loc;
	while( (i <= e_loc) && (is_done == false) ) {
		j = st[i].id;
		cmp = assign_I(symm_algns[j].y.lower+symm_algns[j].yl_diff, symm_algns[j].y.upper-symm_algns[j].yr_diff);
		reg = assign_I(symm_algns[j].x.lower+symm_algns[j].xl_diff, symm_algns[j].x.upper-symm_algns[j].xr_diff);

		if( (symm_algns[j].sign == DELETED) || (symm_algns[j].sign != ortho_algns[id].sign) ) {}
		else if( cmp.lower >= query_x.upper ) {}
		else if( cmp.upper <= query_x.lower ) {}
		else {
			find_overlapping_ends(query_x, cmp, SELF2, symm_algns, j, g, res_b, res_e);
			tmp = assign_I(*res_b, *res_e);
			find_overlapping_ends(cmp, query_x, SELF1, ortho_algns, id, f, res_b, res_e);
			temp_reg = assign_I(*res_b, *res_e);
			if( ((tmp.upper - tmp.lower) < DEL_TH) || ((temp_reg.upper - temp_reg.lower) < DEL_TH) ) {} 
			else if( strict_almost_equal(temp_reg, tmp) == true ) {
				if( cmp.lower < query_x.lower ) b = query_x.lower;
				else b = cmp.lower;

				if( cmp.upper > query_x.upper ) e = query_x.upper;
				else e = cmp.upper;

				list[num_list].id = j;	
				list[num_list].x = assign_I(b, e);	
				list[num_list].y = assign_I(tmp.lower, tmp.upper);	
				num_list++;

				query_x = assign_I(e, query_x.upper);
				if( abs(query_x.upper - e) < DEL_TH ) {
					is_done = true;
				}

				if( ortho_algns[id].sign == 0 ) {
					query_y = assign_I(tmp.upper, query_y.upper);
					if( abs(query_y.upper - tmp.upper) < DEL_TH ) is_done = true;
				}
				else {
					query_y = assign_I(query_y.lower, tmp.lower);
					if( abs(tmp.lower - query_y.lower) < DEL_TH ) is_done = true;
				}
			}	
		}
		i++;
	}

	for( i = 0; i < (num_list-1); i++ ) {
		b = abs(list[i].x.upper - list[i+1].x.lower);
		if( ortho_algns[id].sign == 0 ) {
			e = abs(list[i].y.upper - list[i+1].y.lower);
		}
		else if( ortho_algns[id].sign == 1 ) {
			e = abs(list[i].y.lower - list[i+1].y.upper);
		}
		else {
			fatalf("unexpected sign %d\n", ortho_algns[id].sign);
		}

		if( overlap(list[i].x, list[i+1].x) == true ) b = 0;
		if( overlap(list[i].y, list[i+1].y) == true ) e = 0;

		if( (b <= DEL_TH) || (e <= DEL_TH) ) {
			if( list[i+1].x.upper > list[i].x.lower ) {
				list[i].x = assign_I(list[i].x.lower, list[i+1].x.upper);
			}
			else {
				fatalf("unsorted list [%d,%d] and [%d,%d]\n", list[i].x.lower, list[i].x.upper, list[i+1].x.lower, list[i+1].x.upper);
			}

			if( ortho_algns[id].sign == 0 ) {
				if( list[i+1].y.upper > list[i].y.lower ) {
					list[i].y = assign_I(list[i].y.lower, list[i+1].y.upper);
				}
				else {
					fatalf("unsorted list [%d,%d] and [%d,%d]\n", list[i].y.lower, list[i].y.upper, list[i+1].y.lower, list[i+1].y.upper);
				}
			}
			else {
				if( list[i].y.upper > list[i+1].y.lower ) {
					list[i].y = assign_I(list[i+1].y.lower, list[i].y.upper);
				}
				else {
					fatalf("unsorted list [%d,%d] and [%d,%d]\n", list[i].y.lower, list[i].y.upper, list[i+1].y.lower, list[i+1].y.upper);
				}
			}

			for( j = (i+2); j < num_list; j++) {
				list[j-1].id = list[j].id;
				list[j-1].x = assign_I(list[j].x.lower, list[j].x.upper);
				list[j-1].y = assign_I(list[j].y.lower, list[j].y.upper);
				list[j-1].val = list[j].val;
			}
			num_list--;
		}
	}

	free(res_b);
	free(res_e);
	return(num_list);
}

int toggle_list(struct DotList *ortho_algns, int id, struct short_alist *list, int num_list)
{
	int i = 0, j = 0;
	int b1 = 0, e1 = 0, b2 = 0, e2 = 0;
	int count = 0;
	struct short_alist *temp_list;

	temp_list = (struct short_alist * ) ckalloc((num_list+1) * sizeof(struct short_alist));
	for( i = 0; i <= num_list; i++ ) {
		temp_list[i].id = -1;
		temp_list[i].val = 0;
		temp_list[i].x = assign_I(0, 1);
		temp_list[i].y = assign_I(0, 1);
	}
	
	if( num_list < 0 ) {
		fatalf("negative count: %d", num_list);
	}
	else if( num_list == 0 ) {
		b1 = ortho_algns[id].y.lower+ortho_algns[id].yl_diff;
		e1 = ortho_algns[id].y.upper-ortho_algns[id].yr_diff;	
		b2 = ortho_algns[id].x.lower+ortho_algns[id].xl_diff;
		e2 = ortho_algns[id].x.upper-ortho_algns[id].xr_diff;	
		temp_list[count].y = assign_I(b1, e1);
		temp_list[count].x = assign_I(b2, e2);
		temp_list[count].val = -1;
		count++;
	} 
	else {
		for( i = 0; i < num_list; i++ ) {
			list[i].val = list[i].y.lower;
		}
	
		quick_sort_inc_alist(list, 0, num_list-1);

		b1 = ortho_algns[id].y.lower+ortho_algns[id].yl_diff;
		e1 = list[0].y.lower;	
		if( ortho_algns[id].sign == 0 ) {
			b2 = ortho_algns[id].x.lower+ortho_algns[id].xl_diff;
			e2 = list[0].x.lower;;
		}	
		else {
			b2 = list[0].x.upper;;
			e2 = ortho_algns[id].x.upper-ortho_algns[id].xr_diff;
		}

		if( ((e1-b1) > DEL_TH) && ((e2-b2) > DEL_TH) ) {
			temp_list[count].y = assign_I(b1, e1);
			temp_list[count].x = assign_I(b2, e2);
			temp_list[count].val = ortho_algns[id].index;
			count++;
		}

		i = 1;
		while( i < num_list ) {
			b1 = list[i-1].y.upper;	
			e1 = list[i].y.lower;
			if( ortho_algns[id].sign == 0 ) {
				b2 = list[i-1].x.upper;
				e2 = list[i].x.lower;;
			}	
			else {
				b2 = list[i].x.upper;
				e2 = list[i-1].x.lower;;
			}

			if( ((e1-b1) > DEL_TH) && ((e2-b2) > DEL_TH) ) {
				temp_list[count].y = assign_I(b1, e1);
				temp_list[count].x = assign_I(b2, e2);
				temp_list[count].val = ortho_algns[id].index;
				count++;
			}
			i++;
		}

		b1 = list[num_list-1].y.upper;
		e1 = ortho_algns[id].y.upper-ortho_algns[id].yr_diff;
		if( ortho_algns[id].sign == 0 ) {
			b2 = list[num_list-1].x.upper;
			e2 = ortho_algns[id].x.upper-ortho_algns[id].xr_diff;
		}	
		else {
			b2 = ortho_algns[id].x.lower+ortho_algns[id].xl_diff;
			e2 = list[num_list-1].x.lower;
		}

		if( ((e1-b1) > DEL_TH) && ((e2-b2) > DEL_TH) ) {
			temp_list[count].y = assign_I(b1, e1);
			temp_list[count].x = assign_I(b2, e2);
			temp_list[count].val = ortho_algns[id].index;
			count++;
		}
	}

	j = 0;
	for( i = 0; i < count; i++ ) {
		list[j] = assign_alist(temp_list[i]);
		j++;
	}
	free(temp_list);
	return(count);
}

void adjust_algn_diff(struct DotList *ortho_algns, int num_ortho_algns, struct DotList *algns, int num_algns)
{
	int i = 0, j = 0;
	bool is_in = false;
	int cur_id = -1;

	for( i = 0; i < num_ortho_algns; i++ ) {
		j = 0;
		cur_id = -1;
		is_in = false;
		while( (j < num_algns) && (is_in == false ) ) {
			if( algns[j].indiv_fid == ortho_algns[i].indiv_fid ) {
				is_in = true;				
				cur_id = j;
			}
			j++;
		} 	

		if( is_in == false ) {
			fatalf("%d [%d-%d, %d-%d] not found in initial alignments\n", ortho_algns[i].indiv_fid, ortho_algns[i].x.lower, ortho_algns[i].x.upper, ortho_algns[i].y.lower, ortho_algns[i].y.upper);			
		}
		else {
			if( (ortho_algns[i].indiv_fid != algns[cur_id].indiv_fid) || (ortho_algns[i].sign != algns[cur_id].sign) ) {
				fatal("having different fid or sign in two alignments\n");
			} 
			else {
				ortho_algns[i].index = cur_id;
				ortho_algns[i].fid = algns[cur_id].fid;
				ortho_algns[i].xl_diff = abs(ortho_algns[i].x.lower-algns[cur_id].x.lower);
				ortho_algns[i].xr_diff = abs(ortho_algns[i].x.upper-algns[cur_id].x.upper);
				ortho_algns[i].yl_diff = abs(ortho_algns[i].y.lower-algns[cur_id].y.lower);
				ortho_algns[i].yr_diff = abs(ortho_algns[i].y.upper-algns[cur_id].y.upper);
				ortho_algns[i].x = assign_I(algns[cur_id].x.lower, algns[cur_id].x.upper);
				ortho_algns[i].y = assign_I(algns[cur_id].y.lower, algns[cur_id].y.upper);
			}
		}
	}
}

bool is_pid_higher(struct short_alist cur_list, struct DotList *ortho_algns, int num_ortho, FILE *f, float avg_pid, float cur_pid, int id) 
{
	struct DotList *candi_algns;
	int num_candi = 0;
	float pid = (float)0;
	bool res = false;
	int i = 0;
	float max_pid = (float)0;
	int indiv_fid = -1; 
	struct I x_int, y_int;

	x_int = assign_I(0,1);
	y_int = assign_I(0,1);
	candi_algns = (struct DotList *) ckalloc(num_ortho * sizeof(struct DotList));
	initialize_algns(candi_algns, 0, num_ortho);
	num_candi = search_candi_algns(ortho_algns, num_ortho, cur_list.x, cur_list.y, candi_algns);
	if( id != -1 ) {
		indiv_fid = ortho_algns[id].indiv_fid;
	}

	if( num_candi == 0 ) {
		if( (int)(cur_pid+0.5) >= (int)(avg_pid-SYM_PID_TH-0.5) ) {
			res = true;	
		}
	}
	else {
		max_pid = 0;
		for( i = 0; i < num_candi; i++ ) {
			x_int = assign_I(candi_algns[i].x.lower+candi_algns[i].xl_diff, candi_algns[i].x.upper-candi_algns[i].xr_diff);
			y_int = assign_I(candi_algns[i].y.lower+candi_algns[i].yl_diff, candi_algns[i].y.upper-candi_algns[i].yr_diff);
			if( candi_algns[i].sign == DELETED ) {}
			else if( (x_int.upper - x_int.lower) < 0 ) {}
			else if( (y_int.upper - y_int.lower) < 0 ) {}
			else if( (indiv_fid != -1) && (indiv_fid == candi_algns[i].indiv_fid) ) {}
			else if( (proper_overlap(cur_list.x, x_int) == true) && ((width(intersect(cur_list.x, x_int)) > ERR_SM_TH)) ) {
				pid = cal_pid_part_algn(candi_algns, i, abs(cur_list.x.lower-candi_algns[i].x.lower), abs(cur_list.x.upper-candi_algns[i].x.upper), f, SELF1);
				if( pid > max_pid ) {
					max_pid = pid;
				}
			}
			else if( (proper_overlap(cur_list.y, y_int) == true) && ((width(intersect(cur_list.y, y_int)) > ERR_SM_TH)) ) {
				pid = cal_pid_part_algn(candi_algns, i, abs(cur_list.y.lower-candi_algns[i].y.lower), abs(cur_list.y.upper-candi_algns[i].y.upper), f, SELF2);
				if( pid > max_pid ) {
					max_pid = pid;
				}
			}
		}

		if( max_pid == 0 ) {
			if( (int)(cur_pid+0.5) >= (int)(avg_pid-SYM_PID_TH-0.5) ) {
				res = true;	
			}
		}
		else if( ((int)(cur_pid+0.5) >= (int)(max_pid+0.5)) || ((int)(cur_pid+0.5) >= (int)(avg_pid-SYM_PID_TH-0.5)) ) {
			res = true;
		}
	}

	free(candi_algns);
	return(res);
}

bool is_interval_bigger(struct short_alist cur_list, struct DotList *ortho_algns, int num_ortho, FILE *f, float avg_pid, float cur_pid, int id) 
{
	struct DotList *candi_algns;
	int num_candi = 0;
	float pid = (float)0;
	bool res = true;
	int i = 0;
	float max_pid = (float)0;
	int indiv_fid = ortho_algns[id].indiv_fid;
	struct I x_int, y_int;

	x_int = assign_I(0, 1);
	y_int = assign_I(0, 1);

	candi_algns = (struct DotList *) ckalloc(num_ortho * sizeof(struct DotList));
	initialize_algns(candi_algns, 0, num_ortho);
	num_candi = search_candi_algns(ortho_algns, num_ortho, cur_list.x, cur_list.y, candi_algns);

	if( num_candi == 0 ) {
		if( (int)(cur_pid+0.5) > (int)(avg_pid-5.5) ) {
			res = true;	
		}
		else {
			res = false;
		}
	}
	else {
		max_pid = 0;
		for( i = 0; i < num_candi; i++ ) {
			x_int = assign_I(candi_algns[i].x.lower+candi_algns[i].xl_diff, candi_algns[i].x.upper-candi_algns[i].xr_diff);
			y_int = assign_I(candi_algns[i].y.lower+candi_algns[i].yl_diff, candi_algns[i].y.upper-candi_algns[i].yr_diff);
			if( candi_algns[i].sign == DELETED ) {}
			else if( x_int.upper <= x_int.lower ) {}
			else if( y_int.upper <= y_int.lower ) {}
			else if( (indiv_fid != -1) && (indiv_fid == candi_algns[i].indiv_fid) ) {}
			else if( (proper_overlap(cur_list.x, x_int) == true) && ((width(intersect(cur_list.x, x_int)) > ERR_SM_TH)) ) {
				pid = cal_pid_part_algn(candi_algns, i, abs(cur_list.x.lower-candi_algns[i].x.lower), abs(cur_list.x.upper-candi_algns[i].x.upper), f, SELF1);
				if(strict_almost_equal(cur_list.x, x_int) == true) {}
				else if( strict_subset(cur_list.x, x_int) == true ) {
					res = false;
				}
				if( pid > max_pid ) {
					max_pid = pid;
				}
			}
			else if( (proper_overlap(cur_list.y, y_int) == true) && ((width(intersect(cur_list.y, y_int)) > ERR_SM_TH)) ) {
				pid = cal_pid_part_algn(candi_algns, i, abs(cur_list.y.lower-candi_algns[i].y.lower), abs(cur_list.y.upper-candi_algns[i].y.upper), f, SELF2);
				if(strict_almost_equal(cur_list.y, candi_algns[i].y) == true) {}
				else if( strict_subset(cur_list.y, candi_algns[i].y) == true ) {
					res = false;
				}
				if( pid > max_pid ) {
					max_pid = pid;
				}
			}
		}

		if( res == false ) {}
		else if( max_pid == 0 ) {
			if( (int)(cur_pid+0.5) > (int)(avg_pid-5.5) ) {
				res = true;	
			}
			else {
				res = false;
			}
		}
		else if( ((int)(cur_pid+0.5) >= (int)(max_pid+0.5)) || ((int)(cur_pid+0.5) >= (int)(avg_pid-0.5)) ) {
			res = true;
		}
		else {
			res = false;
		}
	}

	free(candi_algns);
	return(res);
}

int cut_algns(struct DotList **ortho_algns, int id, int num1, struct short_alist *list, int num_list, FILE *f)
{
	int j = 0;

	if( num1 == 0 ) {
	}
	else if( num_list > 1 ) {
		*ortho_algns = (struct DotList *) ckrealloc(*ortho_algns, (num1+num_list-1) * sizeof(struct DotList));
		initialize_algns(*ortho_algns, num1, num1+num_list-1);
	}

	if( num_list > 0 ) {
		for( j = 0; j < (num_list-1); j++ ) {
			assign_algn(*ortho_algns, num1+j, (*ortho_algns)[id]);
			(*ortho_algns)[num1+j].index = num1+j;
		}

		for( j = 0; j < (num_list-1); j++ ) {
			update_algn_diff(*ortho_algns, num1+j, list[j+1].y.lower, f, YL_DIFF);
			update_algn_diff(*ortho_algns, num1+j, list[j+1].y.upper, f, YR_DIFF);
		}
		update_algn_diff(*ortho_algns, id, list[0].y.lower, f, YL_DIFF);
		update_algn_diff(*ortho_algns, id, list[0].y.upper, f, YR_DIFF);
	}

	return(num1);
}

int add_algns(struct short_alist *add_list, int num_add, struct DotList *algns, int num_algns, FILE *f, struct DotList *added_algns, int num_added, int sign)
{
	int count = 0;
	struct DotList *candi_algns; 
	int num_candi = 0;
	int i = 0, j = 0;
	struct I cmp, reg, query_x, query_y, tmp;
	int *res_b, *res_e;

	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));

	count = num_added;
	cmp = assign_I(0,1);
	reg = assign_I(0,1);
	tmp = assign_I(0,1);
	query_x = assign_I(0,1);
	query_y = assign_I(0,1);
	candi_algns = (struct DotList *) ckalloc(num_algns * sizeof(struct DotList));	
	initialize_algns(candi_algns, 0, num_algns);

	for ( i = 0; i < num_add; i++ ) {
		num_candi = search_candi_algns(algns, num_algns, add_list[i].y, add_list[i].x, candi_algns);
		query_x = assign_I(add_list[i].y.lower, add_list[i].y.upper);
		query_y = assign_I(add_list[i].x.lower, add_list[i].x.upper);

		for( j = 0; j < num_candi; j++ ) {	
			cmp = assign_I(candi_algns[j].x.lower, candi_algns[j].x.upper);
			reg = assign_I(candi_algns[j].y.lower, candi_algns[j].y.upper);

			if( (candi_algns[j].sign == DELETED) || (candi_algns[j].sign != sign) ) {}
			else if( cmp.lower >= query_x.upper ) {}
			else if( cmp.upper <= query_x.lower ) {}
			else {
				*res_b = -1;
				*res_e = -1;
				find_overlapping_ends(query_x, cmp, SELF1, candi_algns, j, f, res_b, res_e);
				tmp = assign_I(*res_b, *res_e);
				if( ((*res_e) - (*res_b)) < DEL_TH ) {} 
				else if( strict_almost_equal(query_y, tmp) == true ) {
					assign_algn(added_algns, count, candi_algns[j]);
					update_algn_diff(added_algns, count, query_x.lower, f, XL_DIFF);
					update_algn_diff(added_algns, count, query_x.upper, f, XR_DIFF);
					count++;
				}
/*
				else if( strict_subset(tmp, query_y) == true ) {
					assign_algn(added_algns, count, candi_algns[j]);
					update_algn_diff(added_algns, count, tmp.lower, f, YL_DIFF);
					update_algn_diff(added_algns, count, tmp.upper, f, YR_DIFF);
					count++;
				}
*/
			}
		}
	}

	free(res_b);
	free(res_e);
	free(candi_algns);

	return(count);
}

void chain_back(struct DotList *ortho_algns, int num_ortho)
{
	struct slist *st;
	int i = 0, j = 0;
	int cur_id = 0, tmp_id = 0;
	struct I reg_x1, reg_x2;
	struct I reg_y1, reg_y2;

	reg_x1 = assign_I(0, 1);
	reg_x2 = assign_I(0, 1);
	reg_y1 = assign_I(0, 1);
	reg_y2 = assign_I(0, 1);

	st = (struct slist *) ckalloc(num_ortho * sizeof(struct slist));
	initialize_slist(st, 0, num_ortho);
	sort_init_algns(st, ortho_algns, num_ortho, SELF1);

	for( i = 0; i < num_ortho; i++ ) {
		cur_id = st[i].id;
		if( ortho_algns[cur_id].sign == DELETED ) {}
		else {
			for( j = (i+1); j < num_ortho; j++ ) {
				tmp_id = st[j].id;
				if( ortho_algns[tmp_id].sign == DELETED ) {}
				else if( ortho_algns[cur_id].indiv_fid == ortho_algns[tmp_id].indiv_fid ) {
					reg_x1 = assign_I(ortho_algns[cur_id].x.lower+ortho_algns[cur_id].xl_diff, ortho_algns[cur_id].x.upper-ortho_algns[cur_id].xr_diff);
					reg_x2 = assign_I(ortho_algns[tmp_id].x.lower+ortho_algns[tmp_id].xl_diff, ortho_algns[tmp_id].x.upper-ortho_algns[tmp_id].xr_diff);
					reg_y1 = assign_I(ortho_algns[cur_id].y.lower+ortho_algns[cur_id].yl_diff, ortho_algns[cur_id].y.upper-ortho_algns[cur_id].yr_diff);
					reg_y2 = assign_I(ortho_algns[tmp_id].y.lower+ortho_algns[tmp_id].yl_diff, ortho_algns[tmp_id].y.upper-ortho_algns[tmp_id].yr_diff);

					if( (ortho_algns[cur_id].sign == 0) && (  ( overlap(reg_x1, reg_x2) == true ) || (abs(reg_x2.lower-reg_x1.upper) < DEL_TH) || ( overlap(reg_y1, reg_y2) == true ) || (abs(reg_y2.lower-reg_y1.upper) < DEL_TH) ) )  {
						if( ortho_algns[tmp_id].xl_diff > ortho_algns[cur_id].xl_diff ) {
							ortho_algns[tmp_id].xl_diff = ortho_algns[cur_id].xl_diff;
							ortho_algns[tmp_id].yl_diff = ortho_algns[cur_id].yl_diff;
						}

						ortho_algns[cur_id].sign = DELETED;
					}
					else if( (ortho_algns[cur_id].sign == 1) && (  ( overlap(reg_x1, reg_x2) == true ) || (abs(reg_x2.lower-reg_x1.upper) < DEL_TH) || ( overlap(reg_y1, reg_y2) == true ) || (abs(reg_y1.lower-reg_y2.upper) < DEL_TH) ) )  {
						if( ortho_algns[tmp_id].xl_diff > ortho_algns[cur_id].xl_diff ) {
							ortho_algns[tmp_id].xl_diff = ortho_algns[cur_id].xl_diff;
							ortho_algns[tmp_id].yr_diff = ortho_algns[cur_id].yr_diff;
						}

						ortho_algns[cur_id].sign = DELETED;
					}
				}
			}
		}
	}

	free(st);
}
