// non-human species vs. human

#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "regions.h"
#include "util_input.h"
#include "util_ops.h"
#include "util_gen.h"
#include "find_ps.h"
#include "contigs_op.h"

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

int main(int argc, char **argv)
{
	struct DotList *algns;
	int *num_algns;
	int *size1, *size2;
	int count = 0;
	char species[100], species2[100];
	char temp_name[100];
	FILE *f;
	struct exons_list *exons1, *exons2, *genes1, *genes2; 
	int *num_exons1, *num_genes1;
	int *num_exons2, *num_genes2;
	struct exons_list *candi_ps;
	int *num_candi, *num_alloc;
	int temp_num1 = 0, temp_num2 = 0;
	int i = 0, j = 0;
	int b = 0, e = 0;
	int lo = 0, hi = 1;
  struct n_pair *contigs1, *contigs2;
  int num_contigs1 = 0, num_contigs2 = 0;
  int *len_sum1, *len_sum2;
	int len = 0;
	int ctg_id = -1;


	debug_mode = FALSE;
	if( argc == 7 ) {
		debug_mode = TRUE;
	}
	else if( argc != 6 ) {
		fatal("args: $(species2).human.chained_maf human.codex $(species2).codex ref_sp contigs_list\n");
	}

	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	*size1 = 0;
	*size2 = 0;
	*num_algns = 0;
	count = count_local_algns(argv[1], species, species2); // note species2 is the reference name
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(algns, 0, count);

	 	read_maf(argv[1], PAIR_MODE, algns, num_algns, size1, size2);
	}

	if( count != (*num_algns) ) {
		fatalf("counting error 1: %d vs. %d in %s vs. %s\n", count, *num_algns, species, species2);
	}

	if( (f = ckopen(argv[5], "r")) != NULL ) {
  	num_contigs1 = count_contigs(species, f);
  	num_contigs2 = count_contigs(species2, f);
	}
	else fatalf("%s not exist\n", argv[5]);

  if( num_contigs1 > 0 ) {
    contigs1 = (struct n_pair *) ckalloc(num_contigs1 * sizeof(struct n_pair));
    len_sum1 = (int *) ckalloc(num_contigs1 * sizeof(int));
    read_contigs_file(species, f, contigs1, num_contigs1);
    cal_length_sum_from_ctglist(len_sum1, contigs1, num_contigs1);
  }

  if( num_contigs2 > 0 ) {
    contigs2 = (struct n_pair *) ckalloc(num_contigs2 * sizeof(struct n_pair));
    len_sum2 = (int *) ckalloc(num_contigs2 * sizeof(int));
    read_contigs_file(species2, f, contigs2, num_contigs2);
    cal_length_sum_from_ctglist(len_sum2, contigs2, num_contigs2);
  }
  fclose(f);
	
	if( (*num_algns) > 0 ) {
		adjust_algn_pos(algns, *num_algns, contigs1, num_contigs1, size1, contigs2, num_contigs2, size2, CTG_NOT_ASSIGNED_BUT_LEN);
	}

	num_exons1 = (int *) ckalloc(sizeof(int));
	num_exons2 = (int *) ckalloc(sizeof(int));
	*num_exons1 = 0;
	*num_exons2 = 0;
	num_genes1 = (int *) ckalloc(sizeof(int));
	num_genes2 = (int *) ckalloc(sizeof(int));
	*num_genes1 = 0;
	*num_genes2 = 0;

	count = count_genes(argv[2], species2, species2);
	genes2 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

	if( count > 0 ) {
		initialize_exons_list(genes2, 0, count);
	}

	count = count_exons(argv[2], species2, species2);
	exons2 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

	if( count > 0 ) {
		initialize_exons_list(exons2, 0, count);
		read_only_exons(exons2, num_exons2, genes2, num_genes2, argv[2], species2, species2, contigs2, contigs2, num_contigs2, num_contigs2);
		for( i = 0; i < *num_exons2; i++ ) adjust_pos_exons(exons2, i, len_sum2, num_contigs2);
		for( i = 0; i < *num_genes2; i++ ) adjust_pos_exons(genes2, i, len_sum2, num_contigs2);
	}

	if( count != (*num_exons2) ) {
		fatalf("counting error 2: %d vs. %d in %s\n", count, *num_exons2, species2);
	}

	count = count_genes(argv[3], species, species);
	genes1 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

	if( count > 0 ) {
		initialize_exons_list(genes1, 0, count);
	}

	count = count_exons(argv[3], species, species);
	exons1 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

	if( count > 0 ) {
		initialize_exons_list(exons1, 0, count);
		read_only_exons(exons1, num_exons1, genes1, num_genes1, argv[3], species, species, contigs1, contigs1, num_contigs1, num_contigs1);
		for( i = 0; i < *num_exons1; i++ ) adjust_pos_exons(exons1, i, len_sum1, num_contigs1);
		for( i = 0; i < *num_genes1; i++ ) adjust_pos_exons(genes1, i, len_sum1, num_contigs1);
	}

	if( count != (*num_exons1) ) {
		fatalf("counting error 3: %d vs. %d in %s\n", count, *num_exons1, species);
	}

	f = ckopen(argv[1], "r");
  num_candi = (int *) ckalloc(sizeof(int));
  num_alloc = (int *) ckalloc(sizeof(int));
  *num_alloc = ALLOC_UNIT;
  *num_candi = 0;
  candi_ps = (struct exons_list *) ckalloc(ALLOC_UNIT * sizeof(struct exons_list));
  initialize_exons_list(candi_ps, 0, ALLOC_UNIT);

	if( (*num_algns) > 0 ) {
		find_pseudogenes(genes2, num_genes2, genes1, num_genes1, exons1, num_exons1, algns, *num_algns, f, &candi_ps, num_candi, num_alloc);
	}

	if( (*num_candi) > 0 ) {
		sort_exons(candi_ps, *num_candi);
	}

	i = 0;
	while( i < (*num_candi) ) {
		if( candi_ps[i].sign != 'n' ) {
			j = i + 1;
			while( (j < (*num_candi)) && (candi_ps[j].sign != 'n') && (proper_overlap(candi_ps[i].reg, candi_ps[j].reg) == true) ) {
				if( candi_ps[i].reg.lower < candi_ps[j].reg.lower) lo = candi_ps[i].reg.lower;
				else lo = candi_ps[j].reg.lower;

				if( candi_ps[i].reg.upper < candi_ps[j].reg.upper) hi = candi_ps[j].reg.upper;
				else hi = candi_ps[i].reg.upper;
							
				candi_ps[i].reg = assign_I(lo, hi);
				candi_ps[j].sign = 'n';
				j++;
			}
			i = j;
		}
		else {
			i++;
		}
	}

	temp_num1 = *num_genes2;
	temp_num2 = *num_exons2;
	if( (*num_candi) > 0 ) {
		sort_exons(candi_ps, *num_candi);
		genes2 = ckrealloc(genes2, ((*num_genes2) + (*num_candi)) * sizeof(struct exons_list));
		exons2 = ckrealloc(exons2, ((*num_exons2) + (*num_candi)) * sizeof(struct exons_list));
		j = 1;
		for( i = 0; i < (*num_candi); i++ ) {
			if( candi_ps[i].sign != 'n' ) {
				sprintf(temp_name, "%s_ps", candi_ps[i].name);
				j++;
				genes2[temp_num1] = assign_exons(candi_ps[i]);
				exons2[temp_num2] = assign_exons(candi_ps[i]);
				strcpy(genes2[temp_num1].name, temp_name);
				strcpy(exons2[temp_num2].name, temp_name);
				exons2[temp_num2].fid = temp_num1;
				temp_num1++;
				temp_num2++;
			}
		}
	}

	*num_genes2 = temp_num1;
	*num_exons2 = temp_num2;

	for( i = 0; i < temp_num2; i++ ) {
		j = exons2[i].fid;
		b = i;
		while( (i < (temp_num2-1)) && (exons2[i].fid == exons2[i+1].fid) ) {
			i++;
		}
		e = i;
		genes2[j].cmp_reg = assign_I(b, e); // cmp_reg saves index for start exon and end temporally
	}

	sort_exons(genes2, *num_genes2);
	j = 1;
	for( i = 0; i < (*num_genes2); i++ ) {
		ctg_id = -1;
		if( (strstr(species2, "human") != NULL) || (strcmp(species2, argv[4]) == 0) ) {
			strcpy(temp_name, genes2[i].name);
		}
		else if( strstr(genes2[i].name, "_ps") != NULL ) {
			sprintf(temp_name, "%s%d_ps", species2, j);
			j++;
		}
		else {
			sprintf(temp_name, "%s%d", species2, j);
			j++;
		}
		
		ctg_id = genes2[i].ctg_id;
		if( ctg_id != -1 ) len = len_sum2[ctg_id];
		else len = 0;

		if( genes2[i].sign == '<' ) {
			if( ctg_id == -1 ) {
				printf("%c %d %d %s %s (complement)\n", genes2[i].sign, genes2[i].reg.lower, genes2[i].reg.upper-1, temp_name, species2);
			}
			else {
				printf("%c %d %d %s %s (complement)\n", genes2[i].sign, genes2[i].reg.lower-len, genes2[i].reg.upper-1-len, temp_name, contigs2[ctg_id].name2);
			}
		}
		else if( genes2[i].sign == '>') {
			if( ctg_id == -1 ) {
				printf("%c %d %d %s %s\n", genes2[i].sign, genes2[i].reg.lower, genes2[i].reg.upper-1, temp_name, species2);
			}
			else {
				printf("%c %d %d %s %s\n", genes2[i].sign, genes2[i].reg.lower-len, genes2[i].reg.upper-1-len, temp_name, contigs2[ctg_id].name2);
			}
		}

		if( (genes2[i].sign == '<') || (genes2[i].sign == '>') ) {
			for( b = genes2[i].cmp_reg.lower; b <= genes2[i].cmp_reg.upper; b++ ) {
				printf("%d %d\n", exons2[b].reg.lower-len, exons2[b].reg.upper-1-len);
			}
		}
	}

	fclose(f);
	free(size1);
	free(size2);

	if( (*num_algns) > 0 ) {
		free(algns);
	}

  if( num_contigs1 > 0 ) {
    free(contigs1);
    free(len_sum1);
  }

  if( num_contigs2 > 0 ) {
    free(contigs2);
    free(len_sum2);
  }

	free(num_candi);
	free(num_alloc);
	free(candi_ps);
	free(num_algns);
	free(exons1);
	free(genes1);
	free(num_exons1);
	free(num_genes1);
	free(exons2);
	free(genes2);
	free(num_exons2);
	free(num_genes2);
	return(EXIT_SUCCESS);
}

