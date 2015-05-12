#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_input.h"
#include "util_gen.h"
#include "map_genes.h"
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
	FILE *f;
	struct exons_list *exons1, *genes1; 
	struct exons_list *exons2, *genes2; 
	int *num_exons1, *num_genes1;
	int *num_exons2, *num_genes2;
	float avg_pid = (float)0;
	int i = 0;
	struct n_pair *contigs1, *contigs2; 
	int num_contigs1 = 0, num_contigs2 = 0;
	int *len_sum1, *len_sum2;

	debug_mode = FALSE;
	if( argc == 5 ) {
		debug_mode = TRUE;
	}
	else if( argc != 4 ) {
		fatal("args: $(species1).$(species2).ortho.maf $(species1).codex $(species2).codex contigs.list\n");
	}

	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	*size1 = 0;
	*size2 = 0;
	count = count_local_algns(argv[1], species, species2);
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(algns, 0, count);

	 	read_maf(argv[1], PAIR_MODE, algns, num_algns, size1, size2);
	}

	if( count != (*num_algns) ) {
		fatalf("counting error 1: %d vs. %d\n", count, *num_algns);
	}

	if( (*num_algns) > 0 ) {
		avg_pid = cal_avg_pid(algns, *num_algns);
	}

  f = ckopen(argv[4], "r");
  num_contigs1 = count_contigs(species, f);
  num_contigs2 = count_contigs(species2, f);

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
	*num_exons1 = 0;
	num_genes1 = (int *) ckalloc(sizeof(int));
	*num_genes1 = 0;

	count = count_genes(argv[2], species, species);
	genes1 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));
	if( count > 0 ) {
		initialize_exons_list(genes1, 0, count);
	}

	count = count_exons(argv[2], species, species);
	exons1 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));
	if( count > 0 ) {
		initialize_exons_list(exons1, 0, count);
		read_only_exons(exons1, num_exons1, genes1, num_genes1, argv[2], species, species, contigs1, contigs1, num_contigs1, num_contigs1);
		for( i = 0; i < *num_exons1; i++ ) adjust_pos_exons(exons1, i, len_sum1, num_contigs1);
    for( i = 0; i < *num_genes1; i++ ) adjust_pos_exons(genes1, i, len_sum1, num_contigs1);
	}

	if( count != (*num_exons1) ) {
		fatalf("counting error 3: %d vs. %d\n", count, *num_exons1);
	}

	num_exons2 = (int *) ckalloc(sizeof(int));
	*num_exons2 = 0;
	num_genes2 = (int *) ckalloc(sizeof(int));
	*num_genes2 = 0;

	count = count_genes(argv[3], species2, species2);
	genes2 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));
	if( count > 0 ) {
		initialize_exons_list(genes2, 0, count);
	}

	count = count_exons(argv[3], species2, species2);
	exons2 = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));
	if( count > 0 ) {
		initialize_exons_list(exons2, 0, count);
		read_only_exons(exons2, num_exons2, genes2, num_genes2, argv[3], species2, species2, contigs2, contigs2, num_contigs2, num_contigs2);
		for( i = 0; i < *num_exons2; i++ ) adjust_pos_exons(exons2, i, len_sum2, num_contigs2);
    for( i = 0; i < *num_genes2; i++ ) adjust_pos_exons(genes2, i, len_sum2, num_contigs2);
	}

	f = ckopen(argv[1], "r");

	for( i = 0; i < (*num_genes2); i++ ) {
		genes2[i].val = -1;
	}

	printf("%s", species2);
	if( *num_algns > 0 ) {
//		map_genes(algns, *num_algns, exons1, *num_exons1, genes1, *num_genes1, genes2, *num_genes2, f); 
		map_genes_partition(algns, *num_algns, exons1, *num_exons1, genes1, *num_genes1, genes2, *num_genes2, f); 
	}
	else {
		for( i = 0; i < (*num_genes2); i++ ) {
      if( strstr(genes2[i].name, "_ps") != NULL ) {
        printf("\t\\%d;", -1);
      }
      else {
        printf("\t%d;", -1);
			}
		}
	}
	printf("\n");

	fclose(f);
	free(size1);
	free(size2);

	if( (*num_algns) > 0 ) {
		free(algns);
	}

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

