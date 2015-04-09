#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_input.h"
#include "util_gen.h"
#include "analysis.h"
#include "contigs_op.h"

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

int main(int argc, char **argv)
{
	struct DotList *content_ortho;
	struct DotList *position_ortho;
	int *num_content, *num_position;
	int *size1, *size2;
	int count = 0;
	char species[100], species2[100];
	FILE *f, *g;
	struct exons_list *exons, *genes; 
	int *num_exons, *num_genes;
	float content_avg_pid = (float)0; 
	float position_avg_pid = (float)0; 
	struct n_pair *contigs1, *contigs2;
	int num_contigs1 = 0, num_contigs2 = 0;
	int *len_sum1, *len_sum2;
	int i = 0;

	debug_mode = FALSE;
	if( argc == 5 ) {
		debug_mode = TRUE;
	}
	else if( argc != 4 ) {
		fatal("args: content_ortho_maf position_ortho_maf annot.codex contigs.list\n");
	}

	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num_content = (int *) ckalloc(sizeof(int));
	num_position = (int *) ckalloc(sizeof(int));
	*size1 = 0;
	*size2 = 0;
	*num_content = 0;
	*num_position = 0;
	count = count_local_algns(argv[1], species, species2);
	if( count > 0 ) {
		content_ortho = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(content_ortho, 0, count);

	 	read_maf(argv[1], O_MODE, content_ortho, num_content, size1, size2);
	}

	if( count != (*num_content) ) {
		fatalf("counting error 1: %d vs. %d\n", count, *num_content);
	}

	count = count_local_algns(argv[2], species, species2);
	if( count > 0 ) {
		position_ortho = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(position_ortho, 0, count);

	 	read_maf(argv[2], O_MODE, position_ortho, num_position, size1, size2);
	}

	if( count != (*num_position) ) {
		fatalf("counting error 2: %d vs. %d\n", count, *num_position);
	}

  f = ckopen(argv[4], "r");
  num_contigs1 = count_contigs(species, f);
  num_contigs2 = count_contigs(species2, f);

  if( num_contigs1 > 0 ) {
    contigs1 = (struct n_pair *) ckalloc(num_contigs1 * sizeof(struct n_pair));
    len_sum1 = (int *) ckalloc(num_contigs1 * sizeof(int));
    read_contigs_file(species, f, contigs1, num_contigs1);
    cal_length_sum(len_sum1, contigs1, num_contigs1);
  }

  if( num_contigs2 > 0 ) {
    contigs2 = (struct n_pair *) ckalloc(num_contigs2 * sizeof(struct n_pair));
    len_sum2 = (int *) ckalloc(num_contigs2 * sizeof(int));
    read_contigs_file(species2, f, contigs2, num_contigs2);
    cal_length_sum(len_sum2, contigs2, num_contigs2);
  }
  fclose(f);
	
	if( (*num_content) > 0 ) {
		adjust_algn_pos(content_ortho, *num_content, contigs1, num_contigs1, size1, contigs2, num_contigs2, size2, CTG_NOT_ASSIGNED);
	}

	if( (*num_position) > 0 ) {
		adjust_algn_pos(position_ortho, *num_position, contigs1, num_contigs1, size1, contigs2, num_contigs2, size2, CTG_NOT_ASSIGNED);
	}

	num_exons = (int *) ckalloc(sizeof(int));
	*num_exons = 0;
	num_genes = (int *) ckalloc(sizeof(int));
	*num_genes = 0;

	count = count_genes(argv[3], species, species2);
	genes = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

	count = count_exons(argv[3], species, species2);
	exons = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));
	read_only_exons(exons, num_exons, genes, num_genes, argv[3], species, species2, contigs1, contigs2, num_contigs1, num_contigs2);
  for( i = 0; i < *num_exons; i++ ) {
		if( exons[i].sp_id == SELF1 ) adjust_pos_exons(exons, i, len_sum1, num_contigs1);
		else if( exons[i].sp_id == SELF2 ) adjust_pos_exons(exons, i, len_sum2, num_contigs2);
		else fatalf("species code not assigned: exons %d-%d %s\n", exons[i].reg.lower, exons[i].reg.upper, exons[i].name);
	}

  for( i = 0; i < *num_genes; i++ ) {
		if( genes[i].sp_id == SELF1 ) adjust_pos_exons(genes, i, len_sum1, num_contigs1);
		else if( genes[i].sp_id == SELF2 ) adjust_pos_exons(genes, i, len_sum2, num_contigs2);
		else fatalf("species code not assigned: genes  %d-%d %s\n",  genes[i].reg.lower, genes[i].reg.upper, genes[i].name);
	}

	if( count != (*num_exons) ) {
		fatalf("counting error 3: %d vs. %d\n", count, *num_exons);
	}

	if( (*num_content) > 0 ) {
		content_avg_pid = cal_avg_pid(content_ortho, *num_content);
	}

	if( (*num_position) > 0 ) {
		position_avg_pid = cal_avg_pid(position_ortho, *num_position);
	}
	
	f = ckopen(argv[1], "r");
	g = ckopen(argv[2], "r");

//	get_numbers(content_ortho, *num_content, position_ortho, *num_position, exons, *num_exons, f, g, *size1, *size2); 
	get_numbers(content_ortho, *num_content, position_ortho, *num_position, f, g, *size1, *size2); 

	fclose(g);
	fclose(f);
	free(size1);
	free(size2);

	if( (*num_content) > 0 ) {
		free(content_ortho);
	}

	if( (*num_position) > 0 ) {
		free(position_ortho);
	}

	free(num_content);
	free(num_position);
	free(exons);
	free(genes);
	free(num_exons);
	free(num_genes);
	return(EXIT_SUCCESS);
}

