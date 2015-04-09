#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_input.h"
#include "util_gen.h"
#include "util_algns.h"
#include "util_i.h"
#include "regions.h"
#include "read_algn.h"

#define ORTHO_CONTENT 1
#define ORTHO_POSITION 2
#define SIMPLE_CUT 3

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
	int i = 0;
	FILE *f;
	struct slist *mark1;
	int *num1, *num2;
	int num_bases = 0, num_mbases = 0;

	debug_mode = FALSE;
	if( argc == 2 ) {
		debug_mode = TRUE;
	}
	else if( argc != 1 ) {
		fatal("args: merged_maf_file");
	}

	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	num1 = (int *) ckalloc(sizeof(int));
	num2 = (int *) ckalloc(sizeof(int));
	num_algns = (int *) ckalloc(sizeof(int));
	*size1 = 0;
	*size2 = 0;
	*num_algns = 0;
	count = count_local_algns(argv[1], species, species2);
	if( count > 0 ) {
		algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(algns, 0, count);

		read_maf(argv[1], G_MODE, algns, num_algns, size1, size2);
	}

	if( count != (*num_algns) ) {
		fatalf("counting error 1: %d vs. %d\n", count, *num_algns);
	}

	f = ckopen(argv[1], "r");

	mark1 = (struct slist *) ckalloc(count * sizeof(struct slist));
	initialize_slist(mark1, 0, count);
	for( i = 0; i < (*num_algns); i++ ) {
    mark1[i].val = 0; // the number of nucleotides in a local alignment
    mark1[i].val_red = 0; // the number of nucleotides in common
    mark1[i].sp_state = 0; // the number of aligned nucleotides 
    mark1[i].add_sp_state = 0; // the number of matched nucleotides
    mark1[i].id = algns[i].fid;
    mark1[i].is_x = false;
	}

  for( i = 0; i < (*num_algns); i++ ) {
    *num1 = 0;
    *num2 = 0;
    mark1[i].add_sp_state = count_nucs(algns[i], f, width(algns[i].x), num1, num2);
    mark1[i].val = *num1;
    mark1[i].sp_state = *num2;
  }

	for( i = 0; i < (*num_algns); i++ ) {
		num_bases = num_bases + mark1[i].sp_state;
		num_mbases = num_mbases + mark1[i].add_sp_state;
	}

	printf("%d %d %d ", *num_algns, num_bases, num_mbases);
	fclose(f);
	free(size1);
	free(size2);

	if( (*num_algns) > 0 ) {
		free(algns);
	}

	free(num_algns);
	free(num1);
	free(num2);
	free(mark1);
	return(EXIT_SUCCESS);
}

