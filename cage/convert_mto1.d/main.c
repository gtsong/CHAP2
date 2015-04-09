/*
	converting many-to-many into many-to-one
*/

#include "main.h"
#include "regions.h"
#include "read_maf.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "write_init_maf.h"
#include "map_algns.h"

int debug_mode;
char S[BIG], T[BIG];
char S1[BIG], T1[BIG];

int main(int argc, char **argv)
{
	struct DotList *init_algns; 
	int *num_init_algns; // the number of local alignments in the initial dot plot
  char species[100], species2[100];
	int *size1, *size2;
	FILE *f;
	int count = 0;

	if(argc != 2 )
	{
		fatal("args: dots-file \n");
	}

  num_init_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

	if( (f = fopen(argv[1], "r")) == NULL ) {
		fatalf("file %s does not exist\n", argv[1]);
	}

	while(fgets(S, BIG, f)) {
		if( S[0] == '#' ) {
			while( S[0] == '#' ) {}

			count = 0;	
		}

  	if( S[0] == 'a' ) {
			count++;

			if( count == 1 ) {
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", species) != 1) || (sscanf(T, "%*s %s %*s", species2) != 1)) {}
  		}
		}
  }
	fclose(f);

	init_algns = (struct DotList *) ckalloc(sizeof(struct DotList) * count);
	initialize_algns(init_algns, count);
	read_maf(argv[1], G_MODE, init_algns, num_init_algns, size1, size2);
	write_init_maf(init_algns, *num_init_algns, species, species2, *size1, *size2, f, PAIR);

	map_one_to_one(*num_init_algns, init_algns, f);
	

	free(size2);
	free(size1);
	free(init_algns);
	free(num_init_algns);
	return EXIT_SUCCESS;
}
