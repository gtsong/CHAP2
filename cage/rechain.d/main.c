#include "main.h"
#include "read_maf.h"
#include "rechain.h"
#include "util.h"
#include "util_reduced.h"
#include "regions.h"

char S[BIG], T[BIG];
int debug_mode;

int main(int argc, char **argv)
{
	int i = 0;
	struct DotList *init_algns;
	int *size1, *size2; 
	int *num_algns;
	int count = 0;
	FILE *fp;
	char *status;
	char species[100], species2[100];
	int b1, b2, e1, e2;
  char len1[100], len2[100];
  char strand[100];
	struct slist *st;
	
	debug_mode = FALSE;
	if((argc == 4) && (strcmp(argv[2], "debug-mode") == 0)) {
		debug_mode = TRUE;
	}
	else if(argc != 3 )
	{
		fatal("args: final_ortho_maf chained_maf (and debug-mode)\n");
	}

	num_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

  fp = ckopen(argv[1], "r");
  count = 0;
  if( (fgets(S, BIG, fp) == NULL) || (strncmp(S, "##maf version", 13)))
    fatalf("%s is not a maf file", argv[1]);
  while(S[0] == '#')
    if((status = fgets(S, BIG, fp)) == NULL)
      fatalf("no alignments in %s", argv[1]);

  while((status != NULL) && (strstr(S, "eof") == NULL)) {
    if (S[0] != 'a') {
			fatalf("expecting an a-line in %s, saw %s", argv[1], S);
		}
		else count++;

    if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL))
      fatalf("cannot find alignment in %s", argv[1]);
    if( (sscanf(S, "%*s %s %d %d %*s %s", species, &b1, &e1, len1) != 4) ||
        (sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5))
    {
      fatalf("bad alignment info of %s and %s in %s", S, T, argv[1]);
    }

    if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
      fatalf("bad alignment end in %s", argv[1]);
    status = fgets(S, BIG, fp);
		if( S[0] == '#' ) {
			status = fgets(S, BIG, fp);
		}
  }
  fclose(fp);

	if( count > 0 ) {
 		init_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
 		st = (struct slist *) ckalloc(count * sizeof(struct slist));
 		read_maf(argv[1], B_MODE, init_algns, num_algns, size1, size2);

		for( i = 0; i < (*num_algns); i++ ) {
			st[i].id = i;
		}
		sort_by_bid(st, init_algns, *num_algns);
		merge_overlaps(st, init_algns, *num_algns);
		overwrite_dots(num_algns, init_algns);
		sort_by_bid(st, init_algns, *num_algns);

		fp = ckopen(argv[2], "r");	
		restore_chained_algns(st, init_algns, *num_algns, species, species2, *size1, *size2, fp);
		fclose(fp);
		free(st);
		free(init_algns);
	}
	else {
    printf("##maf version=1 scoring=lastz-pid\n");
	}

	free(num_algns);
	free(size1);
	free(size2);
	return(EXIT_SUCCESS);
}
