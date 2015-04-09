#include "main.h"
#include "regions.h"
#include "read_maf.h"
#include "find_merging.h"
#include "remove_rps.h"
#include "kd_op.h"
#include "util.h"

char S[BIG], T[BIG];
int debug_mode;

int main(int argc, char **argv)
{
	struct DotList *self_algns;
	struct DotList *pair_algns;
	struct kdnode *tree;
	struct perm_pt *p_pts;
	struct kdnode *pair_tree;
	struct perm_pt *pair_pts;

	int *num_algns;
	int *num_pair;
	int size_seq1, size_seq2;
	int *size1, *size2;
	int sp_flag; 
	FILE *fp, *out;
	int count = 0, temp = 0;
	int num_self_algns = 0, num_pair_algns = 0;
	char *status;
	char species[NAME_LEN], species2[NAME_LEN];

	int b1 = 0, e1 = 1, b2 = 0, e2 = 1;
	char len1[NAME_LEN], len2[NAME_LEN];
	char strand[NAME_LEN];

	debug_mode = FALSE;
	if( argc == 7 ) {
		debug_mode = TRUE;
	}	
	else if(argc !=6 )
	{
		fatal("args: self_maf pair_maf new_self new_pair sp_flag\n");
	}

	num_algns = (int *) ckalloc(sizeof(int));
	num_pair = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

	sp_flag = atoi(argv[5]);

	fp = ckopen(argv[1], "r");
	count = 0;
	if( ( (status = fgets(S, BIG, fp)) == NULL) || (strncmp(S, "##maf version", 13)))
		fatalf("%s is not a maf file", argv[1]);
	while((status != NULL) && (S[0] == '#')) { 
		status = fgets(S, BIG, fp);
	}
//			fatalf("no alignments in %s", argv[1]);
			
	while((status != NULL) && (strstr(S, "eof") == NULL)) {
		if (strncmp(S, "a ", 2)) fatalf("expecting an a-line in %s, saw %s", argv[1], S);

		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL))
			fatalf("cannot find alignment in %s", argv[1]);
		if( (sscanf(S, "%*s %s %d %d %*s %s", species, &b1, &e1, len1) != 4) || 
				(sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5))
		{
			fatalf("bad alignment info of %s and %s in %s", S, T, argv[1]);
		}
		else {
			e1 += b1;
			e2 += b2;

			if( strcmp(strand, "-") == 0) {
				temp = b2;
				b2 = atoi(len2) - e2;
				e2 = atoi(len2) - temp;
			}

			b1++;
			b2++;
			e1++;
			e2++;
      
			if( (b1 < b2) && (!((b1 == b2) && (e1 == e2))) ) count++;
		}

  	if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
    	fatalf("bad alignment end in %s", argv[1]);
  	status = fgets(S, BIG, fp);
	}
	fclose(fp);
	
	num_self_algns = count;
	if( count > 0 ) {
		self_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		read_maf(argv[1], D_MODE, self_algns, num_algns, size1, size2);	
	}
	else {
		*num_algns = 0;
	}

	fp = ckopen(argv[2], "r");
	count = 0;
	if( ((status = fgets(S, BIG, fp)) == NULL) || (strncmp(S, "##maf version", 13)))
		fatalf("%s is not a maf file", argv[2]);
	while((status != NULL) && (S[0] == '#')) {
		status = fgets(S, BIG, fp);
	}
//			fatalf("no alignments in %s", argv[2]);
			
	while((status != NULL) && (strstr(S, "eof") == NULL)) {
		if (strncmp(S, "a ", 2)) fatalf("expecting an a-line in %s, saw %s", argv[2], S);
		else count++;

		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL) )
			fatalf("cannot find alignment in %s", argv[2]);
		if( (sscanf(S, "%*s %s %d %d %*s %s", species, &b1, &e1, len1) != 4) || 
				(sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5))
		{		
			fatalf("bad alignment info of %s and %s in %s", S, T, argv[2]);
		}

  	if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
    	fatalf("bad alignment end in %s", argv[2]);
  	status = fgets(S, BIG, fp);
	}
	
	fclose(fp);

	num_pair_algns = count;
	if( count > 0 ) {
		pair_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		read_maf(argv[2], G_MODE, pair_algns, num_pair, size1, size2);	
	}
	else {
		*num_pair = 0;
	}

	size_seq1 = *size1;
	size_seq2 = *size2;

	p_pts = (struct perm_pt *) ckalloc((*num_algns) * sizeof(struct perm_pt));
	pair_pts = (struct perm_pt *) ckalloc((*num_pair) * sizeof(struct perm_pt));

	assign_perm(p_pts, (*num_algns), self_algns, LEFT);
	tree = build_kd(p_pts, 0, (*num_algns)-1);
	assign_perm(pair_pts, (*num_pair), pair_algns, LEFT);
	pair_tree = build_kd(pair_pts, 0, (*num_pair)-1);

	if( (*num_algns) > 0 ) {
		remove_overlapped(self_algns, num_algns, tree, p_pts, pair_tree, pair_pts, size_seq1, size_seq2, pair_algns, num_pair, sp_flag);
	}

	fp = ckopen(argv[1], "r");
	out = ckopen(argv[3], "w");
	fprintf(out, "##maf version=1 scoring=lastz-percentage-identity\n");
	count = 0;
	if( ((status = fgets(S, BIG, fp)) == NULL) || (strncmp(S, "##maf version", 13)))
		fatalf("%s is not a maf file", argv[1]);
	while((status != NULL) && (S[0] == '#')) {
		status = fgets(S, BIG, fp);
	}
//			fatalf("no alignments in %s", argv[1]);
			
	while((status != NULL) && (strstr(S, "eof") == NULL)) {
		if (strncmp(S, "a ", 2)) fatalf("expecting an a-line in %s, saw %s", argv[1], S);

		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL) )
			fatalf("cannot find alignment in %s", argv[1]);
		if( (sscanf(S, "%*s %s %d %d %*s %s", species, &b1, &e1, len1) != 4) || 
				(sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5))
		{
			fatalf("bad alignment info in %s", argv[1]);
		}
		else {
			e1 += b1;
			e2 += b2;

			if( strcmp(strand, "-") == 0 ) {
				temp = b2;
				b2 = atoi(len2) - e2;
				e2 = atoi(len2) - temp;
			}

			b1++;
			b2++;
			e1++;
			e2++;

			if( (b1 < b2) && (!((b1 == b2) && (e1 == e2) )) ) {
				if( self_algns[count].sign == 2 ) {}
				else {
					fprintf(out, "a score=%d\n", self_algns[count].identity); 
					fprintf(out, "%s", S); 
					fprintf(out, "%s", T); 
					fprintf(out, "\n");
				}
				count++;
			}
		}

  	if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
    	fatalf("bad alignment end in %s", argv[2]);
  	status = fgets(S, BIG, fp);
	}
	fclose(out);
	fclose(fp);

	fp = ckopen(argv[2], "r");
	out = ckopen(argv[4], "w");
	fprintf(out, "##maf version=1 scoring=lastz-percetage-identity\n");
	count = 0;
	if( ((status = fgets(S, BIG, fp)) == NULL) || (strncmp(S, "##maf version", 13)))
		fatalf("%s is not a maf file", argv[1]);
	while((status != NULL) && (S[0] == '#')) {
 		status = fgets(S, BIG, fp);
	}
//			fatalf("no alignments in %s", argv[1]);
			
	while((status != NULL) && (strstr(S, "eof") == NULL)) {
		if (strncmp(S, "a ", 2)) fatalf("expecting an a-line in %s, saw %s", argv[1], S);

		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL) )
			fatalf("cannot find alignment in %s", argv[1]);
		if( (sscanf(S, "%*s %s %d %d %*s %s", species, &b1, &e1, len1) != 4) || 
				(sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5))
		{
			fatalf("bad alignment info in %s", argv[1]);
		}
		else {
			if( pair_algns[count].sign == 2 ) {}
			else {
				fprintf(out, "a score=%d\n", pair_algns[count].identity); 
				fprintf(out, "%s", S); 
				fprintf(out, "%s", T); 
				fprintf(out, "\n"); 
			}
			count++;
		}

  	if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
    	fatalf("bad alignment end in %s", argv[2]);
  	status = fgets(S, BIG, fp);
	}
	fclose(out);
	fclose(fp);

	if( num_self_algns > 0 ) {
		free(self_algns);
	}

	if( num_pair_algns > 0 ) {
		free(pair_algns);
	}
	free(size2);
	free(size1);
	free(pair_pts);
	free(p_pts);
	free_kd(pair_tree);
	free_kd(tree);
	free(num_pair);
	free(num_algns);
	return(EXIT_SUCCESS);
}
