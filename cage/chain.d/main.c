#include "main.h"
#include "regions.h"
#include "read_maf.h"
#include "write_maf.h"
#include "adjust_plot.h"
#include "adjust_plot_genes.h"
#include "find_merging.h"
#include "kd_op.h"
#include "util.h"
#include "util_genes.h"
#include "util_gen.h"

char S[BIG], T[BIG];
char S1[BIG], T1[BIG];
int debug_mode;
int self_pair;

int main(int argc, char **argv)
{
	struct DotList *pair_algns;
	struct DotList *init_algns;
	struct kdnode *tree; // k-d tree based on the left point of each alignment
	struct perm_pt *p_pts; // permutation list sorted by k-d tree based on the left point
	int *size1, *size2; 
	int *num_algns;
	int *num_init;
	int *threshold;
	int i;
	int count = 0;
	FILE *fp, *f;
	char *status;
	char species[100], species2[100], name[100];
	int b1, b2, e1, e2;
  char len1[100], len2[100];
  char strand[100];
	struct r_list *rp1, *rp2;
	char buf[500];
	float rate;
	int a, b;
	char type[100];
	int num_rp1 = 0, num_rp2 = 0;
	int num_rp_count1 = 0, num_rp_count2 = 0;
	struct g_list *genes1, *genes2;
	int num_genes1 = 0, num_genes2 = 0;
	
	debug_mode = FALSE;
	if((argc == 6) && (strcmp(argv[5], "debug-mode") == 0)) {
		debug_mode = TRUE;
	}
	else if((argc == 7) || (argc == 8)) {
		if((argc == 8) && (strcmp(argv[7], "debug-mode") == 0)) debug_mode = TRUE;
	  f = fopen(argv[5], "r");
 		while(fgets(buf, 1000, f))
	  {
    	if( buf[0] == '#' ) {}
    	else if((buf[0] == '>') || (buf[0] == '<')) num_genes1++;
  	}
  	genes1 = (struct g_list *) ckalloc(sizeof(struct g_list) * num_genes1);
		fseek(f, 0, SEEK_SET);
		num_genes1 = read_genes(f, genes1, argv[5]);
		fclose(f);

	  f = fopen(argv[6], "r");
 		while(fgets(buf, 1000, f))
	  {
    	if( buf[0] == '#' ) {}
    	else if((buf[0] == '>') || (buf[0] == '<')) num_genes2++;
  	}
  	genes2 = (struct g_list *) ckalloc(sizeof(struct g_list) * num_genes2);
		fseek(f, 0, SEEK_SET);
		num_genes2 = read_genes(f, genes2, argv[6]);
		fclose(f);
	}
	else if(argc !=5 )
	{
		fatal("args: pair_maf output rm-out_file1 rm-out_file2\n");
	}

	num_algns = (int *) ckalloc(sizeof(int));
	num_init = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));
	threshold = (int *) ckalloc(sizeof(int));

//	tree = (struct kdnode *) ckalloc(sizeof(struct kdnode));

	if(strcmp(argv[3], argv[4]) == 0) {
		self_pair = SELF;
	}
	else self_pair = PAIR;

  fp = ckopen(argv[1], "r");
  count = 0;
  if( ((status = fgets(S, BIG, fp)) == NULL) || (strncmp(S, "##maf version", 13)))
    fatalf("%s is not a maf file", argv[1]);
  while((status != NULL) && (S[0] == '#')) {
    status = fgets(S, BIG, fp);
	}
//      fatalf("no alignments in %s", argv[1]);

  while((status != NULL) && (strstr(S, "eof") == NULL)) {
    if (strncmp(S, "a ", 2)) fatalf("expecting an a-line in %s, saw %s", argv[1], S);
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
  }
  fclose(fp);

	if( count > 0 ) {
  	pair_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
  	init_algns = (struct DotList *) ckalloc(count * sizeof(struct DotList));
		initialize_algns(pair_algns, count);
		initialize_algns(init_algns, count);

 		read_maf(argv[1], G_MODE, pair_algns, num_algns, size1, size2);
	  read_maf(argv[1], G_MODE, init_algns, num_init, size1, size2);
		p_pts = (struct perm_pt *) ckalloc((*num_algns) * sizeof(struct perm_pt));
		assign_perm(p_pts, (*num_algns), pair_algns, LEFT);
		tree = build_kd(p_pts, 0, (*num_algns)-1);

	  f = fopen(argv[3], "r");
	  for( i = 0; i < 3; i++ ) fgets(buf, 500, f);
 		num_rp_count1 = 0;
	 	while(fgets(buf, 500, f) != NULL) num_rp_count1++;

		if( num_rp_count1 > 0 ) {
			rp1 = (struct r_list *) ckalloc(num_rp_count1 * sizeof(struct r_list));
			init_rlist(rp1, 0, num_rp_count1-1);

 			fseek(f, 0, SEEK_SET);
 			for( i = 0; i < 3; i++ ) fgets(buf, 500, f);
 			i = 0;
	 		while(fgets(buf, 500, f) != NULL) {
	 	  	sscanf(buf, "%*s %f %*s %*s %s %d %d %*s %*s %*s %s %*s", &rate, name, &a, &b, type);
				if( (strstr(type, "Simple_repeat") != NULL) || (strstr(type, "Low_complexity") != NULL) ) {
					if( (i >= 1) && (strcmp(rp1[i-1].name, name) == 0) && (((b - a) <= 100) || ((rp1[i-1].end - rp1[i-1].start) < 100)) && (rp1[i-1].d_rate == 0)) {
						rp1[i-1].end = b;
					}
					else {
						rp1[i].d_rate = (float)0;
 			 	  	rp1[i].start = a;
 		  			rp1[i].end = b;
						rp1[i].id = i;
						strcpy(rp1[i].name, name);
						i++;
					}
				}
 			  else {
					rp1[i].d_rate = rate;
		 	 	  rp1[i].start = a;
 		  		rp1[i].end = b;
					rp1[i].id = i;
					strcpy(rp1[i].name, name);
					i++;
				}
 		 }
 		 fclose(f);
		}

		num_rp1 = i;
	
 		f = fopen(argv[4], "r");
 		for( i = 0; i < 3; i++ ) fgets(buf, 500, f);
	  num_rp_count2 = 0;
 		while(fgets(buf, 500, f) != NULL) num_rp_count2++;

		if( num_rp_count2 > 0 ) {
		  rp2 = (struct r_list *) ckalloc(num_rp_count2 * sizeof(struct r_list));
			init_rlist(rp2, 0, num_rp_count2-1);

			fseek(f, 0, SEEK_SET);
			for( i = 0; i < 3; i++ ) fgets(buf, 500, f);
			i = 0;
			while(fgets(buf, 500, f) != NULL) {
				sscanf(buf, "%*s %f %*s %*s %s %d %d %*s %*s %*s %s %*s", &rate, name, &a, &b, type);
				if( (strstr(type, "Simple_repeat") != NULL) || (strstr(type, "Low_complexity") != NULL) ) {
					if( (i >= 1) && (strcmp(rp2[i-1].name, name) == 0) && (((b - a) <= 100) || ((rp2[i-1].end - rp2[i-1].start) < 100)) && (rp2[i-1].d_rate == 0)) {
						rp2[i-1].end = b;
					}
					else {
						rp2[i].d_rate = (float)0;
 			 	  	rp2[i].start = a;
 		  			rp2[i].end = b;
						rp2[i].id = i;
						strcpy(rp2[i].name, name);
						i++;
					}
				}
 			  else {
					rp2[i].d_rate = rate;
		 	 	  rp2[i].start = a;
 		  		rp2[i].end = b;
					rp2[i].id = i;
					strcpy(rp2[i].name, name);
					i++;
				}
 			}
 			fclose(f);
		}

		num_rp2 = i;

 		for( i = 0; i < *num_init; i++ ) init_algns[i].c_id = -1;

		f = fopen(argv[1], "r");

		if( (argc == 7) || (argc == 8) ) {
			if( (*num_algns) > 0 ) {
				adjust_plot_pair_genes(pair_algns, num_algns, tree, p_pts, *size1, *size2, rp1, rp2, num_rp1, num_rp2, genes1, genes2, num_genes1, num_genes2, init_algns, f);
			}
		}
		else {
			if( (*num_algns) > 0 ) {
				adjust_plot_pair(pair_algns, num_algns, tree, p_pts, *size1, *size2, rp1, rp2, num_rp1, num_rp2, init_algns, f);
			}
		}

		if( (*num_init) > 0 ) {
			write_maf(argv[2], init_algns, *num_init, rp1, rp2, *size1, *size2, f); 
		}
		else {
			f = fopen(argv[2], "w");
			fprintf(f, "##maf version=1 scoring=lastz-pid\n");
			fclose(f);
		}

		if( (argc == 7) || (argc == 8) ) {
			free(genes1);
			free(genes2);
		}

		fclose(f);
		if( num_rp_count1 > 0 ) {
			free(rp1);
		}

		if( num_rp_count2 > 0 ) {
			free(rp2);
		}
		free(p_pts);
		free(pair_algns);
		free(init_algns);
		free_kd(tree);
	}
	else {
		f = fopen(argv[2], "w");
		fprintf(f, "##maf version=1 scoring=lastz-pid\n");
		fclose(f);
	}

	free(num_init);
	free(num_algns);
	free(size1);
	free(size2);
	free(threshold);
	return(EXIT_SUCCESS);
}
