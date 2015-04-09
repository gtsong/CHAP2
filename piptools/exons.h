/* $Revision: 1.6 $$Date: 2001/03/19 17:58:51 $ */

#ifndef PIPTOOLS_EXONS
#define PIPTOOLS_EXONS

int get_exons (char *filename, int **Beg, int **End, int **Dir,
	char ***Gene_name, char ***Exon_name, int max_errs);

void free_exons (int *beg, int *end, int *dir, char **gene_name,
	char **exon_name, int n);

char *exon_header (char *filename);

#endif
