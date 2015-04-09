#include "main.h"
#include "util.h"
#include "util_input.h"
#include "contigs_op.h"

char S[BIG], T[BIG];
int debug_mode;

void write_only_name(char *gene_name, char *temp_name);

int main(int argc, char **argv)
{
	struct exons_list *exons, *genes;
	int *num_exons, *num_genes;
	int count = 0, i = 0;
	char temp_name[LEN_NAME];
	char temp_num[LEN_NAME];
	int num_ps = 1;
	char *ptr;
	struct n_pair *contigs;
	int num_contigs = 0;
	int *len_sum;
	FILE *f;

	ptr = NULL;
	debug_mode = FALSE;

	if( argc != 4 ) {
	  fatal("args: $(species).codex $(species_name) contigs.list\n");
  }

	num_exons = (int *) ckalloc(sizeof(int));
	num_genes = (int *) ckalloc(sizeof(int));
	*num_exons = 0;
	*num_genes = 0;

  f = ckopen(argv[3], "r");
  num_contigs = count_contigs(argv[2], f);

  if( num_contigs > 0 ) {
    contigs = (struct n_pair *) ckalloc(num_contigs * sizeof(struct n_pair));
    len_sum = (int *) ckalloc(num_contigs * sizeof(int));
    read_contigs_file(argv[2], f, contigs, num_contigs);
    cal_length_sum_from_ctglist(len_sum, contigs, num_contigs);
  }
	fclose(f);

	count = count_genes(argv[1], argv[2], argv[2]);
  genes = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

  if( count > 0 ) {
    initialize_exons_list(genes, 0, count);
  }

  count = count_exons(argv[1], argv[2], argv[2]);
  exons = (struct exons_list *) ckalloc(count * sizeof(struct exons_list));

  if( count > 0 ) {
    initialize_exons_list(exons, 0, count);
    read_only_exons(exons, num_exons, genes, num_genes, argv[1], argv[2], argv[2], contigs, contigs, num_contigs, num_contigs);
    for( i = 0; i < *num_exons; i++ ) adjust_pos_exons(exons, i, len_sum, num_contigs);
    for( i = 0; i < *num_genes; i++ ) adjust_pos_exons(genes, i, len_sum, num_contigs);
  }

  if( count != (*num_exons) ) {
    fatalf("counting error: %d vs. %d\n", count, *num_exons);
  }

	if( (*num_genes) > 0 ) {
		printf("gene_name");
		for( i = 0; i < (*num_genes); i++ ) {
			if( strstr(genes[i].name, "_ps") != NULL ) {
				write_only_name(genes[i].name, temp_name);
				ptr = strstr(temp_name, argv[2]);
				if( ptr != NULL ) {
					ptr = ptr + strlen(argv[2]);
					*ptr = '\0';
				}
				numtostr(num_ps, temp_num);
				strcat(temp_name, temp_num);
				num_ps++;

				printf("\t%s", temp_name);
			}	
			else {
				printf("\t%s", genes[i].name);
			}
		}
		printf("\n");
		printf("%s", argv[2]);
		for( i = 0; i < (*num_genes); i++ ) {
			if( strstr(genes[i].name, "_ps") != NULL ) {
				printf("\t\\%d;", i+1);
			}
			else {
				printf("\t%d;", i+1);
			}
		}
		printf("\n");

		if( (*num_exons) > 0 ) {
			free(exons);
		}

		if( (*num_genes) > 0 ) {
			free(genes);
		}
	}

	free(num_exons);
	free(num_genes);
	return(EXIT_SUCCESS);
}

void write_only_name(char *gene_name, char *temp_name)
{
	int len = 0;

	len = strlen(gene_name);
	strcpy(temp_name, gene_name);

	if( strstr(gene_name, "_ps") != NULL ) {
		temp_name[len-3] = '\0';
	}
}
