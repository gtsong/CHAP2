#ifndef CONTIGS_OP_H
#define CONTIGS_OP_H

#define CTG_ASSIGNED 0
#define CTG_NOT_ASSIGNED 1
#define CTG_NOT_ASSIGNED_BUT_LEN 2

void merge_sep_contigs(struct DotList *algns, int num_algns, struct n_pair **contigs1, int *num_contigs1, int *num_alloc1, int **len_sum1, struct n_pair **contigs2, int *num_contigs2, int *num_alloc2, int **len_sum2);
int add_ctg_name(int id, char *sp_name, char *ctg_name, int len, struct n_pair *pairs, int num_pairs);
int is_ctg_in(char *sp_name, char *ctg_name, struct n_pair *pairs, int num_pairs);
void cal_length_sum_from_ctglist(int *len_sum, struct n_pair *contigs, int num);
void cal_length_sum(int *len_sum1, struct n_pair *contigs1, int num1);
void concat_ctg_name(char *name, char *sp_name, char *ctg_name);
int get_ctg_id(char *sp_name, char *ctg_name, struct n_pair *pairs, int num);
void print_contig_list(FILE *ctg_f, struct n_pair *contigs, int *len_sum, int num_org, int num);
void read_contigs_file(char *name, FILE *f, struct n_pair *contigs, int num_contigs);
int read_one_contig_file(char *fname, FILE *fp, struct n_pair *contigs, int count);
void replace_contigs_len(struct n_pair *contigs, int num);

#endif /* CONTIGS_OP_H */
