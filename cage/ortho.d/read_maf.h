#ifndef READ_MAF_H
#define READ_MAF_H

#define CTG_ASSIGNED 0
#define CTG_NOT_ASSIGNED 1

int cal_pid(char *s, char *t, int ncol);
float cal_pid_maf(char *s, char *t, int ncol);
float cal_pid_maf_beg(char *s, char *t, int beg, int ncol);
void read_maf(char *fname, int mode, struct DotList *algns, int *num_algns, int *size1, int *size2);
void adjust_multi_contig_pos(struct DotList *algns, int num_algns, int *size1, int *size2, struct n_pair *contigs1, int *num_contigs1, struct n_pair *contigs2, int *num_contigs2);
void adjust_algn_pos(struct DotList *algns, int num_algns, struct n_pair *contigs1, int num1, int *size1, struct n_pair *contigs2, int num2, int *size2, int mode);
int add_ctg_name(int id, char *sp_name, char *ctg_name, int len, struct n_pair *pairs, int num_pairs);
int is_ctg_in(char *sp_name, char *ctg_name, struct n_pair *pairs, int num_pairs);
void cal_length_sum(int *len_sum1, struct n_pair *contigs1, int num1);
void concat_ctg_name(char *name, char *sp_name, char *ctg_name);

#endif /* READ_MAF_H */
