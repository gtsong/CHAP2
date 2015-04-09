#ifndef CONTIGS_OP_H
#define CONTIGS_OP_H

#define CTG_ASSIGNED 0
#define CTG_NOT_ASSIGNED 1

int add_ctg_name(int id, char *sp_name, char *ctg_name, int len, struct n_pair *pairs, int num_pairs);
int is_ctg_in(char *sp_name, char *ctg_name, struct n_pair *pairs, int num_pairs);
void cal_length_sum(int *len_sum1, struct n_pair *contigs1, int num1);
void concat_ctg_name(char *name, char *sp_name, char *ctg_name);

#endif /* CONTIGS_OP_H */
