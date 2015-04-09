#ifndef CONTIGS_OP_H
#define CONTIGS_OP_H

#define CTG_ASSIGNED 0
#define CTG_NOT_ASSIGNED 1

#define ALLOC_UNIT 10
#define LEN_NAME 100

struct n_pair{
  char name1[LEN_NAME];
  char name2[LEN_NAME];
  int id;
  int len;
};

int add_ctg_name(int id, char *sp_name, char *ctg_name, int len, struct n_pair *pairs, int num_pairs);
int is_ctg_in(char *sp_name, char *ctg_name, struct n_pair *pairs, int num_pairs);
void cal_length_sum(int *len_sum1, struct n_pair *contigs1, int num1);
void concat_ctg_name(char *name, char *sp_name, char *ctg_name);
void init_pair(struct n_pair *list, int b, int e);

#endif /* CONTIGS_OP_H */
