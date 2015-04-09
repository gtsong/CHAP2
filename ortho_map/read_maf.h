#ifndef READ_MAF_H
#define READ_MAF_H

int cal_pid(char *s, char *t, int ncol);
float cal_pid_maf(char *s, char *t, int ncol);
float cal_pid_maf_beg(char *s, char *t, int beg, int ncol);
void read_maf(char *fname, int mode, struct DotList *algns, int *num_algns, int *size1, int *size2);
int cat_srcblock(char *token);
int cal_match_maf_beg(char *s, char *t, int beg, int ncol);
void adjust_multi_contig_pos(struct DotList *algns, int num_algns, int *size1, int *size2, struct n_pair *contigs1, int *num_contigs1, struct n_pair *contigs2, int *num_contigs2);
void adjust_algn_pos(struct DotList *algns, int num_algns, struct n_pair *contigs1, int num1, int *size1, struct n_pair *contigs2, int num2, int *size2, int mode);

#endif /* READ_MAF_H */
