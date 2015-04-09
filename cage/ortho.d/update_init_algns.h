#ifndef	UPDATE_INIT_ALGNS_H
#define	UPDATE_INIT_ALGNS_H

#define CASE_1 1
#define CASE_2 2
#define CASE_3 3

#define ALL 1 // when rolling-back for all local alignments
#define INS 2 // when rolling-back for an inserted copy
#define DEL 3 // for a deleted copy
#define DUP 4 // for a duplicated copy

void mark_chain(struct DotList *dots, int id1, int id2, struct DotList *init_dots);
void update_init_algn(struct I reg, struct DotList *dots, int id, bool is_x, int cmp_id, struct DotList *init_dots, FILE *fp, int mode);
void adjust_init_algn(struct DotList *init_dots, int id, struct I dup_reg, bool is_x, FILE *fp, int mode, bool dup_is_x);
void rollback_init_dots(struct DotList *dots, int id, bool is_x, int num_algns, struct DotList *init_algns, int num_init_algns, FILE *fp, int num_ops, struct ops_list *ops, int mode, int size);
void write_init_maf(FILE *output_f, struct DotList *init_algns, int num_algns, char *sp1_name, char *sp2_name, int len1, int len2, FILE *fp, int flag);
void write_init_maf_stdout(struct DotList *init_algns, int num_algns, struct n_pair *contigs1, struct n_pair *contigs2, int *len_sum1, int *len_sum2, int num1, int num2, int len1, int len2, FILE *fp, int flag);
bool check_in_cid_list(int index, int *list, int num_list);
void update_init_algn_del(struct DotList *dots, int cmp_id, struct I del_reg, bool is_x, FILE *fp, struct DotList *init_dots);
void update_pid_init_algns(int num_init_algns, struct DotList *init_algns, FILE *f, float avg_pid);

#endif /* UPDATE_INIT_ALGNS_H */
