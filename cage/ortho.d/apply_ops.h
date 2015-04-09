#ifndef APPLY_OPS_H
#define APPLY_OPS_H

void print_conv_on_dup(int avg_pid, struct DotList alg, struct cv_list *cv, int num_cv, int *num_ops, struct ops_list *ops, int run_mode);
int update_conv(struct DotList *algns, int rm_id, bool is_x, struct cv_list *cv, int num_cv);
int update_conv_del(struct I reg, struct cv_list *cv, int num_cv);
int update_exons(struct DotList *algns, int rm_id, bool is_x, struct exons_list *exons, int num_exon);
int update_exons_del(struct I reg, struct exons_list *exons, int num_exons);
int check_exons_bound(struct DotList algn, struct exons_list *exons, int num_exons, struct exons_list *genes, int num_genes);
int cal_cur_pos_ops(int num_ops, struct ops_list *ops, struct ops_list *new_ops, int sp_id, int seq1_len);
void update_cur_pos_ops(struct ops_list *ops, int *num_ops, int size1);

#endif /* APPLY_OPS_H */
