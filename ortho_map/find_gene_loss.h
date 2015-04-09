#ifndef FIND_GENE_LOSS_H
#define FIND_GENE_LOSS_H
#include "kd_op.h"

void check_gene_loss(int *num_list, struct DotList *dots, int sp_code, int rm_sp, int left_sp, int *cur_num, struct ops_list *ops);
int find_opt_gene_loss(struct DotList *pair_alg, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct gap_list *gp, struct DotList *rm_sp_alg, int num_rm_sp);
int check_inc_algns(struct I reg, struct DotList *dots, int num);
void iden_gene_loss(int num_init_algns, struct DotList *init_algns, int num_ops, struct ops_list *ops, FILE *f, struct ops_list **del_ops, int *num_del_ops, int *num_alloc);
int which_side_in_ops_list(struct DotList algn1, struct DotList algn2, int num_ops, struct ops_list *ops);
void update_init_algns_for_loss(struct ops_list *del_ops, int num_del_ops, int num_init_algns, struct DotList *init_algns, int sp_id, FILE *f);

#endif /* FIND_GENE_LOSS_H */
