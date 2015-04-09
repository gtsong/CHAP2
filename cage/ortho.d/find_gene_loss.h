#ifndef FIND_GENE_LOSS_H
#define FIND_GENE_LOSS_H
#include "kd_op.h"

void check_gene_loss(int *num_list, struct DotList *dots, int sp_code, int rm_sp, int left_sp, int *cur_num, struct ops_list *ops);
int find_opt_gene_loss(struct DotList *pair_alg, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct gap_list *gp, struct DotList *rm_sp_alg, int num_rm_sp);
int check_inc_algns(struct I reg, struct DotList *dots, int num);

#endif /* FIND_GENE_LOSS_H */
