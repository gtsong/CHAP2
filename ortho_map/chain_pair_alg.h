#ifndef CHAIN_PAIR_ALG_H
#define CHAIN_PAIR_ALG_H
#include "kd_tree.h"

void chain_pair_alg_update_init(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct DotList *init_algns);
int chain_pair_alg(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct chain_list *clist);
int find_opt_ch_alg(struct DotList *dots, int num_lines, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid);

#endif /* CHAIN_PAIR_ALG_H */
