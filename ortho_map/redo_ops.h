#ifndef REDO_OPS_H
#define REDO_OPS_H

#include "kd_tree.h"

void redo_ops(struct ops_list *ops, int num_ops, int *num_algns, struct DotList *algns, struct cv_list *cv, int *num_cv, int size, struct DotList *init_algns, int num_init_algns, FILE *fp, int size1);
struct slist find_mapping_algn(char ori, struct I src, struct I dst, struct DotList *alg_self, int num_algn_self, int size, struct ID_List *suspend, int *num_suspend_pairs, FILE *fp, struct DotList *init_algns, struct kdnode *tree, struct perm_pt *p_pts);
int search_algns_for_op(char ori, struct I src, struct I dst, struct DotList *algns, int num_algns, bool *is_x);
int make_ops_algns(struct DotList *ops_algns, struct ops_list *ops, int b, int e);

#endif /* REDO_OPS_H */
