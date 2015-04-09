#ifndef ID_INPAR_H
#define ID_INPAR_H

#include "kd_tree.h"

#define BEFORE_SP 1
#define AFTER_SP 2
struct slist iden_inpar(int *num_alg_self, struct DotList *alg_self, int size, int threshold, int *num_dup_copy, struct ID_List *dlist, int mode, struct DotList *init_algns, FILE *fp);
struct slist find_opt_inpar(int num_alg_self, struct DotList *alg_self, int size, int threshold, int *num_dup_copy, struct ID_List *dlist, int mode, struct kdnode *tree, struct perm_pt *p_pts, struct DotList *init_algns, FILE *fp);

#endif /* ID_INPAR_H */
