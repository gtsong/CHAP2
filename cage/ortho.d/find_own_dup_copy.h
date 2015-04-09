#ifndef FIND_OWN_DUP_COPY_H
#define FIND_OWN_DUP_COPY_H

#include "kd_tree.h"

#define DIFF_PID 8 // the difference of percentage identities of two alignments for duplication into the own copy 

void find_own_dup_copy(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size);
void find_opt_own_dup_copy(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int h_fid);

#endif /* FIND_OWN_DUP_COPY_H */
