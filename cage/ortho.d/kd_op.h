#ifndef KD_OP_H
#define KD_OP_H

#include "kd_tree.h"

void print_kd(struct perm_pt *p_pts, struct kdnode *tree);
int find_pred_blk(struct kdnode *p, int x, int y);
int find_successor(struct kdnode *p, int xval, int yval);

#endif /* KD_OP_H */
