#ifndef EXTEND_SLIST_H
#define EXTEND_SLIST_H

#include "kd_tree.h"

void extend_slist(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct ID_List *slist, int num_dup, int res_id, bool is_x, FILE *fp, struct DotList *init_dots);
int find_gap_for_ins(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, int res_id, bool res_is_x, bool *f_is_x, FILE *fp, struct DotList *init_dots);
bool in_slist(int l_id, int r_id, int cur_index, int num_dup, struct ID_List *slist, struct DotList *dots);
int find_alt_ins_id(int id, struct DotList *dots, struct kdnode *tree, struct perm_pt *p_pts, int res_id, bool is_x, bool *f_is_x, int size, FILE *fp, struct DotList *init_dots);

#endif /* EXTEND_SLIST_H */
