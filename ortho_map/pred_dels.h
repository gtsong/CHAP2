#ifndef PRED_DELS_H
#define PRED_DELS_H
#include "kd_op.h"

#define DEL_OVERLAP_TH 80
#define DEL_ID_TH 4

int pred_dels(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct gap_list *del_list, int pid);
int find_opt_del(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct gap_list *gp);
int assign_gap_list(struct gap_list *del_list, int loc, int last_gid, struct gap_list gp);
void check_overlapped_del(struct DotList *dots, struct gap_list *del_list, int num_del);
bool check_inv_same_del(struct gap_list *del_list, int num_del, int cur, int del_id);
void remove_dif_del(struct gap_list *del_list, int num_del, int last_gid, int pid, struct DotList *dots);
int remove_false_dels(struct DotList *dots, struct gap_list *del_list, int num_del, struct ID_List *dlist, int num_dup);
int find_redo_del_list(struct I src, struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct gap_list *del_list);

#endif /* PRED_DELS_H */
