#ifndef FIND_DUP_COPY_H
#define FIND_DUP_COPY_H

#include "kd_tree.h"

#define SP_LEFT 0
#define SP_RIGHT 1
#define SP_BOTH 2
#define NONE 3

int find_dup_copy(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size, struct ID_List *dlist, FILE *fp, struct DotList *init_dots);
int find_opt_du_copy(struct DotList *dots, int num_lines, int id, struct perm_pt *st, struct kdnode *tree, int size, int w_sid, int w_fid, int h_sid, int h_fid, int *cid, bool *x_ins, bool *f_is_x, int *t_ins, FILE *fp, struct DotList *init_dots);
int find_id(struct DotList *dots, struct kdnode *tree, int size, int id, int xval, int yval, int option);
int check_inclusion_close_dup(int id, struct DotList *dots, int num_lines, bool *x_ins, bool *t_ins);
bool check_whole_regions_inclusion(struct DotList *dots, int num_lines, int mid, int left_id, int right_id, bool is_x);
int inserted_dup_copy(int id, struct ID_List *dlist, int *num_dup, struct DotList *dots, int num_lines, int mode, int *add_info, struct kdnode *tree, struct perm_pt *p_pts, int size, FILE *fp, struct DotList *init_dots);
bool check_whole_del_region_inclusion(struct DotList *dots, int num_lines, struct I del_reg);
bool tandem_exist(struct DotList *dots, struct perm_pt *p_pts, struct kdnode *tree, int size, int id1, int id2);
int find_id_len(struct kdnode *tree, int size, int len, int xval, int yval, int option);

#endif /* FIND_DUP_COPY_H */
