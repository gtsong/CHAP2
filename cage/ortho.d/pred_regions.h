#ifndef PRED_REGIONS_H
#define PRED_REGIONS_H

#include "kd_tree.h"

#define DUP_X_TO_Y 1
#define DUP_Y_TO_X 2
#define CON_X_TO_Y 3
#define CON_Y_TO_X 4
#define NEGLECT 5 

#define NO_OVERLAP 0
#define OVERLAP_IN_X 1
#define OVERLAP_IN_Y 2
#define OVERLAP_BOTH 3

#define OVERLAP_LEFT_X 4
#define OVERLAP_RIGHT_X 5
#define OVERLAP_LEFT_Y 6
#define OVERLAP_RIGHT_Y 7

#define TH_ALG_ERR 1

#define TANDEM_CHECK 0
#define GENERAL_CHECK 1

int find_covering_regions(int youngest, bool *is_x, int *regs, int num_list, struct DotList *dots);
int is_left_to_right(int id, int num_list, struct DotList *org);
int is_left_to_right_again(int id, int num_list, struct DotList *dots);
bool check_negligible_line(int id, int num_list, struct DotList *dots);
int is_left_to_right_count(int *num_x, int *num_y, int id, int num_list, struct DotList *org, int t_val, struct ID_List *dlist, int num_dup, int flag, struct kdnode *tree, struct perm_pt *p_pts, int size, FILE *fp, struct DotList *init_algns);
int is_left_to_right_dup_count(int *num_x, int *num_y, int id, int num_list, struct DotList *org, int t_val, int r_id, int l_id);
int is_left_to_right_count_strict(int *num_x, int *num_y, int id, int num_list, struct DotList *org);
void get_elm_list(int *num_id, int *elm_id, int id, bool is_x, int num_list, struct DotList *dots);
bool is_s_list(struct DotList *dots, int ins_id, int cur_id, struct ID_List *dlist, int num_dup, struct kdnode *tree, struct perm_pt *p_pts, int size, FILE *fp, struct DotList *init_dots);

#endif /* PRED_REGIONS_H */
