#ifndef ADJUST_PLOT_H
#define ADJUST_PLOT_H
#include "kd_tree.h"

#define W_SID 0 
#define W_FID 1
#define H_SID 2
#define H_FID 3

void adjust_plot_pair(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct r_list *rm_out1, struct r_list *rm_out2, int num_rp1, int num_rp2, struct DotList *init_dots, FILE *fp);
int find_id(struct DotList *dots, struct kdnode *tree, int size, int id, int xval, int yval, int option);
int find_id_pair(struct DotList *dots, struct kdnode *tree, int size1, int size2, int id, int xval, int yval, int option);
int find_opt_fr(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct r_list *rp1, int num_rp1, struct r_list *rp2, int num_rp2, int *rp1_id, int *rp2_id, FILE *fp);

#endif /* ADJUST_PLOT_H */
