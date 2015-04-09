#ifndef ADJUST_PLOT_GENES_H
#define ADJUST_PLOT_GENES_H
#include "kd_tree.h"

#define W_SID 0 
#define W_FID 1
#define H_SID 2
#define H_FID 3

void adjust_plot_pair_genes(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct r_list *rm_out1, struct r_list *rm_out2, int num_rp1, int num_rp2, struct g_list *genes1, struct g_list *genes2, int num_genes1, int num_genes2, struct DotList *init_dots, FILE *fp);
int find_opt_fr_genes(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct r_list *rp1, int num_rp1, struct r_list *rp2, int num_rp2, struct g_list *genes1, struct g_list *genes2, int num_genes1, int num_genes2, int *rp1_id, int *rp2_id, FILE *fp);
bool is_it_the_same_gene(struct DotList *dots, int id1, int id2, struct g_list *genes1, struct g_list *genes2, int num_genes1, int num_genes2);

#endif /* ADJUST_PLOT_GENES_H */
