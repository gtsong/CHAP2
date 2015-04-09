#ifndef SEC_ROUND_H
#define SEC_ROUND_H

#include "kd_tree.h"

void reorder_dups(struct slist *sorted, int *num_alg_self, struct DotList *alg_self, int size, int threshold, struct DotList *init_algns, int num_init_algns, FILE *fp, float avg_pid, struct cv_list *cv, int *num_cv, struct ops_list *ops, int *num_ops, struct exons_list *exons, int *num_exons, struct exons_list *gene, int *num_genes, int *old_dups, int num_old_dups, int run_mode, int size1);
int choose_candi_algns_pid(struct DotList *candi_algns, struct slist *sorted, struct DotList *algns, int num_algns, struct DotList *pair_algns, int *num_pair);
int choose_candi_algns_cmp_pid(struct DotList *candi_algns, struct slist *sorted, struct DotList *algns, int num_algns, struct DotList *pair_algns, int *num_pair, int *old_dups, int num_old_dups);
bool is_overlap_ortho_algns(struct DotList cur_algn, bool is_x,  struct DotList *pair_algns, int num_pair, bool *is_exist, bool is_in_old_dup);
struct slist latest_dup(struct DotList *candi_algns, int num_candi, int size, int threshold, int *num_suspend_pairs, struct ID_List *suspend, int mode, struct kdnode *tree, struct perm_pt *p_pts, struct DotList *alg_self, int num_algn_self, struct DotList *init_algns, int num_init_algns, FILE *fp, struct DotList *pair_algns, int num_pair, struct exons_list *exons, int num_exons, struct exons_list *genes, int num_genes, float avg_pid, struct ops_list *ops, int num_ops);
float cal_cover_rate(struct I reg, struct slist *st, int b, int e, struct DotList *ortho, int *pid, int sp_id);
int cmp_ortho_mappings(struct DotList cur_algn, struct DotList *ortho, int num_ortho);
bool is_in(int id, int *list, int num_list);
int check_distance_prev_dup(int sign, int id, struct DotList *init_algns, struct ops_list *ops, int num_ops);
void update_algns_sign(struct DotList *alg_self, int id, int num_alg_self, struct DotList *init_algns);
struct slist replace_by_better(struct DotList *alg_self, struct slist res_list, struct slist *op_info, int num_op, int pid_th, int num_overlap);
int check_prev_tandem_dup(struct DotList *init_algns, int num_init_algns, int id, struct ops_list *ops, int num_ops); 
#endif /* SEC_ROUND_H */
