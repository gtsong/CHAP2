#ifndef ID_ORTHO_H
#define ID_ORTHO_H

#define INC_LEN_TH 100 // a threshold value for checking the inclusion of other alignment region, not in the same length 
#define REMOVED -1
#define NON_RM 0

void obtain_ortho_algns(struct DotList *algns, int num_algns, struct DotList *init_algns, int num_init_algns);
void iden_ortho_update_init(int *num_pair, struct DotList *pair_alg, int num_init_algns, struct DotList *init_algns);
void check_loc_for_ortho(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val);
void check_loc_for_ortho_one_side(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val);
void check_loc_one_side(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val);
void check_pid_for_ortho(struct slist *sorted, int *num_pair, struct DotList *pair_alg, int t_val);
bool check_contain(struct I reg1, struct I reg2, int t_val);
void remove_non_ortho(int *num_pair, struct DotList *pair_alg);
int assign_ortho_algns(struct DotList *ortho_algns, struct DotList *init_algns, int num_init_algns, float *temp_val);
void adjust_boundary_init(struct DotList *init_algns, struct DotList *ortho_algns, int id, int num_init_algns, FILE *f, int sp_id);
void redo_dups_for_mtom(int num_ops, struct ops_list *ops, int num_init_algns, struct DotList *init_algns, FILE *f, int sp_id);
void add_ortho_intervals_dups(struct ops_list *ops, int id, struct slist *sorted, struct DotList *ortho_algns, int num_ortho_algns, int cur_sp_id, int *exc_list, FILE *f, float avg_pid);
void add_ortho_intervals(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, int s_new, int e_new, struct DotList *algns, int sp_id, FILE *f, float avg_pid);
int make_exc_list(struct I query, struct ops_list *ops, int cur_loc, int *exc_list);
//bool is_within_bound(struct I temp, struct short_alist *src_list, int num_list, struct I *aj);
bool is_within_bound(struct ops_list ops, struct I temp, struct short_alist *src_list, int num_list, struct I *aj, struct DotList *algns, int id, int sp_id, FILE *f);
void adjust_boundary(int id, struct DotList *temp_algns, struct DotList *algns, int num_algns, FILE *f);

#endif /* ID_ORTHO_H */
