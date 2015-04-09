#ifndef ID_ORTHO_CONV_H
#define ID_ORTHO_CONV_H

extern int debug_mode;

int redo_dups_for_mtom_inc_conv(int num_ops, struct ops_list *ops, int num_init_algns, struct DotList **org_init_algns, FILE *f, int sp_id, int *num_alloc_blocks);
int untag_ortho_intervals_conv(struct ops_list *ops, int id, struct slist *sorted, struct DotList **ortho_algns, int num_ortho_algns, int cur_sp_id, FILE *f, struct DotList **org_init_algns, int num_algns, int *num_alloc_blocks, int *num_blocks, float avg_pid, bool *is_untagged);
int untag_ortho_intervals(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, struct DotList **ortho_algns, int num_ortho, int sp_id, FILE *f, struct DotList **init_algns, int num_init_algns, int *num_alloc_blocks, int *num_blocks, bool *is_untagged);
bool check_source_existence(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, int s_new, int e_new, struct DotList *algns, int sp_id, FILE *f, float avg_pid, bool *skip_untagging);
int chop_match_regions(struct ops_list *ops, int cur_id, struct slist *sorted, int s_loc, int e_loc, struct DotList *algns, int sp_id, FILE *f, struct short_alist *src_list, bool is_source, bool for_untagging);
bool check_later_dup_existence(struct ops_list *ops, int id, int num_ops);
void find_overlapping_ends(struct I ops_reg, struct I cmp, int sp_id, struct DotList *algns, int i, FILE *f, int *res_b, int *res_e);
float cal_ortho_algn_pid(struct I reg, int sp_id, struct slist *st, struct DotList *init_algns, int num_algns, FILE *f);

#endif /* ID_ORTHO_CONV_H */
