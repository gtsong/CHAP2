#ifndef ADJUST_OPS_H
#define ADJUST_OPS_H

bool is_already_in(struct ops_list tmp_ops, struct ops_list *ops, int num_ops);
//int adjust_ops_list(struct ops_list *cur_ops, int num_cur_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list, struct DotList **ortho_algns, int *num_algns, struct slist **sorted);
//bool do_all_species_agree(struct ops_list tmp_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list, struct DotList **ortho_algns, int *num_algns, struct slist **sorted, bool *is_dir_same);
int adjust_ops_list(struct ops_list *cur_ops, int num_cur_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list);
bool do_all_species_agree(struct ops_list tmp_ops, struct ops_list **ops, int *num_ops, int *temp_list, int num_list, bool *is_dir_same);
bool check_one_to_one_mapping(int pid, struct I src, struct I dst, struct DotList *algns, int num_algns, struct slist *sorted);
void adjust_direction(struct ops_list *ops, int num_ops);

#endif /* ADJUST_OPS_H */
