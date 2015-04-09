#ifndef UTIL_OPS_H
#define UTIL_OPS_H

void sort_conv(struct cv_list *a, int num);
void quick_sort_dec_conv(struct cv_list *a, int lo, int hi, int mode);
void quick_sort_inc_conv(struct cv_list *a, int lo, int hi, int mode);
int quick_search_close_conv(struct cv_list *sorted, int i, int j, int query);
struct cv_list assign_conv(struct cv_list a);
struct exons_list assign_exons(struct exons_list a);
struct ops_list assign_ops(struct ops_list a);
void init_conv(struct cv_list *cv, int b, int e);
int count_ops(char *name, FILE *f, int sp_mode);
void read_ops_file(char *name, FILE *f, struct ops_list *ops, int num_ops, int sp_mode, struct n_pair *contigs, int num_contigs, char *ref_name);
int count_final_ops(char *name, FILE *f, int sp_mode);
void read_final_ops_file(char *name, FILE *f, struct ops_list *ops, int num_ops, int sp_mode, struct n_pair *contigs, int num_contigs);
int add_one_ops(struct ops_list *ops, int id, char *buf, struct n_pair *contigs, int num_contigs, int *len_sum, int mode);
void init_ops(struct ops_list *ops, int b, int e);
bool is_on_prev_events(struct I reg, struct ops_list *ops, int from, int to);
void cat_spname(char *name, char *res);

#endif /* UTIL_OPS_H */
