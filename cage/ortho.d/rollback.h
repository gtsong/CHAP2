#ifndef ROLLBACK_H
#define ROLLBACK_H

#define PAIR_1 0
#define PAIR_2 1
#define SELF_1 2
#define SELF_2 3

void shift_regions(int con, struct I young, int *num, struct DotList *dots);
void shift_regions_pair(bool is_x, int con, struct I young, int *num, struct DotList *dots);
int rollback_step_dup_no_overlap(bool is_x, int id, int *num, struct DotList *dots);
int rollback_step_pairwise(bool is_x, struct I reg, int *num, struct DotList *dots);
int rollback_step_dup_overlap(bool is_x, int id, int *num, struct DotList *dots);
int rollback_step_conversion(bool is_x, int id, int *num, struct DotList *dots);
int rollback_step_del(int cur, int num_del, struct gap_list *del_list, int *num, struct DotList *dots, int num_ops, struct ops_list *ops);
struct I redefine_del_list(int *extra, int *cur_count, int *del_loc, int num_del, struct gap_list *del_list, int count, struct DotList *dots);
void chop_off(int mode, struct I seg, int num, struct DotList *dots);
void chop_off_self(bool x_or_y, int id, int num, struct DotList *dots);

#endif /* ROLLBACK_H */
