#ifndef UTIL_ALGNS_H
#define UTIL_ALGNS_H
#include "main.h"

int search_candi_algns(struct DotList *algns, int num_algns, struct I src, struct I dst, struct DotList *candi_algns);
float pick_sim_level(struct I src, struct I dst, struct DotList *algns, int num_algns, FILE *f);
int pick_sim_level_algn(struct DotList *algns, int algn_id, FILE *f, struct exons_list *list1, int num_list1, struct exons_list *list2, int num_list2, int sp_id);
bool is_in_skip_int(int loc, struct I *list, int num_list);
float cal_pid_part_algn(struct DotList *algns, int algn_id, int left_diff, int right_diff, FILE *f, int sp_id);
void cal_algn_beg_end(struct DotList *algns, int id, struct I ov, int sp_id, FILE *f, int *beg1, int *end1);

#endif /* UTIL_ALGNS_H */
