#ifndef UTIL_ALGNS_H
#define UTIL_ALGNS_H
#include "main.h"

int search_candi_algns(struct DotList *algns, int num_algns, struct I src, struct I dst, struct DotList *candi_algns);
float pick_sim_level(struct I src, struct I dst, struct DotList *algns, int num_algns, FILE *f);
bool is_in_skip_int(int loc, struct I *list, int num_list);
float cal_pid_part_algn(struct DotList *algns, int algn_id, int left_diff, int right_diff, FILE *f, int sp_id);

#endif /* UTIL_ALGNS_H */
