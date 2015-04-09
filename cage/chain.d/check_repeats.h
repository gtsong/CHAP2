#ifndef CHECK_REPEATS_H
#define CHECK_REPEATS_H
#include "main.h"

#define Max_Rps 10

bool check_repeats(struct DotList *dots, int loc_id, int comp_id, int dist, bool is_x);
float compute_portion_rps(int from, int to, struct r_list *rp, int n_rps);
int obtain_repeats_list(int from, int to, struct r_list *rp, struct r_list *rp_list, int num_list);
int get_starting_loc(struct DotList *dots, int loc, int comp, bool is_x);
int count_lower(char *buf, int len);
int sum_len_rps(struct r_list *rp, int n_rps);

#endif /* CHECK_REPEATS_H */
