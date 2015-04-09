#ifndef CHECK_COPY_H
#define CHECK_COPY_H

#include "main.h"

int get_score_copy(struct DotList *dots, int num_lines, struct gap_list gps, int *cid, bool *x_ins);
int check_insertion_copy(struct DotList *dots, int num_lines, struct gap_list gp, int *cid, bool *x_ins);
bool check_boundary(struct I x, int from, int to);
bool check_again_equal(struct I cur, struct gap_list gp, struct DotList *dots);
int check_match_gap_aln(struct DotList *dots, struct gap_list gp, int res_id, bool is_x);
bool is_back_to_td(struct DotList *dots, struct gap_list gps);

#endif /* CHECK_COPY_H */
