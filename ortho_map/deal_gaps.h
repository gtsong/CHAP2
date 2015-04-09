#ifndef DEAL_GAPS_H
#define DEAL_GAPS_H
#include "main.h"

struct gap_list define_gap_new_type_inc(struct DotList *dots, int loc_id, int comp_id, bool is_x);
struct gap_list define_gap_new_type(struct DotList *dots, int loc_id, int comp_id, bool is_x);
struct gap_list define_gap(struct DotList *dots, int loc_id, int comp_id, int d, int sd, bool is_x);
int decide_gap(struct r_list *rp1, struct r_list *rp2, int n1, int n2, int sign, struct r_list *rp);
bool check_inclusion_alignments(struct gap_list gp, struct DotList *dots, int num);
int check_into_own(struct DotList *dots, int loc_id, int comp_id);
struct gap_list redefine_for_del(struct DotList *dots, struct gap_list gps);
int get_starting_loc(struct DotList *dots, int loc, int comp, bool is_x);
int get_starting_loc_reg(struct DotList aln_1, struct DotList aln_2, bool is_x, int *from, int *to);

#endif /* DEAL_GAPS_H */
