#ifndef DEAL_GAPS_H
#define DEAL_GAPS_H
#include "main.h"
#include "adjust_plot.h"

struct gap_list define_gap_new_type_inc(struct DotList *dots, int loc_id, int comp_id, bool is_x);
struct gap_list define_gap_new_type(struct DotList *dots, int loc_id, int comp_id, bool is_x);
struct gap_list define_gap(struct DotList *dots, int loc_id, int comp_id, int d, int sd, bool is_x);
bool check_inclusion_alignments(struct gap_list gp, struct DotList *dots, int num);
bool check_insertion_own_aln(int id1, int id2, struct DotList *dots, int num);

#endif /* DEAL_GAPS_H */
