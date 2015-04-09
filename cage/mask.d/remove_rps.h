#ifndef REMOVE_RPS_H
#define REMOVE_RPS_H
#include "main.h"
#include "kd_tree.h"

#define OVERLAP 0
#define NO_OL 1
#define SELF 2

#define X_SIDE 0
#define Y_SIDE 1

void remove_overlapped(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, struct kdnode *pair_tree, struct perm_pt *pair_pts, int size_seq1, int size_seq2, struct DotList *pair_alg, int *num_pair, int sp_flag);
void remove_rps(struct DotList *dots, int *num);
void throw_away_rps(struct DotList *dots, int num_lines, struct IntList self_alg, struct kdnode *tree, struct perm_pt *p_pts, int size, int flag);
void throw_away_rps_pair(struct DotList *dots, int num_lines, struct IntList self_alg, struct kdnode *tree, struct perm_pt *p_pts, int size_seq1, int size_seq2, int sp_flag);
void reduce_rps(struct DotList *dots, int num_lines, struct slist *sorted, int cur, struct kdnode *tree, struct perm_pt *p_pts, int size);

int find_overlapped(struct DotList *dots, int num, struct slist *sorted, int cur, int *clist, int *longest, struct IntList *self_alg);
int find_rps(struct DotList *dots, struct perm_pt *p_pts, int cur, int start, int end, int *clist, int *longest);
int find_clist(struct DotList *dots, int num, struct slist *sorted, int cur, int *clist, int *longest, struct IntList *self_alg, int flag);
int check_longest_one(struct DotList *dots, int *clist, struct IntList *self_alg, int flag, int num_list);
int overlap_len(struct DotList *dots, int i, int id, struct I self);
bool check_contain(int cur, int *clist, int num_list);
void recal_range(struct DotList *dots, int id);
struct I add_intervals(struct DotList *dots, int id, struct I temp);
int find_self_alg(struct DotList *dots, int num_lines, struct slist *sorted, int id, struct IntList *self_alg, struct kdnode *tree, struct perm_pt *p_pts, int size);
bool check_neighbors(struct DotList *dots, int id, struct perm_pt *list, int num, struct I self_alg, int flag);
void remove_odd_algns(struct DotList *dots, int num, struct slist *sorted, int cur_id);
bool is_same_contigs(struct DotList *dots, int id1, int id2);

#endif /* REMOVE_RPS_H */
