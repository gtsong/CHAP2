#ifndef UTIL_GEN_H
#define UTIL_GEN_H
#include "main.h"
#include "util_i.h"
#include "adjust_plot.h"

void adjust_after_merging(struct gap_list *gps, int num_gaps, int pt);
void add_symmetric_points(struct DotList *s_dots, struct DotList *dots, int num);
void sort_list(struct slist *st, struct DotList *dots, int num);
void sort_by_yintercept(struct slist *st, struct DotList *dots, int num);
void quick_sort_dec(struct slist *a, int lo, int hi);
void quick_sort_inc(struct slist *a, int lo, int hi);
void quick_sort_inc_int(struct I *a, int lo, int hi);
struct slist assign_slist(struct slist a);
void quick_sort_plist_x(struct perm_pt *a, int lo, int hi);
void quick_sort_plist_y(struct perm_pt *a, int lo, int hi);
struct perm_pt assign_pm_val(struct perm_pt a);
void print_sorted(struct slist *a, struct DotList *dots, int num);
void print_sorted_plist(struct perm_pt *a, int l, int u, int cutdim);
void initialize_algns(struct DotList *temp, int count);
void sort_init_algns(struct slist *st, struct DotList *dots, int num, int mode);
int search_range_b(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode);
int search_range_e(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode);
void assign_algn(struct DotList *temp, int loc, struct DotList cur);
struct gap_list assign_glist(struct gap_list a);
void init_rlist(struct r_list *a, int b, int e);

#endif /* UTIL_GEN_H */
