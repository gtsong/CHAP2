#ifndef UTIL_GEN_H
#define UTIL_GEN_H
#include "main.h"
#include "kd_tree.h"

void add_symmetric_points(struct DotList *s_dots, struct DotList *dots, int num);
void sort_list(struct slist *st, struct DotList *dots, int num);
void sort_by_yintercept(struct slist *st, struct DotList *dots, int num);
void quick_sort_dec(struct slist *a, int lo, int hi);
void quick_sort_inc(struct slist *a, int lo, int hi);
void quick_sort_inc_int(struct I *a, int lo, int hi);
struct slist assign_slist(struct slist a);
void print_sorted(struct slist *a, struct DotList *dots, int num);
void sort_by_pid(struct slist *st, struct DotList *dots, int num);
void sort_by_width(struct slist *st, struct DotList *dots, int num);
void quick_sort_plist_x(struct perm_pt *a, int lo, int hi);
void quick_sort_plist_y(struct perm_pt *a, int lo, int hi);
struct perm_pt assign_pm_val(struct perm_pt a);

#endif /* UTIL_GEN_H */
