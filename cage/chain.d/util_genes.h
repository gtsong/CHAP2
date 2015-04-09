#ifndef UTIL_GENES_H
#define UTIL_GENES_H

#define POS_BASE 1
#define PID_BASE 2
#define LEN_BASE 3

int read_genes(FILE *f, struct g_list *genes, char *fname);
void quick_sort_dec_genes(struct g_list *a, int lo, int hi, int mode);
void quick_sort_inc_genes(struct g_list *a, int lo, int hi, int mode);
int quick_search_close_genes(struct g_list *sorted, int i, int j, int query);
struct g_list assign_genes(struct g_list a);
bool is_all_digits(char *item);

#endif /* UTIL_GENES_H */
