#ifndef FIND_MERGING_H
#define FIND_MERGING_H

#include "main.h"

#define CHECK_DEL 0
#define CHECK_INS_DUP 1
#define CHECK_GENE_LOSS 2 

bool find_merging(struct slist *st, struct DotList *dots, int num, int loc);
bool check_candi(struct DotList *dots, int id, int comp_id, int flag);
int distance(struct DotList *dots, int loc_id, int comp_id, bool *is_x, int *sd);
void merging_step(struct DotList *dots, int loc, int comp);
int m_val(struct slist *st, struct DotList *dots, int i);
bool is_tandem(struct DotList cur);
int compute_closeness(struct DotList *dots, int loc_id, int comp_id);

#endif /* FIND_MERGING_H */
