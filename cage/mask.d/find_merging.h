#ifndef	FIND_MERGING_H
#define FIND_MERGING_H

#define CHANING 1
#define TRANSITIVE 2

bool find_merging(struct slist *st, struct DotList *dots, int num, int loc);
bool check_candi(struct DotList *dots, int id, int comp_id, bool is_x);
int distance(struct DotList *dots, int loc_id, int comp_id, bool *is_x, int *sd, int flag);
void merging_step(struct DotList *dots, int loc, int comp);
int m_val(struct slist *st, struct DotList *dots, int i);
int compute_closeness(struct DotList *dots, int loc_id, int comp_id);

#endif /* FIND_MERGING_H */
