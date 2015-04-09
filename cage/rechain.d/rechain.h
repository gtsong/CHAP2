#ifndef RECHAIN_H
#define RECHAIN_H

void merge_overlaps(struct slist *st, struct DotList *init_algns, int num_algns);
void restore_chained_algns(struct slist *st, struct DotList *init_algns, int num_algns, char *species, char *species2, int len1, int len2, FILE *fp);

#endif /* RECHAIN_H */
