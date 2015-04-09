#ifndef MAP_ALGNS_H
#define MAP_ALGNS_H

void map_one_to_one(int num_init_algns, struct DotList *init_algns, FILE *f);
void cut_part_algn(struct DotList *init_algns, int cur_id, int cmp_id, int mode, FILE *f);

#endif /* MAP_ALGNS_H */
