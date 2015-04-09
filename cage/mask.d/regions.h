#ifndef REGIONS_H
#define REGIONS_H

bool strict_subset(struct I reg1, struct I reg2);
bool proper_in(int r, struct I reg);
bool proper_overlap(struct I reg1, struct I reg2);
bool almost_subset(struct I reg1, struct I reg2);
bool tight_subset(struct I reg1, struct I reg2);
bool d_tight_subset(struct I reg1, struct I reg2);
bool s_tight_subset(struct I reg1, struct I reg2);
bool too_loosen_subset(struct I reg1, struct I reg2);
bool loosen_subset(struct I reg1, struct I reg2);
void init_array(int *array, int num);
void cut_off(int *num, struct DotList *dots);
void overwrite_dots(int *num, struct DotList *dots);
int compute_distance(struct I x1, struct I y1, struct I x2, struct I y2, int sign);

#endif /* REGIONS_H */
