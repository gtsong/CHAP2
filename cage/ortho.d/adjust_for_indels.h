#ifndef ADJUST_FOR_INDELS_H
#define ADJUST_FOR_INDELS_H

#define NONE -1
#define DEL 0
#define INS 1

void adjust_for_indels(char op, struct I src, struct I dst, int len, int flag);
void reflect_diff(struct I *origin, struct I *dup, struct I *temp_origin, struct I *temp_dup, int num_ops);

#endif /* ADJUST_FOR_INDELS_H */
