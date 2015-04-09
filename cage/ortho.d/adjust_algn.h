#ifndef ADJUST_ALGN_H
#define ADJUST_ALGN_H

void adjust_alignment(struct DotList *dots, int i, struct I to);
int adjust_algn_reg(struct I *reg_x, struct I *reg_y, struct I to, int sign);

#endif /* ADJUST_ALGN_H */
