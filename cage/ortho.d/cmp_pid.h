#ifndef CMP_PID_H
#define CMP_PID_H

#define OP_RATIO_TH 0.1 // threshold for the ratio of the width of overlaps over the entire length

void update_inpar_list(int id, struct slist *alg_id, int num_algns, struct DotList *algns, struct DotList *init_algns, FILE *fp);
int cmp_pid(int cur_id, bool is_x_cur, int cmp_id, bool is_x_cmp, struct DotList *algns, struct DotList *init_algns, FILE *fp);
bool is_cand_cmp(struct I cur_reg, struct I cmp_reg);

#endif /* CMP_PID_H */
