#ifndef PRED_SP_H
#define PRED_SP_H

void predict_sp_op(int sp_code, int rm_sp, int left_sp, int *num_list, struct DotList *dots, int *cur_num, struct ops_list *ops);
void rollback_sp(struct DotList *dots, int *num_dots, struct DotList *init_dots, int num_init_dots, int cur_num, struct ops_list *ops, int sp_id);

#endif /* PRED_SP_H */
