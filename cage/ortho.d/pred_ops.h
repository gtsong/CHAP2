#ifndef PRED_OPS_H
#define PRED_OPS_H

void predict_op(bool is_x, int id, int *num_list, struct DotList *dots, int overlap, int num_ops, struct ops_list *ops);
void pred_dup(int con, char op_ch, int pred_op, bool is_x_to_y, int id, int *num_list, struct DotList *dots, int num_ops, struct ops_list *ops);
void generate_ops(char op, int wide, bool direction, struct I from, struct I to, int flag, int num_ops, struct ops_list *ops, int sp_id);
void generate_ops_del(char op, int location, int wide, int num_ops, struct ops_list *ops);

#endif /* PRED_OPS_H */
