#ifndef HANDLE_TANDEM_DUP_H
#define HANDLE_TANDEM_DUP_H

// #define TD_PID 5

void handle_tandem_dup(struct DotList *dots, int *num, struct DotList *init_dots);
void convert_tandem_region(struct DotList *dots, int num, int id, int *t_list, int num_tandem); 
void conv_td_reg(struct DotList *dots, int num, int id, int *t_list, int num_tandem, struct DotList *init_dots, int flag, int *val1, int *val2, int *val_org);
int find_tandem_list(struct DotList *dots, struct slist *st, int id, int num, int *t_list);
int check_tandem_reg(struct DotList t, struct DotList *dots, int num);
void adjust_init_offset(struct DotList *init_algns, int init_id, struct DotList t1, struct DotList *algns, int cur_id);

#endif /* HANDLE_TANDEM_DUP_H */
