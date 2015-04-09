#ifndef INS_DUP_COPY_H
#define INS_DUP_COPY_H

#define NUM_OF_UNIT 3

#define UNSUSPENDED 0
#define INSERTED_IN_X 1
#define INSERTED_IN_Y 2
#define SP_OVERLAP_IN_X 3
#define SP_OVERLAP_IN_Y 4
#define SUSPENDED 5
#define COVER_INS_SP_IN_X 6
#define COVER_INS_SP_IN_Y 7
#define COVER_INS_SP_BOTH 15
#define INSERTED_IN_BOTH 8
#define INS_IN_X_TAN 9
#define INS_IN_Y_TAN 10
#define SP_OVERLAP_IN_X_L 11
#define SP_OVERLAP_IN_X_R 12
#define SP_OVERLAP_IN_Y_L 13
#define SP_OVERLAP_IN_Y_R 14

int check_ins_dup_copy(int id, struct ID_List *dlist, int num_dup_copy, struct DotList *dots, int mode, int *add_info);
void merge_ins_dup(int dup_copy, int id, struct DotList *dots, struct ID_List *dlist, int num_dup_copy, int *ch_index, int sp1, int sp2);
void rollback_ins_dup(int dup_copy, int id, struct DotList *dots, struct ID_List *dlist, int num_dup_copy, struct DotList *init_dots, FILE *fp);
void merge_two_algns(struct DotList *dots, int left_id, int right_id, bool is_x, struct I cur, int dup_copy);
bool check_later_dup(struct ID_List clist, int id, struct DotList *dots);
int check_cover_ins_sp(struct ID_List clist, int id, struct DotList *dots);

#endif /* INS_DUP_COPY_H */
