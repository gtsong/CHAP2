#ifndef REFINE_ALGNS_H
#define REFINE_ALGNS_H

int adjust_all_algns_diff(struct DotList *ortho_algns, int num_ortho_algns, struct DotList *algns, int num_algns);
void undo_events(struct I reg, struct DotList *algns, int num_algns, int mode, FILE *f);
void remove_false_mtom(struct DotList *algns, int num_algns, struct ops_list *ops, int num_ops, FILE *f);
int where_in_ops(struct I src, struct I dst, struct ops_list *ops, int num_ops);
void remove_algn_overlap(struct DotList *algns, int num_algns, int cur_id, int cmp_id, struct I src, struct I dst, struct ops_list * ops, int num_ops, FILE *f);
bool is_only_algn(struct DotList *algns, int id, struct I reg, int num_algns);
int find_loc_both_chained(struct DotList *algns, int id1, int id2, struct I ov, FILE *f, int sp_id, bool *is_first);
int find_loc_both_chained_sp1(char *S1, char *T1, char *S2, char *T2, int b1, int b2, int e1, int e2, bool *is_first);
int find_loc_both_chained_sp2(char *S1, char *T1, char *S2, char *T2, int b1, int b2, int e1, int e2, int sign1, int sign2, bool *is_first);


#endif /* REFINE_ALGNS_H */
