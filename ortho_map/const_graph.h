#ifndef CONST_GRAPH_H
#define CONST_GRAPH_H
#include "graph_bipar.h"

//void const_graph(int num_algns, struct DotList *algns, FILE *f, struct ops_list *ops, int num_ops, int size);
void const_graph(int num_algns, struct DotList *algns, FILE *f);
void initialize_matrix(struct matrix_ent *init_matrix, int num_left, int num_right);
void initialize_graph(struct bipar_node *graph, int num_nodes);
void initialize_edge_ent(struct edge_ent *edge, int num_edges);
void init_map_in_greedy(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, int num_init_algns, FILE *f, struct bipar_node *graph);
int pick_max_element(struct matrix_ent *init_matrix, int num_left, int num_right, int *id1, int *id2, bool *src_l, bool *src_r, int prev_id1, int prev_id2);
void update_matrix_tmp_val(struct bipar_node *graph, struct matrix_ent *init_matrix, int num_left, int num_right, int l_entry, int r_entry, struct DotList *init_algns, int num_init_algns, FILE *f);
void mark_matched(bool *matched, int num1, int num2, struct matrix_ent *init_matrix, int flag);
void max_bipar_weight_match(struct matrix_ent *init_matrix, int num_left, int num_right);
int shortest_augmenting_path(int v, struct matrix_ent *init_matrix, int num_left, int num_right, int *l_prev, int *r_prev, int *r_dist, bool *r_matched, int *r_weight);  
void pick_max_unvisited(int *l_dist, int *r_dist, int num_left, int num_right, int *id1, int *id2, bool *l_visited, bool *r_visited);
void free_graph(struct bipar_node *graph, int num_nodes);
int cal_inc_weight(struct matrix_ent *init_matrix, int num_left, int num_right, int *l_prev, int *r_prev, int id, int source, struct DotList *init_algns, FILE *f);
int update_init_matrix_val(struct matrix_ent *init_matrix, int l_entry, int r_entry, int num_right,  struct DotList *init_algns, struct DotList *candi_algns, int num_algns, FILE *f);
void update_init_algns_for_ancestral(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, int num_init_algns);
void mark_event_src_nodes(struct bipar_node *graph, int num_left, int num_right, struct ops_list *ops, int num_ops, int size);
void mark_label_on_graph(struct I reg, struct bipar_node *graph, int from, int to, int ops_id);
void init_map_in_dp(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f);
void cal_scores(int **s1, struct point **p1, struct point **temp, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f);
void update_score_in_matrix(int **s1, struct point **p1, struct point **temp, int id1, int id2, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f);
int cal_non_merged_score(struct DotList *algns, int id1, int id2, FILE *f, int score);
void update_score_in_matrix_rev(int **s1, struct point **p1, struct point **temp, int id1, int id2, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f);
void cal_scores_rev(int **s2, struct point **p2, struct point **temp, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f);
int find_max_score(int **s1, struct point **p1, struct point **t1, int b1, int e1, int b2, int e2, int *id1, int *id2, int *end_id1, int *end_id2, int mode);
void combine_scores(int **s1, int **s2, struct point **p1, struct point **p2, struct point **t1, struct point **t2, int **scores, struct point **paths, int num_left, int num_right, int *res_id1, int *res_id2);
void update_matrix_paths(int **s, struct point **p, int start_id1, int start_id2, int prev_id1, int prev_id2, int b1, int e1, int b2, int e2, bool *is_filled1, bool *is_filled2, int num_left, int num_rigth);
bool find_next_bound(struct point **paths, int start_id1, int start_id2, bool *is_filled1, bool *is_filled2, int num_left, int num_right, int *b1, int *e1, int *b2, int *e2, int *prev_id1, int *prev_id2);
int cal_inc_weight_for_path(struct matrix_ent *init_matrix, int num_left, int num_right, int *l_prev, int *r_prev, int id, int source, int run_mode);
void check_unmapped_nodes(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns);
bool is_still_in_the_same_part(int sid1, int sid2, int from, int to, int b1, int e1, int b2, int e2);
bool is_already_visited(int b, int e, struct point **paths, int num_left, int num_right, int s1, int s2, int cur_run);

#endif /* CONST_GRAPH_H */
