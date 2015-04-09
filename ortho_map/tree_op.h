#ifndef TREE_OP_H 
#define TREE_OP_H

#define TREE_PRINT 0
#define INT_PRINT 1

#define NAME 0
#define COLON 1
#define COMMA 2
#define LPAR 3
#define RPAR 4
#define SEMICOLON 5
#define B_LEN 6
#define LBR 7 // left bracket
#define RBR 8 // right bracket

void init_tree(struct p_tree *t);
struct p_tree *read_tree(char *filename);
struct p_tree *read_one_line(char *one_line); 
bool is_sep_char(char ch);
bool is_done_node(struct p_tree *node);
bool is_leaf_node(struct p_tree *node);
void leave_only_taxons(char *one_line);
void free_p_tree(struct p_tree *t);
void mapping_tree(struct p_tree *gt, struct p_tree *sp);
bool is_in_list(int id, int *list, int num_id);
int assign_dmode(struct p_tree *gt, struct p_tree *sp);
struct p_tree *find_lca(int *list, int num_sp, struct p_tree *sp, int *lca);
bool is_covering(int *list, int num_sp, int *cmp, int num_cmp);
bool is_descendant(int nid1, int nid2, struct p_tree *p);
bool is_node_in_tree(int nid, struct p_tree *t);
int get_id_in_list(char *name, struct sp_list *list, int num_sp);
void assign_sp_id(struct p_tree *t, struct sp_list *list, int depth, int num_sp);
void assign_gid(struct p_tree *t, struct g_list *genes, int num_genes);
void get_deepest_sp(struct p_tree *sp, int *sp1, int *sp2, int *out_sp, int *cur_max);
int get_leaf_node(struct p_tree *sp, int dir);
void postorder(struct p_tree *t);
struct p_tree *get_lca_node(struct p_tree *cp, int sp_id);
void make_new_node(struct p_tree *ch_l, struct p_tree *ch_r);
void mark_visited(struct p_tree *c_node);
void print_tree(struct p_tree *t, int mode);
void swap_nodes(struct p_tree *a, struct p_tree *b);
void delete_one_sp(struct p_tree *t, int sp);
void assign_node(struct p_tree *a, struct p_tree *b, int mode);
void write_bk_gname(struct p_tree *gt, struct g_list *genes);
void extract_gtree_3sp(struct p_tree *t, int sp1, int sp2, int out_sp, FILE *f);
bool is_3sp_in_list(int sp1, int sp2, int out_sp, struct p_tree *t);
struct p_tree *add_node(struct p_tree *cur, int id, int num_list, char *name, int sp_code, int sign);
bool is_outgroup(struct p_tree *p, int sp1_code, int sp2_code);
int assign_sp_code(char *line, struct sp_list *sp_code, int num_sp_code);
bool is_within_branch(int out_code, int sp1_code, int sp2_code, struct p_tree *p);

#endif /* TREE_OP_H */
