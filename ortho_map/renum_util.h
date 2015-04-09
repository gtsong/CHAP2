#ifndef RENUM_UTIL_H
#define RENUM_UTIL_H

void read_view_file(char *fname, struct mapping_list *ortho_maps);
int count_id(char *token);
void assign_ids(char *token, int *ids, int num_ids);
void init_ortho_maps(struct mapping_list * maps, int num_list);
void init_o_list(struct o_list * oids, int num_oids);
void init_part_list(struct part_list * pieces, int num_parts);
void print_view_file(struct p_tree *t, struct sp_list *genes, int num_genes, struct mapping_list *ortho_maps, int num_list);
void postorder_print(struct p_tree *t, struct mapping_list *ortho_maps, int num_list);
void assign_sp_id_mapping_list(struct p_tree *t, struct mapping_list *list, int depth, int num_sp);
int get_id_in_mapping_list(char *name, struct mapping_list *list, int num_sp);
int concat_piece(char *line, int loc, char *temp_name);
int count_pieces(char *token);
void flip_tree(struct p_tree *t, int ref_id);

#endif /* RENUM_UTIL_H */
