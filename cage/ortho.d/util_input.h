#ifndef UTIL_INPUT_H
#define UTIL_INPUT_H

void initialize_exons_list(struct exons_list *a, int from, int to);
int count_exons(char *fname, char *species, char *species2);
int count_genes(char *fname, char *species, char *species2);
void read_only_exons(struct exons_list *exons, int *num_exons, struct exons_list
 *genes, int *num_genes, char *fname, char *species, char *species2);
void read_exons(struct exons_list *init_exons, struct exons_list *exons, int *num_exons, struct exons_list *genes, int *num_genes, struct exons_list *skip_reg1, int *num_skip1, struct exons_list *skip_reg2, int *num_skip2, char *fname, struct sp_list *sp_code, int num_sp_code, int sp1_code, int sp2_code);
int count_local_algns(char *fname, char *species, char *species2);
int count_lines(char *fname);
int count_tokens(char *line);
int concat_tokens(char *line, int loc, char *temp_name);
void numtostr(int c, char *str);
bool is_all_digits(char *item);

#endif /* UTIL_INPUT_H */
