#ifndef FIND_PS_H
#define FIND_PS_H

void find_pseudogenes(struct exons_list *genes1, int *num_genes1, struct exons_list *genes2, int *num_genes2, struct exons_list *exons2, int *num_exons2, struct DotList *algns, int num_algns, FILE *f, struct exons_list **candi_ps, int *num_candi, int *num_alloc);
void update_ps_candidates(struct exons_list cur_gene, struct exons_list * exons2, int from, int to, struct exons_list *genes1, int *num_genes1, struct DotList *algns, int num_algns, FILE *f, struct exons_list **candi_ps, int *num_candi, int *num_alloc, struct slist *sorted);
bool is_exons_included(struct I temp_reg, char cur_sign, struct exons_list *exons, int from, int to);
bool is_genes_included(struct I temp_reg, char cur_sign, struct exons_list *exons, int from, int to);
bool is_ps_overlapped(struct I temp_reg, struct I cmp_reg, char cur_sign, struct exons_list *candi_ps, int num);
void remove_redun_candi(struct exons_list *candi_ps, int num, struct exons_list *genes, int num_genes);

#endif /* FIND_PS_H */
