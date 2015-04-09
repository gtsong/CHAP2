#ifndef MAP_GENES_H
#define MAP_GENES_H

void map_genes(struct DotList *algns, int num_algns, struct exons_list *exons1, int num_exons1, struct exons_list *genes1, int num_genes1, struct exons_list *genes2, int num_genes2, FILE *f);
void map_genes_partition(struct DotList *algns, int num_algns, struct exons_list *exons1, int num_exons1, struct exons_list *genes1, int num_genes1, struct exons_list *genes2, int num_genes2, FILE *f);
int get_gene_index(struct I reg, struct exons_list *genes, int num_genes, struct exons_list *exons, int num_exons);

#endif /* MAP_GENES_H */
