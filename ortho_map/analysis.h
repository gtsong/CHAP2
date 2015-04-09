#ifndef ANALYSIS_H
#define ANALYSIS_H

//void get_numbers(struct DotList *content_ortho, int num_content, struct DotList *position_ortho, int num_position, struct exons_list *exons, int num_exons, FILE *f, FILE *g, int size1, int size2);
void get_numbers(struct DotList *content_ortho, int num_content, struct DotList *position_ortho, int num_position, FILE *f, FILE *g, int size1, int size2);
void compute_common_algns(struct DotList *content_ortho, int num_content, struct DotList *position_ortho, int num_position, FILE *f, struct slist *common);

#endif /* ANALYSIS_H */
