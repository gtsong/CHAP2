#ifndef CHECK_GAPS_H
#define CHECK_GAPS_H

#define X_MODE 0
#define Y_MODE 1

int get_score(struct DotList *dots, struct gap_list gps, float *d_rate, struct r_list *rp1, int num_rp1, struct r_list *rp2, int num_rp2, int *id1, int *id2, FILE *fp);
void merge_two_alignments(struct DotList *dots, struct gap_list *gps, int num_gaps, struct r_list *rp, int num_rp);
int check_insertion(struct DotList *dots, struct gap_list gp, float *d_rate, struct r_list *rp_list, int num_list, int flag, int *id);
void check_gap(struct DotList *dots, struct gap_list *gps, int num_gaps, int *temp_list, int num_list);

#endif /* CHECK_GAPS_H */
