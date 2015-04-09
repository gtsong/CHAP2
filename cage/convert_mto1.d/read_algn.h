#ifndef READ_ALGN_H
#define READ_ALGN_H

#define REG 0
#define EXT 1

#define GAP_INC 0
#define NO_GAP_INC 1
#define GAP_INC_IN_Y 2

void get_nth_algn(char *seq1, char *seq2, int id, int b, FILE *fp, struct b_list *a_info, int mode);
char *nucs_b(char *x, int b, int *beg);
char *nucs_b2(char *x, int b, int *count);
bool scan_compare(char *S3, char *T3, char *S2, char *T2, int *avg_diff);
int find_yloc_one(struct DotList algn, FILE *fp, int num_nu, int flag);
int find_end_xloc(struct DotList algn, FILE *fp, int num_nu);
int find_xloc_one(struct DotList algn, FILE *fp, int num_nu, int flag);
void fill_N(int gap_len1, int gap_len2, char *S3, char *T3);
int find_yloc_one_ch(struct DotList *init_dots, struct DotList algn, FILE *fp, int num_nu, int flag);
int count_ncol(struct I reg1, struct I reg2, char *S, int end_loc, int *beg);

#endif /* READ_ALGN_H */
