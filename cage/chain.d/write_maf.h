#ifndef WRITE_MAF_H
#define WRITE_MAF_H

#define REG 0
#define EXT 1

#define GAP_INC 0
#define NO_GAP_INC 1

void write_maf(char *fname, struct DotList *algns, int num_algns, struct r_list *rp1, struct r_list *rp2, int len1, int len2, FILE *fp, char *species, char *species2);
void adjust_two_ch_algns(struct DotList *algns, int id1, int id2, int beg1, int beg2, FILE *fp, int *s_b, int *s_e, bool is_x, struct r_list *rp1, struct r_list *rp2, int rp_id1, int rp_id2);
int get_algn_ch(char *S1, char *T1, struct DotList algn, struct DotList *init_algns, FILE *fp, struct r_list *rp1, struct r_list *rp2, struct b_list *binfo, int *cid, bool *is_written);
void adjust_end_seq(struct DotList *init_algns, int old_id, int cur_id, int beg1, int beg2, FILE *fp, int *fst_end, int *sec_beg);

#endif /* WRITE_MAF_H */
