#ifndef MAIN_H
#define MAIN_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <string.h>

#define LEN_NAME 100

#define DELETED 2 

#define SELF1 3
#define SELF2 4
#define SELF 0
#define PAIR 1
#define INIT_SELF1 5
#define INIT_SELF2 6
#define INIT_PAIR 7

#define TRUE 0
#define FALSE 1

#define WIN_SIZE 100
#define MIN_LEN 500
#define ERR_SM_TH 101
#define ERR_LG_TH 202
#define ERR_TH 305
#define ERR_MAX_TH 400
#define BIG 1000000
#define LG_LEN 3000000
#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

#define MIN_GAP 200 // threshold of the minimum length of a gap to force to chain 
#define C_OP_TH 3000 // Threshold value of the distance between two close overlapped alignments
#define T_OP_TH 100
#define RR_C 0.70 // Threshold of a repeat ratio to decide to be contained
#define RR_NC 0.30 // Threshold of a repeat ratio to decide to be not contained
//#define RR_C 0.65 // Threshold of a repeat ratio to decide to be contained
//#define RR_NC 0.35 // Threshold of a repeat ratio to decide to be not contained
//#define RR_L_C 0.55 // loose threshold of a repeat ratio to decide to be contained
//#define RR_L_NC 0.45 // loose threshold of a repeat ratio to decide to be not contained
#define RR_L_C 0.70 // Threshold of a repeat ratio to decide to be contained
#define RR_L_NC 0.30 // Threshold of a repeat ratio to decide to be not contained
#define RP_BD 200 // Threshold of boundary of a gap and a repeat
#define S_RP_BD 50 // Threshold of boundary of a repeat to ignore
#define RP_TH 0.1 // the threshold for the difference between the boundary of a repeat and a gap
#define RP_P_IDT 86 // 14% is a divergence rate of an interspersed repeat 24MYRs ago, but an interspersed repeat with 14% of a divergence rate in AluY family exists
#define RP_DIFF 5 // the difference of percentage identity of an interspersed repeat and a local alignment
#define P_IDT 95 // the cutoff value of the percentage identity in a dot plot 
#define G_MODE 0 // when reading a dot plot for pairwise alignment between different species
#define D_MODE 1 // when reading a dot plot for self-alignment

#define S_FACTOR 1 // the scaling factor
#define NUM_CHARS_LINE 50 // the number of characters per a line
#define NUM_GAPS 10000 // Maximum of the number of gaps in a dot plot
#define Max_num 5000 // Maximum of the number of total regions in a dot plot 
#define Max_list 4000 // Maximum of the number of total regions in a list 
#define Max_base 500000 // Maximum bases in the dot plot

/* Be careful of setting of threshold values */
#define THRESHOLD 50
#define LOOSEN_T 2000
#define T_DIS_TH 6500 // Threshold for checking the transitivity
#define TIGHT_T_DIS_TH 5000 // Threshold for checking the transitivity
#define DIS_THRESHOLD 1500 // Threshold to use checking the existence of alignments
#define NEW_DIS_TH 3000 // Threshold to check the alignment which looks like almost self alignment
#define NEW_CHECK_TH 50000
#define O_TH 50 // Threshold to be used for checking overlap
#define TIGHT_BASE_TH 5000 // Threshold to be used for checking overlap
#define TIGHT_TH 1500 // Threshold to use for tight_subset()
#define RANGE_TH 9000 // Threshold of the range of alignments
#define S_RANGE_TH 1000 // Threshold of the range of detecting self alignments
#define LEN_TH 2000 // 
#define GAP_THRESHOLD 100000 // Threshold of searching a kd-tree for chaining 
#define MDIS_THRESHOLD 1000000 // Threshold of interspersed regoins to check the mergeability
#define SIZE_OF_ONE_UNIT 5
#define BASIS 3
#define LEN 50

#define DIF_TH 3 // the threshold to decide the difference between two identities
#define C_TH 10
#define M_THRESHOLD 0
#define TD_TH 100 // Threshold to be used for checking overlap
#define M_TH 250 // Threshold in the decision of the occurrence of insertion
#define L_M_TH 500 // Threshold in the decision of the occurrence of insertion
#define S_M_TH 150 // Threshold in the decision of the occurrence of insertion
#define LG_TH 3000 // the standard to discriminate a large alignment and a small alignment
#define T_RATE 0.90
#define RATE_DIF 0.1

#define EFFECTIVE_VALUE 150
// #define EFFECTIVE_VALUE 10
#define NUM_LOOP 5

#define SELF_ID -1

#define Num_Ops 1000
#define ALT_EFFEC_VALUE 10
#define FRAC_P_MATCHES 0.90
#define PERFECT_ID 99
#define RD_MODE 2

enum bool_type {false, true};
typedef enum bool_type bool;

struct I {
  int lower;
  int upper;
};

struct slist{
  int id;
	int val;
	bool is_x;
};

struct DotList{
  struct I x, y;
  struct I x_int, y_int;
  int sign; // if the sign is 0, the direction is the same
  int init_sign; // if the sign is 0, the direction is the same
	          // otherwise, the direction of two regions is reverse
	int fid;
	int identity;
	int l_pid;
	int l_id;
	int m_id;
	int lock;
	int left_diff, right_diff;
	int xl_diff, xr_diff, yl_diff, yr_diff;
	int xl_offset, xr_offset, yl_offset, yr_offset;
	int c_id;
	int rp1_id, rp2_id;
	int pair_self;
	int index;
	struct I m_x, m_y;
	int s_id; // id of a symmetric alignment
	int len1, len2;
	char name1[LEN_NAME], name2[LEN_NAME];
};

struct b_list{
	int b1, e1, len1;
	int b2, e2, len2;
	char strand;
	int pid;
};

struct r_list{
	int id;
	int start;
	int end;
	float d_rate;
	char name[LEN_NAME];
};

struct gap_list{
	int rp_id1, rp_id2;
	int type;
	/*
		0: overlap vs. overlap
		1: overlap vs. no overlap
		2: no overlap vs. overlap
		3: no overlap vs. no overlap
	*/
	int id1;
	int id2;
	int x1, x2; // matched pair - classifier
	int y1, y2; // gap interval
	char name1[LEN_NAME], name2[LEN_NAME];
};

struct ct_list{ // constraints list
	int id1; // precedent alignment 
	bool is_x1;
	int id2; // later alignment
	bool is_x2;
};

struct g_list{
  int gid; // a gene id
  char gname[LEN_NAME];
  char strand; // if the strand is '+', the direction is non-reverse
            // otherwise, the strand is reverse
  int txStart, txEnd; // the location of a gene
	char sp_name[LEN_NAME];
};
#endif /* MAIN_H */
