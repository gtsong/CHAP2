#ifndef MAIN_H
#define MAIN_H
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#define TRUE 0
#define FALSE 1

#define NAME_LEN 50
#define BIG 1000000
#define SP_1 1
#define SP_2 2

#define PERFECT_ID 99 // the basis of the perfect match ID
#define FRAC_P_MATCHES 0.90 // the fraction of the perfect matches

#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

#define MIN_LEN 500
#define ERR_SM_TH 101
#define ERR_LG_TH 202
#define ERR_TH 305

#define RP_BD 200 // Threshold of boundary of a gap and a repeat
#define S_RP_BD 50 // Threshold of boundary of a repeat to ignore
#define RP_TH 0.1 // the threshold for the difference between the boundary of a repeat and a gap
#define RP_P_IDT 86 // 14% is a divergence rate of an interspersed repeat 24MYRs ago, but an interspersed repeat with 14% of a divergence rate in AluY family exists
#define RP_DIFF 3 // the difference of percentage identity of an interspersed repeat and a local alignment
#define P_IDT 50 // the cutoff value of the percentage identity in a dot plot 

#define G_MODE 0 // when reading a dot plot for pairwise alignment between different species
#define D_MODE 1 // when reading a dot plot for self-alignment
#define RD_MODE 2

#define S_FACTOR 1 // the scaling factor
#define NUM_CHARS_LINE 50 // the number of characters per a line
#define NUM_GAPS 10000 // Maximum of the number of gaps in a dot plot
#define Num_Ops 1000
#define Max_num 20000 // Maximum of the number of total regions in a dot plot 
#define Max_list 4000 // Maximum of the number of total regions in a list 
#define Max_base 500000 // Maximum bases in the dot plot

/* Be careful of setting of threshold values */
#define THRESHOLD 50
#define LOOSEN_T 2000
//#define T_DIS_TH 4000 // Threshold for checking of masking out 
// #define T_DIS_TH 2000 // Threshold for checking of masking out 
#define T_DIS_TH 1000 // Threshold for checking of masking out 
// #define SELF_T_DIS_TH 2000 // Threshold for checking of masking out 
#define SELF_T_DIS_TH 500 // Threshold for checking of masking out 
#define LONG_SELF 4000
#define MID_SELF 3000
#define MID2_SELF 1000
#define TIGHT_T_DIS_TH 5000 // Threshold for checking the transitivity
#define DIS_THRESHOLD 1500 // Threshold to use checking the existence of alignments
#define NEW_DIS_TH 1000 // Threshold to check the alignment which looks like almost self alignment
#define NEW_CHECK_TH 50000
#define O_TH 50 // Threshold to be used for checking overlap
#define TIGHT_BASE_TH 5000 // Threshold to classify the threshold of tight_subset()
#define TIGHT_TH 2000 // Threshold to use for tight_subset()
#define RANGE_TH 1000 // Threshold of the range of alignments
#define LEN_TH 2000 // 
#define MDIS_THRESHOLD 1000000 // Threshold of interspersed regoins to check the mergeability
#define SIZE_OF_ONE_UNIT 5
#define BASIS 3
#define LEN 50

#define DIF_TH 2 // the threshold to decide the difference between two identities
#define C_TH 10
#define M_THRESHOLD 0
#define M_TH 250 // Threshold in the decision of the occurrence of insertion
#define T_RATE 0.90
#define RATE_DIF 0.1

// #define S_EF_VALUE 300
#define S_EF_VALUE 10
// #define ALT_EFFEC_VALUE 300
// #define EFFECTIVE_VALUE 10
#define ALT_EFFEC_VALUE 10
#define EFFECTIVE_VALUE 10
#define NUM_LOOP 1

#define SELF_ID -1

enum bool_type {false, true};
typedef enum bool_type bool;

struct I {
  int lower;
  int upper;
};

struct seg_em{
	struct I x, y;
	int pid;
	struct seg_em *next;
};

struct DotList{
  struct I x, y;
  int sign; // if the sign is 0, the direction is the same
	          // otherwise, the direction of two regions is reverse
	struct I m_x, m_y; // id of the starting and end interval of this line in the meaningful split regions - the starting number is 1
	int fid;  // an initial id
	int identity; // the average identity
	int l_pid; // the lowest identity
	int l_id; 
	int lock;
	int left_diff, right_diff;
	int c_id; // alignment id connected in chaining
	int pair_self;
	int rp1_id, rp2_id;
	int size1, size2;
	char name1[NAME_LEN], name2[NAME_LEN];
};

struct breakpoint{
	int pt;
	int reuse_s;
	int reuse_t;
};

struct IntList{
	struct I reg;
	char name1[NAME_LEN];
	char name2[NAME_LEN];
}; //interval list

struct slist{
	int id;
	int val;
	int sp_state;
	bool is_x;
};

struct r_list{
	int start;
	int end;
	int len;
	float d_rate;
	char rn[10];
};

struct gap_list{
	int gid;
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
//	r_list rp[30];
	int num_rp;
};

struct ct_list{ // constraints list
	int id1; // precedent alignment 
	bool is_x1;
	int id2; // later alignment
	bool is_x2;
};

#endif /* MAIN_H */
