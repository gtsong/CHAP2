#ifndef MAIN_H
#define MAIN_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <string.h>

#define WIN_SIZE 100
#define LEN_NAME 100
#define BIG 1000000

#define REG 0
#define EXT 1

#define SELF 0
#define PAIR 1
#define SELF1 2
#define SELF2 3

#define TRUE 0
#define FALSE 1

#define DEL_TH 15
#define PID_TH 1

#define ALT_EFFEC_VALUE 0
#define G_MODE 0 // when reading a dot plot for pairwise alignment between different species
#define D_MODE 1 // when reading a dot plot for self-alignment
#define C_MODE 2 // when reading a dot plot for self-alignment

#define DELETED 2

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
  int sign; // if the sign is 0, the direction is the same
	          // otherwise, the direction of two regions is reverse
	struct I x_int, y_int; // id of the starting and end interval of this line in the meaningful split regions - the starting number is 1
	int fid;
	int identity;
	int l_pid;
	int l_id;
	int lock;
	int left_diff, right_diff;
	int c_id;
	int rp1_id, rp2_id;
	int pair_self;
	struct I m_x, m_y;
	int s_id; // id of a symmetric alignment
};

struct b_list{
	int b1, e1, len1;
	int b2, e2, len2;
	char strand;
	int pid;
};

#endif /* MAIN_H */
