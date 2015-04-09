#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"

extern char S[BIG], T[BIG];
extern int self_pair;

char *nucs(char *x) {
  int i;

  if ((*x != 's') || (*++x != ' '))
    fatal("expecting an s-line");
  for (i = 0; i < 5; ++i) {
    while (*++x != ' ')
      ;
    while (*++x == ' ')
      ;
  }
  if (!(strchr("ACGTN", toupper(*x)) || strchr("-", *x)))
    fatalf("expecting nucleotides in %s", x);
  return x;
}

int cal_pid(char *s, char *t, int ncol) {
  int i, match = 0, aligned = 0;
  int pid;

  for( i = 0; i < ncol; i++ )
  {
    if( (s[i] != '-') && (t[i] != '-') )
    {
      aligned++;
      if( toupper(s[i]) == toupper(t[i]) )
      {
        match++;
      }
    }
  }

  pid = (int)((((float)(match*100))/((float)aligned)) + 0.5);
  return(pid);
}

void read_maf(char *fname, int mode, struct DotList *algns, int *num_algns, int *size1, int *size2) {
	FILE *fp;
	char *status;
	int i = 0;
	int count = 0;
	int temp;
	int a_pid;
	int b1, e1, b2, e2;
	char strand[100], len1[100], len2[100];
	char species1[100], species2[100];
	char *s, *t;
	int ncol = 0;
	int max_len1 = 0, max_len2 = 0;

	fp = ckopen(fname, "r");
	if ((fgets(S, BIG, fp) == NULL) || strncmp(S, "##maf version", 13))
		fatalf("%s is not a maf file", fname);
	while (S[0] == '#')
		if ((status = fgets(S, BIG, fp)) == NULL) {}
//			fatalf("no alignments in %s", fname);

	while ((status != NULL) && (strstr(S, "eof") == NULL)) {
		if (strncmp(S, "a ", 2))
			fatalf("expecting an a-line in %s, saw %s",
			  fname, S);
		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL))
			fatalf("cannot find alignment in %s", fname);
		if ((sscanf(S, "%*s %s %d %d %*s %s", species1, &b1, &e1, len1) != 4) ||
		    (sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5))
		{
			fatalf("bad alignment info in %s", fname);
		}
		// aligned interval given as base-0 start and length
		e1 += b1;
		e2 += b2;

		if( strcmp(strand, "-") == 0) {
			temp = b2;
			b2 = atoi(len2) - e2;
			e2 = atoi(len2) - temp;	
		}			

		b1++;
		b2++;
		e1++;
		e2++;

		s = nucs(S);
		t = nucs(T);
		ncol = strlen(s);
//		if( ncol >= 600 ) {
//			a_pid = (int) (cal_pid_maf_beg(s, t, 150, strlen(s)-150) + 0.5);
//		}
//		else {
			a_pid = cal_pid(s, t, strlen(s)-1);
//		}

		if( (mode == D_MODE) && ( (b1 >= b2) || ((b1 == b2) && (e1 == e2)) ) ) {}
		else  {
			algns[count].x = assign_I(b1, e1);
			if( b2 < e2 ) algns[count].y = assign_I(b2, e2);
			else algns[count].y = assign_I(e2, b2);
			algns[count].identity = a_pid;

			if( strcmp(strand, "+") == 0 ) algns[count].sign = 0;
			else if( strcmp(strand, "-") == 0 ) algns[count].sign = 1;
			algns[count].init_sign = algns[count].sign;
			algns[count].fid = i; // ith alignment
			algns[count].index = count; // ith alignment
      algns[count].c_id = -1; // not chained alignment
      algns[count].rp1_id = -1; // the inserted repeat id of the chained alignment in first seq
      algns[count].rp2_id = -1; // the inserted repeat id of the chained alignment in second seq 
      algns[count].l_id = -1;
      algns[count].lock = -1;  
      algns[count].m_x = assign_I(algns[count].x.lower, algns[count].x.upper);
      algns[count].m_y = assign_I(algns[count].y.lower, algns[count].y.upper);
      algns[count].left_diff = 0;
      algns[count].right_diff = 0;
      algns[count].pair_self = -1;
      algns[count].l_pid = -1;
      algns[count].len1 = atoi(len1);
      algns[count].len2 = atoi(len2);
			strcpy(algns[count].name1, species1);
			strcpy(algns[count].name2, species2);

			if( algns[count].len1 > max_len1 ) max_len1 = algns[count].len1;
			if( algns[count].len2 > max_len2 ) max_len2 = algns[count].len2;
			
			count++;
		}

		if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
			fatalf("bad alignment end in %s", fname);
		status = fgets(S, BIG, fp);
		i++; // ith alignment 
	}

	*size1 = max_len1;
	*size2 = max_len2;
	*num_algns = count;
	fclose(fp);
}

float cal_pid_maf(char *s, char *t, int ncol) {
  int i = 0, match = 0, aligned = 0;
  float pid = (float) 0;

  for( i = 0; i < ncol; i++ )
  {
    if( (strchr("ACGT", toupper(s[i]))) && (strchr("ACGT", toupper(t[i]))))
    {
      aligned++;
      if( toupper(s[i]) == toupper(t[i]) )
      {
        match++;
      }
    }
  }

  pid = ((float)(match*100))/((float)aligned);
  return(pid);
}

float cal_pid_maf_beg(char *s, char *t, int beg, int ncol) {
  int i, match = 0, aligned = 0;
  float pid;

  for( i = beg; i < ncol; i++ )
  {
    if( (strchr("ACGT", toupper(s[i]))) && (strchr("ACGT", toupper(t[i]))))
    {
      aligned++;
      if( toupper(s[i]) == toupper(t[i]) )
      {
        match++;
      }
    }
  }

  pid = ((float)(match*100))/((float)aligned);
  return(pid);
}
