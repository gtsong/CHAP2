#include "main.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"
#include "contigs_op.h"

extern char S[BIG], T[BIG];
extern int debug_mode;

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
  float pid;

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
	char strand[LEN_NAME], len1[LEN_NAME], len2[LEN_NAME];
	char name1[LEN_NAME], name2[LEN_NAME];
	char *s, *t;
	int algn_type = (SELF1 - 1);
	int j = 0;

	fp = ckopen(fname, "r");
	if (((status = fgets(S, BIG, fp)) == NULL) || strncmp(S, "##maf version", 13))
		fatalf("%s is not a maf file", fname);
	while ((status != NULL) && (S[0] == '#')) {
		if( (mode == C_MODE) && (strncmp(S, "##maf", 5) == 0) ) {
			algn_type++;
			j = 0;
		}
		status = fgets(S, BIG, fp);
	}

	while ((status != NULL) && (strstr(S, "eof") == NULL)) {
		if(S[0] == '#') {
			if(mode == C_MODE) {
				while(S[0] == '#') {
					if ((status = fgets(S, BIG, fp)) == NULL)
						fatalf("no alignments in %s", fname);
					if( strncmp(S, "##maf", 5) == 0 ) algn_type++;
				}	
				algn_type++;
				if( algn_type > PAIR ) fatal("too many alignments are combined\n");
			}
			else {
				fatalf("not supported maf format\n");
			}
			j = 0;
		}

		if (S[0] != 'a')
			fatalf("expecting an a-line in %s, saw %s",
			  fname, S);
		if ((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL))
			fatalf("cannot find alignment in %s", fname);
		if ((sscanf(S, "%*s %s %d %d %*s %s", name1, &b1, &e1, len1) != 4) || (sscanf(T, "%*s %s %d %d %s %s", name2, &b2, &e2, strand, len2) != 5))
		{
			fatalf("bad alignment info of 2 in %s", fname);
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
		a_pid = cal_pid(s, t, strlen(s)-1);

		if( ((mode == D_MODE) || ((mode == C_MODE) && (algn_type <= PAIR))) && (( (algn_type != PAIR) && (b1 >= b2)) || ((algn_type != PAIR) && (abs(b1-b2) <= DEL_TH) && (abs(e1-e2) <=DEL_TH)) || ((e1-b1) < ALT_EFFEC_VALUE) || (a_pid <= PID_TH) )) {}
		else  {
			algns[count].x = assign_I(b1, e1);
			if( b2 < e2 ) algns[count].y = assign_I(b2, e2);
			else algns[count].y = assign_I(e2, b2);
			algns[count].identity = a_pid;

			if( strcmp(strand, "+") == 0 ) {
				algns[count].sign = 0;
				algns[count].init_sign = 0;
			}	
			else if( strcmp(strand, "-") == 0 ) {
				algns[count].sign = 1;
				algns[count].init_sign = 1;
			}
			else {
				algns[count].sign = DELETED;
				algns[count].init_sign = DELETED;
			}

			algns[count].fid = i; // ith alignment
			algns[count].indiv_fid = j; // j alignment
			algns[count].index = count; // ith alignment
      algns[count].c_id = -1; // not chained alignment
      algns[count].m_id = -1; // not chained alignment
      algns[count].rp1_id = -1; // the inserted repeat id of the chained alignment in first seq
      algns[count].rp2_id = -1; // the inserted repeat id of the chained alignment in second seq 
      algns[count].l_id = -1;
      algns[count].lock = -1;  
      algns[count].m_x = assign_I(0,1);
      algns[count].m_y = assign_I(0,1);
      algns[count].xl_diff = 0; // the offset of the left end
      algns[count].yl_diff = 0; // the offset of the left end
      algns[count].xr_diff = 0; // the offset of the right end
      algns[count].yr_diff = 0; // the offset of the right end
      algns[count].pair_self = -1;
      algns[count].l_pid = -1;
			algns[count].sp_id = algn_type; // SELF1 for first self-alignment, SELF2 for second self-alignment and PAIR for pairwise alignment
      algns[count].xl_offset = 0; // the offset of low of x
      algns[count].yl_offset = 0; // the offset of up of x
      algns[count].xr_offset = 0; // the offset of low of y 
			if( algn_type == PAIR ) algns[count].pair_self = PAIR;
			else algns[count].pair_self = SELF;
			strcpy(algns[count].name1, name1);
			strcpy(algns[count].name2, name2);
			algns[count].len1 = atoi(len1);
			algns[count].len2 = atoi(len2);
			algns[count].ctg_id1 = -1;
			algns[count].ctg_id2 = -1;

			count++;
		}

		if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
			fatalf("bad alignment end in %s", fname);
		status = fgets(S, BIG, fp);
		i++; // ith alignment 
		j++;
	}

	*size1 = atoi(len1);
	*size2 = atoi(len2);
	*num_algns = count;
	fclose(fp);
}

float cal_pid_maf(char *s, char *t, int ncol) {
  int i, match = 0, aligned = 0;
  float pid;

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

void adjust_multi_contig_pos(struct DotList *algns, int num_algns, int *size1, int *size2, struct n_pair *contigs1, int *num_contigs1, struct n_pair *contigs2, int *num_contigs2)
{
	int i = 0, id = -1;
	int num1 = 0, num2 = 0;
	char sp_name1[LEN_NAME], ctg_name1[LEN_NAME];
	char sp_name2[LEN_NAME], ctg_name2[LEN_NAME];
	char name1[LEN_NAME], name2[LEN_NAME];

	strcpy(sp_name1, "");
	strcpy(sp_name2, "");
	strcpy(ctg_name1, "");
	strcpy(ctg_name2, "");

	for( i = 0; i < num_algns; i++ ) {
		strcpy(name1, algns[i].name1);
		concat_ctg_name(name1, sp_name1, ctg_name1);
		strcpy(name2, algns[i].name2);
		concat_ctg_name(name2, sp_name2, ctg_name2);

		if( algns[i].sp_id == SELF1 ) {
			id = add_ctg_name(algns[i].index, sp_name1, ctg_name1, algns[i].len1, contigs1, num1); 
			if( id == num1 ) num1++;
			algns[i].ctg_id1 = id;

			id = add_ctg_name(algns[i].index, sp_name2, ctg_name2, algns[i].len2, contigs1, num1); 
			if( id == num1 ) num1++;
			algns[i].ctg_id2 = id;
		}
		else if( algns[i].sp_id == SELF2 ) {
			id = add_ctg_name(algns[i].index, sp_name1, ctg_name1, algns[i].len1, contigs2, num2); 
			if( id == num2 ) num2++;
			algns[i].ctg_id1 = id;

			id = add_ctg_name(algns[i].index, sp_name2, ctg_name2, algns[i].len2, contigs2, num2); 
			if( id == num2 ) num2++;
			algns[i].ctg_id2 = id;
		}
		else if( algns[i].sp_id == PAIR ) {
			id = add_ctg_name(algns[i].index, sp_name1, ctg_name1, algns[i].len1, contigs1, num1); 
			if( id == num1 ) num1++;
			algns[i].ctg_id1 = id;

			id = add_ctg_name(algns[i].index, sp_name2, ctg_name2, algns[i].len2, contigs2, num2); 
			if( id == num2 ) num2++;
			algns[i].ctg_id2 = id;
		}
		else {
			fatalf("wrong alignment type in %s:%d-%d, %s:%d-%d\n", algns[i].name1, algns[i].x.lower, algns[i].x.upper, algns[i].name2, algns[i].y.lower, algns[i].y.upper);
		}
	}

	adjust_algn_pos(algns, num_algns, contigs1, num1, size1, contigs2, num2, size2, CTG_ASSIGNED);

	*num_contigs1 = num1;
	*num_contigs2 = num2;
}

void adjust_algn_pos(struct DotList *algns, int num_algns, struct n_pair *contigs1, int num1, int *size1, struct n_pair *contigs2, int num2, int *size2, int mode)
{
	int *len_sum1, *len_sum2;
	int i = 0;
	int id1 = 0, id2 = 0;
	char name[LEN_NAME], sp_name[LEN_NAME], ctg_name[LEN_NAME];
	int ctg_id = -1;
	
	if( num1 > 0 ) len_sum1 = (int *) ckalloc(sizeof(int) * num1);
	if( num2 > 0 ) len_sum2 = (int *) ckalloc(sizeof(int) * num2);
	
	for( i = 0; i < num1; i++ ) len_sum1[i] = 0;
	for( i = 0; i < num2; i++ ) len_sum2[i] = 0;

  cal_length_sum(len_sum1, contigs1, num1);
  cal_length_sum(len_sum2, contigs2, num2);

	for( i = 0; i < num_algns; i++ ) {
		if( mode == CTG_NOT_ASSIGNED ) {
			strcpy(name, algns[i].name1);
			if( algns[i].sp_id == SELF2 ) {
				concat_ctg_name(name, sp_name, ctg_name);
				ctg_id = is_ctg_in(sp_name, ctg_name, contigs2, num2);
			}
			else {
				concat_ctg_name(name, sp_name, ctg_name);
				ctg_id = is_ctg_in(sp_name, ctg_name, contigs1, num1);
			}
			if( ctg_id == -1 ) {
				fatalf("Contig %s not assigned in the list\n", ctg_name);
			}
			else {
				algns[i].ctg_id1 = ctg_id;
			}

			strcpy(name, algns[i].name2);
			if( algns[i].sp_id == SELF1 ) {
				concat_ctg_name(name, sp_name, ctg_name);
				ctg_id = is_ctg_in(sp_name, ctg_name, contigs1, num1);
			}
			else {
				concat_ctg_name(name, sp_name, ctg_name);
				ctg_id = is_ctg_in(sp_name, ctg_name, contigs2, num2);
			}

			if( ctg_id == -1 ) {
				fatalf("Contig %s not assigned in the list\n", ctg_name);
			}
			else {
				algns[i].ctg_id2 = ctg_id;
			}
		}
		
		if( algns[i].sp_id == SELF1 ) {
			id1 = algns[i].ctg_id1;
			if( id1 >= num1 ) fatalf("%d: not valid, larger than %d\n", id1, num1);
			if( (id1 == -1)	&& (num1 > 0) ) {
				fatalf("wrong contig num assigned: %s - %s in read_maf.c\n", algns[i].name1, algns[i].name2);
			}

			if( id1 != -1 ) algns[i].x = assign_I(algns[i].x.lower + len_sum1[id1], algns[i].x.upper + len_sum1[id1]);

			id2 = algns[i].ctg_id2;
			if( id2 >= num1 ) {
				fatalf("%d: not valid, larger than %d\n", id2, num1);
			}	
			if( (id2 == -1)	&& (num1 > 0) ) {
				fatalf("wrong contig num assigned: %s - %s in read_maf.c\n", algns[i].name1, algns[i].name2);
			}

			if( id2 != -1 ) algns[i].y = assign_I(algns[i].y.lower + len_sum1[id2], algns[i].y.upper + len_sum1[id2]);
		}
		else if( algns[i].sp_id == SELF2 ) {
			id1 = algns[i].ctg_id1;
			if( id1 >= num2 ) {
				fatalf("%d: not valid, larger than %d\n", id1, num2);
			}	
			if( (id1 == -1)	&& (num2 > 0) ) {
				fatalf("wrong contig num assigned: %s - %s in read_maf.c\n", algns[i].name1, algns[i].name2);
			}

			if( id1 != -1 ) algns[i].x = assign_I(algns[i].x.lower + len_sum2[id1], algns[i].x.upper + len_sum2[id1]);

			id2 = algns[i].ctg_id2;
			if( id2 >= num2 ) {
				fatalf("%d: not valid, larger than %d\n", id2, num2);
			}	
			if( (id2 == -1)	&& (num2 > 0) ) {
				fatalf("wrong contig num assigned: %s - %s in read_maf.c\n", algns[i].name1, algns[i].name2);
			}

			if( id2 != -1 ) algns[i].y = assign_I(algns[i].y.lower + len_sum2[id2], algns[i].y.upper + len_sum2[id2]);
		}
		else if( algns[i].sp_id == PAIR ) {
			id1 = algns[i].ctg_id1;
			if( id1 >= num1 ) {
				fatalf("%d: not valid, larger than %d\n", id1, num1);
			}	
			if( (id1 == -1)	&& (num1 > 0) ) {
				fatalf("wrong contig num assigned: %s - %s in read_maf.c\n", algns[i].name1, algns[i].name2);
			}

			if( id1 != -1 ) algns[i].x = assign_I(algns[i].x.lower + len_sum1[id1], algns[i].x.upper + len_sum1[id1]);

			id2 = algns[i].ctg_id2;
			if( id2 >= num2 ) {
				fatalf("%d: not valid, larger than %d\n", id2, num2);
			}	
			if( (id2 == -1)	&& (num2 > 0) ) {
				fatalf("wrong contig num assigned: %s - %s in read_maf.c\n", algns[i].name1, algns[i].name2);
			}

			if( id2 != -1 ) algns[i].y = assign_I(algns[i].y.lower + len_sum2[id2], algns[i].y.upper + len_sum2[id2]);
		}
	}

	if( num1 > 0 ) {
		*size1 = len_sum1[num1-1] + contigs1[num1-1].len;
		free(len_sum1);
	}
	if( num2 > 0 ) {
		*size2 = len_sum2[num2-1] + contigs2[num2-1].len;
		free(len_sum2);
	}
}

