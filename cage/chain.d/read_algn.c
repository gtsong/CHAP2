#include "main.h"
#include "regions.h"
#include "read_algn.h"
#include "util.h"
#include "util_gen.h"
#include "read_maf.h"

extern char S[BIG], T[BIG];
extern char S1[BIG], T1[BIG];
extern int debug_mode;

// In REG mode, b is a starting position based on the number of nucleotides and in EXT mode, b is based on counts including the gap character '-'
// a_info is a data structure to return the length of the concatenated sequence  and the starting position for '+' strand or (the end position + 1) for '-' strand
void get_nth_algn(char *seq1, char *seq2, int id, int b, FILE *fp, struct b_list *a_info, int mode) 
{
	int i = 0;
	char *status;
	int b1, e1, b2, e2;
	char strand[100], len1[100], len2[100];
	int *beg, *count;
	bool is_found = false;
	int pid;

	beg = (int *) ckalloc(sizeof(int));
	count = (int *) ckalloc(sizeof(int));
	*beg = 0;
	*count = 0;
	fseek(fp, 0, SEEK_SET);
  if ((fgets(S, BIG, fp) == NULL) || strncmp(S, "##maf version", 13))
    fatalf("not a maf file in %s\n", S);
  while (S[0] == '#')
    if ((status = fgets(S, BIG, fp)) == NULL)
      fatal("no alignments");

  while ((is_found == false) && (status != NULL) && (strstr(S, "eof") == NULL)) 
	{
  	while ((S[0] == '#') || (S[0] =='\n')) fgets(S, BIG, fp);

    if(S[0] != 'a')
      fatalf("expecting an a-line in %s", S);
		else {
			sscanf(S+8, "%d", &pid);
		}

    if((fgets(S, BIG, fp) == NULL) || (fgets(T, BIG, fp) == NULL))
      fatal("cannot find alignment");

    if ((sscanf(S, "%*s %*s %d %d %*s %s", &b1, &e1, len1) != 3) || (sscanf(T, "%*s %*s %d %d %s %s", &b2, &e2, strand, len2) != 4))
		{
			fatal("bad alignment info");
		}    
		else if( i == id ) {
			is_found = true;
			if( mode == REG ) {
				if( (strcmp( S, "" ) != 0) && (strcmp( T, "" ) != 0) )
				{
					strcpy(seq1, nucs_b(S, b, beg));
					strcpy(seq2, nucs_b2(T, *beg, count));
					(*a_info).b1 = b1+b+1;
					(*a_info).len1 = e1-b;
					if( strcmp(strand, "+") == 0 ) {
						(*a_info).b2 = b2+(*count)+1;
						(*a_info).len2 = e2-(*count);
						(*a_info).strand = '+';
						(*a_info).pid = pid;
					}
					else {
						(*a_info).b2 = (atoi(len2)-b2)-(*count)+1;
						(*a_info).len2 = e2-(*count);
						(*a_info).strand = '-';
						(*a_info).pid = pid;
					}
				}
				else {
					fatalf("empty sequence in %d\n", id);
				}
			}
			else if( mode == EXT ) {
				if( strcmp( seq1, "" ) != 0 ) 
				{
					strcat(seq1, nucs_b2(S, b, beg));
					strcat(seq2, nucs_b2(T, b, count));
				}
				else {
					strcpy(seq1, nucs_b2(S, b, beg));
					strcpy(seq2, nucs_b2(T, b, count));
				}
				(*a_info).b1 = b1+(*beg)+1;
				(*a_info).len1 = e1-(*beg);
				(*a_info).pid = pid;
				if( strcmp(strand, "+") == 0 ) {
					(*a_info).b2 = b2+(*count)+1;
					(*a_info).len2 = e2-(*count);
					(*a_info).strand = '+';
				}
				else {
					(*a_info).b2 = (atoi(len2)-b2)-(*count)+1;
					(*a_info).len2 = e2-(*count);
					(*a_info).strand = '-';
				}
			}
		}

    if ((fgets(S, BIG, fp) == NULL) || (S[0] != '\n'))
      fatal("bad alignment end");
    status = fgets(S, BIG, fp);
    i++; // ith alignment 
	}
	
	seq1[BIG-1]='\0';
	seq2[BIG-1]='\0';
	free(count);
	free(beg);
	if( is_found == false ) fatalf("%d th alignment does not exist\n", id);
}

char *nucs_b(char *x, int b, int *beg) {
  int i;
	char *temp;
	int count = 0;
	int j = 1;

	temp = x;

  if ( (*x != 's') || (*++x != ' ') )
    fatal("expecting an s-line");

  for (i = 0; i < 5; ++i) {
    while (*++x != ' ') j++;
    while (*++x == ' ') j++;
		j++;
  }
	
	for( i = 0; i < b; ++i ) {
		while (*++x == '-') {
			count++;
			j++;
			if( j >= (BIG-2) ) {
				strcpy(x, "\n");
				*beg = count;
				return x;
			}
		}
		count++;
		j++;
		if( j >= (BIG-2) ) {
			strcpy(x, "\n");
			*beg = count;
			return x;
		}
	}

	if( strlen(x) >= BIG ) {
		x[count+BIG-2] = '\n';
		x[count+BIG-1] = '\0';
	}

	if( (*x == '\n') || (*x == '\0') ) strcpy(x, "\n");
  else if (!(strchr("ACGTN", toupper(*x)) || strchr("-", *x)))
    fatalf("expecting nucleotides in %s", temp);

	*beg = count;
  return x;
}

char *nucs_b2(char *x, int b, int *beg) {
  int i;
	int count = 0;
	int j = 1;

  if ( (*x != 's') || (*++x != ' ') )
    fatal("expecting an s-line");
  for (i = 0; i < 5; ++i) {
    while (*++x != ' ') j++;
    while (*++x == ' ') j++;
		j++;
  }

	for( i = 0; i < b; ++i ) {
		if(*x != '-') count++;
		++x;
		j++;
		if( j >= (BIG-2) ) {
			strcpy(x, "\n");
			*beg = count;
			return x;
		}
	}

	if( strlen(x) >= BIG ) {
		x[count+BIG-2] = '\n';
		x[count+BIG-1] = '\0';
	}

	if( (*x == '\n') || (*x == '\0') ) strcpy(x, "\n");
  else if (!(strchr("ACGTN", toupper(*x)) || strchr("-", *x)))
    fatalf("expecting nucleotides in %s", x);
	*beg = count;
  return x;
}

bool scan_compare(char *S3, char *T3, char *S2, char *T2, int *avg_diff)
{
	int i, j;
	int num_nu = 0;
	int b, next_b;
	int diff;
	int num_win;
	float pid1, pid2;
	bool res = true;
	float diff_pid = 0, pre_dpid = 0;
	int temp_diff = 0;
	int len;

	len = strlen(S3) - 1;

	for( i = 0; i < len; i++ ) {
		if( strchr("ACGT", toupper(S3[i])) && strchr("ACGT", toupper(T3[i])) && strchr("ACGT", toupper(S2[i])) && strchr("ACGT", toupper(T2[i])) ) num_nu++;
	}

	if( num_nu <= WIN_SIZE ) {
		b = 0;
		num_win = 1;
	}
	else {
		num_win = (int)(num_nu/WIN_SIZE);
		diff = (num_nu - (num_win * WIN_SIZE))/2;
		i = 0;
		j = 0;
		b = 0;
		while( i < diff ) {
			if( strchr("ACGT", toupper(S3[j])) && strchr("ACGT", toupper(T3[j])) && strchr("ACGT", toupper(S2[j])) && strchr("ACGT", toupper(T2[j])) ) i++;
			j++;	
		}
		b = j;
	}

	for( i = 0; i < num_win; i++ ) {
		next_b = b;
		j = 0;
		while( j < WIN_SIZE ) {
			if( strchr("ACGT", toupper(S3[next_b])) && strchr("ACGT", toupper(T3[next_b])) && strchr("ACGT", toupper(S2[next_b])) && strchr("ACGT", toupper(T2[next_b])) ) j++;
			next_b++;	
		}
		pre_dpid = diff_pid;
		pid1 = cal_pid_maf_beg(S3, T3, b, next_b);
		pid2 = cal_pid_maf_beg(S2, T2, b, next_b);
		diff_pid = pid1 - pid2;
		b = next_b;
		if( (pre_dpid > 0) && (diff_pid < 0) ) {
			printf("change of alignment with higher similarity %f -> %f\n", pre_dpid, diff_pid);
			res = false;
		}
		else if( (pre_dpid < 0) && (diff_pid > 0) ) {
			printf("change of alignment with higher similarity %f -> %f\n", pre_dpid, diff_pid);
			res = false;
		}
		temp_diff = temp_diff + diff_pid;
	}
	
	temp_diff = (temp_diff+0.5)/num_win;
	*avg_diff = temp_diff;
	return(res);
}

int find_yloc_one(struct DotList algn, FILE *fp, int num_nu, int flag)
{
	int i = 0, j = 0, k = 0;
	struct b_list *a_info;
	int yloc = -1;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;

	strcpy(S1, "");
	strcpy(T1, "");
	get_nth_algn(S1, T1, algn.fid, 0, fp, a_info, REG);
	free(a_info);

	while( (i < BIG) && (((flag == NO_GAP_INC) && (j < num_nu) && (S1[i] != '\0') && (S1[i] != '\n')) || ((flag == GAP_INC) && (i < num_nu) && (S1[i] != '\0') && (S1[i] != '\n')) || ((flag == GAP_INC_IN_Y) && (k < num_nu) && (T1[i] != '\0') && (T1[i] != '\n')))) {
		if( strchr("ACGTN", toupper(S1[i])) ) j++;
		if( strchr("ACGTN", toupper(T1[i]))) k++;
		i++;
	}
	
	if( algn.init_sign == 0 ) {
		yloc = algn.y.lower + k;
//		if( algn.yl_offset != 0 ) {
//			yloc = yloc - algn.yl_offset;
//		}
	}
	else if( algn.init_sign == 1 ) {
		yloc = algn.y.upper - k;
//		if( algn.yr_offset != 0 ) {
//			yloc = yloc - algn.yr_offset;
//		}
	}
	else fatalf("find_yloc_one: unexpected value %d for algn.init_sign\n", algn.init_sign);

	if( flag == GAP_INC_IN_Y ) return(i);
	else return(yloc);
}

int find_end_xloc(struct DotList algn, FILE *fp, int num_nu) // only NO_GAP_COUNT case needs this function
{
	int j = 0, k = 0;
	struct b_list *a_info;
//	char S1[BIG], T1[BIG];

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;

	get_nth_algn(S1, T1, algn.fid, 0, fp, a_info, REG);
	free(a_info);

	while( (j < num_nu) && (S1[k] != '\0') && (S1[k] != '\n')) {
		if( strchr("ACGTN", toupper(S1[k])) ) j++;
		k++;
	}

	return(k);
}

int find_xloc_one(struct DotList algn, FILE *fp, int num_nu, int flag) // only NO_GAP_COUNT case needs this function
{
	int j = 0, k = 0;
	struct b_list *a_info;
	int xloc = -1;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;

	get_nth_algn(S1, T1, algn.fid, 0, fp, a_info, REG);
	free(a_info);

	while( (k < BIG) && (((flag == NO_GAP_INC) && (j < num_nu) && (S1[k] != '\0') && (S1[k] != '\n')) || ((flag == GAP_INC) && (k < num_nu) && (S1[k] != '\0') && (S1[k] != '\n')) )) {
		if( strchr("ACGTN", toupper(S1[k])) ) j++;
		k++;
	}

	if( flag == GAP_INC ) {
		xloc = algn.x.lower + j;
		if( algn.xl_offset != 0 ) {
			xloc = xloc - algn.xl_offset;
		}
	}
	else if( flag == NO_GAP_INC ) xloc = k;
	else fatalf("find_xloc_one: unexpected value %d for algn.init_sign\n", algn.init_sign);

	return(xloc);
}

void fill_N(int gap_len1, int gap_len2, char *S3, char *T3)
{
	int min_len;
	bool is_x;
	int i;
	int cur_len;

	if( S3[strlen(S3)-1] == '\n' ) cur_len = strlen(S3)-1;
	else cur_len = strlen(S3);

	if( gap_len1 > gap_len2 ) {
		min_len = gap_len2;
		is_x = true; // gap in x 
	}
	else {
		min_len = gap_len1;
		is_x = false; // gap in y
	}

	for( i = 0; i < gap_len1; i++ ) {
		S3[cur_len+i] = 'N';
		T3[cur_len+i] = '-';
	}

	for( i = 0; i < gap_len2; i++ ) {
		S3[cur_len+gap_len1+i] = '-';
		T3[cur_len+gap_len1+i] = 'N';
	}

	S3[cur_len+gap_len1+gap_len2] = '\0';
	T3[cur_len+gap_len1+gap_len2] = '\0';
}

int find_yloc_one_ch(struct DotList *init_dots, struct DotList algn, FILE *fp, int offset, int flag) // offset is different from one in find_yloc_one. it is an offset from the end location of x region.
{
	int res; 
	int b, e;
	int old_id, cid, mid;
	int cur_offset;
	int old_e;
	int index;

	index = algn.index;

	if( init_dots[index].c_id == -1 ) {
//		res = algn.y.lower + find_yloc_one(init_dots[index], fp, offset, flag) - init_dots[index].y.lower;	
		cur_offset = algn.x.upper - algn.x.lower - offset;
		res = find_yloc_one(init_dots[index], fp, offset, flag);	
	}	
	else {
		cid = init_dots[index].c_id;
		old_id = index;

		while( cid != -1 ) {
			old_id = cid;
			cid = init_dots[cid].c_id;
		}

		b = init_dots[old_id].x.lower + init_dots[old_id].xl_diff; 
		e = init_dots[old_id].x.upper - init_dots[old_id].xr_diff;

		if( (e - b) <= 0 ) {
			return(-1);
		}

		cur_offset = offset;		
		mid = old_id;
		while( (mid != -1) && ((e-b) < cur_offset) && (cur_offset >= 0)) {
			old_e = e;
			b = init_dots[mid].x.lower + init_dots[mid].xl_diff; 
			e = init_dots[mid].x.upper - init_dots[mid].xr_diff;
			cur_offset = cur_offset - (e - old_e);
			old_id = mid;
			mid = init_dots[mid].m_id;				
		}

		if( mid != -1 ) {
			cur_offset = init_dots[mid].x.upper - init_dots[mid].xr_offset - init_dots[mid].x.lower + init_dots[mid].xl_offset - offset;
			res = find_yloc_one(init_dots[mid], fp, cur_offset, flag);
		}	
		else {
			cur_offset = init_dots[old_id].x.upper - init_dots[mid].xr_offset - init_dots[old_id].x.lower + init_dots[old_id].xl_offset - offset;
			res = find_yloc_one(init_dots[old_id], fp, cur_offset, flag);
		}
	}

	return(res);
}

// return the column location of the end of the overlap between two intervals reg1 and reg2 including gaps
int count_ncol(struct I reg1, struct I reg2, char *seq, int end_loc, int *beg)
{
	int k = 0, count = 0;
	int ncol = 0;

	if( reg1.lower >= reg2.lower ) *beg = 0;
	else {
		k = 0;
		count = 0;
		while( (k < end_loc) && (count < abs(reg1.lower - reg2.lower)) ) {
			while( (k < end_loc) && (seq[k] == '-') ) k++;
				count++;
				k++;
		}
  	*beg = k;
	}

  if( reg1.upper <= reg2.upper ) ncol = end_loc;
	else {
		k = end_loc - 1;
		count = 0;
		while( (k >= 0) && (count < abs(reg1.upper - reg2.upper)) ) {
			while( (k >= 0) && (seq[k] == '-') ) k--;
			count++;
			k--;
		}
		ncol = k+1;
	}
	return(ncol);
}

void update_algn_diff(struct DotList *algns, int id, int point, FILE *f, int mode)
{
	int cur_diff = 0;
	int cur = 0, old = 0, loc = 0;

	if( mode == XL_DIFF ) {
    cur_diff = algns[id].xl_diff;
    if( cur_diff < abs(point-algns[id].x.lower) ) {
      algns[id].xl_diff = abs(point - algns[id].x.lower);
      cur = find_yloc_one(algns[id], f, abs(point - algns[id].x.lower), NO_GAP_INC);
      if( algns[id].init_sign == 0 ) {
        old = algns[id].y.lower;
        if( cur > old ) algns[id].yl_diff = cur - old;
      }
      else if( algns[id].init_sign == 1 ) {
        old = algns[id].y.upper;
        if( cur < old ) algns[id].yr_diff = old - cur;
      }
    }
	}
	else if( mode == XR_DIFF ) {
    cur_diff = algns[id].xr_diff;
    if( cur_diff < (algns[id].x.upper - point) ) {
      algns[id].xr_diff = algns[id].x.upper - point;
      cur = find_yloc_one(algns[id], f, abs(point - algns[id].x.lower), NO_GAP_INC);
      if( algns[id].init_sign == 0 ) {
        old = algns[id].y.upper;
        if( cur < old ) algns[id].yr_diff = old - cur;
      }
      else if( algns[id].init_sign == 1 ) {
        old = algns[id].y.lower;
        if( cur > old ) algns[id].yl_diff = cur - old;
      }
    }
  }
  else if( mode == YL_DIFF ) {
    cur_diff = algns[id].yl_diff;
    if( cur_diff < abs(point - algns[id].y.lower) ) {
      algns[id].yl_diff = abs(point - algns[id].y.lower);
      if( algns[id].init_sign == 0 ) {
        loc = find_yloc_one(algns[id], f, abs(point - algns[id].y.lower), GAP_INC_IN_Y);
        cur = find_xloc_one(algns[id], f, loc, GAP_INC);
        old = algns[id].x.lower;
        if( cur > old ) algns[id].xl_diff = cur - old;
      }
      else if( algns[id].init_sign == 1 ) {
        loc = find_yloc_one(algns[id], f, abs(algns[id].y.upper - point), GAP_INC_IN_Y);
        cur = find_xloc_one(algns[id], f, loc, GAP_INC);
        old = algns[id].x.upper;
        if( cur < old ) algns[id].xr_diff = old - cur;
      }
    }
	}
	else if( mode == YR_DIFF ) {
    cur_diff = algns[id].yr_diff;
    if( cur_diff < (algns[id].y.upper - point) ) {
      algns[id].yr_diff = algns[id].y.upper - point;
      if( algns[id].init_sign == 0 ) {
        loc = find_yloc_one(algns[id], f, abs(point - algns[id].y.lower), GAP_INC_IN_Y);
        cur = find_xloc_one(algns[id], f, loc, GAP_INC);
        old = algns[id].x.upper;
        if( cur < old ) algns[id].xr_diff = old - cur;
      }
      else if( algns[id].init_sign == 1 ) {
        loc = find_yloc_one(algns[id], f, abs(algns[id].y.upper - point), GAP_INC_IN_Y);
        cur = find_xloc_one(algns[id], f, loc, GAP_INC);
        old = algns[id].x.lower;
        if( cur > old ) algns[id].xl_diff = cur - old;
      }
    }
	}
}

int count_nucs(struct DotList algn, FILE *fp, int num_nu, int *num1, int *num2)
{
	int i = 0, j = 0;
	struct b_list *a_info;
	int count1 = 0, count2 = 0, count3 = 0;

	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;
  a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;

	strcpy(S1, "");
	strcpy(T1, "");
	get_nth_algn(S1, T1, algn.fid, 0, fp, a_info, REG);
	free(a_info);

	while( (i < BIG) && (j < num_nu) && (S1[i] != '\0') && (S1[i] != '\n')) {
		if( strchr("ACGTN", toupper(S1[i])) ) j++;

		if( strchr("ACGT", toupper(S1[i])) ) count1++;
		if( strchr("ACGT", toupper(S1[i])) && strchr("ACGT", toupper(T1[i]) ) ) {
			count2++;
			if( S1[i] == T1[i] ) {
				count3++;
			}
		}
		i++;
	}
	
	*num1 = count1;
	*num2 = count2;
	return(count3);
}
