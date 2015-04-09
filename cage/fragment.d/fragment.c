#include "main.h"
#include "fragment.h"
#include "util.h"
#include "read_maf.h"
#include "read_algn.h"

extern int debug_mode;

void fragment_algns(struct DotList *init_algns, int num_algns, char *sp1_name, char *sp2_name, int len1, int len2, FILE *fp, int cut_th)
{
	int e1 = 0, e2 = 0;
	int next_e1 = 0, next_e2 = 0;
	int cur_b = 0, next_b = 0;
	int cur_col = 0, next_col = 0;
	int i = 0, j = 0;
	struct b_list *a_info;
	int cur_len = 0, next_len = 0;
	int pid1 = 0, pid2 = 0;
	int num_gaps = 0;
	int gap1 = 0, gap2 = 0;
	char S1[BIG], T1[BIG];
	char S2[BIG], T2[BIG];
	int beg = 0, end = 0;
	int len = 0;
	bool is_N = false;

  a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
  a_info->b1 = 0;
  a_info->e1 = 1;
  a_info->len1 = 0;  
	a_info->b2 = 0;
  a_info->e2 = 1;
  a_info->len2 = 0;
  a_info->strand = '+';
  a_info->pid = 0;
		
	printf("##maf version=1 scoring=lastz-pid\n");

	for( i = 0; i < num_algns; i++ ) 
	{
    get_nth_algn(S2, T2, init_algns[i].fid, beg, fp, a_info, REG);
		end = strlen(S2);
    (*a_info).len1 = len1;
    (*a_info).len2 = len2;
		(*a_info).b1 = init_algns[i].x.lower - 1;

    if( init_algns[i].sign == 0 ) {
      (*a_info).b2 = init_algns[i].y.lower-1;
      (*a_info).strand = '+';
    }
    else if( init_algns[i].sign == 1 ) {
      (*a_info).b2 = len2 - init_algns[i].y.upper + 1;
      (*a_info).strand = '-';
    }

		cur_len = 0;
		cur_b = 0;
		e1 = 0;
		e2 = 0;
		cur_col = 0;
    for( j = 0; (j < end) && ((S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n') && (S2[j] != 'N') && (T2[j] != 'N')) && ((cur_len <= (DEL_TH * cut_th)) || ((S2[j] != '-') && (T2[j] != '-'))); j++ ) {
	    if( strchr("ACGTN", toupper(S2[j])) ) e1++;
      if( strchr("ACGTN", toupper(T2[j])) ) e2++;
			if( strchr("ACGT", toupper(S2[j])) && strchr("ACGT", toupper(T2[j])) ) cur_len++;
    }
		cur_col = j;

		while( (j < end) && (S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n') ) 
		{
			num_gaps = 0;
			gap1 = 0;
			gap2 = 0;
			is_N = false;
    	while( (j < end) && ((S2[j] == '-') || (S2[j] == 'N') || (T2[j] == '-') || (T2[j] == 'N'))) {
	    	if( strchr("ACGTN", toupper(S2[j])) ) gap1++;
      	if( strchr("ACGTN", toupper(T2[j])) ) gap2++;
				if( (toupper(S2[j]) == 'N') || (toupper(T2[j]) == 'N') ) {
					is_N = true;
				}
				j++;
				num_gaps++;
			}
			next_b = j;
			next_e1 = 0;
			next_e2 = 0;
			next_len = 0;	
    	while( (j < end) && ((S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n')) && (((S2[j] != 'N') && (T2[j] != 'N')) && ((next_len <= (DEL_TH * cut_th)) || ((S2[j] != '-') && (T2[j] != '-') ))) ) {
   	  	if( strchr("ACGTN", toupper(S2[j])) ) next_e1++;
   	  	if( strchr("ACGTN", toupper(T2[j])) ) next_e2++;
				if( strchr("ACGT", toupper(S2[j])) && strchr("ACGT", toupper(T2[j])) ) next_len++;
				j++;
    	}

			next_col = j;
			pid1 = (int)(cal_pid_maf_beg(S2, T2, cur_b, cur_col));
			pid2 = (int)(cal_pid_maf_beg(S2, T2, next_b, next_col));

			if( (e1 > 0) && (e2 > 0) && ((is_N == true) || (abs(pid1 - pid2) > (PID_TH*cut_th))) ) {
      	(*a_info).e1 = e1;
      	(*a_info).e2 = e2;

				strcpy(S1,S2);
				strcpy(T1,T2);
				S1[cur_col]='\0';	
				T1[cur_col]='\0';	
      	if( ((a_info->b1)-1+(a_info->e1)) > (a_info->len1) ) {
        	fatalf("Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
      	}

      	printf("a srcblocks=%d\n", init_algns[i].fid);
      	printf("s %s %d %d + %d %s", sp1_name, (a_info->b1), a_info->e1, a_info->len1, S1+cur_b);
      	len = strlen(S1);
      	if( S1[len-1] != '\n' ) printf("\n");

      	if( ((a_info->len2) - (a_info->b2) - (a_info->e2)) < 0 ) {
        	fatalf("Bad coordinates %d %d len %d\n", a_info->b2-1, a_info->e2, a_info->len2);
      	}

      	printf("s %s %d %d %c %d %s", sp2_name, a_info->b2, a_info->e2, a_info->strand, a_info->len2, T1+cur_b);
      	len = strlen(T1);
      	if( T1[len-1] != '\n' ) printf("\n");
      	printf("\n");

				e1 = next_e1;
				e2 = next_e2;
				(*a_info).b1 = (*a_info).b1 + (*a_info).e1 + gap1;
				if( init_algns[i].sign == 0 ) {
					(*a_info).b2 = (*a_info).b2 + (*a_info).e2 + gap2;
				}
				else if( init_algns[i].sign == 1 ) {
					(*a_info).b2 = (*a_info).b2 - (*a_info).e2 - gap2;
				}

				cur_b = next_b;
				cur_col = next_col;
				cur_len = next_len;
			}
			else {
				cur_col = next_col; 
				e1 = e1 + gap1 + next_e1;
				e2 = e2 + gap2 + next_e2;
				cur_len = cur_len + next_len;
			}
		}

		j = cur_col;
		next_e1 = 0;
		next_e2 = 0;
  	while ((j < end) && ((S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n')) ) {
   		if( strchr("ACGTN", toupper(S2[j])) ) next_e1++;
   	 	if( strchr("ACGTN", toupper(T2[j])) ) next_e2++;
			j++;
    }

		e1 = e1 + next_e1;
		e2 = e2 + next_e2;
		cur_col = j;
		if( (e1 > 0) && (e2 > 0) ) {
 			(*a_info).e1 = e1;
 			(*a_info).e2 = e2;

			strcpy(S1,S2);
			strcpy(T1,T2);
			S1[cur_col]='\0';	
			T1[cur_col]='\0';	
 			if( ((a_info->b1)-1+(a_info->e1)) > (a_info->len1) ) {
 		  	fatalf("Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
 			}

 			printf("a srcblocks=%d\n", init_algns[i].fid);
 			printf("s %s %d %d + %d %s", sp1_name, (a_info->b1), a_info->e1, a_info->len1, S1+cur_b);
 			len = strlen(S1);
 			if( S1[len-1] != '\n' ) printf("\n");

 			if( ((a_info->len2) - (a_info->b2) - (a_info->e2)) < 0 ) {
 		  	fatalf("Bad coordinates %d %d len %d\n", a_info->b2-1, a_info->e2, a_info->len2);
 			}

 			printf("s %s %d %d %c %d %s", sp2_name, a_info->b2, a_info->e2, a_info->strand, a_info->len2, T1+cur_b);
 			len = strlen(T1);
		 	if( T1[len-1] != '\n' ) printf("\n");
 			printf("\n");
		}
	}
	free(a_info);
}
