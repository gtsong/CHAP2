#include "main.h"
#include "rechain.h"
#include "util.h"
#include "regions.h"
#include "util_i.h"
#include "read_maf.h"
#include "read_algn.h"

extern int debug_mode;

void restore_chained_algns(struct slist *st, struct DotList *init_algns, int num_algns, char *sp1_name, char *sp2_name, int len1, int len2, FILE *fp)
{
	int e1 = 0, e2 = 0;
	int old_e1 = 0, old_e2 = 0;
	int cur_b = 0;
	int cur_col = 0;
	int i = 0, j = 0;
	int cur_id = 0;
	struct b_list *a_info;
	int num_gaps = 0;
	int gap1 = 0, gap2 = 0;
	char S1[BIG], T1[BIG];
	char S2[BIG], T2[BIG];
	int beg = 0, end = 0;
	int len = 0;
	int cur_bid = 0;
	int b1 = 0, b2 = 0;

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

	i = 0;
	while( i < num_algns ) 
	{
		cur_id = st[i].id;
    get_nth_algn(S2, T2, init_algns[cur_id].fid, beg, fp, a_info, REG);
		end = strlen(S2);
    (*a_info).len1 = len1;
    (*a_info).len2 = len2;
		b1 = (*a_info).b1;
		b2 = (*a_info).b2;
		(*a_info).b1 = init_algns[cur_id].x.lower - 1;

    if( init_algns[cur_id].sign == 0 ) {
      (*a_info).b2 = init_algns[cur_id].y.lower - 1;
      (*a_info).strand = '+';
    }
    else if( init_algns[cur_id].sign == 1 ) {
      (*a_info).b2 = len2 - init_algns[cur_id].y.upper + 1;
      (*a_info).strand = '-';
    }
		cur_bid = init_algns[cur_id].fid;

		cur_b = 0;
		cur_col = 0;
		e1 = 0;
		e2 = 0;
		j = 0;
    while( (i < num_algns) && ((*a_info).b1 >= (b1+e1)) && (j < end) && (S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n') )
		{
	    if( strchr("ACGTN", toupper(S2[j])) ) e1++;
      if( strchr("ACGTN", toupper(T2[j])) ) e2++;
			j++;
		}

		old_e1 = e1;
		old_e2 = e2;
		cur_b = j;
    while( (i < num_algns) && ((b1+e1) < init_algns[st[i].id].x.upper) && (j < end) && (S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n') )
		{
	    if( strchr("ACGTN", toupper(S2[j])) ) e1++;
      if( strchr("ACGTN", toupper(T2[j])) ) e2++;
			j++;
			
			if( ((i+1) < num_algns) &&  ((b1+e1) >= init_algns[st[i].id].x.upper) && (cur_bid == init_algns[st[i+1].id].fid) ) {
				if( i < num_algns ) {
					gap1 = 0;
					gap2 = 0;
					num_gaps = 0;
    			while( (j < end) && ((S2[j] == '-') || (S2[j] == 'N') || (T2[j] == '-') || (T2[j] == 'N'))) {
	    			if( strchr("ACGTN", toupper(S2[j])) ) gap1++;
      			if( strchr("ACGTN", toupper(T2[j])) ) gap2++;
						j++;
						num_gaps++;
					}
					
					if( (b1 + e1 + gap1) >= init_algns[st[i+1].id].x.lower ) {
						e1 = e1 + gap1;
						e2 = e2 + gap2;
						i++;
					}
					else {
						cur_col = j - num_gaps;

      			(*a_info).e1 = e1 - old_e1;
 			     	(*a_info).e2 = e2 - old_e2;

						strcpy(S1,S2);
						strcpy(T1,T2);
						S1[cur_col]='\0';	
						T1[cur_col]='\0';	
  		    	if( ((a_info->b1)-1+(a_info->e1)) > (a_info->len1) ) {
    		    	fatalf("Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
     			 	}

  		    	printf("\na srcblock=%d\n", init_algns[cur_id].fid);
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

						i++;
						cur_id = st[i].id; // i is already increased by 1
						b1 = b1 + e1 + gap1;
						(*a_info).b1 = init_algns[cur_id].x.lower - 1;
    				if( init_algns[cur_id].sign == 0 ) {
      				(*a_info).b2 = init_algns[cur_id].y.lower - 1;
    				}
    				else if( init_algns[cur_id].sign == 1 ) {
      				(*a_info).b2 = len2 - init_algns[cur_id].y.upper + 1;
						}
						e1 = 0;
						e2 = 0;
    				while( (i < num_algns) && ((*a_info).b1 >= (b1+e1)) && (j < end) && (S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n') )
						{
	    				if( strchr("ACGTN", toupper(S2[j])) ) e1++;
      				if( strchr("ACGTN", toupper(T2[j])) ) e2++;
							j++;
						}
						cur_b = j;
						old_e1 = e1;
						old_e2 = e2;
					}
				}
			}
		}

		cur_col = j;
    (*a_info).e1 = e1-old_e1;
 		(*a_info).e2 = e2-old_e2;

		strcpy(S1,S2);
		strcpy(T1,T2);
		S1[cur_col]='\0';	
		T1[cur_col]='\0';	
  	if( ((a_info->b1)-1+(a_info->e1)) > (a_info->len1) ) {
     	fatalf("Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
    }

  	printf("\na srcblock=%d\n", init_algns[cur_id].fid);
    printf("s %s %d %d + %d %s", sp1_name, (a_info->b1), a_info->e1, a_info->len1, S1+cur_b);
    len = strlen(S1);
    if( S1[len-1] != '\n' ) printf("\n");

   	if( ((a_info->len2) - (a_info->b2) - (a_info->e2)) < 0 ) {
    	fatalf("Bad coordinates %d %d len %d\n", a_info->b2-1, a_info->e2, a_info->len2);
    }

   	printf("s %s %d %d %c %d %s", sp2_name, a_info->b2, a_info->e2, a_info->strand, a_info->len2, T1+cur_b);
   	len = strlen(T1);
   	if( T1[len-1] != '\n' ) printf("\n");

		i++;
	}

	free(a_info);
}

void merge_overlaps(struct slist *st, struct DotList *init_algns, int num_algns)
{
	int i = 0, j = 0;
	int cur_id = 0, cmp_id = 0;
	int lo = 0, hi = 0;

	while( i < num_algns ) {
		cur_id = st[i].id;
		j = 1;
		cmp_id = st[i+j].id;
		while( ((i+j) < num_algns) && (init_algns[cur_id].fid == init_algns[cmp_id].fid) && (proper_overlap(init_algns[cur_id].x, init_algns[cmp_id].x) == true) ) {
			if( init_algns[cur_id].x.lower < init_algns[cmp_id].x.lower ) {
				lo = init_algns[cur_id].x.lower;
			}
			else lo = init_algns[cmp_id].x.lower;

			if( init_algns[cur_id].x.upper < init_algns[cmp_id].x.upper ) {
				hi = init_algns[cmp_id].x.upper;
			}
			else hi = init_algns[cur_id].x.upper;
			init_algns[cur_id].x = assign_I(lo, hi);

			if( init_algns[cur_id].y.lower < init_algns[cmp_id].y.lower ) {
				lo = init_algns[cur_id].y.lower;
			}
			else lo = init_algns[cmp_id].y.lower;

			if( init_algns[cur_id].y.upper < init_algns[cmp_id].y.upper ) {
				hi = init_algns[cmp_id].y.upper;
			}
			else hi = init_algns[cur_id].y.upper;
			init_algns[cur_id].y = assign_I(lo, hi);
			init_algns[cmp_id].sign = DELETED;

			j++;
			if( (i+j) < num_algns ) {
				cmp_id = st[i+j].id;
			}
		}

		i = i + j;
	}
}
