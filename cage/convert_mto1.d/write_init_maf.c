#include "main.h"
#include "regions.h"
#include "write_init_maf.h"
#include "read_algn.h"
#include "util_gen.h"
#include "read_maf.h"
#include "util.h"
#include "util_i.h"

// flag is SELF1, SELF2 or PAIR
void write_init_maf(struct DotList *init_algns, int num_algns, char *sp1_name, char *sp2_name, int len1, int len2, FILE *fp, int flag)
{
	char S2[BIG], T2[BIG];
	struct b_list *a_info;
	int i, j;
	int beg, end;
	int pid;
	int y_b, y_e;
	int e1, e2;
	int gap1, gap2;
	int skip_nu;
	int len;

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
		if( (init_algns[i].sp_id == flag) && (init_algns[i].sign != DELETED)) 
		{
			if( init_algns[i].rp1_id != -1 ) {
				beg = init_algns[i].xl_diff+init_algns[i].xl_offset;
				skip_nu = find_xloc_one(init_algns[i], fp, init_algns[i].xl_diff+init_algns[i].xl_offset, NO_GAP_INC);
				end = find_end_xloc(init_algns[i], fp, init_algns[i].x.upper - init_algns[i].x.lower - init_algns[i].xr_diff + init_algns[i].xr_offset);
			}
			else {
				beg = init_algns[i].xl_diff;
				skip_nu = find_xloc_one(init_algns[i], fp, init_algns[i].xl_diff, NO_GAP_INC);
				end = find_end_xloc(init_algns[i], fp, init_algns[i].x.upper - init_algns[i].x.lower - init_algns[i].xr_diff);
			}
			(*a_info).b1 = init_algns[i].x.lower + init_algns[i].xl_diff - 1;
			y_b = init_algns[i].y.lower + init_algns[i].yl_diff;
			y_e = init_algns[i].y.upper - init_algns[i].yr_diff;

			get_nth_algn(S2, T2, init_algns[i].fid, beg, fp, a_info, REG);
			(*a_info).len1 = len1;
			(*a_info).len2 = len2;
			end = end - skip_nu;
			S2[end] = '\0';
			T2[end] = '\0';

			if( init_algns[i].sign == 0 ) {
				(*a_info).b2 = y_b-1;
				(*a_info).strand = '+';
			}
			else if( init_algns[i].sign == 1 ) {
				(*a_info).b2 = (a_info->len2) - y_e + 1;
				(*a_info).strand = '-';
			}

		  e1 = 0;
		  e2 = 0;
		  gap1 = 0;
		  gap2 = 0;
		  for( j = 0; (j < end) && ((S2[j] != '\0') && (S2[j] != '\n') && (T2[j] != '\0') && (T2[j] != '\n')); j++ ) {
    		if( strchr("ACGTN", toupper(S2[j])) ) e1++;
				if( S2[j] == '-' ) gap1++;
    		if( strchr("ACGTN", toupper(T2[j])) ) e2++;
				if( T2[j] == '-' ) gap2++;
  		}

			if( (e1 >= MIN_GAP) && (e2 >= MIN_GAP) ) {
				(*a_info).e1 = e1;
				(*a_info).e2 = e2;
			
   	 	  pid = (int)(cal_pid_maf(S2, T2, end) + 0.5);

				if( ((a_info->b1)-1+(a_info->e1)) > (a_info->len1) ) {
					fatalf("Bad coordinates %d %d len %d\n", a_info->b1-1, a_info->e1, a_info->len1);
				}

     		printf("a srcblock=%d\n", init_algns[i].indiv_fid);
      	printf("s %s %d %d + %d %s", sp1_name, (a_info->b1)-1, a_info->e1, a_info->len1, S2);
      	len = strlen(S2);
      	if( S2[len-1] != '\n' ) printf("\n");

				if( ((a_info->len2) - (a_info->b2) - (a_info->e2)) < 0 ) {
					fatalf("Bad coordinates %d %d len %d\n", a_info->b2-1, a_info->e2, a_info->len2);
				}

      	printf("s %s %d %d %c %d %s", sp2_name, a_info->b2, a_info->e2, a_info->strand, a_info->len2, T2);
      	len = strlen(T2);
      	if( T2[len-1] != '\n' ) printf("\n");
      	printf("\n");
			}
		}
	}

	free(a_info);
}
