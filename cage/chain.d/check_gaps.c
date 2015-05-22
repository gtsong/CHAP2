#include "main.h"
#include "regions.h"
#include "check_gaps.h"
#include "find_merging.h"
#include "util_gen.h"
#include "check_repeats.h"
#include "util_i.h"
#include "util.h"
#include "util_algns.h"

int get_score(struct DotList *dots, struct gap_list gps, float *d_rate, struct r_list *rp1, int num_rp1, struct r_list *rp2, int num_rp2, int *id1, int *id2, FILE *fp)
{
	int score = -1;
	int temp_score = -1;
	int pid1 = 0, pid2 = 0;
	int lid = 0, rid = 0;
	struct gap_list temp_gps = {0, 0, 0, 0, 0, 0, 0, 0, 0, "", ""}; 

	temp_gps = assign_glist(gps);
	*id1 = -1;
	*id2 = -1;
	if( temp_gps.type == 0 ) 
	{
		if( dots[temp_gps.id1].x.lower <= dots[temp_gps.id2].x.lower ) 
		{
			lid = temp_gps.id1;
			rid = temp_gps.id2;
		}
		else {
			rid = temp_gps.id1;
			lid = temp_gps.id2;
		}

		if( width(dots[lid].x) >= 500 ) {
			pid1 = (int)(cal_pid_part_algn(dots, lid, width(dots[lid].x)-400, 100, fp, SELF1)+0.5);
		}
		else {
			pid1 = dots[lid].identity;
		}

		if( width(dots[rid].x) >= 500 ) {
			pid2 = (int)(cal_pid_part_algn(dots, rid, 100, width(dots[rid].x)-400, fp, SELF1)+0.5);
		}
		else {
			pid2 = dots[rid].identity;
		}

		score = abs(pid1 - pid2);
		if( score <= DIF_TH ) return(score);
		else return(-1);
	}
	else if( temp_gps.type == -1 ) return(-1);
	else if( temp_gps.type == 3)
	{

		if( (strcmp(dots[temp_gps.id1].name1, dots[temp_gps.id2].name1) != 0) || ((strcmp(dots[temp_gps.id1].name2, dots[temp_gps.id2].name2) != 0)) ) {
			printf("error: two alignments from different contigs e.g. %s-%s and %s-%s in check_gaps.c:55\n", dots[temp_gps.id1].name1, dots[temp_gps.id2].name1, dots[temp_gps.id1].name2, dots[temp_gps.id2].name2);
		}

		strcpy(temp_gps.name1, dots[temp_gps.id1].name1);
		strcpy(temp_gps.name2, dots[temp_gps.id1].name2);
		if( (score = check_insertion(dots, temp_gps, d_rate, rp2, num_rp2, Y_MODE, id2)) != -1 ) 
		{
			temp_score = score;
			if( (score = check_insertion(dots, temp_gps, d_rate, rp1, num_rp1, X_MODE, id1)) != -1 ) 
			{
				if( temp_score < score ) return(temp_score);
				else return(score);
			}
			else return(-1);
		}
		else return(-1);
	}
	else if( (temp_gps.type == 1) || (temp_gps.type == 21) ) {
		if( (strcmp(dots[temp_gps.id1].name1, dots[temp_gps.id2].name1) != 0) || ((strcmp(dots[temp_gps.id1].name2, dots[temp_gps.id2].name2) != 0)) ) {
			printf("error: two alignments from different contigs e.g. %s-%s and %s-%s in check_gaps.c:74\n", dots[temp_gps.id1].name1, dots[temp_gps.id2].name1, dots[temp_gps.id1].name2, dots[temp_gps.id2].name2);
		}

		strcpy(temp_gps.name1, dots[temp_gps.id1].name1);
		strcpy(temp_gps.name2, dots[temp_gps.id1].name2);
		if( (score = check_insertion(dots, temp_gps, d_rate, rp2, num_rp2, Y_MODE, id2)) != -1 ) return(score);
		else return(-1);
	}
	else if( (temp_gps.type == 2) || (temp_gps.type == 22) ) {
		if( (strcmp(dots[temp_gps.id1].name1, dots[temp_gps.id2].name1) != 0) || ((strcmp(dots[temp_gps.id1].name2, dots[temp_gps.id2].name2) != 0)) ) {
			printf("error: two alignments from different contigs e.g. %s-%s and %s-%s in check_gaps.c:84\n", dots[temp_gps.id1].name1, dots[temp_gps.id2].name1, dots[temp_gps.id1].name2, dots[temp_gps.id2].name2);
		}

		strcpy(temp_gps.name1, dots[temp_gps.id1].name2);
		strcpy(temp_gps.name2, dots[temp_gps.id1].name1);
		if( (score = check_insertion(dots, temp_gps, d_rate, rp1, num_rp1, Y_MODE, id1)) != -1 ) return(score);
		else return(-1);
	}
	else return(-1);
}

void merge_two_alignments(struct DotList *dots, struct gap_list *gps, int num_gaps, struct r_list *rp, int num_rp)
{
	int i = 0;
	int *temp_val;
	float *d_rate;
	int *id;

	temp_val = (int *) ckalloc(sizeof(int));
	id = (int *) ckalloc(sizeof(int));
	d_rate = (float *) ckalloc(sizeof(float));

	*d_rate = 100;

	for( i = 0; i < num_gaps; i++ )
	{
		if( gps[i].type == 0 ) 
		{
			if( (abs(dots[gps[i].id1].identity - dots[gps[i].id2].identity)) <= DIF_TH )
			{
				merging_step(dots, gps[i].id1, gps[i].id2);
				adjust_after_merging(gps, num_gaps, i);
			}
		}
		else if( gps[i].type == -1 ) 
		{
		}
		else
		{
			if( check_insertion(dots, gps[i], d_rate, rp, num_rp, Y_MODE, id) != -1 )
			{
				merging_step(dots, gps[i].id1, gps[i].id2);
				adjust_after_merging(gps, num_gaps, i);
			}
		}
	}

	free(d_rate);
	free(id);
	free(temp_val);
}

int check_insertion(struct DotList *dots, struct gap_list gp, float *d_rate, struct r_list *rp_list, int num_list, int flag, int *id)
{
	int num_rp = -1;
	struct r_list *rp = NULL, *temp_rp = NULL;
	int num_temp = 0;
	int pid_rp = 100, pid_algns = -1;
	int score = -1;
	struct I x1 = {0, 1}, x2 = {0, 1};
	float cr = -1;
	int width_gap = -1;
	int new_len = -1;
	int from = -1, to = -1;
	int temp_from = -1, temp_to = -1;
	int cur_id = -1, start = -1, end = 0;
	float cur_rate = 100;
	int i = 0, j = 0;
	char name[LEN_NAME] = "";
	
	*d_rate = 100; 

	strcpy(name, "");
	if( flag == Y_MODE ) {
		strcpy(name, gp.name2);
	}
	else strcpy(name, "NONE");
	if( flag == X_MODE ) {
		strcpy(name, gp.name1);
	}
	else strcpy(name, "NONE");

	x1 = assign_I(0,1);
	x2 = assign_I(0,1);

	if( num_list > 0 ) {
		rp = (struct r_list *) ckalloc(sizeof(struct r_list) * num_list);
		temp_rp = (struct r_list *) ckalloc(sizeof(struct r_list) * num_list);
	}

	j = 0;
	for( i = 0; i < num_list; i++ ) {
		rp[i].id = -1;
		rp[i].start = -1;
		rp[i].end = 0;
		rp[i].d_rate = 100;
		strcpy(rp[i].name, "");

		if( strcmp(rp_list[i].name, name) == 0 ) {
			temp_rp[j].id = rp_list[i].id;
			temp_rp[j].start = rp_list[i].start;
			temp_rp[j].end = rp_list[i].end;
			temp_rp[j].d_rate = rp_list[i].d_rate;
			strcpy(temp_rp[j].name, rp_list[i].name);
			j++;
		}
	}

	num_temp = j;

	if( flag == Y_MODE ) 
	{
		num_rp = obtain_repeats_list(gp.y1, gp.y2, rp, temp_rp, num_temp);	
		width_gap = gp.y2 - gp.y1;
		to = gp.y2;
		from = gp.y1;
	}
	else if( flag == X_MODE )
	{
		num_rp = obtain_repeats_list(gp.x1, gp.x2, rp, temp_rp, num_temp);	
		width_gap = gp.x2 - gp.x1;
		to = gp.x2;
		from = gp.x1;
	}

	if( (num_rp == -1) || (num_rp == 0) || (num_rp > 1)) {
    if( (loosen_subset(dots[gp.id1].x, dots[gp.id2].x) == true) || (loosen_subset(dots[gp.id2].x, dots[gp.id1].x) == true) || (loosen_subset(dots[gp.id1].y, dots[gp.id2].y) == true) || (loosen_subset(dots[gp.id2].y, dots[gp.id1].y) == true) ) {}
    else if( proper_overlap(dots[gp.id1].x, dots[gp.id2].x) == true ) {
      new_len = width(intersect(dots[gp.id1].x, dots[gp.id2].x));
    }
    else if( proper_overlap(dots[gp.id1].y, dots[gp.id2].y) == true ) {
      new_len = width(intersect(dots[gp.id1].y, dots[gp.id2].y));
		}

		if( (gp.type == 21) || (gp.type == 22) ) {
      from = gp.y1 - abs(gp.y2 - gp.y1);
      to = gp.y1;
      num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
      if( (num_rp == -1) || (num_rp == 0) || (num_rp > 1)) {
        from = gp.y1 - (abs(gp.y2 - gp.y1)/2);
        to = gp.y1 + (abs(gp.y2 - gp.y1)/2);
        num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
      }
    }
		else if( new_len > 0 ) {
      from = gp.y1;
      to = gp.y2 + new_len;
      num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
      if( (num_rp == -1) || (num_rp == 0) || (num_rp > 1)) {
        from = gp.y1 - new_len;
        to = gp.y2;
        num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
        if( (num_rp == -1) || (num_rp == 0) || (num_rp > 1)) {
          from = gp.y1 - (new_len/2);
          to = gp.y2 + (new_len/2);
          num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
        }
      }
    }
	}
	else if( (num_rp == 1) && ( (gp.type == 21) || (gp.type == 22) ) )
	{
		cur_id = rp[0].id;
		start = rp[0].start;
		end = rp[0].end;
		cur_rate = rp[0].d_rate;
		temp_from = from;
		temp_to = to;

    from = gp.y1 - abs(gp.y2 - gp.y1);
    to = gp.y1;
    num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
		if( (num_rp == 0) || (num_rp == -1) || (num_rp > 2) || ((num_rp == 1) && (rp[0].d_rate > cur_rate))) {
			rp[0].id = cur_id;
			rp[0].start = start;
			rp[0].end = end;
			rp[0].d_rate = cur_rate;
			from = temp_from;
			to = temp_to;
			num_rp = 1;
		}

		cur_id = rp[0].id;
		start = rp[0].start;
		end = rp[0].end;
		cur_rate = rp[0].d_rate;
		temp_from = from;
		temp_to = to;

    from = gp.y1 - (abs(gp.y2 - gp.y1)/2);
    to = gp.y1 + (abs(gp.y2 - gp.y1)/2);
    num_rp = obtain_repeats_list(from, to, rp, temp_rp, num_temp);
		if( (num_rp == 0) || (num_rp == -1) || (num_rp > 2) || ((num_rp == 1) && (rp[0].d_rate > cur_rate))) {
			rp[0].id = cur_id;
			rp[0].start = start;
			rp[0].end = end;
			rp[0].d_rate = cur_rate;
			from = temp_from;
			to = temp_to;
			num_rp = 1;
		}
	}

	if( (num_rp == -1) || (num_rp == 0) )
	{
		*id = -1;
		if( num_list > 0 ) {
			free(temp_rp);
			free(rp);
		}
		return(-1);
	}
	else if( num_rp > 1 )
	{
		*id = -1;
		if( num_list > 0 ) {
			free(temp_rp);
			free(rp);
		}
		return(-1);
	}
	else if( num_rp == 1 )
	{
		if( abs(rp[0].end - rp[0].start) > abs(to - from) ) width_gap = abs(rp[0].end - rp[0].start);
		else width_gap = abs(to - from);

		cr = ((float)abs(to - rp[0].end) + (float)(from - rp[0].start))/((float)width_gap);
		x1 = assign_I(dots[gp.id1].x.lower, dots[gp.id1].x.upper);
		x2 = assign_I(dots[gp.id2].x.lower, dots[gp.id2].x.upper);

		pid_algns = ((dots[gp.id1].identity * width(x1)) + (dots[gp.id2].identity * width(x2))) / (width(x1) + width(x2));

		score = abs((dots[gp.id1].identity) - (dots[gp.id2].identity));
		pid_rp = 100 - ((int) rp[0].d_rate);
		if( (*d_rate) > rp[0].d_rate ) *d_rate = rp[0].d_rate;

		if((cr <= RP_TH) && (pid_rp >= RP_P_IDT))
		{
			*id = rp[0].id;
			if( num_list > 0 ) {
				free(temp_rp);
				free(rp);
			}
			return( score );
		}
		else if( pid_rp >= (pid_algns - RP_DIFF) ) 
		{
			*id = rp[0].id;
			if( num_list > 0 ) {
				free(temp_rp);
				free(rp);
			}
			return( score );
		}
		else 
		{
			*id = -1;
			if( num_list > 0 ) {
				free(temp_rp);
				free(rp);
			}
			return(-1);
		}
	}
	if( num_list > 0 ) {
		free(temp_rp);
		free(rp);
	}
	return(score);
}

void check_gap(struct DotList *dots, struct gap_list *gps, int num_gaps, int *temp_list, int num_list)
{
	int i, j;
	int pt1, pt2;
	int num_merge = 0;
	int type = 0;

	for( i = 0; i < num_list; i++ )
	{
		pt1 = temp_list[i];
		if( gps[pt1].type != -1 )
		{
			if( gps[pt1].type == 0 )
			{
				merging_step(dots, gps[pt1].id1, gps[pt1].id2);
				adjust_after_merging(gps, num_gaps, pt1);
			}
			else if((gps[pt1].type == 1) || (gps[pt1].type == 2))
			{
				for( j = (i+1); j < num_list; j++ )
				{
					pt2 = temp_list[j];
					if( gps[pt2].type != -1 )
					{
						type = 1;
						if( (type == 1) || (type == 2) )
						{
							merging_step(dots, gps[pt2].id1, gps[pt2].id2);
							adjust_after_merging(gps, num_gaps, pt1);
							num_merge++;
						}
					}
				}
				if( num_merge > 0 ) 
				{
					merging_step(dots, gps[pt1].id1, gps[pt1].id2);
					adjust_after_merging(gps, num_gaps, pt1);
					num_merge = 0;
				}	
			}
		}
	}
}
