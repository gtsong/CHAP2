#include "main.h"
#include "regions.h"
#include "chain_pair_alg.h"
#include "find_merging.h"
#include "deal_gaps.h"
#include "util_gen.h"
#include "kd_op.h"
#include "find_dup_copy.h"
#include "util.h"
#include "update_init_algns.h"

extern int debug_mode;

void chain_pair_alg_update_init(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct DotList *init_algns)
{
	struct slist *sorted;
	int i = 0;
	int num_lines = 0;
	int xval = 0, yval = 0;
	int w_sid = 0, w_fid = 0, h_sid = 0, h_fid = 0;
	int opt_id = -1;
	bool *is_x;

	is_x = (bool *) malloc(sizeof(bool));
	num_lines = *num;

	sorted = (struct slist *) ckalloc((*num) * sizeof(struct slist));

	if( num_lines > 0 ) sort_list(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		opt_id = -1;
		if( dots[sorted[i].id].sign == DELETED )
		{
		}
		else
		{
			xval = dots[sorted[i].id].x.lower;
			if( dots[sorted[i].id].sign == 0 )
			{
				yval = dots[sorted[i].id].y.lower;
			}
			else 
			{
				yval = dots[sorted[i].id].y.upper;
			}	
			
			w_sid = find_id(dots, tree, size1, sorted[i].id, xval, yval, W_SID);
			w_fid = find_id(dots, tree, size1, sorted[i].id, xval, yval, W_FID);
			h_sid = find_id(dots, tree, size2, sorted[i].id, xval, yval, H_SID);
			h_fid = find_id(dots, tree, size2, sorted[i].id, xval, yval, H_FID);

			opt_id = find_opt_ch_alg(dots, num_lines, sorted[i].id, p_pts, w_sid, w_fid, h_sid, h_fid);
		}

		if( opt_id != -1 ) 
		{
			if( debug_mode == TRUE ) {
				printf("Combined: %d-%d,%d-%d\n", dots[opt_id].x.lower, dots[opt_id].x.upper, dots[opt_id].y.lower, dots[opt_id].y.upper);
				printf("%d-%d,%d-%d\n", dots[sorted[i].id].x.lower, dots[sorted[i].id].x.upper, dots[sorted[i].id].y.lower, dots[sorted[i].id].y.upper); 
			}

			mark_chain(dots, opt_id, sorted[i].id, init_algns);
			merging_step(dots, opt_id, sorted[i].id);
		}	
	}
	overwrite_dots(num, dots);

	free(sorted);
	free(is_x);
}

int chain_pair_alg(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct chain_list *clist)
{
	int num_clist = 0;
	struct slist *sorted;
	int i = 0;
	int num_lines;
	int xval, yval;
	int w_sid, w_fid, h_sid, h_fid;
	int opt_id;
	bool *is_x;
	int *sd;
	int cur_num;
	int j, cur_id, num_list;

	is_x = (bool *) malloc(sizeof(bool));
	sd = (int *) malloc(sizeof(int));
	num_lines = *num;

	sorted = (struct slist *) ckalloc(num_lines * sizeof(struct slist));
	sort_list(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
		opt_id = -1;
		if( dots[sorted[i].id].sign == 2 )
		{
		}
		else
		{
			xval = dots[sorted[i].id].x.lower;
			if( dots[sorted[i].id].sign == 0 )
			{
				yval = dots[sorted[i].id].y.lower;
			}
			else 
			{
				yval = dots[sorted[i].id].y.upper;
			}	
			
			w_sid = find_id(dots, tree, size1, sorted[i].id, xval, yval, W_SID);
			w_fid = find_id(dots, tree, size1, sorted[i].id, xval, yval, W_FID);
			h_sid = find_id(dots, tree, size2, sorted[i].id, xval, yval, H_SID);
			h_fid = find_id(dots, tree, size2, sorted[i].id, xval, yval, H_FID);

			opt_id = find_opt_ch_alg(dots, num_lines, sorted[i].id, p_pts, w_sid, w_fid, h_sid, h_fid);
		}

		if( opt_id != -1 ) 
		{
			if((dots[opt_id].c_id == -1) && (dots[sorted[i].id].c_id == -1))
			{
				dots[opt_id].c_id = num_clist;
				clist[num_clist].c_id = num_clist;
				cur_num = 0;
				clist[num_clist].fid[cur_num] = dots[sorted[i].id].fid;
				clist[num_clist].num_algs = 1;
				num_clist++;
			}
			else if(dots[opt_id].c_id == -1)
			{
				dots[opt_id].c_id = dots[sorted[i].id].c_id;
				cur_num = clist[dots[opt_id].c_id].num_algs;
				clist[dots[opt_id].c_id].fid[cur_num] = dots[sorted[i].id].fid;
				clist[dots[opt_id].c_id].num_algs = clist[dots[opt_id].c_id].num_algs + 1;
			}
			else if(dots[sorted[i].id].c_id == -1 )
			{
				dots[opt_id].c_id = dots[opt_id].c_id;
				cur_num = clist[dots[opt_id].c_id].num_algs;
				clist[dots[opt_id].c_id].fid[cur_num] = dots[sorted[i].id].fid;
				clist[dots[opt_id].c_id].num_algs = clist[dots[opt_id].c_id].num_algs + 1;
			}
			else
			{
				cur_num = clist[dots[opt_id].c_id].num_algs;
				cur_id = dots[sorted[i].id].c_id;
				num_list = clist[cur_id].num_algs;
				for( j = 0; j < num_list; j++ )
				{
					clist[dots[opt_id].c_id].fid[cur_num] = clist[cur_id].fid[j];
					cur_num++;
				}
				clist[dots[opt_id].c_id].num_algs = cur_num;
				clist[cur_id].c_id = -1;
				clist[cur_id].num_algs = 0;
			}

			printf("Combined: %d-%d,%d-%d\n", dots[opt_id].x.lower, dots[opt_id].x.upper, dots[opt_id].y.lower, dots[opt_id].y.upper);
			printf("%d-%d,%d-%d\n", dots[sorted[i].id].x.lower, dots[sorted[i].id].x.upper, dots[sorted[i].id].y.lower, dots[sorted[i].id].y.upper); 

			merging_step(dots, opt_id, sorted[i].id);
		}	
	}
	overwrite_dots(num, dots);

	free(sorted);
	free(sd);
	free(is_x);
	return(num_clist);
}

int find_opt_ch_alg(struct DotList *dots, int num_lines, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid)
{
	int i;
	int min_score = 1000;
	int max_id = -1;
	bool *is_x;
	int *sd;
	int d;
	struct gap_list gps;
	int temp_score;
	int start, mid1 = -1, mid2 = -1, end;
	float *d_rate;

	is_x = (bool *) malloc(sizeof(bool));
	sd = (int *) malloc(sizeof(int));
	d_rate = (float *) malloc(sizeof(float));

	if( w_sid < h_sid )
	{
		start = w_sid;
		if( w_fid < h_sid )
		{
			mid1 = w_fid;
			mid2 = h_sid;
			end = h_fid;
		}
		else 
		{
			if( w_fid < h_fid ) end = h_fid;
			else end = w_fid;
		}
	}
	else 
	{
		start = h_sid;
		if( h_fid < w_sid )
		{
			mid1 = h_fid;
			mid2 = w_sid;
			end = w_fid;
		}
		else
		{
			if( h_fid < w_fid ) end = w_fid;
			else end = h_fid;
		}
	}

	if( end >= num_lines) end = num_lines - 1;

	for( i = start; i <= end; i++ )
	{
		if( st[i].id == id ) {}	
		else if( dots[st[i].id].x.lower > dots[id].x.lower ) {}
		else if( (dots[st[i].id].sign == 0) && (dots[st[i].id].y.lower > dots[id].y.lower )) {}
		else if( (dots[st[i].id].sign == 1) && (dots[st[i].id].y.lower < dots[id].y.lower )) {}
		else if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
			if((dots[st[i].id].sign != 2) && (dots[st[i].id].sign == dots[id].sign) && ((d = distance(dots, st[i].id, id, is_x, sd)) <= CLOSE_TH))
			{
				if( d <= M_THRESHOLD )
				{
					gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
					if(gps.type == -1) // this gap is meaningless
					{
					}
					else 
					{
						gps.gid = 0;
						temp_score = abs(gps.y1 - gps.y2);
						if( temp_score != -1 )
						{
							if( min_score >= temp_score )
							{
								min_score = temp_score;
								max_id = st[i].id;
							}
						}
					}
				}
				else if( d > M_THRESHOLD )
				{
					if( ((*sd) <= RANGE_TH) && (d <= RANGE_TH) )
					{
						if( check_candi(dots, id, st[i].id, *is_x) == false ) {}
						else
						{
							gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
							if(gps.type == -1) // this gap is meaningless
							{
							}
							else
							{
								gps.gid = 0;
								temp_score = abs(gps.y1 - gps.y2);
								if( temp_score != -1 )
								{
									if( min_score > temp_score )
									{
										min_score = temp_score;
										max_id = st[i].id;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if( max_id == -1 ) 
	{
		free(d_rate);
		free(sd);
		free(is_x);
		return(-1);
	}
	else 
	{
		free(d_rate);
		free(sd);
		free(is_x);
		return(max_id);
	}
}
