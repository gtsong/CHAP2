#include "main.h"
#include "regions.h"
#include "adjust_plot.h"
#include "adjust_plot_genes.h"
#include "find_merging.h"
#include "check_repeats.h"
#include "deal_gaps.h"
#include "util_gen.h"
#include "check_gaps.h"
#include "kd_op.h"
#include "write_maf.h"
#include "read_algn.h"
#include "util_i.h"
#include "util.h"

void adjust_plot_pair_genes(struct DotList *dots, int *num, struct kdnode *tree, struct perm_pt *p_pts, int size1, int size2, struct r_list *rp1, struct r_list *rp2, int num_rp1, int num_rp2, struct g_list *genes1, struct g_list *genes2, int num_genes1, int num_genes2, struct DotList *init_dots, FILE *fp)
{
	struct slist *sorted;
	int i = 0;
	int num_lines;
	int xval, yval;
	int w_sid, w_fid, h_sid, h_fid;
	int opt_id;
	bool *is_x;
	int *sd;
	int *rp1_id, *rp2_id;

	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));
	num_lines = *num;
	sorted = (struct slist *) ckalloc(num_lines * (sizeof(struct slist)));
	rp1_id = (int *) ckalloc(sizeof(int));
	rp2_id = (int *) ckalloc(sizeof(int));

	sort_list(sorted, dots, num_lines);

	for( i = 0; i < num_lines; i++ )
	{
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
			
			w_sid = find_id_pair(dots, tree, size1, size2, sorted[i].id, xval, yval, W_SID);
			w_fid = find_id_pair(dots, tree, size1, size2, sorted[i].id, xval, yval, W_FID);
			h_sid = find_id_pair(dots, tree, size1, size2, sorted[i].id, xval, yval, H_SID);
			h_fid = find_id_pair(dots, tree, size1, size2, sorted[i].id, xval, yval, H_FID);

			opt_id = find_opt_fr_genes(dots, sorted[i].id, p_pts, w_sid, w_fid, h_sid, h_fid, rp1, num_rp1, rp2, num_rp2, genes1, genes2, num_genes1, num_genes2, rp1_id, rp2_id, fp);
		}

		if( opt_id != -1 ) 
		{
			init_dots[sorted[i].id].c_id = init_dots[opt_id].fid;
			init_dots[sorted[i].id].rp1_id = *rp1_id;
			init_dots[sorted[i].id].rp2_id = *rp2_id;
			merging_step(dots, opt_id, sorted[i].id);
		}	
	}
	overwrite_dots(num, dots);

	free(rp1_id);
	free(rp2_id);
	free(sd);
	free(is_x);
	free(sorted);
}

int find_opt_fr_genes(struct DotList *dots, int id, struct perm_pt *st, int w_sid, int w_fid, int h_sid, int h_fid, struct r_list *rp1, int num_rp1, struct r_list *rp2, int num_rp2, struct g_list *genes1, struct g_list *genes2, int num_genes1, int num_genes2, int *rp1_id, int *rp2_id, FILE *fp)
{
	int i = 0;
	int min_score = 1000;
	int max_id = -1;
	bool *is_x;
	int *sd;
	int d = 0;
	struct gap_list gps;
	int temp_score = -1;
	int start, mid1 = -1, mid2 = -1, end;
	float *d_rate;
	float min_rate = 100;
	int len1 = 0, len2 = 0, len = 0, m_th = 0;
	int op_len = 0, op_len_x = 0, op_len_y = 0;
	int closeness = 0;
	struct I temp;
	int y_cur = 0, y_old = 0;
	int *id1, *id2;
	
	is_x = (bool *) ckalloc(sizeof(bool));
	sd = (int *) ckalloc(sizeof(int));
	d_rate = (float *) ckalloc(sizeof(float));
	id1 = (int *) ckalloc(sizeof(int));
	id2 = (int *) ckalloc(sizeof(int));
	*id1 = -1;
	*id2 = -1;
	*rp1_id = -1;
	*rp2_id = -1;
	temp = assign_I(0, 1);

  gps.type = -1;
  gps.id1 = -1;
  gps.id2 = -1;
  gps.x1 = 0;
  gps.x2 = 1;
  gps.y1 = 0;
  gps.y2 = 1;
  strcpy(gps.name1, "");
  strcpy(gps.name2, "");

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

	// m_x and m_y save the coordinated of the initial alignment before getting chained
	for( i = start; i <= end; i++ )
	{
		if( st[i].id == id ) {}	
    else if( (strcmp(dots[st[i].id].name1, dots[id].name1) != 0) || (strcmp(dots[st[i].id].name2, dots[id].name2) != 0) ) {}
		else if( dots[st[i].id].x.lower > dots[id].x.lower ) {}
		else if( dots[st[i].id].x.lower > dots[id].m_x.lower ) {}
		else if( (dots[st[i].id].sign == 0) && (dots[st[i].id].y.lower > dots[id].y.lower )) {}
		else if( (dots[st[i].id].sign == 1) && (dots[st[i].id].y.lower < dots[id].y.lower )) {}
		else if( subset(dots[st[i].id].m_x, dots[id].x) || subset(dots[st[i].id].m_y, dots[id].y) || subset(dots[id].m_x, dots[st[i].id].x) || subset(dots[id].m_y, dots[st[i].id].y) ) {}
		else if( (mid1 != -1) && (i > mid1) && (i < mid2)) {}
		else
		{
// is_x of 'distance' function is true if x region is larger 
			if((dots[st[i].id].sign != 2) && (dots[st[i].id].sign == dots[id].sign) && ((d = distance(dots, st[i].id, id, is_x, sd)) <= MDIS_THRESHOLD))
			{
				len1 = width(dots[st[i].id].x);
				len2 = width(dots[id].x);

				if( len1 > len2 ) len = len2;
				else len = len1;

				if( (len1 >= LG_TH) && (len2 >= LG_TH ) ) m_th = L_M_TH;
				else m_th = M_TH;

				gps.type = -1;
				temp_score = -1;
				if((*sd) <= len)
				{
					op_len = 0;
					op_len_x = 0;
					op_len_y = 0;

       		if( proper_overlap(dots[st[i].id].x, dots[id].x) == true )
          {
            op_len_x = width(intersect(dots[st[i].id].x, dots[id].x));
            op_len = op_len_x;
          }

          if( proper_overlap(dots[st[i].id].y, dots[id].y) == true )
          {
            op_len_y = width(intersect(dots[st[i].id].y, dots[id].y));
            if( op_len < op_len_y ) op_len = op_len_y;
					}

          if( ((*sd) > m_th) || (op_len > m_th) )
          {
            if( (strict_almost_equal(dots[st[i].id].x, dots[id].x) == true) || (strict_almost_equal(dots[st[i].id].y, dots[id].y) == true) )
            {
              gps.type = -1;
            }
						else if( ((closeness = compute_closeness(dots, st[i].id, id)) <= C_OP_TH) && ((*sd) > m_th) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true )) {
							gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
						}
						else if( ((*sd) > m_th) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true )) {
							gps.type = -1;
						}
						else if((proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true ))
						{
							temp = intersect(dots[st[i].id].x, dots[id].x);
							if( dots[id].sign == 0 ) {
								y_cur = dots[id].y.lower;
								y_old = find_yloc_one(dots[st[i].id], fp, temp.lower-dots[st[i].id].x.lower, NO_GAP_INC);
							}
							else if( dots[id].sign == 1 ) {
								y_cur = find_yloc_one(dots[st[i].id], fp, temp.lower-dots[st[i].id].x.lower, NO_GAP_INC);
								y_old = dots[id].y.upper;
							}

							if( y_old >= y_cur ) {
								gps = define_gap_new_type(dots, st[i].id, id, false);
								if( gps.type == -1 ) gps = define_gap_new_type(dots, st[i].id, id, true);
							}
							else if( y_old < y_cur ) {
								gps = define_gap_new_type(dots, st[i].id, id, true);
								if( gps.type == -1 ) gps = define_gap_new_type(dots, st[i].id, id, false);
							}
						}
            else if( proper_overlap(dots[st[i].id].x, dots[id].x) == true )
            {
              gps = define_gap_new_type(dots, st[i].id, id, true);
            }
            else if( proper_overlap(dots[st[i].id].y, dots[id].y) == true )
            {
              gps = define_gap_new_type(dots, st[i].id, id, false);
            }
						else gps.type = -1;
					}
					else if( ((closeness = compute_closeness(dots, st[i].id, id)) > C_OP_TH) && (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true )) 
					{
						closeness = compute_closeness(dots, st[i].id, id);
						temp = intersect(dots[st[i].id].x, dots[id].x);
						if( dots[id].sign == 0 ) {
							y_cur = dots[id].y.lower;
							y_old = find_yloc_one(dots[st[i].id], fp, temp.lower-dots[st[i].id].x.lower, NO_GAP_INC);
						}
						else if( dots[id].sign == 1 ) {
							y_cur = find_yloc_one(dots[st[i].id], fp, temp.lower-dots[st[i].id].x.lower, NO_GAP_INC);
							y_old = dots[id].y.upper;
						}
						if( y_old >= y_cur ) {
							gps = define_gap_new_type(dots, st[i].id, id, false);
							if( gps.type == -1 ) gps = define_gap_new_type(dots, st[i].id, id, true);
						}
						else if( y_old < y_cur ) {
							gps = define_gap_new_type(dots, st[i].id, id, true);
							if( gps.type == -1 ) gps = define_gap_new_type(dots, st[i].id, id, false);
						}
					}
					else if( (check_candi(dots, id, st[i].id, *is_x) == false) && ( (*sd) > TD_TH)) gps.type = -1;
			  	else
					{
            if( (proper_overlap(dots[st[i].id].x, dots[id].x) == true) && (proper_overlap(dots[st[i].id].y, dots[id].y) == true ) ) {
							gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
						}
            else if( proper_overlap(dots[st[i].id].x, dots[id].x) == true )
            {
              gps = define_gap_new_type(dots, st[i].id, id, true);
            }
            else if( proper_overlap(dots[st[i].id].y, dots[id].y) == true )
            {

              gps = define_gap_new_type(dots, st[i].id, id, false);
            }
						else gps = define_gap(dots, st[i].id, id, d, *sd, *is_x);
					}

					if((gps.type == -1) || (gps.type == 3)) // this gap is meaningless
					{
					}
					else
					{
						gps.rp_id1 = -1;
						gps.rp_id2 = -1;
            if( abs(gps.y2 - gps.y1) < ERR_LG_TH ) {
              gps.type = 0;
            }
						temp_score = get_score(dots, gps, d_rate, rp1, num_rp1, rp2, num_rp2, id1, id2, fp);

						if( temp_score == -1 ) {
							if( (gps.type == 21) || (gps.type == 22) ) {
								if( gps.type == 21 ) {
									gps = define_gap_new_type(dots, st[i].id, id, false);
								}
								else if( gps.type == 22 ) {
									gps = define_gap_new_type(dots, st[i].id, id, true);
								}
					
								if((gps.type == -1) || (gps.type == 3)) {}// this gap is meaningless
								else temp_score = get_score(dots, gps, d_rate, rp1, num_rp1, rp2, num_rp2, id1, id2, fp);
							}
						}

					  if( temp_score != -1 )
						{
							if( min_score > temp_score )
							{
								min_score = temp_score;
								max_id = st[i].id;
								*rp1_id = *id1;
								*rp2_id = *id2;
							}
							else if( min_score == temp_score )
							{
								if( (*d_rate) <= min_rate ) 
								{
									min_rate = (*d_rate);
									min_score = temp_score;
									max_id = st[i].id;
									*rp1_id = *id1;
									*rp2_id = *id2;
								}
							}	
						}
					}

					if( (temp_score == -1) && (d <= 1000) && ((*sd) <= 1000) && (is_it_the_same_gene(dots, st[i].id, id, genes1, genes2, num_genes1, num_genes2) == true) )
					{
          	temp_score = abs(dots[st[i].id].identity - dots[id].identity);    
         		if( temp_score != -1 )
          	{
            	if( min_score > temp_score )
            	{
              	min_score = temp_score;
              	max_id = st[i].id;
              	*rp1_id = *id1;
              	*rp2_id = *id2;
            	}
            	else if( min_score == temp_score )
            	{
              	if( (*d_rate) <= min_rate )
              	{
                	min_rate = (*d_rate);
                	min_score = temp_score;
                	max_id = st[i].id;
                	*rp1_id = *id1;
                	*rp2_id = *id2;
              	}
            	}
						}
					}
				}
				else if( (d <= 1000) && ((*sd) <= 1000) && is_it_the_same_gene(dots, st[i].id, id, genes1, genes2, num_genes1, num_genes2) == true ) {
					temp_score = abs(dots[st[i].id].identity - dots[id].identity);

          if( temp_score != -1 )
          {
            if( min_score > temp_score )
            {
              min_score = temp_score;
              max_id = st[i].id;
              *rp1_id = *id1;
              *rp2_id = *id2;
            }
            else if( min_score == temp_score )
            {
              if( (*d_rate) <= min_rate )
              {
                min_rate = (*d_rate);
                min_score = temp_score;
                max_id = st[i].id;
                *rp1_id = *id1;
                *rp2_id = *id2;
              }
            }
          }
				}
			}
		}
	}

	free(id1);
	free(id2);
	free(d_rate);
	free(sd);
	free(is_x);
	if( max_id == -1 ) 
	{
		return(-1);
	}
	else 
	{
		return(max_id);
	}
}

bool is_it_the_same_gene(struct DotList *dots, int id1, int id2, struct g_list *genes1, struct g_list *genes2, int num_genes1, int num_genes2) 
{
	bool is_sp1 = false, is_sp2 = false;
	int i = 0;
	struct I cur_gene;
	struct I cur_reg;
	int low = -1, hi = 0;

	cur_gene = assign_I(-1,0);
	cur_reg = assign_I(-1,0);
	
	if( dots[id1].x.lower < dots[id2].x.lower ) low = dots[id1].x.lower;
	else low = dots[id2].x.lower;

	if( dots[id1].x.upper < dots[id2].x.upper ) hi = dots[id2].x.upper;
	else hi = dots[id1].x.upper;

	cur_reg = assign_I(low, hi);
	while( (i < num_genes1) && (is_sp1 == false) ) {
		cur_gene = assign_I(genes1[i].txStart, genes1[i].txEnd);
		if( ((strcmp(genes1[i].sp_name, "") == 0) || (strcmp(genes1[i].sp_name, dots[id1].name1) == 0) ) && (subset(cur_gene, dots[id1].x) == false) && (proper_overlap(dots[id1].x, cur_gene) == true)	&& (subset(cur_gene, dots[id2].x) == false) && (proper_overlap(dots[id2].x, cur_gene) == true ) && (subset(cur_gene, cur_reg) == true) )
		{
			is_sp1 = true;
		}
		i++;
	}

	if( dots[id1].y.lower < dots[id2].y.lower ) low = dots[id1].y.lower;
	else low = dots[id2].y.lower;

	if( dots[id1].y.upper < dots[id2].y.upper ) hi = dots[id2].y.upper;
	else hi = dots[id1].y.upper;

	cur_reg = assign_I(low, hi);
	i = 0;
	while( (i < num_genes2) && (is_sp2 == false) && (cur_reg.lower != -1) && (cur_gene.lower != -1)) {
		cur_gene = assign_I(genes2[i].txStart, genes2[i].txEnd);
		if( ((strcmp(genes2[i].sp_name, "") == 0) || (strcmp(genes2[i].sp_name, dots[id1].name2) == 0)) && (subset(cur_gene, dots[id1].y) == false) && (proper_overlap(dots[id1].y, cur_gene) == true)	&& (subset(cur_gene, dots[id2].y) == false) && (proper_overlap(dots[id2].y, cur_gene) == true ) && (subset(cur_gene, cur_reg) == true) )
		{
			is_sp2 = true;
		}
		i++;
	}

	if( (is_sp1 == true) && (is_sp2 == true) ) {
		return(true);
	}
	else return(false);
}
