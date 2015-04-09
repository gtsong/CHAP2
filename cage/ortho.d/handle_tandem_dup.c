#include "main.h"
#include "regions.h"
#include "handle_tandem_dup.h"
#include "find_dup_copy.h"
#include "find_merging.h"
#include "deal_gaps.h"
#include "util_gen.h"
#include "check_copy.h"
#include "util.h"
#include "util_i.h"

void handle_tandem_dup(struct DotList *dots, int *num, struct DotList *init_dots)
{
	struct slist *sorted;
	int i = 0;
	int cur_id = 0;
	struct DotList *self;
	int count = 0;
	int j = 0;
	int temp = 0;
	int num_lines;
	int *t_list; // a list of tandem dups
	int *cur_tlist;
	int *val1, *val2, *val_org;
	int num_tandem = 0;

	for( i = 0; i < *num; i++ )
	{
		if( dots[i].pair_self == SELF )
			count++;
	}

	if( count > 0 ) {
		self = (struct DotList *) ckalloc(count * (sizeof(struct DotList)));
		sorted = (struct slist *) ckalloc(count * (sizeof(struct slist)));
		t_list = (int *) ckalloc(count * (sizeof(int)));
		cur_tlist = (int *) ckalloc(count * (sizeof(int)));
		val1 = (int *) ckalloc(count * (sizeof(int)));
		val2 = (int *) ckalloc(count * (sizeof(int)));
		val_org = (int *) ckalloc(count * (sizeof(int)));

		initialize_algns(self, count);

		j = 0;
		for( i = 0; i < *num; i++ )
		{
			if( dots[i].pair_self == SELF ) 
			{
				assign_algn(self, j, dots[i]);	
				self[j].c_id = i;
				j++;
			}
		}
	
		count = j;
		num_lines = *num;
		
		for( i = 0; i < count; i++ ) {
			sorted[i].id = i;
			t_list[i] = 0;
			cur_tlist[i] = 0;
			val1[i] = -1;
			val2[i] = -1;
			val_org[i] = -1;
		}
		sort_by_width(sorted, self, count);
	
		for( i = 0; i < count; i++ )
		{
			cur_id = sorted[i].id;
			if( (self[cur_id].sign == 2) || (self[cur_id].pair_self == PAIR) ) {}
			else if( proper_overlap( self[cur_id].x, self[cur_id].y ) == true )
			{
				num_tandem = 0;
				num_tandem = find_tandem_list(self, sorted, i, count, t_list);
				if( num_tandem > 0 ) 
				{
					for( j = 0; j < num_tandem; j++ )
					{
						temp = t_list[j];
						cur_tlist[j] = self[temp].c_id;
					}
					
					conv_td_reg(dots, num_lines, self[cur_id].c_id, cur_tlist, num_tandem, init_dots, FIRST_RUN, val1, val2, val_org);
					conv_td_reg(self, count, cur_id, t_list, num_tandem, init_dots, SECOND_RUN, val1, val2, val_org);
				}
			}
		}

		free(val_org);
		free(val1);
		free(val2);
		free(cur_tlist);
		free(t_list);
		free(sorted);
		free(self);
	}
}

int find_tandem_list(struct DotList *dots, struct slist *st, int id, int num, int *t_list)
{
	int i = 0;
	int cur_id, temp = 0;
	int j = 0;
	
	cur_id = st[id].id;

	for( i = (id + 1); i < num; i++ )
	{
		temp = st[i].id;
		if( (abs(dots[temp].identity - dots[cur_id].identity) <= TD_PID ) && (dots[temp].ctg_id1 == dots[cur_id].ctg_id1) && ( strict_subset( dots[temp].x, dots[cur_id].x ) == true ) && (dots[temp].ctg_id2 == dots[cur_id].ctg_id2) && ( strict_subset( dots[temp].y, dots[cur_id].y ) == true ) )
		{
			t_list[j] = temp;
			j++;
		}
	}

	return(j);
}

void convert_tandem_region(struct DotList *dots, int num, int id, int *t_list, int num_tandem) 
{
	int i;
	int cur_id, cmp_id;
	struct DotList t1, t2;
	int len_x, len_y;
	int cur_len = 0;
	int val_t1, val_t2;
	
	for( i = 0; i < num_tandem; i++ )
	{
		val_t1 = -1;
		val_t2 = -1;
		t1.x = assign_I(-1, 0);
		t2.x = assign_I(-1, 0);
		t1.y = assign_I(-1, 0);
		t2.y = assign_I(-1, 0);
		cmp_id = t_list[i];
		if( i == 0 ) cur_id = id;
		else cur_id = t_list[i-1];

		if( ( strict_almost_equal( dots[cmp_id].x, dots[cur_id].x ) == true ) || ( strict_almost_equal( dots[cmp_id].y, dots[cur_id].y) == true ) ) {}
		else if( ( strict_subset( dots[cmp_id].x, dots[cur_id].x ) == true ) && ( strict_subset( dots[cmp_id].y, dots[cur_id].y ) == true ) )
		{
			if( abs(dots[cur_id].x.upper - dots[cmp_id].x.upper) > abs(dots[cur_id].x.lower - dots[cmp_id].x.lower)	)
			{
				if( ( dots[cur_id].x.upper - dots[cmp_id].x.upper ) <= 0 ) t1.x = assign_I(-1, 0);
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t1.x = assign_I(dots[cmp_id].x.upper, dots[cur_id].x.upper);
					cur_len = (int)(((float)(width(t1.x)) * ((float)len_y)/(float)len_x));
					if( cur_len < DEL_TH ) {
						t1.x = assign_I(-1, 0);
						t1.y = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 0 ) {
						if( dots[cur_id].y.upper > (dots[cur_id].y.upper - cur_len)) t1.y = assign_I(dots[cur_id].y.upper - cur_len, dots[cur_id].y.upper);
						else t1.x = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 1 ) {
						if( (dots[cur_id].y.lower + cur_len) > dots[cur_id].y.lower ) t1.y = assign_I(dots[cur_id].y.lower, dots[cur_id].y.lower + cur_len);
						else t1.x = assign_I(-1, 0);
					}
				}
			}
			else
			{
				if( ( dots[cur_id].x.lower - dots[cmp_id].x.lower ) >= 0 ) t1.x = assign_I(-1, 0);
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t1.x = assign_I(dots[cur_id].x.lower, dots[cmp_id].x.lower);
					cur_len = (int)(((float)(width(t1.x)) * ((float)len_y)/(float)len_x));
					if( cur_len < DEL_TH ) {
						t1.x = assign_I(-1, 0);
						t1.y = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 0 ) {
						if( (dots[cur_id].y.lower + cur_len) > dots[cur_id].y.lower ) t1.y = assign_I(dots[cur_id].y.lower, dots[cur_id].y.lower + cur_len);
						else t1.x = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 1 ) {
						if( dots[cur_id].y.upper > dots[cur_id].y.upper ) t1.y = assign_I(dots[cur_id].y.upper - cur_len, dots[cur_id].y.upper);
						else t1.x = assign_I(-1, 0);
					}
				}
			}

			if( abs(dots[cmp_id].y.lower - dots[cur_id].y.lower) > abs(dots[cur_id].y.upper - dots[cmp_id].y.upper)	)
			{
				if( ( dots[cmp_id].y.lower - dots[cur_id].y.lower ) <= 0 ) t2.x = assign_I(-1, 0); 
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t2.y = assign_I(dots[cur_id].y.lower, dots[cmp_id].y.lower);
					cur_len = (int)(((float)(width(t2.y)) * ((float)len_x)/(float)len_y));
					if( cur_len < DEL_TH ) {
						t2.x = assign_I(-1, 0);
						t2.y = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 0 ) {
						if( (dots[cur_id].x.lower + cur_len) > dots[cur_id].x.lower ) t2.x = assign_I(dots[cur_id].x.lower, dots[cur_id].x.lower + cur_len);
						else t2.x = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 1 ) {
						if( dots[cur_id].x.upper > (dots[cur_id].x.upper - cur_len) ) t2.x = assign_I(dots[cur_id].x.upper - cur_len, dots[cur_id].x.upper);
						else t2.x = assign_I(-1, 0);
					}
				}
			}
			else
			{
				if( ( dots[cur_id].y.upper - dots[cmp_id].y.upper ) <= 0 ) t2.x = assign_I(-1, 0);
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t2.y = assign_I(dots[cmp_id].y.upper, dots[cur_id].y.upper);
					cur_len = (int)(((float)(width(t2.y)) * ((float)len_x)/(float)len_y));
					if( cur_len < DEL_TH ) {
						t2.x = assign_I(-1, 0);
						t2.y = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 0 ) {
						if( dots[cur_id].x.upper > (dots[cur_id].x.upper - cur_len) ) t2.x = assign_I(dots[cur_id].x.upper - cur_len, dots[cur_id].x.upper);
						else t2.x = assign_I(-1, 0);
					}
					else if( dots[cur_id].sign == 1 ) {
						if( (dots[cur_id].x.lower + cur_len) > dots[cur_id].x.lower ) t2.x = assign_I(dots[cur_id].x.lower, dots[cur_id].x.lower + cur_len);
						else t2.x = assign_I(-1, 0);
					}
				}
			}
		}

		if( (t1.x.lower >= 0) && (t1.y.lower >= 0)) val_t1 = check_tandem_reg( t1, dots, num );
		else val_t1 = -1;

		if( (t2.x.lower >= 0) && (t2.y.lower >= 0) ) val_t2 = check_tandem_reg( t2, dots, num );
		else val_t2 = -1;

		if( (val_t1 != -1) && (val_t2 != -1) ) 
		{
			if( val_t1 <= val_t2 )
			{
				dots[cur_id].x = assign_I(t1.x.lower, t1.x.upper);
				dots[cur_id].y = assign_I(t1.y.lower, t1.y.upper);
				dots[cur_id].rp1_id = 0;
			}
			else 
			{
				dots[cur_id].x = assign_I(t2.x.lower, t2.x.upper);
				dots[cur_id].y = assign_I(t2.y.lower, t2.y.upper);
				dots[cur_id].rp1_id = 0;
			}
		}
		else if( val_t1 != -1 )
		{
			dots[cur_id].x = assign_I(t1.x.lower, t1.x.upper);
			dots[cur_id].y = assign_I(t1.y.lower, t1.y.upper);
			dots[cur_id].rp1_id = 0;
		}
		else if( val_t2 != -1 )
		{
			dots[cur_id].x = assign_I(t2.x.lower, t2.x.upper);
			dots[cur_id].y = assign_I(t2.y.lower, t2.y.upper);
			dots[cur_id].rp1_id = 0;
		}
	}
}

void adjust_init_offset(struct DotList *init_algns, int init_id, struct DotList t1, struct DotList *algns, int cur_id)
{

	if( ((init_algns[init_id].x.upper + t1.x.upper - algns[cur_id].x.upper) > (init_algns[init_id].x.lower + t1.x.lower - algns[cur_id].x.lower)) && ((init_algns[init_id].y.lower + t1.y.lower - algns[cur_id].y.lower) < (init_algns[init_id].y.upper + t1.y.upper - algns[cur_id].y.upper)) ) 
	{
		init_algns[init_id].xl_offset = init_algns[init_id].xl_offset + t1.x.lower - algns[cur_id].x.lower;
		init_algns[init_id].xr_offset = init_algns[init_id].xr_offset + t1.x.upper - algns[cur_id].x.upper;
		init_algns[init_id].yl_offset = init_algns[init_id].yl_offset + t1.y.lower - algns[cur_id].y.lower;
		init_algns[init_id].yr_offset = init_algns[init_id].yr_offset + t1.y.upper - algns[cur_id].y.upper;
		init_algns[init_id].x = assign_I(init_algns[init_id].x.lower + t1.x.lower - algns[cur_id].x.lower, init_algns[init_id].x.upper + t1.x.upper - algns[cur_id].x.upper);
		init_algns[init_id].y = assign_I(init_algns[init_id].y.lower + t1.y.lower - algns[cur_id].y.lower, init_algns[init_id].y.upper + t1.y.upper - algns[cur_id].y.upper);
		init_algns[init_id].rp1_id = 0;
	}
}

void conv_td_reg(struct DotList *dots, int num, int id, int *t_list, int num_tandem, struct DotList *init_dots, int flag, int *val1, int *val2, int *val_org)
{
	int i;
	int cur_id, cmp_id;
	struct DotList t1, t2;
	struct DotList *cur_t;
	int len_x, len_y;
	int cur_len = 0;
	int val_t1, val_t2, val_org_reg;
	int init_id;

	cur_t = (struct DotList *) ckalloc(sizeof(struct DotList));
	
	for( i = 0; i < num_tandem; i++ )
	{

		if( flag == FIRST_RUN ) {
			val_org_reg = -1;
			val_t1 = -1;
			val_t2 = -1;
		}
		else {
			val_org_reg = val_org[i];
			val_t1 = val1[i];
			val_t2 = val2[i];
		}

		t1.x = assign_I(-1, 0);
		t2.x = assign_I(-1, 0);
		t1.y = assign_I(-1, 0);
		t2.y = assign_I(-1, 0);
		cmp_id = t_list[i];
		if( i == 0 ) cur_id = id;
		else cur_id = t_list[i-1];

		if( dots[cmp_id].ctg_id1 != dots[cur_id].ctg_id1 ) {
			fatalf("error: handling alignments from different contigs %s vs %s in handling_tandem_duplications.c\n", dots[cmp_id].name1, dots[cur_id].name1);
		}
		
		if( dots[cmp_id].ctg_id2 != dots[cur_id].ctg_id2 ) {
			fatalf("error: handling alignments from different contigs %s vs %s in handling_tandem_duplications.c\n", dots[cmp_id].name2, dots[cur_id].name2);
		}

		if( ( strict_almost_equal( dots[cmp_id].x, dots[cur_id].x ) == true ) || ( strict_almost_equal( dots[cmp_id].y, dots[cur_id].y) == true ) ) {}
		else if( ( strict_subset( dots[cmp_id].x, dots[cur_id].x ) == true ) && ( strict_subset( dots[cmp_id].y, dots[cur_id].y ) == true ) )
		{
			if( abs(dots[cur_id].x.upper - dots[cmp_id].x.upper) > abs(dots[cur_id].x.lower - dots[cmp_id].x.lower)	)
			{
				if( ( dots[cur_id].x.upper - dots[cmp_id].x.upper ) <= 0 ) t1.x = assign_I(-1, 0);
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t1.x = assign_I(dots[cmp_id].x.upper, dots[cur_id].x.upper);
					cur_len = (int)(((float)(width(t1.x)) * ((float)len_y)/(float)len_x));
					t1.y = assign_I(dots[cur_id].x.upper, dots[cur_id].x.upper + cur_len);
				}
			}
			else
			{
				if( ( dots[cur_id].x.lower - dots[cmp_id].x.lower ) >= 0 ) t1.x = assign_I(-1, 0);
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t1.x = assign_I(dots[cur_id].x.lower, dots[cmp_id].x.lower);
					cur_len = (int)(((float)(width(t1.x)) * ((float)len_y)/(float)len_x));
					t1.y = assign_I(dots[cmp_id].x.lower, dots[cmp_id].x.lower + cur_len);
				}
			}

			if( abs(dots[cmp_id].y.lower - dots[cur_id].y.lower) > abs(dots[cur_id].y.upper - dots[cmp_id].y.upper)	)
			{
				if( ( dots[cmp_id].y.lower - dots[cur_id].y.lower ) <= 0 ) t2.x = assign_I(-1, 0); 
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t2.y = assign_I(dots[cur_id].y.lower, dots[cmp_id].y.lower);
					cur_len = (int)(((float)(width(t2.y)) * ((float)len_x)/(float)len_y));
					t2.x = assign_I(dots[cur_id].y.lower - cur_len, dots[cur_id].y.lower);
				}
			}
			else
			{
				if( ( dots[cur_id].y.upper - dots[cmp_id].y.upper ) <= 0 ) t2.x = assign_I(-1, 0);
				else
				{
					len_x = width(dots[cur_id].x);
					len_y = width(dots[cur_id].y);

					t2.y = assign_I(dots[cmp_id].y.upper, dots[cur_id].y.upper);
					cur_len = (int)(((float)(width(t2.y)) * ((float)len_x)/(float)len_y));
					t2.x = assign_I(dots[cmp_id].y.upper - cur_len, dots[cmp_id].y.upper);
				}
			}
		}

		val_org_reg = -1;
		if( !proper_overlap(dots[cur_id].x, dots[cur_id].y) ) {
			val_org_reg = STRICT;
			val_org_reg = check_tandem_reg( dots[cur_id], dots, num );
		}

		if( flag == FIRST_RUN ) {
			if( (t1.x.lower >= 0) && (t1.y.lower >= 0) ) {
				val_t1 = check_tandem_reg( t1, dots, num );
			}
			else val_t1 = -1;

			if( (t2.x.lower >= 0) && (t2.y.lower >= 0)) {
				val_t2 = check_tandem_reg( t2, dots, num );
			}
			else val_t2 = -1;

			if( (val_t1 == -1) && (val_t2 == -1) ) {
				if( t1.x.lower >= 0 ) val_t1 = LOOSE;
				else if( t2.x.lower >= 0 ) val_t2 = LOOSE;
			}

			val_org[i] = val_org_reg;
			val1[i] = val_t1;
			val2[i] = val_t2;
		}

		if( val_org_reg != -1 ) {}
		else if( (val_t1 != -1) && (val_t2 != -1) && (t1.x.lower >= 0) && (t1.y.lower >= 0) && (t2.x.lower >= 0) && (t2.y.lower >= 0)) 
		{
			if( val_t1 <= val_t2 )
			{
				init_id = dots[cur_id].index;
				if( (flag == FIRST_RUN) && (init_dots[init_id].c_id == -1) && (init_dots[init_id].m_id == -1) ) {
// in order to get the original boundaries, offsets defined here should be just substrated.
					adjust_init_offset(init_dots, init_id, t1, dots, cur_id);
				}

				dots[cur_id].x = assign_I(t1.x.lower, t1.x.upper);
				dots[cur_id].y = assign_I(t1.y.lower, t1.y.upper);
				dots[cur_id].rp1_id = 0;
			}
			else 
			{
				init_id = dots[cur_id].index;
				if( (flag == FIRST_RUN) && (init_dots[init_id].c_id == -1) && (init_dots[init_id].m_id == -1) ) {
					adjust_init_offset(init_dots, init_id, t2, dots, cur_id);
				}
				dots[cur_id].x = assign_I(t2.x.lower, t2.x.upper);
				dots[cur_id].y = assign_I(t2.y.lower, t2.y.upper);
				dots[cur_id].rp1_id = 0;
				init_id = dots[cur_id].index;
			}
		}
		else if( (val_t1 != -1) && (t1.x.lower >= 0) && (t1.y.lower >= 0))
		{
			init_id = dots[cur_id].index;
			if( (flag == FIRST_RUN) && (init_dots[init_id].c_id == -1) && (init_dots[init_id].m_id == -1) ) {
// in order to reflect the change of the boundaries, offsets defined here should be just added.
				adjust_init_offset(init_dots, init_id, t1, dots, cur_id);
			}

			dots[cur_id].x = assign_I(t1.x.lower, t1.x.upper);
			dots[cur_id].y = assign_I(t1.y.lower, t1.y.upper);
			dots[cur_id].rp1_id = 0;
			init_id = dots[cur_id].index;
		}
		else if( (val_t2 != -1) && (t2.x.lower >= 0) && (t2.y.lower >= 0))
		{
			init_id = dots[cur_id].index;
			if( (flag == FIRST_RUN) && (init_dots[init_id].c_id == -1) && (init_dots[init_id].m_id == -1) ) {
				adjust_init_offset(init_dots, init_id, t2, dots, cur_id);
			}

			dots[cur_id].x = assign_I(t2.x.lower, t2.x.upper);
			dots[cur_id].y = assign_I(t2.y.lower, t2.y.upper);
			dots[cur_id].rp1_id = 0;
			init_id = dots[cur_id].index;
		}
	}

	val_org_reg = -1;
	cmp_id = t_list[num_tandem-1];
	len_x = width(dots[cmp_id].x);
	len_y = width(dots[cmp_id].y);
	if( proper_overlap(dots[cmp_id].x, dots[cmp_id].y) ) {
		t1.x = assign_I(dots[cmp_id].x.lower, dots[cmp_id].x.lower + (dots[cmp_id].y.upper - dots[cmp_id].x.lower)/2);
		t1.y = assign_I(dots[cmp_id].x.lower, dots[cmp_id].x.lower + (dots[cmp_id].y.upper - dots[cmp_id].x.lower)/2);
	}
	else {
		val_org_reg = STRICT;
		t1.x = assign_I(dots[cmp_id].x.lower, dots[cmp_id].x.upper);
		t1.y = assign_I(dots[cmp_id].y.lower, dots[cmp_id].y.upper);
	}

	cur_len = (int)(((float)(width(t1.x)) * ((float)len_y)/(float)len_x));
	t1.y = assign_I(t1.x.upper, t1.x.upper + cur_len);
	if( t2.y.lower != -1 ) {
		cur_len = (int)(((float)(width(t2.y)) * ((float)len_x)/(float)len_y));
		t2.x = assign_I(t2.y.lower - cur_len, t2.y.lower);
	}
	else t2.x = assign_I(-1,0);

	if( flag == FIRST_RUN ) {
		if( val_org_reg != -1 ) val_org_reg = check_tandem_reg(dots[cmp_id], dots, num);
		if( (t1.x.lower >= 0) && (t1.y.lower >= 0) ) val_t1 = check_tandem_reg(t1, dots, num);
		else val_t1 = -1;

		if( (t2.x.lower < 0) || (t2.y.lower < 0) ) val_t2 = -1;
		else val_t2 = check_tandem_reg(t2, dots, num);
		val_org[num_tandem] = val_org_reg;
		val1[num_tandem] = val_t1;
		val2[num_tandem] = val_t2;
	}
	else {
		val_org_reg = val_org[num_tandem];
		val_t1 = val1[num_tandem];
		val_t2 = val2[num_tandem];
	}

	if( (t1.x.lower < 0) && (t1.y.lower < 0) ) val_t1 = -1;
	if( (t2.x.lower < 0) && (t2.y.lower < 0) ) val_t2 = -1;

	if( val_org_reg != -1 ) {}
	else if( (val_t1 != -1) && (val_t2 != -1) ) {
		if( val_t1 < val_t2 ) {
			assign_algn(cur_t, 0, t1);
		}
		else assign_algn(cur_t, 0, t2);
	}		
	else if( val_t1 != -1 ) assign_algn(cur_t, 0, t1);
	else if( val_t2 != -1 ) assign_algn(cur_t, 0, t2);

	if( val_org_reg != -1 ) {}
	else if( (val_t1 != -1) || (val_t2 != -1) ) {
		init_id = dots[cmp_id].index;
		if( (flag == FIRST_RUN) && (init_dots[init_id].c_id == -1) && (init_dots[init_id].m_id == -1) ) {
// in order to reflect the change of the boundaries, offsets defined here should be just added.
			adjust_init_offset(init_dots, init_id, *cur_t, dots, cmp_id);
		}

		dots[cmp_id].x = assign_I((*cur_t).x.lower, (*cur_t).x.upper);
		dots[cmp_id].y = assign_I((*cur_t).y.lower, (*cur_t).y.upper);
		dots[cmp_id].rp1_id = 0;
	}

	free(cur_t);
}

int check_tandem_reg(struct DotList t, struct DotList *dots, int num )
{
	int i = 0;
	bool is_end = false;
	bool is_x = true, is_y = true;
	int t_val = STRICT;

	if( (t.x.upper < t.x.lower) || (t.y.upper < t.y.lower) ) {
		fatalf("empty interval %d-%d, %d-%d\n", t.x.lower, t.x.upper, t.y.lower, t.y.upper);
	}

	while((i < num) && (is_end == false))
	{
		if( (dots[i].sign == 2) || (dots[i].pair_self == SELF) ) {}
		else 
		{
			if( ((f_loose_subset(dots[i].x, t.x, t_val) == false) && (f_loose_subset(t.x, dots[i].x, t_val) == false) && (f_loose_overlap(dots[i].x, t.x, t_val) == true)) || ((f_loose_subset(dots[i].y, t.x, t_val) == false) && (f_loose_subset(t.x, dots[i].y, t_val) == false) && (f_loose_overlap(dots[i].y, t.x, t_val) == true))) is_x = false;

			if( ((f_loose_subset(dots[i].x, t.y, t_val) == false) && (f_loose_subset(t.y, dots[i].x, t_val) == false) && (f_loose_overlap(dots[i].x, t.y, t_val) == true)) || ((f_loose_subset(dots[i].y, t.y, t_val) == false) && (f_loose_subset(t.y, dots[i].y, t_val) == false) && (f_loose_overlap(dots[i].y, t.y, t_val) == true))) is_y = false;
		}

		if( (is_x == false) && (is_y == false) )
		{
			if( t_val == LOOSE )
			{
				is_end = true;
			}
			else
			{
				i = 0;
				is_x = true;
				is_y = true;
				t_val = LOOSE;
			}
		}
		else if( i < (num - 1) ) i++;
		else
		{
			is_end = true;
		}
	}

	if( (is_x == true ) || (is_y == true) ) {}
	else t_val = -1;

	return( t_val );
}
