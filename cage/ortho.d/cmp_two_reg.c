#include "main.h"
#include "regions.h"
#include "cmp_two_reg.h"
#include "pred_regions.h"
#include "util_gen.h"
#include "util_i.h"

int det_dup_reg_in_self(int num_list, struct DotList *self, struct I check_x, struct I check_y, int sign)
{
	int i = 0, j = 0;
	struct I temp_1, temp_2;
	struct I cur_1, cur_2, cmp_1, cmp_2;
	bool is_end = false;
	int count_left = 0, count_right = 0;
	int left_pid = 0, right_pid = 0;
	int cut_len_left = 0, cut_len_right = 0;
	int res = TIE;
	int distance = 0;
	bool is_candi = false;
	bool is_assigned = false;

	temp_1 = assign_I(0, 1);
	temp_2 = assign_I(0, 1);

	for( i = 0; i < num_list; i++ )
	{
		cur_1 = assign_I(0, 1);
		cur_2 = assign_I(0, 1);
		cmp_1 = assign_I(0, 1);
		cmp_2 = assign_I(0, 1);
		is_candi = false;
		left_pid = self[i].identity;
		right_pid = self[i].identity;

		if( ( is_assigned == true ) && ((width(temp_1) <= MIN_INTERVAL) || (width(temp_2) <= MIN_INTERVAL)) ) {
			is_end = true;
		}
		else if( (width(self[i].x) <= MIN_INTERVAL) || (width(self[i].y) <= MIN_INTERVAL) )
		{
			is_candi = false;
		}
		else if( (almost_subset(self[i].x, check_x) == true) || ((f_loose_overlap(self[i].x, check_x, SECOND_RUN) == true) && (width(intersect(self[i].x, check_x)) >= MIN_LEN) ))
		{
			temp_1 = assign_I(self[i].x.lower, self[i].x.upper);
			temp_2 = assign_I(self[i].y.lower, self[i].y.upper);
			is_assigned = true;
			is_candi = true;
		}
		else if( (almost_subset(self[i].y, check_x) == true) || (f_loose_overlap(self[i].y, check_x, SECOND_RUN) == true) )
		{
			temp_1 = assign_I(self[i].y.lower, self[i].y.upper);
			temp_2 = assign_I(self[i].x.lower, self[i].x.upper);
			is_assigned = true;
			is_candi = true;
		}
		else {
			is_candi = false;
		}

		if( (is_end == false) && (is_candi == true) )
		{
			if( (temp_1.lower < check_x.lower) && (temp_1.upper > check_x.upper) )
			{
				if( width(check_x) > 0 ) {
					cur_1 = assign_I(0, width(check_x));

					cut_len_left = check_x.lower - temp_1.lower;
					cut_len_right = temp_1.upper - check_x.upper; 
					if( (temp_2.upper - cut_len_right) <= (temp_2.lower + cut_len_left) )
					{
						is_end = true;
					}
					else
					{
						cur_2 = assign_I(temp_2.lower + cut_len_left, temp_2.upper - cut_len_right);
					}
				}
				else {
					is_end = true;
				}
			}
			else
			{
				if( temp_1.lower < check_x.lower )
				{
					if( temp_1.upper > check_x.lower ) {
						cur_1 = assign_I( 0, temp_1.upper - check_x.lower );

						cut_len_left = check_x.lower - temp_1.lower;
						if( temp_2.upper <= (temp_2.lower + cut_len_left) )
						{
							is_end = true;
						}
						else
						{
							cur_2 = assign_I( temp_2.lower + cut_len_left, temp_2.upper );
						}
					}
					else {
						is_end = true;
					}	
				}
				else if( temp_1.upper > check_x.upper )
				{
					if( width(check_x) > (temp_1.lower - check_x.lower) ) {
						cur_1 = assign_I( temp_1.lower - check_x.lower, width(check_x) );
						cut_len_right = temp_1.upper - check_x.upper;
						if( (temp_2.upper - cut_len_right) <= temp_2.lower )
						{
							is_end = true;
						}
						else
						{
							cur_2 = assign_I( temp_2.lower, temp_2.upper - cut_len_right );
						}			
					}
					else {
						is_end = true;
					}
				}
				else
				{
					if( ( temp_1.upper - check_x.lower ) <= (temp_1.lower - check_x.lower) ) {
						is_end = true;
					}
					else {
						cur_1 = assign_I( temp_1.lower - check_x.lower, temp_1.upper - check_x.lower);
						cur_2 = assign_I( temp_2.lower, temp_2.upper );
					}
				}
			}

			if( is_end == false )
			{
				if( sign == 0 )
				{
					if( (cur_1.upper + check_y.lower) > (cur_1.lower + check_y.lower) ) {
						cmp_1 = assign_I(cur_1.lower + check_y.lower, cur_1.upper + check_y.lower);
						cmp_2 = assign_I(cur_2.lower, cur_2.upper);
					}
					else is_end = true;
				}
				else if( sign == 1 )
				{
					if( (check_y.upper - cur_1.lower) > (check_y.upper - cur_1.upper) ) {
						cmp_1 = assign_I(check_y.upper - cur_1.upper, check_y.upper - cur_1.lower);
						cmp_2 = assign_I(cur_2.lower, cur_2.upper);
					}
					else is_end = true;
				}

				if( is_end == false ) {
					if( (width(cmp_1) <= MIN_INTERVAL) || (width(cmp_2) <= MIN_INTERVAL) ) is_end = true;
				}
			}

			j = 0;
			while( (is_end == false) && (j < num_list) )
			{
				if( j == i ) {}
				else
				{
					if( (strict_almost_equal(self[j].x, cmp_1) == true) && ( strict_almost_equal(self[j].y, cmp_2) == true ) )
					{
						right_pid = self[j].identity;
						is_end = true;
						if( left_pid > right_pid ) count_left++;
						else if( left_pid < right_pid ) count_right++;
					}	
					else if( (loose_subset(cmp_1, self[j].x) == true) && (loose_subset(cmp_2, self[j].y) == true ) )
					{
						distance = compute_distance(cmp_1, cmp_2, self[j].x, self[j].y, sign);
						if( distance <= DIS_THRESHOLD )
						{
							right_pid = self[j].identity;
							is_end = true;
							if( left_pid > right_pid ) count_left++;
							else if( left_pid < right_pid ) count_right++;
						}
					}
				}
				j++;
			}
			is_end = false;
		}
	}	

	if( count_left > count_right )
	{
		res = LEFT_SIDE;
	}
	else if( count_left < count_right )
	{
		res = RIGHT_SIDE;
	}
	else
	{
		res = TIE;
	}

	return(res);
}
