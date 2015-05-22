#include "regions.h"
#include "find_merging.h"
#include "deal_gaps.h"
#include "util.h"
#include "util_i.h"

int distance(struct DotList *dots, int loc_id, int comp_id, bool *is_x, int *sd)
{
	int x = 0, y = 1;
	int t_x1 = 0, t_x2 = 0, t_y1 = 1, t_y2 = 1;
	int res = 0;
	int left_id = 0, right_id = 0;

	if( dots[loc_id].x.lower < dots[comp_id].x.lower )
	{
		left_id = loc_id;
		right_id = comp_id;
	}
	else if( dots[loc_id].x.lower > dots[comp_id].x.lower )
	{
		left_id = comp_id;
		right_id = loc_id;
	}
	else
	{
		if( dots[loc_id].x.upper <= dots[comp_id].x.upper )
		{
			left_id = loc_id;
			right_id = comp_id;
		}
		else
		{
			left_id = comp_id;
			right_id = loc_id;
		}
	}

	t_x1 = dots[left_id].x.upper;
	t_x2 = dots[right_id].x.lower;
	if( dots[left_id].sign == 0 )
	{
		t_y1 = dots[left_id].y.upper;
		t_y2 = dots[right_id].y.lower;
	}	
	else if( dots[left_id].sign == 1 )
	{
		t_y1 = dots[left_id].y.lower;
		t_y2 = dots[right_id].y.upper;
	}
	
 	x = abs(t_x1 - t_x2);
	y = abs(t_y1 - t_y2);

	if( x >= y ) 
	{
		*is_x = true;
		res = x;
		*sd = y;
	}	
	else 
	{
		*is_x = false;
		res = y;
		*sd = x;
	}
	return(res);
}

int m_val(struct slist *st, struct DotList *dots, int i)
{
	if(st[i].is_x) return(dots[st[i].id].x.lower);
	else return(dots[st[i].id].y.lower);
}

bool check_candi(struct DotList *dots, int id, int comp_id, int flag)
{
	int th_val;

	if( (flag == CHECK_INS_DUP) || (flag == CHECK_GENE_LOSS) ) th_val = 3*RP_BD;
	else if( flag == CHECK_DEL ) th_val = S_RP_BD;

	if( proper_overlap(dots[id].x, dots[comp_id].x) == true)
	{
		if( (strict_almost_equal(dots[id].x, dots[comp_id].x) == true)  || (width(intersect(dots[id].x, dots[comp_id].x)) >= th_val)) return(false);
	}

	if( flag != CHECK_GENE_LOSS )
	{
		if( proper_overlap(dots[id].y, dots[comp_id].y) == true)
		{
			if( (strict_almost_equal(dots[id].y, dots[comp_id].y) == true ) || (width(intersect(dots[id].y, dots[comp_id].y)) >= th_val)) return(false);
		}
	}

	return(true);
}

bool find_merging(struct slist *st, struct DotList *dots, int num, int loc)
{
  int i;
	int d;
  int *sd;
  bool *is_x;
  struct gap_list *gps;
  int num_gaps = 0;

  is_x = (bool *) ckalloc(sizeof(bool));
  sd = (int *) ckalloc(sizeof(int));
	gps = (struct gap_list *) ckalloc(sizeof(struct gap_list) * num);
  i = loc + 1;
  while( i < num )
  {
    if( (dots[st[i].id].sign != 2) && (dots[st[loc].id].sign == dots[st[i].id].sign) && ((d = distance(dots, st[loc].id, st[i].id, is_x, sd)) <= MDIS_THRESHOLD) )
    {
      gps[num_gaps] = define_gap(dots, st[loc].id, st[i].id, d, *sd, *is_x);
      if(gps[num_gaps].type == -1) // this gap is meaningless
      {
      }
      else num_gaps++;
    }
    i++;
  }

	free(gps);
	free(sd);
  free(is_x);
  return(false);
}

void merging_step(struct DotList *dots, int loc, int comp)
{
  struct I t_x, t_y;  
	struct I x, y;  
	int from, to;
	      
	x = assign_I(dots[loc].x.lower, dots[loc].x.upper);
	y = assign_I(dots[loc].y.lower, dots[loc].y.upper);
	t_x = assign_I(dots[comp].x.lower, dots[comp].x.upper);
	t_y = assign_I(dots[comp].y.lower, dots[comp].y.upper);

  dots[loc].identity = ((dots[loc].identity * width(x)) + (dots[comp].identity * width(t_x))) / (width(x) + width(t_x));

  if( x.lower <= t_x.lower ) from = x.lower;
  else from = t_x.lower;
  if( t_x.upper >= x.upper ) to = t_x.upper;
  else to = x.upper;
													  
  dots[loc].x = assign_I(from, to);
															  
	if( y.lower <= t_y.lower ) from = y.lower;
  else from = t_y.lower;

	if( t_y.upper >= y.upper ) to = t_y.upper;
  else to = y.upper;
  dots[loc].y = assign_I(from, to);
																					    
  dots[comp].sign = DELETED;
}

bool is_tandem(struct DotList cur)
{
	bool res;

	if( (cur.x.upper <= cur.x.lower) || (cur.y.upper <= cur.y.lower) ) {
		fatalf("%d-%d, %d-%d\n", cur.x.lower, cur.y.upper, cur.y.lower, cur.y.upper);
	}

	if( proper_overlap(cur.x, cur.y) == true )
	{
		res = true;
	}
	else if( cur.rp1_id == 0 ) res = true;
	else
	{
		if( abs(cur.x.upper - cur.y.lower ) <= O_TH )
		{
			res = true;
		}
		else
		{
			res = false;
		}
	}

	return(res);
}

int compute_closeness(struct DotList *dots, int loc_id, int comp_id)
{
  int res;
  struct I self_alg;

  self_alg = assign_I(1, 100000);

  if( comp_id == SELF_ID )
  {
    res = compute_distance(dots[loc_id].x, dots[loc_id].y, self_alg, self_alg, dots[loc_id].sign);
  }
  else
  {
    res = compute_distance(dots[loc_id].x, dots[loc_id].y, dots[comp_id].x, dots[comp_id].y, dots[loc_id].sign);
  }

  return(res);
}
