#include "main.h"
#include "regions.h"
#include "adjust_plot.h"
#include "find_merging.h"
#include "util_gen.h"

void adjust_after_merging(struct gap_list *gps, int num_gaps, int pt)
{
	int i;

	for( i = 0; i < num_gaps; i++ )
	{
		if( i != pt )
		{
			if(gps[i].id1 == gps[pt].id2) 
			{
				gps[i].id1 = gps[pt].id1;	
			}
			else if(gps[i].id2 == gps[pt].id2)
			{
				gps[i].id2 = gps[pt].id1;	
			}

			if( gps[i].id1 == gps[i].id2 )
				gps[i].type = -1;
		}
	}
	gps[pt].type = -1;
}

void sort_by_yintercept(struct slist *st, struct DotList *dots, int num)
{
	int i, j;
	int cur_val;

	for( i = 0; i < num; i++ )
	{
		st[i].id = i;
		if( dots[i].sign == 0 ) 
		{
			st[i].val = dots[i].y.lower - dots[i].x.lower;
		}
		else if( dots[i].sign == 1 )
		{
			st[i].val = dots[i].x.upper + dots[i].y.lower;
		}
	}
	quick_sort_inc(st, 0, num-1);

	i = 0;
	while(i < num)
	{
		j = 0;
		cur_val = st[i].val;
		while(((i+j) < num) && (cur_val == st[i+j].val))
		{
			if( dots[st[i+j].id].sign == 0 ) st[i+j].val = dots[st[i+j].id].x.lower;
			else if( dots[st[i+j].id].sign == 1 ) st[i+j].val = dots[st[i+j].id].x.lower;
			j++;
		}
		quick_sort_inc(st, i, i+j-1);
		i = i+j;
	}
}

void sort_list(struct slist *st, struct DotList *dots, int num)
{
	int i, j;
	int cur_val;

	for( i = 0; i < num; i++ )
	{
		st[i].id = i;
		st[i].val = dots[i].x.lower;
	}
	quick_sort_dec(st, 0, num-1);

	i = 0;
	while(i < num)
	{
		j = 0;
		cur_val = st[i].val;
		while(((i+j) < num) && (cur_val == st[i+j].val))
		{
			if( dots[st[i+j].id].sign == 0 ) st[i+j].val = dots[st[i+j].id].y.lower;
			else if( dots[st[i+j].id].sign == 1 ) st[i+j].val = dots[st[i+j].id].y.upper;
			j++;
		}
		quick_sort_inc(st, i, i+j-1);
		i = i+j;
	}
}

void quick_sort_dec(struct slist *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct slist h;
	int x=a[(lo+hi)/2].val;

//  partition
	do
	{    
		while (a[i].val>x) i++; 
		while (a[j].val<x) j--;
		if (i<=j)
		{
			h = assign_slist(a[i]);
			a[i] = assign_slist(a[j]);
			a[j] = assign_slist(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec(a, lo, j);
	if (i < hi) quick_sort_dec(a, i, hi);
}

void quick_sort_inc(struct slist *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct slist h;
	int x=a[(lo+hi)/2].val;

//  partition
	do
	{    
		while (a[i].val<x) i++; 
		while (a[j].val>x) j--;
		if (i<=j)
		{
			h = assign_slist(a[i]);
			a[i] = assign_slist(a[j]);
			a[j] = assign_slist(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc(a, lo, j);
	if (i < hi) quick_sort_inc(a, i, hi);
}

struct slist assign_slist(struct slist a)
{
	struct slist res;

	res.id = a.id;
	res.val = a.val;
	res.is_x = a.is_x;

	return(res);
}

void quick_sort_inc_int(struct I *a, int lo, int hi)
{	
	int i=lo, j=hi;
	struct I h;
	int x=a[(lo+hi)/2].lower;
	do
	{    
		while (a[i].lower<x) i++; 
		while (a[j].lower>x) j--;
		if (i<=j)
		{
			h = assign_I(a[i].lower, a[i].upper);
			a[i] = assign_I(a[j].lower, a[j].upper);
			a[j] = assign_I(h.lower, h.upper);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_int(a, lo, j);
	if (i < hi) quick_sort_inc_int(a, i, hi);
}

struct perm_pt assign_pm_val(struct perm_pt a)
{
	struct perm_pt res;

	res.id = a.id;
	res.x_pt = a.x_pt;
	res.y_pt = a.y_pt;
	res.sign = a.sign;

	return(res);
}

void quick_sort_plist_x(struct perm_pt *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct perm_pt h;
	int x;
	
	x = a[(lo+hi)/2].x_pt;

//  partition
	do
	{    
		while (a[i].x_pt<x) 
		{
			i++; 
		}
		while (a[j].x_pt>x) 
		{
			j--;
		}

		if (i<=j)
		{
			h = assign_pm_val(a[i]);
			a[i] = assign_pm_val(a[j]);
			a[j] = assign_pm_val(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_plist_x(a, lo, j);
	if (i < hi) quick_sort_plist_x(a, i, hi);
}

void quick_sort_plist_y(struct perm_pt *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct perm_pt h;
	int x;
	
	x = a[(lo+hi)/2].y_pt;

//  partition
	do
	{    
		while (a[i].y_pt<x) 
		{
			i++; 
		}
		while (a[j].y_pt>x) 
		{
			j--;
		}

		if (i<=j)
		{
			h = assign_pm_val(a[i]);
			a[i] = assign_pm_val(a[j]);
			a[j] = assign_pm_val(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_plist_y(a, lo, j);
	if (i < hi) quick_sort_plist_y(a, i, hi);
}

void print_sorted(struct slist *a, struct DotList *dots, int num)
{
	int i;

	for( i = 0; i < num; i++ )
	{
		printf("%d-%d %d-%d\n", dots[a[i].id].x.lower, dots[a[i].id].x.upper, dots[a[i].id].y.lower, dots[a[i].id].y.upper);
	}
}

void print_sorted_plist(struct perm_pt *a, int l, int u, int cutdim)
{
	int i;

	for( i = l; i <= u; i++ )
	{
		if( cutdim == 1) printf("%d, ", a[i].x_pt);
		else if( cutdim == 2) printf("%d, ", a[i].y_pt);
		else if( cutdim == 0) 
		{
			printf("%d-%d, ", a[i].x_pt, a[i].y_pt);
		}
	}
	printf("\n");
}

void add_symmetric_points(struct DotList *s_dots, struct DotList *dots, int num)
{
	int i;

	for( i = 0; i < num; i++ )
	{
		s_dots[i].x = assign_I(dots[i].y.lower, dots[i].y.upper);
		s_dots[i].y = assign_I(dots[i].x.lower, dots[i].x.upper);
		s_dots[i].x_int = assign_I(dots[i].y_int.lower, dots[i].y_int.upper);
		s_dots[i].y_int = assign_I(dots[i].x_int.lower, dots[i].x_int.upper);
		s_dots[i].sign = dots[i].sign;
		s_dots[i].identity = dots[i].identity;
	}
}

void initialize_algns(struct DotList *temp, int count)
{
  int i = 0;

  for( i = 0; i < count; i++ ) {
    temp[i].x = assign_I(0, 1);
    temp[i].y = assign_I(0, 1);
    temp[i].identity = 0;
    temp[i].l_id = -1;
    temp[i].m_id = -1;
    temp[i].l_pid = -1;
    temp[i].fid = 0;
    temp[i].lock = -1;
    temp[i].left_diff = 0;
    temp[i].right_diff = 0;
    temp[i].sign = DELETED;
    temp[i].init_sign = DELETED;
    temp[i].c_id = -1;
    temp[i].s_id = -1;
    temp[i].m_x = assign_I(0,1);
    temp[i].x_int = assign_I(0,1);
    temp[i].m_y = assign_I(0,1);
    temp[i].y_int = assign_I(0,1);
    temp[i].pair_self = SELF;
    temp[i].rp1_id = -1;
    temp[i].rp2_id = -1;
		temp[i].xr_diff = 0; 
		temp[i].xl_diff = 0; 
		temp[i].yr_diff = 0; 
		temp[i].yl_diff = 0; 
		temp[i].xr_offset = 0; 
		temp[i].xl_offset = 0; 
		temp[i].yr_offset = 0; 
		temp[i].yl_offset = 0; 
		temp[i].len1 = 0; 
		temp[i].len2 = 0; 
		strcpy(temp[i].name1, "");
		strcpy(temp[i].name2, "");
  }
}
void sort_init_algns(struct slist *st, struct DotList *dots, int num, int mode)
{
  int i, j;
  int cur_val;

  for( i = 0; i < num; i++ )
  {
    st[i].id = i;
    if( mode == SELF1 ) st[i].val = dots[i].x.lower + dots[i].xl_diff;
    else if( mode == SELF2 ) st[i].val = dots[i].y.lower + dots[i].yl_diff;
    else if( mode == INIT_SELF1) st[i].val = dots[i].x.lower;
    else if( mode == INIT_SELF2 ) st[i].val = dots[i].y.lower;
    else if( mode == INIT_PAIR ) st[i].val = dots[i].x.lower + dots[i].xl_diff;
  }
  quick_sort_inc(st, 0, num-1);

  i = 0;
  while(i < num)
  {
    j = 0;
    cur_val = st[i].val;
    while(((i+j) < num) && (cur_val == st[i+j].val))
    {
      if( mode == SELF1 ) st[i+j].val = width(dots[st[i+j].id].x) - dots[st[i+j].id].xl_diff - dots[st[i+j].id].xr_diff;
      else if( mode == SELF2 ) st[i+j].val = width(dots[st[i+j].id].y) - dots[st[i+j].id].yl_diff - dots[st[i+j].id].yr_diff;
      else if( mode == INIT_SELF1) st[i+j].val = width(dots[st[i+j].id].x);
      else if( mode == INIT_SELF2 ) st[i+j].val = width(dots[st[i+j].id].y);
      if( mode == INIT_PAIR ) st[i+j].val = width(dots[st[i+j].id].x) - dots[st[i+j].id].xl_diff - dots[st[i+j].id].xr_diff;
      j++;
    }

    if( mode == INIT_PAIR ) quick_sort_inc(st, i, i+j-1);
    else quick_sort_dec(st, i, i+j-1);
    i = i+j;
  }
}

int search_range_b(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode)
{
  struct I val;
  int res;
  int i, cur;

  i = 0;
  cur = sorted[i].id;

  if(mode == SELF1) val = assign_I(algns[cur].x.lower, algns[cur].x.upper);
  else if(mode == SELF2) val = assign_I(algns[cur].y.lower, algns[cur].y.upper);

  while( ((i+1) < num_algns) && (val.upper < query) ) {
    i++;
    cur = sorted[i].id;
    if(mode == SELF1) val = assign_I(algns[cur].x.lower, algns[cur].x.upper);
    else if(mode == SELF2) val = assign_I(algns[cur].y.lower, algns[cur].y.upper);
  }

  if( i >= (num_algns-1)) res = num_algns-1;
  else res = i;

  return(res);
}

int search_range_e(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode)
{
  struct I val;
  int res;
  int i, cur;

  i = 0;
  cur = sorted[i].id;

  if(mode == SELF1) val = assign_I(algns[cur].x.lower, algns[cur].x.upper);
  else if(mode == SELF2) val = assign_I(algns[cur].y.lower, algns[cur].y.upper);

  while( ((i+1) < num_algns) && (val.lower <= query) ) {
    i++;
    cur = sorted[i].id;
    if(mode == SELF1) val = assign_I(algns[cur].x.lower, algns[cur].x.upper);
    else if(mode == SELF2) val = assign_I(algns[cur].y.lower, algns[cur].y.upper);
  }

  if( i >= (num_algns-1)) res = num_algns-1;
  else res = i;

  return(res);
}

void assign_algn(struct DotList *temp, int loc, struct DotList cur)
{
  temp[loc].x = assign_I(cur.x.lower, cur.x.upper);
  temp[loc].y = assign_I(cur.y.lower, cur.y.upper);
  temp[loc].identity = cur.identity;
  temp[loc].fid = cur.fid;
  temp[loc].sign = cur.sign;
  temp[loc].l_id = cur.l_id;
  temp[loc].l_pid = cur.l_pid;
  temp[loc].lock = cur.lock;
  temp[loc].left_diff = cur.left_diff;
  temp[loc].right_diff = cur.right_diff;
  temp[loc].c_id = cur.c_id;
  temp[loc].m_id = cur.m_id;
  temp[loc].m_x = assign_I(cur.m_x.lower, cur.m_x.upper);
  temp[loc].m_y = assign_I(cur.m_y.lower, cur.m_y.upper);
  temp[loc].pair_self = cur.pair_self;
  temp[loc].rp1_id = cur.rp1_id;
  temp[loc].rp2_id = cur.rp2_id;
  temp[loc].xl_diff = cur.xl_diff;
  temp[loc].xr_diff = cur.xr_diff;
  temp[loc].yl_diff = cur.yl_diff;
  temp[loc].yr_diff = cur.yr_diff;
  temp[loc].xl_offset = cur.xl_offset;
  temp[loc].xr_offset = cur.xr_offset;
  temp[loc].yl_offset = cur.yl_offset;
  temp[loc].yr_offset = cur.yr_offset;
  temp[loc].index = cur.index;
  temp[loc].init_sign = cur.init_sign;
}

struct gap_list assign_glist(struct gap_list a)
{
  struct gap_list res;

  res.type = a.type;
  res.id1 = a.id1;
  res.id2 = a.id2;
  res.x1 = a.x1;
  res.x2 = a.x2;
  res.y1 = a.y1;
  res.y2 = a.y2;
	strcpy(res.name1, a.name1);
	strcpy(res.name2, a.name2);

  return(res);
}

void init_rlist(struct r_list *a, int b, int e)
{
	int i = 0;
	
	for( i = b; i <=e; i++ ) {
		a[i].id = -1;
		a[i].start = 0;
		a[i].end = 1;
		a[i].d_rate = (float)100;
		strcpy(a[i].name, "");
	}
}
