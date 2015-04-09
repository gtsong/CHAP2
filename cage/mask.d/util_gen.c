#include "main.h"
#include "regions.h"
#include "util_gen.h"
#include "util_i.h"

void sort_by_width(struct slist *st, struct DotList *dots, int num)
{
	int i, j;
	int cur_val;

	for( i = 0; i < num; i++ )
	{
		st[i].val = width(dots[st[i].id].x); 
	}
	quick_sort_dec(st, 0, num-1);

	i = 0;
	while(i < num)
	{
		j = 0;
		cur_val = st[i].val;
		while(((i+j) < num) && (cur_val == st[i+j].val))
		{
			st[i+j].val = dots[st[i+j].id].identity; // This part could be changed
			j++;
		}
		quick_sort_dec(st, i, i+j-1);
		i = i+j;
	}
}

void sort_by_pid(struct slist *st, struct DotList *dots, int num)
{
	int i, j;
	int cur_val;

	for( i = 0; i < num; i++ )
	{
		st[i].val = dots[st[i].id].identity;
	}
	quick_sort_dec(st, 0, num-1);

	i = 0;
	while(i < num)
	{
		j = 0;
		cur_val = st[i].val;
		while(((i+j) < num) && (cur_val == st[i+j].val))
		{
			st[i+j].val = width(dots[st[i+j].id].x); // This part could be changed
			j++;
		}
		quick_sort_dec(st, i, i+j-1);
		i = i+j;
	}
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
		while( ((i+j) < num) && (cur_val == st[i+j].val))
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
	res.sp_state = a.sp_state;
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

void print_sorted(struct slist *a, struct DotList *dots, int num)
{
	int i;

	for( i = 0; i < num; i++ )
	{
		printf("%d-%d %d-%d\n", dots[a[i].id].x.lower, dots[a[i].id].x.upper, dots[a[i].id].y.lower, dots[a[i].id].y.upper);
	}
}

void add_symmetric_points(struct DotList *s_dots, struct DotList *dots, int num)
{
	int i;

	for( i = 0; i < num; i++ )
	{
		s_dots[i].x = assign_I(dots[i].y.lower, dots[i].y.upper);
		s_dots[i].y = assign_I(dots[i].x.lower, dots[i].x.upper);
		s_dots[i].m_x = assign_I(dots[i].m_y.lower, dots[i].m_y.upper);
		s_dots[i].m_y = assign_I(dots[i].m_x.lower, dots[i].m_x.upper);
		s_dots[i].sign = dots[i].sign;
		s_dots[i].identity = dots[i].identity;
		s_dots[i].l_id = dots[i].l_id;
		s_dots[i].lock = dots[i].lock;
	}
}

void quick_sort_plist_x(struct perm_pt *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
  int i=lo, j=hi;
	struct perm_pt h;
  int x;
	//  int a_val, b_val;

  x = a[(lo+hi)/2].x_pt;
	//  a_val = a[lo].x_pt;
	//  b_val = a[hi].x_pt;

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
//  int a_val, b_val;

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

struct perm_pt assign_pm_val(struct perm_pt a)
{
  struct perm_pt res;

  res.id = a.id;
  res.x_pt = a.x_pt;
  res.y_pt = a.y_pt;
  res.sign = a.sign;

  return(res);
}
