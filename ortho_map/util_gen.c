#include "main.h"
#include "regions.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "util_ops.h"

extern bool debug_mode;

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

void sort_by_width_y(struct slist *st, struct DotList *dots, int num)
{
	int i, j;
	int cur_val;

	for( i = 0; i < num; i++ )
	{
		st[i].val = width(dots[st[i].id].y); 
	}
	quick_sort_dec(st, 0, num-1);

	i = 0;
	while(i < num)
	{
		j = 0;
		cur_val = st[i].val;
		while( ((i+j) < num) && (cur_val == st[i+j].val))
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
			st[i+j].val = width(dots[st[i+j].id].y); // This part could be changed
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

void sort_by_loc_short_alist(struct short_alist *list, int num)
{
  int i, j;
  int cur_val;

  quick_sort_inc_alist(list, 0, num-1);

  i = 0;
  while(i < num)
  {
    j = 0;
    cur_val = list[i].val;
    while(((i+j) < num) && (cur_val == list[i+j].val)) {
			list[i+j].val = abs(list[i+j].y.upper - list[i+j].y.lower);
			j++;
		}
    quick_sort_dec_alist(list, i, i+j-1);
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

void sort_init_algns(struct slist *st, struct DotList *dots, int num, int mode)
{
 	int i = 0, j = 0, k = 0, n = 0;
	int cur_val = 0;

 	for( i = 0; i < num; i++ )
	{
		st[i].id = i;
		if( mode == SELF1 ) st[i].val = dots[i].ctg_id1;
		else if( mode == SELF2 ) st[i].val = dots[i].ctg_id2;
    else if( mode == INIT_SELF1) st[i].val = dots[i].ctg_id1;
    else if( mode == INIT_SELF2 ) st[i].val = dots[i].ctg_id2;
    else if( mode == INIT_PAIR ) st[i].val = dots[i].ctg_id1;
  }

  quick_sort_inc(st, 0, num-1);

  i = 0;
  while( i < num ) 
  {
    j = 0;
    cur_val = st[i].val;

    while(((i+j) < num) && (cur_val == st[i+j].val) ) {
      if( mode == SELF1 ) st[i+j].val = dots[st[i+j].id].x.lower + dots[st[i+j].id].xl_diff;
      else if( mode == SELF2 ) st[i+j].val = dots[st[i+j].id].y.lower + dots[st[i+j].id].yl_diff;
      else if( mode == INIT_SELF1) st[i+j].val = dots[st[i+j].id].x.lower;
      else if( mode == INIT_SELF2 ) st[i+j].val = dots[st[i+j].id].y.lower;
      else if( mode == INIT_PAIR ) st[i+j].val = dots[st[i+j].id].x.lower + dots[st[i+j].id].xl_diff;
      j++;
    }
    quick_sort_inc(st, i, i+j-1);

    k = i;
    while( k < (i+j) ) {
      cur_val = st[k].val;
      n = 0;
      while(((k+n) < (i+j)) && (cur_val == st[n+k].val))
      {
        if( mode == SELF1 ) st[n+k].val = width(dots[st[n+k].id].x) - dots[st[n+k].id].xl_diff - dots[st[n+k].id].xr_diff;
        else if( mode == SELF2 ) st[n+k].val = width(dots[st[n+k].id].y) - dots[st[n+k].id].yl_diff - dots[st[n+k].id].yr_diff;
        else if( mode == INIT_SELF1) st[n+k].val = width(dots[st[n+k].id].x);
        else if( mode == INIT_SELF2 ) st[n+k].val = width(dots[st[n+k].id].y);
        else if( mode == INIT_PAIR ) st[n+k].val = width(dots[st[n+k].id].x) - dots[st[n+k].id].xl_diff - dots[st[n+k].id].xr_diff;
        n++;
      }

      if( mode == INIT_PAIR ) quick_sort_inc(st, k, n+k-1);
      else quick_sort_dec(st, k, n+k-1);
      k = k+n;
    }

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
		while ((i <= hi) && (a[i].val>x)) i++; 
		while ((j >= lo) && (a[j].val<x)) j--;
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
		while ((i <= hi) && (a[i].val<x)) i++; 
		while ((j >= lo) && (a[j].val>x)) j--;
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

void quick_sort_dec_alist(struct short_alist *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct short_alist h;
	int x=a[(lo+hi)/2].val;

//  partition
	do
	{    
		while ((i <= hi) && (a[i].val>x)) i++; 
		while ((j >= lo) && (a[j].val<x)) j--;
		if (i<=j)
		{
			h = assign_alist(a[i]);
			a[i] = assign_alist(a[j]);
			a[j] = assign_alist(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_alist(a, lo, j);
	if (i < hi) quick_sort_dec_alist(a, i, hi);
}

void quick_sort_inc_alist(struct short_alist *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct short_alist h;
	int x=a[(lo+hi)/2].val;

//  partition
	do
	{    
		while ((i <= hi) && (a[i].val<x)) i++; 
		while ((j >= lo) && (a[j].val>x)) j--;
		if (i<=j)
		{
			h = assign_alist(a[i]);
			a[i] = assign_alist(a[j]);
			a[j] = assign_alist(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_alist(a, lo, j);
	if (i < hi) quick_sort_inc_alist(a, i, hi);
}

struct short_alist assign_alist(struct short_alist a)
{
	struct short_alist res;

	res.id = a.id;
	res.x = assign_I(a.x.lower, a.x.upper);
	res.y = assign_I(a.y.lower, a.y.upper);
	res.val = a.val;

	return(res);
}

void initialize_alist(struct short_alist *a, int from, int to)
{
	int i = 0;
	for( i = from; i < to; i++ ) {
		a[i].id = -1;
		a[i].x = assign_I(0,1);
		a[i].y = assign_I(0,1);
		a[i].val = 0;
	}
}

void initialize_slist(struct slist *a, int from, int to)
{
	int i = 0;
	for( i = from; i < to; i++ ) {
		a[i].id = i;
		a[i].val = 0;
		a[i].val_red = 0;
		a[i].sp_state = 0;
		a[i].add_sp_state = 0;
		a[i].is_x = true;
	}
}

struct slist assign_slist(struct slist a)
{
	struct slist res;

	res.id = a.id;
	res.val = a.val;
	res.sp_state = a.sp_state;
	res.add_sp_state = a.add_sp_state;
	res.is_x = a.is_x;
	res.val_red = a.val_red;

	return(res);
}

struct gap_list assign_glist(struct gap_list a)
{
	struct gap_list res;

	res.gid = a.gid;
	res.type = a.type;
	res.id1 = a.id1;
	res.id2 = a.id2;
	res.x1 = a.x1;
	res.x2 = a.x2;
	res.y1 = a.y1;
	res.y2 = a.y2;
	res.offset = a.offset;

	return(res);
}

void sort_int(struct I *a, int num)
{
	int i, j;

  quick_sort_inc_int(a, 0, num-1, POS_BASE);  
	i = 0;  
	while(i < num)  
	{    
		j = 0;    
		while(((i+j) < num) && (a[i].lower == a[i+j].lower)) j++;
    quick_sort_dec_int(a, i, i+j-1, LEN_BASE);
    i = i+j;
  }
}

void sort_exons(struct exons_list *a, int num)
{
	int i, j;

	quick_sort_inc_exons(a, 0, num-1, CTG_BASE);
	i = 0;
	while(i < num) {
		j = 0;
		while(((i+j) < num) && (a[i].ctg_id == a[i+j].ctg_id)) j++;
    quick_sort_inc_exons(a, i, i+j-1, POS_BASE);
    i = i+j;
	}

  quick_sort_inc_exons(a, 0, num-1, POS_BASE);  
	i = 0;  
	while(i < num)  
	{    
		j = 0;    
		while(((i+j) < num) && (a[i].ctg_id == a[i+j].ctg_id) && (a[i].reg.lower == a[i+j].reg.lower)) j++;
    quick_sort_dec_exons(a, i, i+j-1, LEN_BASE);
    i = i+j;
  }
}

void quick_sort_inc_exons(struct exons_list *a, int lo, int hi, int mode)
{	
	int i=lo, j=hi;
	struct exons_list h;
	int x;
	if( mode == POS_BASE ) x=a[(lo+hi)/2].reg.lower;
	else if( mode == LEN_BASE ) x=width(a[(lo+hi)/2].reg);
	else if( mode == CTG_BASE ) x=a[(lo+hi)/2].ctg_id;
	else {
		fatalf("not supported mode: %d in quick_sort_inc_exons\n", mode);
	}

	do
	{    
		if( mode == POS_BASE ) {
			while ((i <= hi) && (a[i].reg.lower<x)) i++; 
			while ((j >= lo) && (a[j].reg.lower>x)) j--;
		}
		else if( mode == LEN_BASE ) {
			while ((i <= hi) && (width(a[i].reg)<x)) i++; 
			while ((j >= lo) && (width(a[j].reg)>x)) j--;
		}
		else if( mode == CTG_BASE ) {
			while ((i <= hi) && (a[i].ctg_id<x)) i++; 
			while ((j >= lo) && (a[j].ctg_id>x)) j--;
		}

		if (i<=j)
		{
			h = assign_exons(a[i]);
			strcpy(h.name, a[i].name);
			if( i != j) {
				a[i] = assign_exons(a[j]);
				strcpy(a[i].name, a[j].name);
			}
			a[j] = assign_exons(h);
			strcpy(a[j].name, h.name);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_exons(a, lo, j, mode);
	if (i < hi) quick_sort_inc_exons(a, i, hi, mode);
}

void quick_sort_dec_exons(struct exons_list *a, int lo, int hi, int mode)
{	
	int i=lo, j=hi;
	struct exons_list h;
	int x;

	if( mode == POS_BASE ) x=a[(lo+hi)/2].reg.lower;
	else if( mode == LEN_BASE ) x=width(a[(lo+hi)/2].reg);
	do
	{    
		if( mode == POS_BASE ) {
			while ((i <= hi) && (a[i].reg.lower>x)) i++; 
			while ((j >= lo) && (a[j].reg.lower<x)) j--;
		}
		else if( mode == LEN_BASE ) {
			while ((i <= hi) && (width(a[i].reg)>x)) i++; 
			while ((j >= lo) && (width(a[j].reg)<x)) j--;
		}
		if (i<=j)
		{
			h = assign_exons(a[i]);
			strcpy(h.name, a[i].name);

			if( i != j ) {
				a[i] = assign_exons(a[j]);
				strcpy(a[i].name, a[j].name);
			}
			a[j] = assign_exons(h);
			strcpy(a[j].name, h.name);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_exons(a, lo, j, mode);
	if (i < hi) quick_sort_dec_exons(a, i, hi, mode);
}

void quick_sort_num_list_inc(int *a, int lo, int hi)
{	
	int i=lo, j=hi;
	int h;
	int x;
	x=a[(lo+hi)/2];

	do
	{    
		while ((i <= hi) && (a[i]<x)) i++; 
		while ((j >= lo) && (a[j]>x)) j--;

		if (i<=j)
		{
			h = a[i];
			a[i] = a[j];
			a[j] = h;
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_num_list_inc(a, lo, j);
	if (i < hi) quick_sort_num_list_inc(a, i, hi);
}

void quick_sort_inc_int(struct I *a, int lo, int hi, int mode)
{	
	int i=lo, j=hi;
	struct I h;
	int x;
	if( mode == POS_BASE ) x=a[(lo+hi)/2].lower;
	else if( mode == LEN_BASE ) x=width(a[(lo+hi)/2]);

	do
	{    
		if( mode == POS_BASE ) {
			while ((i <= hi) && (a[i].lower<x)) i++; 
			while ((j >= lo) && (a[j].lower>x)) j--;
		}
		else if( mode == LEN_BASE ) {
			while ((i <= hi) && (width(a[i])<x)) i++; 
			while ((j >= lo) && (width(a[j])>x)) j--;
		}

		if (i<=j)
		{
			h = assign_I(a[i].lower, a[i].upper);
			a[i] = assign_I(a[j].lower, a[j].upper);
			a[j] = assign_I(h.lower, h.upper);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_int(a, lo, j, mode);
	if (i < hi) quick_sort_inc_int(a, i, hi, mode);
}

//  lo is the lower index, hi is the upper index
void quick_sort_dec_int(struct I *a, int lo, int hi, int mode)
{	
	int i=lo, j=hi;
	struct I h;
	int x;

	if( mode == POS_BASE ) x=a[(lo+hi)/2].lower;
	else if( mode == LEN_BASE ) x=width(a[(lo+hi)/2]);
	do
	{    
		if( mode == POS_BASE ) {
			while ((i >= hi) && (a[i].lower>x)) i++; 
			while ((j <= lo) && (a[j].lower<x)) j--;
		}
		else if( mode == LEN_BASE ) {
			while ((i >= hi) && (width(a[i])>x)) i++; 
			while ((j <= lo) && (width(a[j])<x)) j--;
		}
		if (i<=j)
		{
			h = assign_I(a[i].lower, a[i].upper);
			a[i] = assign_I(a[j].lower, a[j].upper);
			a[j] = assign_I(h.lower, h.upper);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_int(a, lo, j, mode);
	if (i < hi) quick_sort_dec_int(a, i, hi, mode);
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
		s_dots[i].pair_self = dots[i].pair_self;
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
    while ((i <= hi) && (a[i].x_pt<x))
    {
      i++;
    }
    while ((j >= lo) && (a[j].x_pt>x))
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
    while ((i <= hi) && (a[i].y_pt<x))
    {
      i++;
    }
    while ((j >= lo) && (a[j].y_pt>x))
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

void make_into_one(int num, struct DotList *dots, struct DotList *s_dots, struct DotList *t_dots)
{
  int i;
  for(i = 0; i < num; i++)
  {
    assign_item(t_dots, 2*i, dots[i].x, dots[i].y, dots[i].identity, dots[i].sign);
    assign_item(t_dots, 2*i+1, s_dots[i].x, s_dots[i].y, s_dots[i].identity, s_dots[i].sign);
  }
}

int increase_count(struct DotList *org, struct DotList *dots, int id, int c, bool is_x)
{
  struct I cur_reg;
  struct I temp;

  if( is_x == true ) cur_reg = assign_I(org[id].x.lower, org[id].x.upper);
  else cur_reg = assign_I(org[id].y.lower, org[id].y.upper);

  if( dots[c].sign != 2 )
  {
  	temp = assign_I(dots[c].y.lower, dots[c].y.upper);

    if( loose_subset(temp, cur_reg) == true )
    {
      return 1;
    }
    else return 0;
  }
  else return 0;
}

void assign_item(struct DotList *dots, int loc, struct I x, struct I y, int id, int s)
{
  dots[loc].x = assign_I(x.lower, x.upper);
  dots[loc].y = assign_I(y.lower, y.upper);
  dots[loc].identity = id;
  dots[loc].sign = s;
  dots[loc].l_id = -1;
  dots[loc].lock = -1;
	dots[loc].c_id = -1;
}

void assign_algn(struct DotList *temp, int loc, struct DotList cur)
{
  temp[loc].x = assign_I(cur.x.lower, cur.x.upper);
  temp[loc].y = assign_I(cur.y.lower, cur.y.upper);
  temp[loc].identity = cur.identity;
	temp[loc].fid = cur.fid;
	temp[loc].indiv_fid = cur.indiv_fid;
	temp[loc].m_pid = cur.m_pid;
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
	temp[loc].sp_id = cur.sp_id;
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
	temp[loc].len1 = cur.len1;
  temp[loc].len2 = cur.len2;
  strcpy(temp[loc].name1, cur.name1);
  strcpy(temp[loc].name2, cur.name2);
	temp[loc].ctg_id1 = cur.ctg_id1;
  temp[loc].ctg_id2 = cur.ctg_id2;
}

void initialize_algns(struct DotList *temp, int b, int e)
{
	int i = 0;

	for( i = b; i < e; i++ ) {
	  temp[i].x = assign_I(0, 1);
	  temp[i].y = assign_I(0, 1);
 		temp[i].identity = 0;
		temp[i].fid = -1;
		temp[i].m_pid = 0;
	  temp[i].sign = DELETED;
 		temp[i].l_id = -1;
 		temp[i].l_pid = -1;
		temp[i].lock = -1;
		temp[i].left_diff = 0;
		temp[i].right_diff = 0;
		temp[i].c_id = -1;
		temp[i].m_id = -1;
 		temp[i].m_x = assign_I(0,1);
 		temp[i].m_y = assign_I(0,1);
		temp[i].pair_self = SELF;
		temp[i].sp_id = SELF1;
		temp[i].rp1_id = -1;
		temp[i].rp2_id = -1;
		temp[i].xl_diff = 0;
		temp[i].xr_diff = 0;
		temp[i].yl_diff = 0;
		temp[i].yr_diff = 0;
		temp[i].xl_offset = 0;
		temp[i].xr_offset = 0;
		temp[i].yl_offset = 0;
		temp[i].yr_offset = 0;
		temp[i].index = -1;
		temp[i].init_sign = DELETED;
		temp[i].indiv_fid = -1;
		temp[i].len1 = 0;
	  temp[i].len2 = 0;
 		strcpy(temp[i].name1, "");
 		strcpy(temp[i].name2, "");
		temp[i].ctg_id1 = -1;
	  temp[i].ctg_id2 = -1;
	}
}

void mark_match_algn(struct DotList cur, struct DotList *algns, int num_algns, int sp_id, int mode)
{
	int i = 0, j = 0;
	int num_candi = 0;
	struct DotList *cur_algns;
	struct slist *st;
	int num_count = 0;
	bool is_x;
	bool is_merged;
	struct DotList *temp;
	int cur_id;

	temp = (struct DotList *) ckalloc(sizeof(struct DotList));

	for( i = 0; i < num_algns; i++ )
	{
		if( algns[i].sp_id == sp_id )
		{
			if( ((is_same_sign(cur.sign, algns[i].sign)) && proper_overlap(cur.x, algns[i].x) && (width(intersect(cur.x, algns[i].x)) >= 50) ) && ( is_same_sign(cur.sign, algns[i].sign) && proper_overlap(cur.y, algns[i].y) && (width(intersect(cur.y, algns[i].y)) >= 50) )) num_candi++;
		}
	}

	cur_algns = (struct DotList *) ckalloc(num_candi * sizeof(struct DotList));
	st = (struct slist *) ckalloc(num_candi * sizeof(struct slist));

	for( i = 0; i < num_algns; i++ )
	{
		if( algns[i].sp_id == sp_id )
		{
			if( ( (is_same_sign(cur.sign, algns[i].sign)) && proper_overlap(cur.x, algns[i].x) && (width(intersect(cur.x, algns[i].x)) >= 50) ) && ( is_same_sign(cur.sign, algns[i].sign) && proper_overlap(cur.y, algns[i].y) && (width(intersect(cur.y, algns[i].y)) >= 50) )) 
			{
				assign_algn(cur_algns, j, algns[i]);	
				cur_algns[j].fid = i;
				j++;
			}
		}
	}

	for( i = 0; i < num_candi; i++ ) st[i].id = i;

	if( width(cur.x) >= width(cur.y) ) {
		sort_by_width_y(st, cur_algns, num_candi);
		is_x = false;
	}
	else {
		sort_by_width(st, cur_algns, num_candi);
		is_x = true;
	}

	for( i = 0; i < num_candi; i++ )
	{
		if( cur_algns[i].sign == DELETED ) {}
		else {
			for( j = (i+1); j < num_candi; j++ )
			{
				if( cur_algns[i].sign == DELETED ) break;
				else if( cur_algns[j].sign == DELETED ) {}
				else if( !is_x ){
					if( strict_almost_equal(cur_algns[st[i].id].y, cur_algns[st[j].id].y) && (cur_algns[st[i].id].identity <= cur_algns[st[j].id].identity) ) cur_algns[st[i].id].sign = DELETED;	
					else if( strict_subset(cur_algns[st[j].id].y, cur_algns[st[i].id].y) ) 
					{
						if(cur_algns[st[j].id].identity < (cur_algns[st[i].id].identity + 2)) 
							cur_algns[st[j].id].sign = DELETED;
						else cur_algns[st[i].id].identity = DELETED;
					}
					else if( proper_overlap(cur_algns[st[i].id].y, cur_algns[st[j].id].y) && (width(intersect(cur_algns[st[i].id].y, cur_algns[st[j].id].y)) >= ERR_TH)) {
						if(cur_algns[st[i].id].identity < (cur_algns[st[j].id].identity - 1)) cur_algns[st[i].id].sign = DELETED;
						else cur_algns[st[j].id].sign = DELETED;
					}
				}
				else {
					if( strict_almost_equal(cur_algns[st[i].id].x, cur_algns[st[j].id].x) && (cur_algns[st[i].id].identity <= cur_algns[st[j].id].identity) ) cur_algns[st[i].id].sign = DELETED;	
					else if( strict_subset(cur_algns[st[j].id].x, cur_algns[st[i].id].x) ) 
					{
						if(cur_algns[st[j].id].identity < (cur_algns[st[i].id].identity + 2)) 
							cur_algns[st[j].id].sign = DELETED;
						else cur_algns[st[i].id].identity = DELETED;
					}
					else if( proper_overlap(cur_algns[st[i].id].x, cur_algns[st[j].id].x) && (width(intersect(cur_algns[st[i].id].x, cur_algns[st[j].id].x)) >= ERR_TH)) {
						if(cur_algns[st[i].id].identity < (cur_algns[st[j].id].identity - 1)) cur_algns[st[i].id].sign = DELETED;
						else cur_algns[st[j].id].sign = DELETED;
					}
				}
			}

			if( cur_algns[st[i].id].sign != DELETED )
			{
				if( num_count == 0 ) {
					assign_algn(temp, 0, cur_algns[st[i].id]);
					is_merged = true;
				}
				else {
					is_merged = false;
					is_merged = merge_into_one_algn(temp, cur_algns[st[i].id]);
				}

				if( is_merged ) {
					cur_id = cur_algns[st[i].id].fid;
					algns[cur_id].sign = assign_sign( algns[cur_id].sign, mode);
					if( debug_mode ) printf("%d-%d, %d-%d changes a sign into %d\n", algns[cur_id].x.lower, algns[cur_id].x.upper, algns[cur_id].y.lower, algns[cur_id].y.upper, algns[cur_id].sign);
				}
				num_count++;
			}
		}
	}

	free(temp);
	free(cur_algns);
	free(st);
}

int assign_sign(int sign, int mode)
{
	int res = sign;

	if( (sign == ORTHO) || (sign == ORTHO_COMP) ) res = sign;
	else if( mode == CANDI_MODE ) 
  {
	  if( sign == 0 ) res = TEMP_HIDDEN;
    else if( sign == 1 ) res = TEMP_HIDDEN_COMP;
  }
  else if( mode == ORTHO_MODE ) 
  {
    if( (sign == TEMP_HIDDEN) || (sign == 0) ) res = ORTHO;
		else if( (sign == TEMP_HIDDEN_COMP) || (sign == 1) ) res = ORTHO_COMP;
	}
	else if( mode == BACK_MODE )
	{
    if( sign == TEMP_HIDDEN) res  = 0;
		else if(sign == TEMP_HIDDEN_COMP) res = 1;
	}
	else res = sign;

	return(res);
}

bool is_same_sign(int sign1, int sign2)
{
	bool res = false;

	if( sign1 == 0 ) {
		if( (sign2 == 0) || (sign2 == TEMP_HIDDEN) || (sign2 == ORTHO) ) res = true;
	}
	else if( sign1 == 1 ) {
		if( (sign2 == 1) || (sign2 == TEMP_HIDDEN_COMP) || (sign2 == ORTHO_COMP) ) res = true;
	}

	return(res);
}

int find_match_algn(struct DotList cur, struct DotList *algns, int num_algns, struct DotList *temp, int sp_id)
{
	int i = 0, j = 0;
	int num_candi = 0;
	struct DotList *cur_algns;
	struct slist *st;
	int num_count = 0;
	bool is_x;
	bool is_merged;

	for( i = 0; i < num_algns; i++ )
	{
		if( algns[i].sp_id == sp_id )
		{
			if( ((is_same_sign(cur.sign, algns[i].sign)) && proper_overlap(cur.x, algns[i].x) && (width(intersect(cur.x, algns[i].x)) >= 50) ) && ( is_same_sign(cur.sign, algns[i].sign) && proper_overlap(cur.y, algns[i].y) && (width(intersect(cur.y, algns[i].y)) >= 50) )) num_candi++;
		}
	}

	cur_algns = (struct DotList *) ckalloc(num_candi * sizeof(struct DotList));
	st = (struct slist *) ckalloc(num_candi * sizeof(struct slist));

	for( i = 0; i < num_algns; i++ )
	{
		if( algns[i].sp_id == sp_id )
		{
			if( ( (is_same_sign(cur.sign, algns[i].sign)) && proper_overlap(cur.x, algns[i].x) && (width(intersect(cur.x, algns[i].x)) >= 50) ) && ( (is_same_sign(cur.sign, algns[i].sign)) && proper_overlap(cur.y, algns[i].y) && (width(intersect(cur.y, algns[i].y)) >= 50) )) 
			{
				assign_algn(cur_algns, j, algns[i]);	
				cur_algns[j].sign = cur.sign;
				j++;
			}
		}
	}

	for( i = 0; i < num_candi; i++ ) st[i].id = i;

	if( width(cur.x) >= width(cur.y) ) {
		sort_by_width_y(st, cur_algns, num_candi);
		is_x = false;
	}
	else {
		sort_by_width(st, cur_algns, num_candi);
		is_x = true;
	}

	for( i = 0; i < num_candi; i++ )
	{
		if( cur_algns[i].sign == DELETED ) {}
		else {
			for( j = (i+1); j < num_candi; j++ )
			{
				if( cur_algns[i].sign == DELETED ) break;
				else if( cur_algns[j].sign == DELETED ) {}
				else if( !is_x ){
					if( strict_almost_equal(cur_algns[st[i].id].y, cur_algns[st[j].id].y) && (cur_algns[st[i].id].identity <= cur_algns[st[j].id].identity) ) cur_algns[st[i].id].sign = DELETED;	
					else if( strict_subset(cur_algns[st[j].id].y, cur_algns[st[i].id].y) ) 
					{
						if(cur_algns[st[j].id].identity < (cur_algns[st[i].id].identity + 2)) 
							cur_algns[st[j].id].sign = DELETED;
						else cur_algns[st[i].id].identity = DELETED;
					}
					else if( proper_overlap(cur_algns[st[i].id].y, cur_algns[st[j].id].y) && (width(intersect(cur_algns[st[i].id].y, cur_algns[st[j].id].y)) >= ERR_TH)) {
						if(cur_algns[st[i].id].identity < (cur_algns[st[j].id].identity - 1)) cur_algns[st[i].id].sign = DELETED;
						else cur_algns[st[j].id].sign = DELETED;
					}
				}
				else {
					if( strict_almost_equal(cur_algns[st[i].id].x, cur_algns[st[j].id].x) && (cur_algns[st[i].id].identity <= cur_algns[st[j].id].identity) ) cur_algns[st[i].id].sign = DELETED;	
					else if( strict_subset(cur_algns[st[j].id].x, cur_algns[st[i].id].x) ) 
					{
						if(cur_algns[st[j].id].identity < (cur_algns[st[i].id].identity + 2)) 
							cur_algns[st[j].id].sign = DELETED;
						else cur_algns[st[i].id].identity = DELETED;
					}
					else if( proper_overlap(cur_algns[st[i].id].x, cur_algns[st[j].id].x) && (width(intersect(cur_algns[st[i].id].x, cur_algns[st[j].id].x)) >= ERR_TH)) {
						if(cur_algns[st[i].id].identity < (cur_algns[st[j].id].identity - 1)) cur_algns[st[i].id].sign = DELETED;
						else cur_algns[st[j].id].sign = DELETED;
					}
				}
			}

			if( cur_algns[st[i].id].sign != DELETED )
			{
				if( num_count == 0 ) {
					assign_algn(temp, 0, cur_algns[st[i].id]);
					is_merged = true;
				}
				else {
					is_merged = false;
					is_merged = merge_into_one_algn(temp, cur_algns[st[i].id]);
				}
				if( is_merged ) num_count++;
			}
		}
	}

	free(cur_algns);
	free(st);

	return(num_count);
}

bool is_match_algn(struct DotList cur, struct DotList cmp)
{
	struct I t_x, t_y;
	struct I temp;
	int b, e;
	int sign;

	if( (cmp.sign == ORTHO) || (cmp.sign == TEMP_HIDDEN) ) sign = 0;
	else if( (cmp.sign == ORTHO_COMP) || (cmp.sign == TEMP_HIDDEN_COMP) ) sign = 1;
	else sign = cmp.sign;

	if( (cur.sign == sign) && proper_overlap(cur.x, cmp.x) && proper_overlap(cur.y, cmp.y) ) {
		t_x = intersect(cur.x, cmp.x);
		t_y = intersect(cur.y, cmp.y);

		if( cur.sign == 0 ) {
			b = cur.y.lower + abs(cur.x.lower - t_x.lower); 	
			e = b + width(t_x);
			temp = assign_I(b, e);
			if( strict_almost_equal(t_y, temp) ) return(true);
			else return(false);
		}
		else if( cur.sign == 1 ) {
			e = cur.y.upper - abs(cur.x.lower - t_x.lower);
			b = e - width(t_x);
			temp = assign_I(b, e);
			if( strict_almost_equal(t_y, temp) ) return(true);
			else return(false);
		}
		else return(false);
	}
	else return(false);
}

bool merge_into_one_algn(struct DotList *temp, struct DotList cur)
{
	int b, e;
	bool is_merged = false;

	if(temp->sign == -1) {
		assign_algn(temp, 0, cur);
		is_merged = true;
	}
	else if( almost_subset(cur.x, (*temp).x) || almost_subset(cur.y, (*temp).y) || f_loose_subset(cur.x, (*temp).x, STRICT) || f_loose_subset(cur.y, (*temp).y, STRICT)) is_merged = false;
	else if( (temp->sign == cur.sign) && (temp->sp_id == cur.sp_id) ) {
		is_merged = true;
	 (*temp).identity = (((*temp).identity * width((*temp).x)) + (cur.identity * width(cur.x))) / (width((*temp).x) + width(cur.x));

		if( (*temp).x.lower < cur.x.lower ) b = (*temp).x.lower;
		else b = cur.x.lower;

		if( (*temp).x.upper > cur.x.upper ) e = (*temp).x.upper;
		else e = cur.x.upper;

		if( b < e ) (*temp).x = assign_I(b, e);
		else {
			printf("merged in a wrong way\n");
			exit(EXIT_FAILURE);
		}

		if( (*temp).y.lower < cur.y.lower ) b = (*temp).y.lower;
		else b = cur.y.lower;

		if( (*temp).y.upper > cur.y.upper ) e = (*temp).y.upper;
		else e = cur.y.upper;

		if( b < e ) (*temp).y = assign_I(b, e);
		else {
			printf("merged in a wrong way\n");
			exit(EXIT_FAILURE);
		}
	}
	return(is_merged);
}

void sort_rev_init_algns(struct slist *st, struct DotList *dots, int num, int mode)
{
	int i, j;
	int cur_val;

	for( i = 0; i < num; i++ )
	{
		st[i].id = i;
		if( mode == SELF1 ) st[i].val = dots[i].x.upper - dots[i].xr_diff;
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
			else if(mode == SELF2 ) st[i+j].val = width(dots[st[i+j].id].y) - dots[st[i+j].id].yl_diff - dots[st[i+j].id].yr_diff;
			j++;
		}
		if( mode == SELF1 ) quick_sort_dec(st, i, i+j-1);
		else if( mode == SELF2 ) quick_sort_inc(st, i, i+j-1);
		i = i+j;
	}
}

int search_range_b(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode)
{
	struct I val;
	int res;
	int i, cur;

	if( num_algns > 0 ) {
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
	}
	else res = 0;
		
	return(res);
}

int search_range_e(struct slist *sorted, struct DotList *algns, int num_algns,  int query, int mode)
{
	struct I val;
	int res;
	int i, cur;

	if( num_algns > 0 ) {
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
	}
	else res = 0;
	return(res);
}

int quick_search_range_b(struct slist *sorted, struct DotList *algns, int i, int j, int query, int mode)
{
	int mid;
	struct I val;
	int res;

	mid = (i+j)/2;
	if(mode == SELF1) val = assign_I(algns[sorted[mid].id].x.lower, algns[sorted[mid].id].x.upper);
	else if(mode == SELF2) val = assign_I(algns[sorted[mid].id].y.lower, algns[sorted[mid].id].y.upper);

	if((val.lower > query) || ((val.lower <= query) && (val.upper >= query))) {
		if( i >= (mid-1) ) return (mid-1);	
		else res = quick_search_range_b(sorted, algns, i, mid-1, query, mode);
	} 
	else {
		if( j <= (mid+1) ) return j;
		else res = quick_search_range_b(sorted, algns, mid+1, j, query, mode);
	}

	return(res);
}

int quick_search_range_e(struct slist *sorted, struct DotList *algns, int i, int j, int query, int mode)
{
	int mid;
	struct I val;
	int res;

	mid = (i+j)/2;
	if(mode == SELF1) val = assign_I(algns[sorted[mid].id].x.lower, algns[sorted[mid].id].x.upper);
	else if(mode == SELF2) val = assign_I(algns[sorted[mid].id].y.lower, algns[sorted[mid].id].y.upper);
	if((val.upper < query) || ((val.lower <= query) && (val.upper >= query))) {
		if( j <= (mid+1) ) return (mid+1);	
		else res = quick_search_range_e(sorted, algns, mid+1, j, query, mode);
	} 
	else {
		if( i >= (mid-1) ) return i;
		else res = quick_search_range_e(sorted, algns, i, mid, query, mode);
	}

	return(res);
}

float cal_avg_pid(struct DotList *algns, int num_algns)
{
  int i = 0, j = 0;
  float avg_pid = (float) 0;
  int cur_len = 0;
  int sum_int = 0;

  j = 0;
  avg_pid = (float)0;
  sum_int = 0;
  for( i = 0; i < num_algns; i++ ) {
    if( algns[i].sign != DELETED ) {
      cur_len = (algns[i].x.upper - algns[i].xr_diff) - (algns[i].x.lower + algns[i].xl_diff);
      if( cur_len > 0 ) {
        avg_pid = avg_pid + (float)(((float)algns[i].identity)/((float)100)) * ((float) cur_len);
        sum_int = sum_int + cur_len;
      }
    }
    j++;
  }

  if( sum_int > 0 )
  {
    avg_pid = (avg_pid / ((float)(sum_int))) * ((float)100);
  }

  return(avg_pid);
}

void initialize_I(struct I *regs, int b, int e)
{
	int i = 0;

	for( i = b; i <= e; i++ ) {
		regs[i] = assign_I(0, 1);
	}
}

void initialize_int_list(int *list, int b, int e)
{
	int i = 0;

	for( i = b; i <= e; i++ ) {
		list[i] = 0;
	}
}

void init_bool_list(bool *list, int num_list)
{
	int i = 0;

	for( i = 0; i < num_list; i++ ) {
		list[i] = false;
	}
}

bool is_all_full(bool *list, int num_list)
{
	bool res = true;
	int i = 0;

	while( (i < num_list) && (res == true) ) {
		if( list[i] == false ) res = false;
		i++;
	}

	return(res);
}

void init_n_pair(struct n_pair *list, int b, int e)
{
	int i = 0;

	for( i = b; i <= e; i++ ) {
		strcpy(list[i].name1, "");
		strcpy(list[i].name2, "");
		list[i].id = -1;
		list[i].len = 0;
	}
}

void switch_xy_coordinates(struct DotList *algns, int num_algns)
{
	int i = 0, val = 0;
	struct I temp;
	char temp_name[LEN_NAME];

	temp = assign_I(0, 1);
	strcpy(temp_name, "");

	for( i = 0; i < num_algns; i++ ) {
		if( (algns[i].pair_self != PAIR) && (algns[i].x.lower > algns[i].y.lower) ) 
		{
  		temp = assign_I(algns[i].x.lower, algns[i].x.upper);
			algns[i].x = assign_I(algns[i].y.lower, algns[i].y.upper);
			algns[i].y = assign_I(temp.lower, temp.upper);
			val = algns[i].len1;
			algns[i].len1 = algns[i].len2;
			algns[i].len2 = val;
  		strcpy(temp_name, algns[i].name1);
  		strcpy(algns[i].name1, algns[i].name2);
  		strcpy(algns[i].name2, temp_name);
			val = algns[i].ctg_id1;
			algns[i].ctg_id1 = algns[i].ctg_id2;
			algns[i].ctg_id2 = val;
		}
	}
}
