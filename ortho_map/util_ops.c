#include "main.h"
#include "regions.h"
#include "util_gen.h"
#include "util_ops.h"
#include "util.h"
#include "util_i.h"
#include "util_input.h"
#include "contigs_op.h"

extern int debug_mode;

void sort_conv(struct cv_list *a, int num)
{
	int i, j;
	quick_sort_inc_conv(a, 0, num-1, POS_BASE);
	i = 0;
	while( i < num )
	{
		j = 0;
		while( ((i+j) < num) && (a[i].a1==a[i+j].a1) ) j++;
		quick_sort_dec_conv(a, i, i+j-1, LEN_BASE);
		i = i+j;
	}
}

void quick_sort_dec_conv(struct cv_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct cv_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].a1));
	else if( mode == PID_BASE ) x = a[(lo+hi)/2].conv_pid;
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].a2 - a[(lo+hi)/2].a1);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while ((i <= hi) && (((float)(a[i].a1))>x)) i++; 
			while ((j >= lo) && (((float)(a[j].a1))<x)) j--;
		}
		else if( mode == PID_BASE ) {
			while ((i <= hi) && (a[i].conv_pid>x)) i++; 
			while ((j >= lo) && (a[j].conv_pid<x)) j--;
		}
		else if( mode == LEN_BASE ) {
			while ((i <= hi) && (abs(a[i].a2-a[i].a1)>x)) i++; 
			while ((j >= lo) && (abs(a[j].a2-a[j].a1)<x)) j--;
		}

		if (i<=j)
		{
			h = assign_conv(a[i]);
			a[i] = assign_conv(a[j]);
			a[j] = assign_conv(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_conv(a, lo, j, mode);
	if (i < hi) quick_sort_dec_conv(a, i, hi, mode);
}

void quick_sort_inc_conv(struct cv_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct cv_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].a1));
	else if( mode == PID_BASE ) x = a[(lo+hi)/2].conv_pid;
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].a2 - a[(lo+hi)/2].a1);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while ((i <= hi) && (((float)(a[i].a1))<x)) i++; 
			while ((j >= lo) && (((float)(a[j].a1))>x)) j--;
		}
		else if( mode == PID_BASE ) {
			while ((i <= hi) && (a[i].conv_pid<x)) i++; 
			while ((j >= lo) && (a[j].conv_pid>x)) j--;
		}
		else if( mode == LEN_BASE ) {
			while ((i <= hi) && (abs(a[i].a2-a[i].a1)<x)) i++; 
			while ((j >= lo) && (abs(a[j].a2-a[j].a1)>x)) j--;
		}

		if (i<=j)
		{
			h = assign_conv(a[i]);
			a[i] = assign_conv(a[j]);
			a[j] = assign_conv(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_conv(a, lo, j, mode);
	if (i < hi) quick_sort_inc_conv(a, i, hi, mode);
}

int quick_search_close_conv(struct cv_list *sorted, int i, int j, int query)
{
	int mid;
	int val;
	int res;

	mid = (i+j)/2;
	val = sorted[mid].a1;

	if(val > query) {
		if( j <= (mid+1) ) return (mid+1);	
		else res = quick_search_close_conv(sorted, mid+1, j, query);
	} 
	else if(val < query) {
		if( i >= (mid-1) ) return i;
		else res = quick_search_close_conv(sorted, i, mid-1, query);
	}
	else return mid;

	return(res);
}

struct cv_list assign_conv(struct cv_list a)
{
  struct cv_list res;

  res.fid = a.fid;
  res.s1 = a.s1;
  res.s2 = a.s2;
  res.t1 = a.t1;
  res.t2 = a.t2;
  res.a1 = a.a1;
  res.a2 = a.a2;
  res.b1 = a.b1;
  res.b2 = a.b2;
  res.ori = a.ori;
  res.algn_pid = a.algn_pid;
  res.conv_pid = a.conv_pid;
  res.sp_id = a.sp_id;
  res.dir = a.dir;
	res.ctg_id1 = a.ctg_id1;
	res.ctg_id2 = a.ctg_id2;
	res.sp_code = a.sp_code;
	res.out_code = a.out_code;

  return(res);
}

void init_conv(struct cv_list *cv, int b, int e)
{
	int i = 0;

	for( i = b; i <= e; i++ ) {
		cv[i].fid = -1;
 		cv[i].s1 = 0;
 		cv[i].s2 = 0;
 		cv[i].t1 = 0;
 		cv[i].t2 = 0;
 		cv[i].a1 = 0;
 		cv[i].a2 = 0;
 		cv[i].b1 = 0;
 		cv[i].b2 = 0;
 		cv[i].ori = 'd';
 		cv[i].algn_pid = (float)0;
 		cv[i].conv_pid = (float)0;
 		cv[i].sp_id = 0;
 		cv[i].dir = 0;
 		cv[i].ctg_id1 = -1;
 		cv[i].ctg_id2 = -1;
 		cv[i].sp_code = -1;
 		cv[i].out_code = -1;
	}
}

struct exons_list assign_exons(struct exons_list a)
{
  struct exons_list res;

  res.fid = a.fid; // id in the inital list
  res.reg = assign_I(a.reg.lower, a.reg.upper);
  res.cmp_reg = assign_I(a.cmp_reg.lower, a.cmp_reg.upper);
  res.sp_id = a.sp_id;
  res.val = a.val;
  res.sign = a.sign; // '<' or '>'
  res.ctg_id = a.ctg_id;

  return(res);
}

struct ops_list assign_ops(struct ops_list a)
{
  struct ops_list res;

  res.sign = a.sign;
  res.id = a.id;
  res.srcStart = a.srcStart;
  res.srcEnd = a.srcEnd;
  res.dstStart = a.dstStart;
  res.dstEnd = a.dstEnd;
  res.src_b = a.src_b;
  res.src_e = a.src_e;
  res.dst_b = a.dst_b;
  res.dst_e = a.dst_e;
  res.pid = a.pid;
  res.sp_id = a.sp_id;
  res.ctg_id1 = a.ctg_id1;
  res.ctg_id2 = a.ctg_id2;
  res.dir = a.dir;

  return(res);
}

int count_ops(char *name, FILE *f, int sp_mode)
{
	char buf[1000];
	bool is_in = false;
	char temp_name[LEN_NAME];
	int count = 0;
	int res = 0;
	int mode = -1;

	fseek(f, 0, SEEK_SET);
	while(fgets(buf, 500, f)) {
		if( buf[0] == '#' ) {
			if( is_in == true ) {
				res = count;
			}
			is_in = false;
			sscanf(buf+2, "%*s %s", temp_name); 
			if( strcmp(temp_name, name) == 0 ) {
				is_in = true;	
			}
			count = 0;
		}
		else {
			if( is_in == true ) {
				sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %d %*s", &mode); 
				if( mode == sp_mode ) {
					count++;
				}
			}
		}
	}	

	if( is_in == true ) {
		res = count;
	}

	return(res);
}

int count_final_ops(char *name, FILE *f, int sp_mode)
{
	char buf[1000];
	bool is_in = false;
	char temp_name[LEN_NAME];
	int count = 0;
	int res = 0;
	int mode = -1;
	int i = 0;
	int len = 0;

	fseek(f, 0, SEEK_SET);
	while((is_in == false) && (fgets(buf, 500, f))) {
		if( buf[0] == '#' ) {
			i = 2;
			len = strlen(buf);
			while( (i < len) && (buf[i] != '\0') && (buf[i] != '\n') ) {
				i = concat_tokens(buf, i, temp_name); 
				if( strcmp(temp_name, name) == 0 ) {
					is_in = true;
				}
			}
		}
		else {
			if( is_in == false ) {
				sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %d %*s", &mode);
				if( mode == sp_mode ) {
					count++;
				}
			}
		}
	}	

	res = count;
	return(res);
}

void read_final_ops_file(char *name, FILE *f, struct ops_list *ops, int num_ops, int sp_mode, struct n_pair *contigs, int num_contigs)
{
	char buf[1000];
	bool is_in = false;
	char temp_name[LEN_NAME];
	int count = 0;
	int i = 0;
	int len = 0;
	int *len_sum;

	if( num_contigs > 0 ) {
		len_sum = (int *) ckalloc(sizeof(int) * num_contigs);
//		cal_length_sum(len_sum, contigs, num_contigs);
		for( i = 0; i < num_contigs; i++ ) {
			len_sum[i] = contigs[i].len;
		}
	}

	fseek(f, 0, SEEK_SET);
	while((is_in == false) && (fgets(buf, 500, f))) {
		if( buf[0] == '#' ) {
			i = 2;
			len = strlen(buf);
			while( (i < len) && (buf[i] != '\0') && (buf[i] != '\n') ) {
				i = concat_tokens(buf, i, temp_name); 
				if( strcmp(temp_name, name) == 0 ) {
					is_in = true;
				}
			}
		}
		else {
			if( is_in == false ) {
				count = add_one_ops(ops, count, buf, contigs, num_contigs, len_sum, sp_mode);
			}
			if( count > num_ops ) {
				fatalf("counting events mismatches: %d in %s\n", num_ops, name);
			}
		}
	}
	if( num_contigs > 0 ) free(len_sum);
}

int add_one_ops(struct ops_list *ops, int id, char *buf, struct n_pair *contigs, int num_contigs, int *len_sum, int sp_mode)
{
	int ctg_id = -1;
	char sp_name[LEN_NAME], ctg_name[LEN_NAME];
	char name1[LEN_NAME], name2[LEN_NAME];
	int count = id;
	int mode = -1;

	sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %d %*s", &mode);
	if( mode == sp_mode ) {
		sscanf(buf, "%c %s %d %d %s %d %d %d %d %f", &ops[count].sign, name1, &ops[count].srcStart, &ops[count].srcEnd, name2, &ops[count].dstStart, &ops[count].dstEnd, &ops[count].sp_id, &ops[count].dir, &ops[count].pid);
					
		ops[count].ctg_id1 = -1;
		ops[count].ctg_id2 = -1;
		if( num_contigs > 0 ) {
			concat_ctg_name(name1, sp_name, ctg_name);
			ctg_id = get_ctg_id(sp_name, ctg_name, contigs, num_contigs);

			if( ctg_id >= num_contigs ) {
				fatalf("contig id %d exceeds %d in util_ops.c\n", ctg_id, num_contigs);
			}
			else if( ctg_id != -1 ) {
				ops[count].srcStart = ops[count].srcStart + len_sum[ctg_id];
				ops[count].srcEnd = ops[count].srcEnd + len_sum[ctg_id];
				ops[count].ctg_id1 = ctg_id;
			}

			if( ops[count].sign != 'd' ) {
				concat_ctg_name(name2, sp_name, ctg_name);
				ctg_id = get_ctg_id(sp_name, ctg_name, contigs, num_contigs);
			}
			else {
				ctg_id = -1;
			}

			if( (ops[count].sign != 'd') && ((ctg_id >= num_contigs) || (ctg_id < 0) ) ) {
				fatalf("contig id %d exceeds %d in util_ops.c\n", ctg_id, num_contigs);
			}
			else if( ctg_id != -1 ) {
				ops[count].dstStart = ops[count].dstStart + len_sum[ctg_id];
				ops[count].dstEnd = ops[count].dstEnd + len_sum[ctg_id];
				ops[count].ctg_id2 = ctg_id;
			}
			else {
				ops[count].ctg_id2 = ctg_id;
			}
		}
		count = count + 1;
	}
	return(count);
}

void read_ops_file(char *name, FILE *f, struct ops_list *ops, int num_ops, int sp_mode, struct n_pair *contigs, int num_contigs, char *ref_name)
{
	char buf[1000];
	bool is_in = false;
	char temp_name1[LEN_NAME], temp_name2[LEN_NAME];
	int *len_sum;
	int count = 0;

	if( num_contigs > 0 ) {
		len_sum = (int *) ckalloc(sizeof(int) * num_contigs);
		cal_length_sum(len_sum, contigs, num_contigs);
	}

	strcpy(temp_name1, "");
	strcpy(temp_name2, "");
	fseek(f, 0, SEEK_SET);
	while(fgets(buf, 500, f)) {
		if( buf[0] == '#' ) {
			is_in = false;
			sscanf(buf, "%*s %s %s", temp_name1, temp_name2); 
			if( (strcmp(temp_name1, ref_name) == 0) && (strcmp(temp_name2, name) == 0 )) {
				is_in = true;	
			}
			count = 0;
		}
		else {
			if( is_in == true ) {
				if( sp_mode == SELF1 ) {
					count = add_one_ops(ops, count, buf, contigs, num_contigs, len_sum, sp_mode);
				}
				else if( sp_mode == SELF2 ) {
					count = add_one_ops(ops, count, buf, contigs, num_contigs, len_sum, sp_mode);
				}
				else {
					fatalf("%d: mode must be %d or %d in util_ops.c\n", sp_mode, SELF1, SELF2);
				}
			}

			if( count > num_ops ) {
					fatalf("counting events mismatches: %d in %s\n", num_ops, name);
			}
		}
	}	

	if( num_contigs > 0 ) free(len_sum);
}

void init_ops(struct ops_list *ops, int b, int e)
{
	int i = 0;

	for( i = b; i < e; i++ ) {
 		ops[i].sign = 'n';
	  ops[i].id = -1;
 		ops[i].srcStart = 0;
	 	ops[i].srcEnd = 1;
	 	ops[i].dstStart = 0;
	 	ops[i].dstEnd = 1;
	 	ops[i].src_b = 0;
	 	ops[i].src_e = 1;
	 	ops[i].dst_b = 0;
	 	ops[i].dst_e = 1;
	 	ops[i].pid = 0;
	 	ops[i].sp_id = 0;
	 	ops[i].ctg_id1 = -1;
	 	ops[i].ctg_id2 = -1;
	 	ops[i].dir = 0;
	}
}

bool is_on_prev_events(struct I reg, struct ops_list *ops, int from, int to)
{
	int i = 0;
	struct I src, dst;
	bool res = false;

	i = from;
	while( (i <= to) && (res == false) ) {
		if( (ops[i].sign == '+') || (ops[i].sign == '-') ) {
			src = assign_I(ops[i].srcStart, ops[i].srcEnd);
			dst = assign_I(ops[i].dstStart, ops[i].dstEnd);
			if( (f_loose_subset(reg, src, STRICT) == true ) || (f_loose_subset(reg, dst, STRICT) == true) )
			{
				res = true;
			}
		}
		i++;
	}

	return(res);
}
