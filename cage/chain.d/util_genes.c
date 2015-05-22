#include "main.h"
#include "util_genes.h"
#include "util.h"
#include "util_i.h"

extern int debug_mode;

void quick_sort_dec_genes(struct g_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct g_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].txStart));
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].txEnd - a[(lo+hi)/2].txStart);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while (((float)(a[i].txStart))>x) i++; 
			while (((float)(a[j].txStart))<x) j--;
		}
		else if( mode == LEN_BASE ) {
			while (abs(a[i].txEnd-a[i].txStart)>x) i++; 
			while (abs(a[j].txEnd-a[j].txStart)<x) j--;
		}

		if (i<=j)
		{
			h = assign_genes(a[i]);
			a[i] = assign_genes(a[j]);
			a[j] = assign_genes(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_dec_genes(a, lo, j, mode);
	if (i < hi) quick_sort_dec_genes(a, i, hi, mode);
}

void quick_sort_inc_genes(struct g_list *a, int lo, int hi, int mode)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi;
	struct g_list h;
	float x;
	
	if( mode == POS_BASE ) x = ((float) (a[(lo+hi)/2].txStart));
	else if( mode == LEN_BASE ) x = abs(a[(lo+hi)/2].txEnd - a[(lo+hi)/2].txStart);

//  partition
	do
	{    
		if( mode == POS_BASE ) {
			while (((float)(a[i].txStart))<x) i++; 
			while (((float)(a[j].txStart))>x) j--;
		}
		else if( mode == LEN_BASE ) {
			while (abs(a[i].txEnd-a[i].txStart)<x) i++; 
			while (abs(a[j].txEnd-a[j].txStart)>x) j--;
		}

		if (i<=j)
		{
			h = assign_genes(a[i]);
			a[i] = assign_genes(a[j]);
			a[j] = assign_genes(h);
			i++; 
			j--;
		}
	} while (i <= j);
	if (lo < j) quick_sort_inc_genes(a, lo, j, mode);
	if (i < hi) quick_sort_inc_genes(a, i, hi, mode);
}

int quick_search_close_genes(struct g_list *sorted, int i, int j, int query)
{
	int mid;
	int val;
	int res;

	mid = (i+j)/2;
	val = sorted[mid].txStart;

	if(val > query) {
		if( j <= (mid+1) ) return (mid+1);	
		else res = quick_search_close_genes(sorted, mid+1, j, query);
	} 
	else if(val < query) {
		if( i >= (mid-1) ) return i;
		else res = quick_search_close_genes(sorted, i, mid-1, query);
	}
	else return mid;

	return(res);
}

struct g_list assign_genes(struct g_list a)
{
  struct g_list res;

  res.gid = a.gid;
  res.strand = a.strand;
  res.txStart = a.txStart;
  res.txEnd = a.txEnd;
	strcpy(res.gname, a.gname);

  return(res);
}

int read_genes(FILE *f, struct g_list *genes, char *fname)
{
  int num_genes = 0;
  char buf[1000] = "";
  int i = 0;
	int b = 0, e = 0;
	char name[LEN_NAME] = "";
	char item1[1000] = "", item2[1000] = "", item3[1000] = "", item4[1000] = "";

  i = 0;
  while(fgets(buf, 1000, f))
  {
    if( buf[0] == '#' ) {}
		else if( (buf[0] == '>') || (buf[0] == '<')) {
	    if(buf[0] == '>') genes[i].strand = '+';
	    else if(buf[0] == '<' ) genes[i].strand = '-';

			if( sscanf(buf+1, "%s %s %s %s %*s", item1, item2, item3, item4) == 4 ) {
				strcpy(name, item3);
				if( item4[0] != '#' ) {
					strcpy(genes[i].sp_name, item4);
				}
			}
      else if( sscanf(buf+1, "%s %s %s", item1, item2, item3) != 3 ) {
        fatalf("%s: unsupported format in the %s file\n", buf, fname);
      }

      if( (is_all_digits(item1) == true) && (is_all_digits(item2) == true) )
      {
        b = atoi(item1);
        e = atoi(item2);
        if( b >= e ) {
          fatalf("%s includes an incorrect interval in the %s file\n", buf, fname);
        }
      }
      else {
        fatalf("%s: unsupported format in the %s codex file\n", buf, fname);
      }
			strcpy(genes[i].sp_name, "");

//			sscanf(buf, "%*s %d %d %s %*s", &b, &e, name);
	    genes[i].gid = i;
	   	genes[i].txStart = b;
 	  	genes[i].txEnd = e;
  	  strcpy(genes[i].gname, name);
    	i++;
		}
		else {
      if( sscanf(buf, "%s %s", item1, item2) != 2 ) {
        fatalf("%s: unsupported format in the %s file\n", buf, fname);
      }
      else {
        if( (is_all_digits(item1) == true) && (is_all_digits(item2) == true) )
        {
          b = atoi(item1);
          e = atoi(item2);
          if( b >= e ) {
            fatalf("%s includes an incorrect interval in the %s codex file\n", buf, fname);
          }
        }
        else {
          fatalf("%s: unsupported format in the %s codex file\n", buf, fname);
        }
      }
		}
  }
	num_genes = i;

	return(num_genes);
}

bool is_all_digits(char *item)
{
  int len = strlen(item);
  int i = 0;
  bool res = true;

  while( isspace(item[len-1]) != 0 ) {
    len--;
  }

  for(i = 0; i < len; i++) {
    if( isdigit(item[i]) == 0 ) res = false;
  }

  return(res);
}

