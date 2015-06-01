#include "main.h"
#include "util.h"
#include "util_input.h"
#include "util_gen.h"
#include "contigs_op.h"

void merge_sep_contigs(struct DotList *algns, int num_algns, struct n_pair **contigs1, int *num_contigs1, int *num_alloc1, int **len_sum1, struct n_pair **contigs2, int *num_contigs2, int *num_alloc2, int **len_sum2)
{
	int i = 0, j = 0, id = -1;
	char *name, *sp_name, *ctg_name;
	int num1 = 0, num2 = 0;

	num1 = *num_contigs1;
	num2 = *num_contigs2;

	name = (char *) ckalloc(LEN_NAME * sizeof(char));
	sp_name = (char *) ckalloc(LEN_NAME * sizeof(char));
	ctg_name = (char *) ckalloc(LEN_NAME * sizeof(char));
	for( i = 0; i < num_algns; i++ ) {
		strcpy(name, "");
		strcpy(sp_name, "");
		strcpy(ctg_name, "");
		strcpy(name, algns[i].name1);
		concat_ctg_name(name, sp_name, ctg_name);
		id = add_ctg_name(algns[i].index, sp_name, ctg_name, algns[i].len1, *contigs1, num1);
		
		algns[i].ctg_id1 = id;
		if( id == num1 ) {
			num1++;
			(*len_sum1)[num1-1] = (*len_sum1)[num1-2] + (*contigs1)[num1-2].len;
		}

		if( num1 >=  (*num_alloc1) ) {
			*num_alloc1 = *num_alloc1 + ALLOC_UNIT;
			*contigs1 = (struct n_pair *) ckrealloc(*contigs1, sizeof(struct n_pair) * (*num_alloc1));
			*len_sum1 = (int *) ckrealloc(*len_sum1, (*num_alloc1) * sizeof(int));	
			for( j = ((*num_alloc1)-ALLOC_UNIT); j < (*num_alloc1); j++ ) {
				(*len_sum1)[j] = 0;
			}
			init_n_pair(*contigs1, (*num_alloc1) - ALLOC_UNIT, (*num_alloc1)-1);
		}

		strcpy(name, algns[i].name2);
		concat_ctg_name(name, sp_name, ctg_name);
		id = add_ctg_name(algns[i].index, sp_name, ctg_name, algns[i].len2, *contigs2, num2);
		
		algns[i].ctg_id2 = id;
		if( id == num2 ) {
			num2++;
			(*len_sum2)[num2-1] = (*len_sum2)[num2-2] + (*contigs2)[num2-2].len;
		}

		if( num2 >= (*num_alloc2) ) {
			*num_alloc2 = *num_alloc2 + ALLOC_UNIT;
			*len_sum2 = (int *) ckrealloc(*len_sum2, (*num_alloc2) * sizeof(int));	
			for( j = ((*num_alloc2)-ALLOC_UNIT); j < (*num_alloc2); j++ ) {
				(*len_sum2)[j] = 0;
			}
			*contigs2 = (struct n_pair *) ckrealloc(*contigs2, sizeof(struct n_pair) * (*num_alloc2));
			init_n_pair(*contigs2, (*num_alloc2) - ALLOC_UNIT, (*num_alloc2)-1);
		}
	}

	*num_contigs1 = num1;
	*num_contigs2 = num2;
	free(name);
	free(ctg_name);
	free(sp_name);
}

int add_ctg_name(int id, char *sp_name, char *ctg_name, int len, struct n_pair *pairs, int num_pairs)
{
	int num = 0;
	int ctg_id = -1;

	num = num_pairs;

	ctg_id = is_ctg_in(sp_name, ctg_name, pairs, num_pairs);

	if( ctg_id == -1 ) { // if it's a new contig
		strcpy(pairs[num].name1, sp_name);
		strcpy(pairs[num].name2, ctg_name);
		pairs[num].len = len;
		pairs[num].id = id;	
		ctg_id = num;
	}

	if( ctg_id == -1 ) {
		fatalf("contig id not assigned yet: %s.%s %d\n", sp_name, ctg_name, len); 
	}
	return(ctg_id);
}

int is_ctg_in(char *sp_name, char *ctg_name, struct n_pair *pairs, int num_pairs)
{
	bool is_in = false;
	int i = 0;
	int res = -1;

	while( (is_in == false) && (i < num_pairs) ) {
		if( (strcmp(sp_name, pairs[i].name1) == 0) && (strcmp(ctg_name, pairs[i].name2) == 0 ) ) {
			is_in = true;
			res = i;
		}
		else if( strcmp(sp_name, pairs[i].name1) != 0 ) {
			fatalf("sequence names not match: %s.%s vs. %s.%s in contigs_op.c\n", sp_name, ctg_name, pairs[i].name1, pairs[i].name2);
		}
		i++;
	}

	return(res);
}

void cal_length_sum_from_ctglist(int *len_sum, struct n_pair *contigs, int num)
{
	int i = 0;

	for( i = 0; i < num; i++ ) {
		len_sum[i] = contigs[i].len;
	}
}

void cal_length_sum(int *len_sum1, struct n_pair *contigs1, int num1)
{
  int i = 0, j = 0;

  for( i = 0; i < num1; i++ ) len_sum1[i] = 0;

  for( i = 0; i < num1; i++ ) {
    for( j = (i+1); j < num1; j++ ) {
      len_sum1[j] = len_sum1[j] + contigs1[i].len;
    }
  }
}

void concat_ctg_name(char *name, char *sp_name, char *ctg_name)
{
	char *token = NULL;

	if( strstr(name, ".") != NULL ) {
		token = strtok(name, ".");	
		strcpy(sp_name, token);
		token = strtok(NULL, ".");	
		strcpy(ctg_name, token);
	}
	else {
		strcpy(sp_name, name);
		strcpy(ctg_name, name);
	}
}

int count_contigs(char *name, FILE *f)
{
  char buf[1000];
  bool is_in = false;
  char temp_name[LEN_NAME];
  int count = 0;
  int res = 0;

	strcpy(temp_name, "");

  fseek(f, 0, SEEK_SET);
  while(fgets(buf, 500, f)) {
    if( buf[0] == '#' ) {
      if( is_in == true ) {
        res = count;
      }
      is_in = false;
      sscanf(buf+1, "%s %*s", temp_name);
      if( strcmp(temp_name, name) == 0 ) {
        is_in = true;
      }
      count = 0;
    }
    else {
      if( is_in == true ) {
        count++;
      }
    }
  }

  if( is_in == true ) {
    res = count;
  }

  return(res);
}

void read_contigs_file(char *name, FILE *f, struct n_pair *contigs, int num_contigs)
{
  char buf[1000];
  bool is_in = false;
  char temp_name[LEN_NAME], sp_name[LEN_NAME], ctg_name[LEN_NAME], len_str[LEN_NAME], len_sum[LEN_NAME];
  int count = 0;
  int res = 0;

	strcpy(temp_name, "");
	strcpy(sp_name, "");
	strcpy(ctg_name, "");
	strcpy(len_str, "");

  fseek(f, 0, SEEK_SET);
  while(fgets(buf, 500, f)) {
    if( buf[0] == '#' ) {
      if( is_in == true ) {
        res = count;
      }
      is_in = false;
      sscanf(buf+1, "%s %*s", temp_name);
      if( strcmp(temp_name, name) == 0 ) {
        is_in = true;
      }
      count = 0;
    }
    else {
      if( is_in == true ) {
				if( sscanf(buf, "%s %s %s %s", sp_name, ctg_name, len_str, len_sum) != 4 ) {
					fatalf("unexpected line %s\n", buf);
				}
				else {
					if( count >= num_contigs ) {
						fatalf("contig counting error: %d exceeds %d\n", count, num_contigs);
					}
					else {
						strcpy( contigs[count].name1, sp_name);
						strcpy( contigs[count].name2, ctg_name);
						contigs[count].id = atoi(len_str);
						contigs[count].len = atoi(len_sum);
					}
        	count++;
				}
      }
    }
  }
}

int get_ctg_id(char *sp_name, char *ctg_name, struct n_pair *pairs, int num)
{
  bool is_in = false;
  int i = 0;
  int res = -1;

  while( (is_in == false) && (i < num) ) {
    if( (strcmp(sp_name, pairs[i].name1) == 0) && (strcmp(ctg_name, pairs[i].name2) == 0 ) ) {
      is_in = true;
      res = i;
    }
    else if( strcmp(sp_name, pairs[i].name1) != 0 ) {
      fatalf("sequence names not match: %s.%s vs. %s.%s in contigs_op.c\n", sp_name, ctg_name, pairs[i].name1, pairs[i].name2);
    }
    i++;
  }

	if( (res < 0) || (res >= num) ) {
    fatalf("contig %s.%s not exists (contigs_op.c)\n", sp_name, ctg_name);
	}

  return(res);
}

void print_contig_list(FILE *ctg_f, struct n_pair *contigs, int *len_sum, int num_org, int num)
{
  int i = 0;

  for( i = num_org; i < num; i++ )
  {
    fprintf(ctg_f, "%s %s %d %d\n", contigs[i].name1, contigs[i].name2, contigs[i].len, len_sum[i]);
  }
}

int read_one_contig_file(char *fname, FILE *fp, struct n_pair *contigs, int count)
{
	int i = 0, len1 = 0, len2 = 0;
	int num_contigs = 0;
	char buf[1000];
	char name1[LEN_NAME], name2[LEN_NAME];
	
	i = 0;
  fseek(fp, 0, SEEK_SET);
  while( fgets(buf, 1000, fp) ) {
   if( i >= count ) {
     fatalf("counting error in %s\n", fname);
   }
   sscanf(buf, "%s %s %d %d", name1, name2, &len1, &len2);
   strcpy(contigs[i].name1, name1);
   strcpy(contigs[i].name2, name2);
   contigs[i].len = len2;
   contigs[i].id = len1;
   i++;
  }

	num_contigs = i;
	return(num_contigs);
}

void replace_contigs_len(struct n_pair *contigs, int num)
{
	int i = 0;

	for( i = 0; i < num; i++ ) contigs[i].len = contigs[i].id;
}
