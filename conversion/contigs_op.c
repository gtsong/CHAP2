#include "util.h"
#include "contigs_op.h"

int add_ctg_name(int id, char *sp_name, char *ctg_name, int len, struct n_pair *pairs, int num_pairs)
{
	int num = 0;
	int ctg_id = -1;

	num = num_pairs;

	if( num_pairs > 0 ) {
		ctg_id = is_ctg_in(sp_name, ctg_name, pairs, num_pairs);
	}

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
	char *token;

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

void init_pair(struct n_pair *list, int b, int e)
{
  int i = 0;

  for( i = b; i <= e; i++ ) {
    strcpy(list[i].name1, "");
    strcpy(list[i].name2, "");
    list[i].id = -1;
    list[i].len = 0;
  }
}
