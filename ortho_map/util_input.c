#include "main.h"
#include "util_input.h"
#include "util.h"
#include "util_i.h"
#include "util_gen.h"
#include "read_maf.h"
#include "contigs_op.h"

extern char S[BIG], T[BIG];

void initialize_exons_list(struct exons_list *a, int from, int to)
{
  int i = 0;
  for( i = from; i < to; i++ ) {
    a[i].fid = -1;
    a[i].reg = assign_I(0,1);
    a[i].cmp_reg = assign_I(0,1);
    a[i].sp_id = -1;
    a[i].ctg_id = -1;
    a[i].val = 0;
		a[i].sign = '<';
		strcpy(a[i].name, "");
  }
}

void adjust_pos_conv(struct cv_list *cv, int id, int *len_sum, int num_contigs)
{
	int cur_id = -1;
	
	cur_id = cv[id].ctg_id1;
	if( cur_id != -1 ) {
		if( (cur_id < 0) || (cur_id >= num_contigs) ) fatalf("contigs %d exceeds %d in util_input.c\n", cur_id, num_contigs);
		cv[id].a1 = cv[id].a1 + len_sum[cur_id];
		cv[id].a2 = cv[id].a2 + len_sum[cur_id];
		cv[id].s1 = cv[id].s1 + len_sum[cur_id];
		cv[id].s2 = cv[id].s2 + len_sum[cur_id];
	}

	cur_id = cv[id].ctg_id2;
	if( cur_id != -1 ) {
		if( (cur_id < 0) || (cur_id >= num_contigs) ) fatalf("contigs %d exceeds %d in util_input.c\n", cur_id, num_contigs);
		cv[id].b1 = cv[id].b1 + len_sum[cur_id];
		cv[id].b2 = cv[id].b2 + len_sum[cur_id];
		cv[id].t1 = cv[id].t1 + len_sum[cur_id];
		cv[id].t2 = cv[id].t2 + len_sum[cur_id];
	}
}

void adjust_pos_exons(struct exons_list *exons, int id, int *len_sum, int num_contigs)
{
	int cur_id = -1;
	
	cur_id = exons[id].ctg_id;
	if( cur_id != -1 ) {
		if( (cur_id < 0) || (cur_id >= num_contigs) ) fatalf("contigs %d exceeds %d in util_input.c\n", cur_id, num_contigs);
		exons[id].reg = assign_I(exons[id].reg.lower + len_sum[cur_id], exons[id].reg.upper + len_sum[cur_id]);
		exons[id].cmp_reg = assign_I(exons[id].cmp_reg.lower + len_sum[cur_id], exons[id].cmp_reg.upper + len_sum[cur_id]);
	}
}

int count_exons(char *fname, char *species, char *species2)
{
	char buf[1000];
	char temp_name[50];
	int count = 0;
	FILE *f;
	int temp_code = 0;
	int sp1_code = 1, sp2_code = 2;
	int b = 0, e = 0;

	strcpy(temp_name, "");
	f = ckopen(fname, "r");

	if( strcmp(species, species2) == 0 ) { // codex file includes annotation information for only one species
		temp_code = sp1_code;
	}

  fseek(f, 0, SEEK_SET);
  while( fgets(buf, 1000, f) ) {
		if( ((isspace(buf[0]) != 0) && ((sscanf(buf, "%s %*s", temp_name) !=1))) || (buf[0] == '#') ) {}
    else if( buf[0] == '@' ) {
      sscanf(buf, "%*s %s", temp_name);
			if( strcmp(temp_name, species) == 0 ) {
				temp_code = sp1_code;
			}
			else if( strcmp(temp_name, species2) == 0 ) {
				temp_code = sp2_code;
			}
    }
    else if( (buf[0] == '>') || (buf[0] == '<') ) {}
		else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && (sscanf(buf, "%d %d", &b, &e) == 2)) {
			count++;
		}
	}

	fclose(f);
	return(count);
}

int count_genes(char *fname, char *species, char *species2)
{
	char buf[1000];
	char temp_name[50];
	int count = 0;
	FILE *f;
	int temp_code = 0;
	int sp1_code = 1, sp2_code = 2;

	f = ckopen(fname, "r");

	if( strcmp(species, species2) == 0 ) { // codex file includes annotation information for only one species
		temp_code = sp1_code;
	}
  fseek(f, 0, SEEK_SET);
  while( fgets(buf, 1000, f) ) {
		if( ((isspace(buf[0]) != 0) && (sscanf(buf, "%s %*s", temp_name) !=1)) || (buf[0] == '#') ) {}
    else if( buf[0] == '@' ) {
      sscanf(buf, "%*s %s", temp_name);
			if( strcmp(temp_name, species) == 0 ) {
				temp_code = sp1_code;
			}
			else if( strcmp(temp_name, species2) == 0 ) {
				temp_code = sp2_code;
			}
    }
    else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((buf[0] == '>') || (buf[0] == '<')) ) {
			count++;
		}
	}

	fclose(f);
	return(count);
}

void read_only_exons(struct exons_list *exons, int *num_exons, struct exons_list *genes, int *num_genes, char *fname, char *species, char *species2, struct n_pair *contigs1, struct n_pair *contigs2, int num_contigs1, int num_contigs2)
{
	FILE *f;
	int i = 0, j = 0;
	char buf[1000] = "", temp_name[100] = "";
	int b = 0, e = 0;
	int temp_code = -1;
	int sp1_code = 1, sp2_code = 2;
	char item1[LEN_NAME] = "", item2[LEN_NAME] = "", item3[LEN_NAME] = "", item4[LEN_NAME] = "";
	int ctg_id = -1;

	f = ckopen(fname, "r");

	if( strcmp(species, species2) == 0 ) { // codex file includes annotation information for only one species
		temp_code = sp1_code;
		strcpy(temp_name, species);
	}
  fseek(f, 0, SEEK_SET);
  while( fgets(buf, 1000, f) ) {
		if( ((isspace(buf[0]) != 0) && (sscanf(buf, "%s %*s", temp_name) != 1)) || (buf[0] == '#')) {}
    else if( buf[0] == '@' ) {
      sscanf(buf, "%*s %s", temp_name);
			if( strcmp(temp_name, species) == 0 ) {
				temp_code = sp1_code;
			}
			else if( strcmp(temp_name, species2) == 0 ) {
				temp_code = sp2_code;
			}
			else temp_code = -1;
    }
    else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((buf[0] == '>') || (buf[0] == '<')) ) {
			ctg_id = -1;
      if( sscanf(buf+1, "%s %s %s %s %*s", item1, item2, item3, item4) == 4 ) {
				if( item4[0] == '#' ) {
				}
        else if( temp_code == sp1_code ) {
          ctg_id = get_ctg_id(temp_name, item4, contigs1, num_contigs1);
        }
        else if( temp_code == sp2_code ) {
          ctg_id = get_ctg_id(temp_name, item4, contigs2, num_contigs2);
        }
      }
      else if( sscanf(buf+1, "%s %s %s", item1, item2, item3) != 3 ) {
        fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
      }

      if( (is_all_digits(item1) == true) && (is_all_digits(item2) == true) )
      {
        b = atoi(item1);
        e = atoi(item2);
        if( b >= e ) {
          fatalf("%s includes an incorrect interval in the %s codex file\n", buf, temp_name);
        }
      }
      else {
      	fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
      }

      genes[j].reg = assign_I(b, e+1);
      if( temp_code == sp1_code ) genes[j].sp_id = SELF1;
      else if( temp_code == sp2_code ) genes[j].sp_id = SELF2;
      genes[j].val = 100; // 100% of the gene boundary remains
      genes[j].fid = j;
			genes[j].sign = buf[0];
			genes[j].ctg_id = ctg_id;
			strcpy(genes[j].name, item3);
      j++;
    }
    else if( (temp_code == sp1_code) || (temp_code == sp2_code) ){
      if( sscanf(buf, "%s %s %*s", item1, item2) != 2 ) {
        fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
			}

      if( (is_all_digits(item1) == true) && (is_all_digits(item2) == true) )
      {
        b = atoi(item1);
        e = atoi(item2);
        if( b >= e ) {
          fatalf("%s includes an incorrect interval in the %s codex file\n", buf, temp_name);
        }
      }
      else {
      	fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
      }
      exons[i].reg = assign_I(b, e+1);
      exons[i].val = 100; // 100 % of the exon boundary remains
      exons[i].fid = j-1;
      exons[i].sign = genes[j-1].sign;
			exons[i].ctg_id = ctg_id;
			strcpy(exons[i].name, item3);
      if( temp_code == sp1_code ) exons[i].sp_id = SELF1;
      else if( temp_code == sp2_code ) exons[i].sp_id = SELF2;
      i++;
    }
 	}  
	fclose(f);
	*num_genes = j;
	*num_exons = i;
}

void read_exons(struct exons_list *init_exons, struct exons_list *exons, int *num_exons, struct exons_list *genes, int *num_genes, struct exons_list *skip_reg1, int *num_skip1, struct exons_list *skip_reg2, int *num_skip2, char *fname, struct sp_list *sp_code, int num_sp_code, int sp1_code, int sp2_code, struct n_pair *contigs1, int num_contigs1, int *len_sum1, struct n_pair *contigs2, int num_contigs2, int *len_sum2)
{
	FILE *f;
	int i = 0, j = 0, k = 0;
	char buf[1000], temp_name[100];
	int b = 0, e = 0;
	int temp_code = -1;
	int num_reg1 = 0, num_reg2 = 0;
	char item1[1000], item2[1000], item3[1000], item4[1000];
	int ctg_id = -1;

	num_reg1 = *num_skip1;
	num_reg2 = *num_skip2;

	f = ckopen(fname, "r");

  fseek(f, 0, SEEK_SET);
  while( fgets(buf, 1000, f) ) {
		if( ((isspace(buf[0]) != 0) && ((sscanf(buf, "%s %*s", temp_name) != 1))) || (buf[0] == '#')) {}
    else if( buf[0] == '@' ) {
			ctg_id = -1;
      sscanf(buf, "%*s %s", temp_name);
      for( k = 0; k < num_sp_code; k++ ) {
        if( strcmp(sp_code[k].name, temp_name) == 0 ) temp_code = sp_code[k].id;
      }
    }
    else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((buf[0] == '>') || (buf[0] == '<')) ) {
			ctg_id = -1;
			if( sscanf(buf+1, "%s %s %s %s %*s", item1, item2, item3, item4) == 4 ) {
				if( item4[0] == '#') {
					ctg_id = 0;
				}
				else if( temp_code == sp1_code ) {
					ctg_id = get_ctg_id(temp_name, item4, contigs1, num_contigs1);
				}
				else if( temp_code == sp2_code ) {
					ctg_id = get_ctg_id(temp_name, item4, contigs2, num_contigs2);
				}
				else {
					fatalf("%s: unexpected species name %s in util_input.c\n", sp_code[k].name);
				}
			}
			else if( sscanf(buf+1, "%s %s %s", item1, item2, item3) != 3 ) {
				fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
			}
			else {
				ctg_id = 0; // only one sequence is given and counted as one contig 
			}

			if( (is_all_digits(item1) == true) && (is_all_digits(item2) == true) ) 
			{
				b = atoi(item1);
				e = atoi(item2);
				if( b >= e ) {
					fatalf("%s includes an incorrect interval in the %s codex file\n", buf, temp_name);
				}
			}
			else {
				fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
			}
      genes[j].reg = assign_I(b, e+1);
      if( temp_code == sp1_code ) genes[j].sp_id = SELF1;
      else if( temp_code == sp2_code ) genes[j].sp_id = SELF2;
      genes[j].val = 100; // 100% of the gene boundary remains
      genes[j].fid = j;
			genes[j].ctg_id = ctg_id;
      j++;
    }
    else if( (temp_code == sp1_code) || (temp_code == sp2_code) ){
			if( sscanf(buf, "%s %s", item1, item2) != 2 ) {
				fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
			}
			else {
				if( (is_all_digits(item1) == true) && (is_all_digits(item2) == true) ) 
				{
					b = atoi(item1);
					e = atoi(item2);
					if( temp_code == sp1_code ) {
						if( ( ctg_id < 0 ) || ( ctg_id >= num_contigs1 ) ) {
							fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
						}
						else {
							b = b + len_sum1[ctg_id];
							e = e + len_sum1[ctg_id];
						}
					}
					else if( temp_code == sp2_code ) {
						if( ( ctg_id < 0 ) || ( ctg_id >= num_contigs2 ) ) {
							fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
						}
						else {
							b = b + len_sum2[ctg_id];
							e = e + len_sum2[ctg_id];
						}
					}
					else {

					}

					if( b >= e ) {
						fatalf("%s includes an incorrect interval in the %s codex file\n", buf, temp_name);
					}
				}
				else {
					fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
				}
			}
      exons[i].reg = assign_I(b, e);
      exons[i].val = 100; // 100 % of the exon boundary remains
      exons[i].fid = i;
			exons[i].ctg_id = ctg_id;
      init_exons[i].reg = assign_I(b, e);
      if( temp_code == sp1_code ) exons[i].sp_id = SELF1;
      else if( temp_code == sp2_code ) exons[i].sp_id = SELF2;
      init_exons[i].sp_id = exons[i].sp_id;
      if( init_exons[i].sp_id == SELF1 ) {
        skip_reg1[num_reg1].reg = assign_I(init_exons[i].reg.lower, init_exons[i].reg.upper);
        skip_reg1[num_reg1].val = 100; // exons are always ignored
        num_reg1++;
      }
      else if( init_exons[i].sp_id == SELF2 ) {
        skip_reg2[num_reg2].reg = assign_I(init_exons[i].reg.lower, init_exons[i].reg.upper); 
				skip_reg2[num_reg2].val = 100;
        num_reg2++;
      }
      init_exons[i].fid = i;
      exons[i].fid = i;
      i++;
    }
 	}  
	fclose(f);
	*num_genes = j;
	*num_exons = i;
	*num_skip1 = num_reg1;
	*num_skip2 = num_reg2;
}

int count_local_algns(char *fname, char *species, char *species2)
{
	int count = 0;
	FILE *f = NULL;
	int b1 = 0, e1 = 0, b2 = 0, e2 = 0, temp = 0;
	char len1[LEN_NAME] = "", len2[LEN_NAME] = "", strand[LEN_NAME] = "";
	char name1[LEN_NAME] = "", name2[LEN_NAME] = "";
	bool is_read = false;

	if((f = fopen(fname, "r")) == NULL) {
		fatalf("%s not exist\n", fname);
	}
	else {
    while(fgets(S, BIG, f)) {
      if( S[0] == '#' ) {
        while( (S[0] == '#') && fgets(S, BIG, f) ) {
          if( strncmp(S, "##maf", 5) == 0 ) {
            count = 0;
          }
        }
        count = 0;
      }

      if( S[0] == 'a' ) {
				strcpy(strand, "+");
				strcpy(len1, "0");
				strcpy(len2, "0");
				strcpy(name1, "");
				strcpy(name2, "");
        if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))
          fatalf("cannot find alignment in %s", fname);
        if( (sscanf(S, "%*s %s %d %d %*s %s", name1, &b1, &e1, len1) != 4) || (sscanf(T, "%*s %s %d %d %s %s", name2, &b2, &e2, strand, len2) != 5)) {
			   fatalf("bad alignment info of 2 in %s", fname);
				}
	 		  e1 += b1;
 		  	e2 += b2;

  	  	if( strcmp(strand, "-") == 0) {
   	   		temp = b2;
   	   		b2 = atoi(len2) - e2;
   	   		e2 = atoi(len2) - temp;
   	 		}

   	 		b1++;
   			b2++;
   			e1++;
   			e2++;

		 		if( (abs(e1-b1) <= ERR_SM_TH) || (abs(e2-b2) <= ERR_SM_TH) ) {}
		 		else {
     			count++;
				}
			}
    }
	}

	if( count > 0 ) {
		concat_ctg_name(name1, species, len1);	
		concat_ctg_name(name2, species2, len1);	
	}
	else {
		fseek(f, 0, SEEK_SET);
    while((is_read == false) && fgets(S, BIG, f)) {
      if( S[0] == '#' ) {
         if( strncmp(S, "##maf", 5) == 0 ) {
					is_read = true;
        	if( sscanf(S, "%*s %*s %*s %s %s", species, species2) != 2) {
			   		fatalf("species info expected in %s (util_input.c)\n", S);
					}
				}
			}
		}
	}

  fclose(f);
	return(count);
}

int count_lines(FILE *f)
{
	char buf[1000] = "";
	int count = 0;

	fseek(f, 0, SEEK_SET);
	while(fgets(buf, 1000, f)) count++;	
	return(count);
}

int count_tokens(char *line)
{
	int len = 0, i = 0;
	int num_tokens = 0;

	len = strlen(line);

	while((line[i] != '\0') && (line[i] != '\n'))
	{
		if( isspace(line[i]) != 0 ) {
			num_tokens++;
			while((line[i] != '\0') && (line[i] != '\n') && (isspace(line[i]) != 0)) i++;
		}
		else i++;	
	}
	num_tokens++;
	return(num_tokens);
}

int concat_tokens(char *line, int loc, char *temp_name)
{
	int len = 0, i = 0, j = 0;

	len = strlen(line);
	i = loc;
	if( isspace(line[i]) != 0 ) {
		while((line[i] != '\0') && (line[i] != '\n') && (isspace(line[i]) != 0)) i++;
	}

	while((line[i] != '\0') && (line[i] != '\n') && (isspace(line[i]) == 0)) 
	{
		temp_name[j] = line[i];
		j++;
		if( j >= LEN_NAME ) {
			fatalf("Over the max length for a gene name ( 100 characters ) in %s\n", line);
		}
		i++;
	}
	temp_name[j] = '\0';

	return(i);
}

void numtostr(int c, char *str)
{
  int i = 0;
  int digit = 0;
  char id[10];

  while( c/10 >= 1 )
  {
    id[digit] = (c % 10) + '0';
    digit++;
    c = c/10;
  }
  id[digit] = c + '0';

  for( i = 0; i <= digit; i++ )
  {
    str[i] = id[digit-i];
  }
  str[digit+1] = '\0';
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
