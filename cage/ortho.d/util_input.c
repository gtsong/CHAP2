#include "main.h"
#include "util_input.h"
#include "util.h"
#include "util_i.h"
#include "util_gen.h"

extern char S[BIG], T[BIG];

void initialize_exons_list(struct exons_list *a, int from, int to)
{
  int i = 0;
  for( i = from; i < to; i++ ) {
    a[i].fid = -1;
    a[i].reg = assign_I(0,1);
    a[i].cmp_reg = assign_I(0,1);
    a[i].sp_id = -1;
    a[i].val = 0;
		a[i].sign = '<';
		strcpy(a[i].name, "");
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
		if( isspace(buf[0]) != 0 ) {}
    else if( buf[0] == '#' ) {
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
		if( isspace(buf[0]) != 0 ) {}
    else if( buf[0] == '#' ) {
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

void read_only_exons(struct exons_list *exons, int *num_exons, struct exons_list *genes, int *num_genes, char *fname, char *species, char *species2)
{
	FILE *f;
	int i = 0, j = 0;
	char buf[1000], temp_name[100], name[100];
	int b = 0, e = 0;
	int temp_code = -1;
	int sp1_code = 1, sp2_code = 2;

	f = ckopen(fname, "r");

	if( strcmp(species, species2) == 0 ) { // codex file includes annotation information for only one species
		temp_code = sp1_code;
	}
  fseek(f, 0, SEEK_SET);
  while( fgets(buf, 1000, f) ) {
		if( isspace(buf[0]) != 0 ) {}
    else if( buf[0] == '#' ) {
      sscanf(buf, "%*s %s", temp_name);
			if( strcmp(temp_name, species) == 0 ) {
				temp_code = sp1_code;
			}
			else if( strcmp(temp_name, species2) == 0 ) {
				temp_code = sp2_code;
			}
    }
    else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((buf[0] == '>') || (buf[0] == '<')) ) {
      sscanf(buf+1, "%d %d %s %*s", &b, &e, name);
      genes[j].reg = assign_I(b, e+1);
      if( temp_code == sp1_code ) genes[j].sp_id = SELF1;
      else if( temp_code == sp2_code ) genes[j].sp_id = SELF2;
      genes[j].val = 100; // 100% of the gene boundary remains
      genes[j].fid = j;
			genes[j].sign = buf[0];
			strcpy(genes[j].name, name);
      j++;
    }
    else if( (temp_code == sp1_code) || (temp_code == sp2_code) ){
      sscanf(buf, "%d %d", &b, &e);
      exons[i].reg = assign_I(b, e+1);
      exons[i].val = 100; // 100 % of the exon boundary remains
      exons[i].fid = j-1;
      exons[i].sign = genes[j-1].sign;
			strcpy(exons[i].name, name);
      if( temp_code == sp1_code ) exons[i].sp_id = SELF1;
      else if( temp_code == sp2_code ) exons[i].sp_id = SELF2;
      i++;
    }
 	}  
	fclose(f);
	*num_genes = j;
	*num_exons = i;
}

void read_exons(struct exons_list *init_exons, struct exons_list *exons, int *num_exons, struct exons_list *genes, int *num_genes, struct exons_list *skip_reg1, int *num_skip1, struct exons_list *skip_reg2, int *num_skip2, char *fname, struct sp_list *sp_code, int num_sp_code, int sp1_code, int sp2_code)
{
	FILE *f;
	int i = 0, j = 0, k = 0;
	char buf[1000], temp_name[100];
	int b = 0, e = 0;
	int temp_code = -1;
	int num_reg1 = 0, num_reg2 = 0;
	char item1[1000], item2[1000], item3[1000];

	num_reg1 = *num_skip1;
	num_reg2 = *num_skip2;

	f = ckopen(fname, "r");

  fseek(f, 0, SEEK_SET);
  while( fgets(buf, 1000, f) ) {
		if( isspace(buf[0]) != 0 ) {}
    else if( buf[0] == '#' ) {
      sscanf(buf, "%*s %s", temp_name);
      for( k = 0; k < num_sp_code; k++ ) {
        if( strcmp(sp_code[k].name, temp_name) == 0 ) temp_code = sp_code[k].id;
      }
    }
    else if( ((temp_code == sp1_code) || (temp_code == sp2_code)) && ((buf[0] == '>') || (buf[0] == '<')) ) {
			if( sscanf(buf+1, "%s %s %s", item1, item2, item3) != 3 ) {
				fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
			}
			else {
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
			}
      genes[j].reg = assign_I(b, e+1);
      if( temp_code == sp1_code ) genes[j].sp_id = SELF1;
      else if( temp_code == sp2_code ) genes[j].sp_id = SELF2;
      genes[j].val = 100; // 100% of the gene boundary remains
      genes[j].fid = j;
      j++;
    }
    else if( (temp_code == sp1_code) || (temp_code == sp2_code) ){
			if( sscanf(buf+1, "%s %s", item1, item2) != 2 ) {
				fatalf("%s: unsupported format in the %s codex file\n", buf, temp_name);
			}
			else {
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
			}
      exons[i].reg = assign_I(b, e);
      exons[i].val = 100; // 100 % of the exon boundary remains
      exons[i].fid = i;
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
	FILE *f;
	int b1, e1, b2, e2, temp;
	char len1[100], len2[100], strand[100];

	if((f = fopen(fname, "r")) == NULL) {
	}
	else {
    while(fgets(S, BIG, f)) {
      if( S[0] == '#' ) {
        while( S[0] == '#' ) {
          fgets(S, BIG, f);
          if( strncmp(S, "##maf", 5) == 0 ) {
            count = 0;
          }
        }
        count = 0;
      }

      if( S[0] == 'a' ) {
        if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))
          fatalf("cannot find alignment in %s", fname);
        if( (sscanf(S, "%*s %s %d %d %*s %s", species, &b1, &e1, len1) != 4) || (sscanf(T, "%*s %s %d %d %s %s", species2, &b2, &e2, strand, len2) != 5)) {
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
  fclose(f);
	return(count);
}

int count_lines(char *fname)
{
	FILE *f;
	char buf[1000];
	int count = 0;

	if((f = ckopen(fname, "r")) != NULL ) {
		while(fgets(buf, 1000, f)) count++;	

		fclose(f);
	}
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
