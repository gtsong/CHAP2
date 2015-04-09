/* FDR_cluster - Bonferroni correction is applied. The smallest p-value for each triplet test is multiplied by the number of test */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SPECIES_NUM 1000

double average_len;

struct Pvalue {
	double pvalue;
	double GC_len;
	int max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2;
};

struct Pvalue *pvalue[10000];
int pvalue_num;

static int compar(const void *a, const void *b) {
  if((*(struct Pvalue **)a)->pvalue < (*(struct Pvalue **)b)->pvalue) {
	  return 1;
	}
  else if((*(struct Pvalue **)a)->pvalue == (*(struct Pvalue **)b)->pvalue)
	{
  	if((*(struct Pvalue **)a)->GC_len > (*(struct Pvalue **)b)->GC_len) {
			return 1;
		}
    else if( (*(struct Pvalue **)a)->GC_len == (*(struct Pvalue **)b)->GC_len) 
    {
      if((((*(struct Pvalue **)a)->max_end1) - ((*(struct Pvalue **)a)->max_beg1)) > (((*(struct Pvalue **)b)->max_end1) - ((*(struct Pvalue **)b)->max_beg1)) )			{
				return 1;
			}
      else if((((*(struct Pvalue **)a)->max_end1) - ((*(struct Pvalue **)a)->max_beg1)) == (((*(struct Pvalue **)b)->max_end1) - ((*(struct Pvalue **)b)->max_beg1)) )			
			{
  			if((*(struct Pvalue **)a)->max_beg1 > (*(struct Pvalue **)b)->max_beg1){
					return 1;
				}
				else return -1;
			}
			else return -1;
    }
		else {
			return -1;
		}
	}
	else {
	  return -1;
	}
}

int main(int argc, char *argv[]) {
	FILE *fp[MAX_SPECIES_NUM], *fp_index;
	char buf[5000], chr1[100], chr2[100], orient;
	int index, next_index, b1, e1, b2, e2, length, m1, n1, k1, m2, n2, k2, m3, n3, k3, i, done = 0, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2;
	char identity[100], identity1[100], pvalue1[100], pvalue2[100], pvalue3[100], d_conf[100];
	int b1_1, e1_1, b2_1, e2_1, length_1, identity_num, first_index;
	double identity_1;
	int species_num;

	species_num = argc - 1;
	
	for(i=0; i<species_num; i++)
		fp[i] = fopen(argv[i+1], "r");
	fp_index = fopen(argv[1], "r");
	fgets(buf, 5000, fp_index);
	sscanf(buf, "%d", &next_index);
	while(!done) {
		strcpy(identity1, "0.000000000");
		pvalue_num = 0;
		first_index = 1;
		while(1) {
			for(i=0; i<species_num; i++) {
				if(fgets(buf, 5000, fp[i]) == NULL) {
					done = 1;
					break;
				}
				sscanf(buf, "%d %s %d %d %s %d %d %c %d %d %d %d %s %d %d %d %s %d %d %d %s %s %d %d %d %d %d %d %d %d %d %d %s\n", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, &m1, &n1, &k1, pvalue1, &m2, &n2, &k2, pvalue2, &m3, &n3, &k3, pvalue3, identity, &min_len, &max_len, &max_beg1, &min_beg1, &min_end1, &max_end1, &max_beg2, &min_beg2, &min_end2, &max_end2, d_conf);
				if(strcmp(chr1, "chrY") == 0 || strcmp(chr2, "chrY") == 0)
					continue;
				if(atof(identity) > atof(identity1))
					strcpy(identity1, identity);
				if(strcmp(pvalue3, "unknown") != 0 && !(m3==0 && n3==0)) {
					pvalue[pvalue_num] = (struct Pvalue *)malloc(sizeof(struct Pvalue));
					pvalue[pvalue_num]->pvalue = atof(pvalue3);
					pvalue[pvalue_num]->GC_len = (float)(min_len + max_len) / 2.0;
					pvalue[pvalue_num]->max_beg1 = b1 + max_beg1;
					pvalue[pvalue_num]->min_beg1 = b1 + min_beg1;
					pvalue[pvalue_num]->min_end1 = b1 + min_end1;
					pvalue[pvalue_num]->max_end1 = b1 + max_end1;
					if(orient == '+') {
						pvalue[pvalue_num]->max_beg2 = b2 + max_beg2;
						pvalue[pvalue_num]->min_beg2 = b2 + min_beg2;
						pvalue[pvalue_num]->min_end2 = b2 + min_end2;
						pvalue[pvalue_num]->max_end2 = b2 + max_end2;
					}
					else {
						pvalue[pvalue_num]->max_beg2 = e2 - max_end2;
						pvalue[pvalue_num]->min_beg2 = e2 - min_end2;
						pvalue[pvalue_num]->min_end2 = e2 - min_beg2;
						pvalue[pvalue_num]->max_end2 = e2 - max_beg2;
					}
					pvalue_num++;
				}
			}
			
			if(first_index == 1) {
				b1_1 = b1;
				e1_1 = e1;
				b2_1 = b2;
				e2_1 = e2;
				length_1 = length;
				identity_1 = atof(identity1);
				if(identity_1 > 0.0)
					identity_num = 1;
				else
					identity_num = 0;
			}
			else {
				e1_1 = e1;
				if(orient == '+')
					e2_1 = e2;
				else
					b2_1 = b2;
				length_1 = e1_1 - b1_1 + 1;
				if(atof(identity) > 0.0) {
					identity_1 += atof(identity1);
					identity_num++;
				}
			}
			
			if(fgets(buf, 5000, fp_index)) {
				sscanf(buf, "%d", &next_index);
				if(next_index != index)
					break;
				else
					first_index = 0;
			}
			else
				break;
			
		}
		
		if(!done) {
			qsort((void **)pvalue, pvalue_num, sizeof(struct Pvalue *), compar);
			if(pvalue_num > 0) {
					printf("%f\n", pvalue[pvalue_num-1]->pvalue*pvalue_num);
			}
		}
	}
	for(i=0; i<species_num; i++)
		fclose(fp[i]);
		
	fclose(fp_index);

	return 0;
}

