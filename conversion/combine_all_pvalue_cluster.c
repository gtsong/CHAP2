#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SPECIES_NUM 100
struct Pvalue {
	double pvalue;
	double GC_len;
	char outgroup[100];
	int max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2;
	float d_conf;
	char outgroup_name[2][100], outgroup_orient[2], orthologs_block[2][100];
	int outgroup_start[2], outgroup_end[2], triplet_status;
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
			return(strcmp((*(struct Pvalue **)a)->outgroup, (*(struct Pvalue **)b)->outgroup));
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
	FILE *fp[MAX_SPECIES_NUM], *fp_index, *fp_out;
	char buf[5000], chr1[100], chr2[100], orient;
	int index, next_index, b1, e1, b2, e2, length, m1, n1, k1, m2, n2, k2, m3, n3, k3, i, j, done = 0, status1, status2, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2;
	char identity[100], identity1[100], pvalue1[100], pvalue2[100], pvalue3[100], d_conf[100];
	int b1_1 = 0, e1_1, b2_1, e2_1, length_1, identity_num = 0, first_index, direction;
	double identity_1 = 0.0;
	int species_num;
	char *outgroup[MAX_SPECIES_NUM];
	char all_outgroup[1000];
	float max_d_conf;
	char outgroup_name[2][100], outgroup_orient[2], check_outgroup[100], orthologs_block[2][100];
	int total_pair = 0, total_GC = 0, outgroup_start[2], outgroup_end[2], long_alignment, triplet_status;

	species_num = argc - 3;

	fp_index = fopen(argv[1], "r");
	for(i=0; i<species_num; i++) {
		fp[i] = fopen(argv[i+1], "r");
		strtok(argv[i+1], "/");
		outgroup[i] = strtok(NULL, ".");
	}
	fp_out = fopen(argv[argc-1], "w");
	fgets(buf, 5000, fp_index);
	sscanf(buf, "%d", &next_index);
	while(!done) {
		strcpy(identity1, "0.000000000");
		pvalue_num = 0;
		//printf("%d\n", next_index);
		first_index = 1;
		while(1) {
			for(i=0; i<species_num; i++) {
				if(fgets(buf, 5000, fp[i]) == NULL) {
					done = 1;
					break;
				}
				sscanf(buf, "%d %s %d %d %s %d %d %c %d %d %d %d %s %d %d %d %s %d %d %d %s %s %d %d %d %d %d %d %d %d %d %d %s %s %d %d %c %s %d %d %c %s %s %d", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, &m1, &n1, &k1, pvalue1, &m2, &n2, &k2, pvalue2, &m3, &n3, &k3, pvalue3, identity, &min_len, &max_len, &max_beg1, &min_beg1, &min_end1, &max_end1, &max_beg2, &min_beg2, &min_end2, &max_end2, d_conf, outgroup_name[0], &outgroup_start[0], &outgroup_end[0], &outgroup_orient[0], outgroup_name[1], &outgroup_start[1], &outgroup_end[1], &outgroup_orient[1], orthologs_block[0], orthologs_block[1], &triplet_status);
				if(atof(identity) > atof(identity1))
					strcpy(identity1, identity);
				//if(triplet_status < 3)
				//	continue;
				if(strcmp(pvalue3, "unknown") != 0 && !(m3==0 && n3==0)) {
					pvalue[pvalue_num] = (struct Pvalue *)malloc(sizeof(struct Pvalue));
					pvalue[pvalue_num]->pvalue = atof(pvalue3);
					pvalue[pvalue_num]->GC_len = (float)min_len;// + max_len) / 2.0;
					pvalue[pvalue_num]->max_beg1 = max_beg1;
					pvalue[pvalue_num]->min_beg1 = min_beg1;
					pvalue[pvalue_num]->min_end1 = min_end1;
					pvalue[pvalue_num]->max_end1 = max_end1;
					pvalue[pvalue_num]->max_beg2 = max_beg2;
					pvalue[pvalue_num]->min_beg2 = min_beg2;
					pvalue[pvalue_num]->min_end2 = min_end2;
					pvalue[pvalue_num]->max_end2 = max_end2;
					pvalue[pvalue_num]->d_conf = atof(d_conf);
					strcpy(pvalue[pvalue_num]->outgroup, outgroup[i]);
					strcpy(pvalue[pvalue_num]->outgroup_name[0], outgroup_name[0]);
					strcpy(pvalue[pvalue_num]->outgroup_name[1], outgroup_name[1]);
					pvalue[pvalue_num]->outgroup_start[0] = outgroup_start[0];
					pvalue[pvalue_num]->outgroup_start[1] = outgroup_start[1];
					pvalue[pvalue_num]->outgroup_end[0] = outgroup_end[0];
                                        pvalue[pvalue_num]->outgroup_end[1] = outgroup_end[1];
					pvalue[pvalue_num]->outgroup_orient[0] = outgroup_orient[0];
                                        pvalue[pvalue_num]->outgroup_orient[1] = outgroup_orient[1];
					strcpy(pvalue[pvalue_num]->orthologs_block[0], orthologs_block[0]);
					strcpy(pvalue[pvalue_num]->orthologs_block[1], orthologs_block[1]);
					pvalue[pvalue_num]->triplet_status = triplet_status;
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
			status1 = status2 = 0;
			direction = 0;
			total_pair++;
			if(pvalue_num > 0) {
				if((pvalue[pvalue_num-1]->pvalue * pvalue_num) < atof(argv[argc-2])) {
					status1 = status2 = 2;
					total_GC++;
				}
				else
					status1 = status2 = 1;
				max_d_conf = 0.0;
				strcpy(all_outgroup, "no");
				for(i=pvalue_num-1; i>=0; i--) {
					if((pvalue[i]->pvalue * pvalue_num) < atof(argv[argc-2])) {
						if(fabs(pvalue[i]->d_conf) > fabs(max_d_conf))
							max_d_conf = pvalue[i]->d_conf;
						
						if(pvalue[i]->d_conf > 0.0)
							direction = 1;
						else if(pvalue[i]->d_conf < 0.0)
							direction = 2;
						else
							direction = 0;
							
						long_alignment = 0;
						if(strcmp(pvalue[i]->outgroup_name[0], "NAN") == 0)
							strcpy(check_outgroup, pvalue[i]->outgroup_name[1]);
						else
							strcpy(check_outgroup, pvalue[i]->outgroup_name[0]);
						for(j=pvalue_num-1; j>i; j--) {
							if(strcmp(check_outgroup, pvalue[j]->outgroup_name[0]) == 0)
								long_alignment = 1;
							if(strcmp(check_outgroup, pvalue[j]->outgroup_name[1]) == 0)
								long_alignment = 1;
						}
						
						if(long_alignment == 0) {
							if(orient == '+')
								fprintf(fp_out, "%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%f\t%d\t%d\t%f\t%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%s\t%d\n", index, chr1, b1, e1, chr2, b2, e2, orient, length_1, identity_1/(double)identity_num, status1, status2, (float)(pvalue[i]->min_end1 - pvalue[i]->min_beg1 + 1), pvalue[i]->pvalue, pvalue[i]->max_beg1, pvalue[i]->min_beg1, pvalue[i]->min_end1, pvalue[i]->max_end1, pvalue[i]->max_beg2, pvalue[i]->min_beg2, pvalue[i]->min_end2, pvalue[i]->max_end2, direction, pvalue[i]->outgroup_name[0], pvalue[i]->outgroup_start[0], pvalue[i]->outgroup_end[0], pvalue[i]->outgroup_orient[0], pvalue[i]->outgroup_name[1], pvalue[i]->outgroup_start[1], pvalue[i]->outgroup_end[1], pvalue[i]->outgroup_orient[1], pvalue[i]->orthologs_block[0], pvalue[i]->orthologs_block[1], pvalue[i]->triplet_status);
							else
								fprintf(fp_out, "%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%f\t%d\t%d\t%f\t%e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%s\t%d\n", index, chr1, b1, e1, chr2, b2, e2, orient, length_1, identity_1/(double)identity_num, status1, status2, (float)(pvalue[i]->min_end1 - pvalue[i]->min_beg1 + 1), pvalue[i]->pvalue, pvalue[i]->max_beg1, pvalue[i]->min_beg1, pvalue[i]->min_end1, pvalue[i]->max_end1, pvalue[i]->max_end2, pvalue[i]->min_end2, pvalue[i]->min_beg2, pvalue[i]->max_beg2, direction, pvalue[i]->outgroup_name[0], pvalue[i]->outgroup_start[0], pvalue[i]->outgroup_end[0], pvalue[i]->outgroup_orient[0], pvalue[i]->outgroup_name[1], pvalue[i]->outgroup_start[1], pvalue[i]->outgroup_end[1], pvalue[i]->outgroup_orient[1], pvalue[i]->orthologs_block[0], pvalue[i]->orthologs_block[1], pvalue[i]->triplet_status);
						}
						
						if(strcmp(all_outgroup, "no") !=  0)
							sprintf(all_outgroup, "%s&%s", all_outgroup, pvalue[i]->outgroup);
						else
							strcpy(all_outgroup, pvalue[i]->outgroup);
					}
				}
							
				if(max_d_conf > 0.0)
					direction = 1;
				else if(max_d_conf < 0.0)
					direction = 2;
				else
					direction = 0;

			}
		}
	}
	for(i=0; i<species_num; i++)
		fclose(fp[i]);
		
//	printf("%s", chr1);
//	for(i=strlen(chr1); i<30; i++)
//		printf(" ");
//	printf("%d \t\t   %d (%2.1f%%)\n", total_pair, total_GC, (float)total_GC*100.0/(float)total_pair);

	fclose(fp_index);
	fclose(fp_out);

	return 0;
}

