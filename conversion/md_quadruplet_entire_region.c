#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "maf.h"
#include "util.h"
#include "contigs_op.h"

#define BIG 1000000
#define MAX_GAP 5000
#define MIN_COVERAGE 0.5
#define MIN_MATCH 0.8	
#define MAX_GC_RATIO 0.8
#define ORTHOLOGS_OVERLAP_RATIO 0.5

typedef struct Net {
	int start, end, start2, end2, len, block_index;
	char *net1, *net2, chr1[100], contig1[100], chr2[100], contig2[100], ctg_id1, ctg_id2, orient;
} net_t;

typedef struct Net_level {
	int chain_id, end, chr_size;
	char chr2[100];
} net_level_t;

typedef struct Aln {
	int start1, len1, src_len1, start2, len2, src_len2;
	int ctg_id1, ctg_id2;
	char text1[BIG], text2[BIG];
	struct Aln *next;
} aln_t;

typedef struct Chain {
	char chr1[100], chr2[100], orient;
	char contig1[100], contig2[100];
	int ctg_id1, ctg_id2;
	int start1, end1, start2, end2;
	int len1, len2;
	aln_t *aln;
} chain_t;

net_level_t net_level[1000];
int net_level_num;
net_t net[2][BIG];
int rescored[1][BIG];
int net_num[2];
char buf[BIG];
char seqs[2][3][BIG];
int seq2_pos[2][BIG];
int IS[3][BIG], IS_num[3], IS_sites[4][BIG], IS_ref[BIG];
int IS_quadruplet[BIG], IS_quadruplet_num, IS_quadruplet_sites[2][BIG];
int d_beg_index[10000], d_end_index[10000], d_index_num;
int sites_num[2];
int m[2], n[2], max_descent[3], beg_index = 0, end_index = 0;
double identity;
int max_len, min_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, start[2], end[2];
int offset[3];
int site_beg[2], site_end[2];
int mismatch_net[100000];
int mismatch_net_num;
int long_seqs;
//score parameters
int if_rescoring = 0;
int fill_score = -100, T = 1, Z = 1, O = 400, E = 30, X = 910, Y = 9400, K = 3000, L = 3000; 
int score[4][4] = {{91, -114, -31, -123}, {-114, 100, -125, -31}, {-31, -125, 100, -114}, {-123, -31, -114, 91}};
char outgroup_name[2][100], outgroup_orient[2];
int outgroup_start[2], outgroup_end[2], orthologs_block[2][100], orthologs_block_num[2];
double max_gc_ratio;

void read_score(char *filename) {
	FILE *fp;
	char buf[10000], *temp_str;
	
	fp = fopen(filename, "r");
	while(fgets(buf, 10000, fp)) {
		if(buf[0] == '#')
			continue;
		if(strncmp(buf, "fill_score", 10) == 0) {
			strtok(buf, "=");
			temp_str = strtok(NULL, "=");
			sscanf(temp_str, "%d", &fill_score);
		}
		else if (buf[1] == '=') {
			if(buf[0] == 'T') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &T);
			}
			else if(buf[0] == 'O') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &O);
			}
			else if(buf[0] == 'E') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &E);
			}
			else if(buf[0] == 'X') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &X);
			}
			else if(buf[0] == 'Y') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &Y);
			}
			else if(buf[0] == 'K') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &K);
			}
			else if(buf[0] == 'L') {
				strtok(buf, "=");
				temp_str = strtok(NULL, "=");
				sscanf(temp_str, "%d", &L);
			}
		}
		else if(buf[0] == 'A')
			sscanf(buf, "A %d %d %d %d", &score[0][0], &score[0][1], &score[0][2], &score[0][3]);
		else if(buf[0] == 'C')
			sscanf(buf, "C %d %d %d %d", &score[1][0], &score[1][1], &score[1][2], &score[1][3]);
		else if(buf[0] == 'G')
			sscanf(buf, "G %d %d %d %d", &score[2][0], &score[2][1], &score[2][2], &score[2][3]);
		else if(buf[0] == 'T')
			sscanf(buf, "T %d %d %d %d", &score[3][0], &score[3][1], &score[3][2], &score[3][3]);
	}
	/*
	printf("fill_score = %d; T = %d; O = %d; E = %d; X = %d; Y = %d; K = %d; L = %d\n", fill_score, T, O, E, X, Y, K, L);
	for(i=0; i<4; i++) {
		for(j=0; j<4; j++)
			printf("%d\t", score[i][j]);
		printf("\n");
	}
	*/
}

int cal_score(char A, char B, int gap_opend) {
	if(A == '-' || B == '-') {
		if(gap_opend == 1)
			return -E;
		else
			return -O;
	}
	if(toupper(A) == 'A') {
		if(toupper(B) == 'A')
			return score[0][0];
		else if(toupper(B) == 'C')
			return score[0][1];
		else if(toupper(B) == 'G')
			return score[0][2];
		else if(toupper(B) == 'T')
			return score[0][3];
		else
			return fill_score;
	}
	else if(toupper(A) == 'C') {
		if(toupper(B) == 'A')
			return score[1][0];
		else if(toupper(B) == 'C')
			return score[1][1];
		else if(toupper(B) == 'G')
			return score[1][2];
		else if(toupper(B) == 'T')
			return score[1][3];
		else
			return fill_score;
	}
	else if(toupper(A) == 'G') {
		if(toupper(B) == 'A')
			return score[2][0];
		else if(toupper(B) == 'C')
			return score[2][1];
		else if(toupper(B) == 'G')
			return score[2][2];
		else if(toupper(B) == 'T')
			return score[2][3];
		else
			return fill_score;
	}
	else if(toupper(A) == 'T') {
		if(toupper(B) == 'A')
			return score[3][0];
		else if(toupper(B) == 'C')
			return score[3][1];
		else if(toupper(B) == 'G')
			return score[3][2];
		else if(toupper(B) == 'T')
			return score[3][3];
		else
			return fill_score;
	}
	else
		return fill_score;
}

void rescoring(int i) {
	int j, k, max_score, cur_score, gap_opend, aln_beg, aln_max, gap_free, net_len;

	max_score = cur_score = gap_opend = aln_beg = aln_max = 0;
	gap_free = 1;
	net_len = strlen(net[0][i].net1);
	for(j=0; j<net_len; j++) {
		cur_score += cal_score(net[0][i].net1[j], net[0][i].net2[j], gap_opend);	
		if(net[0][i].net1[j] == '-' || net[0][i].net2[j] == '-') {
			gap_opend = 1;
			gap_free = 0;
		}
		else
			gap_opend = 0;
				
		if(cur_score >= max_score) {
			max_score = cur_score;
			aln_max = j;
		}

		//if(cur_score < 0 || (gap_free == 1 && (max_score - cur_score) > X) || (gap_free == 0 && (max_score - cur_score) > Y) || j == net_len - 1) { 
		if(cur_score < 0 || (max_score - cur_score) > Y || j == net_len - 1) {
			if(max_score < L) //throw away whole segment
				for(k=aln_beg; k<=j; k++) {
					if(net[0][i].net1[k] != '-')
						net[0][i].net1[k] = 'X';
					if(net[0][i].net2[k] != '-')
						net[0][i].net2[k] = 'X';
				}
			else {	//throw away the dropping region
				for(k=aln_max+1; k<=j; k++) {
					if(net[0][i].net1[k] != '-')
						net[0][i].net1[k] = 'X';
					if(net[0][i].net2[k] != '-')
						net[0][i].net2[k] = 'X';
				}
			}
			aln_beg = aln_max = j + 1;
			max_score = cur_score = 0;
			gap_free = 1;
		}
	}
	
	return;
}

int sort_net(void* a, void* b) {
  if( (((net_t*)a)->start - ((net_t*)b)->start) != 0 )
  {
    return(((net_t*)a)->start - ((net_t*)b)->start);
  }
  else if( (((net_t*)a)->start2 - ((net_t*)b)->start2) != 0 )
  {
    return(((net_t*)a)->start2 - ((net_t*)b)->start2);
  }
  else if( ( (((net_t*)a)->end - ((net_t*)a)->start) - (((net_t*)b)->end - ((net_t*)b)->start) ) != 0 ) {
    return( (((net_t*)a)->end - ((net_t*)a)->start) - (((net_t*)b)->end - ((net_t*)b)->start) );
  }
  else {
    return( (((net_t*)a)->end2 - ((net_t*)a)->start2) - (((net_t*)b)->end2 - ((net_t*)b)->start2) );
  }
}

int read_orthologs(char filename[1000], int index, struct n_pair **contigs1, int *num_contigs1, int *num_blocks1, int **len_sum1, struct n_pair **contigs2, int *num_contigs2, int *num_blocks2, int **len_sum2) {
	struct mafFile *mf;
	struct mafAli *ali;
	struct mafComp *mc1, *mc2;
	int num1 = 0, num2 = 0, i = 0, id = -1;
	int org_num1 = 0, org_num2 = 0;

	num1 = *num_contigs1;
	num2 = *num_contigs2;
	org_num1 = num1;
	org_num2 = num2;
	
	net_num[index] = 0;
	mf = mafOpen(filename, 0);
	while((ali = mafNext(mf)) != NULL) {
		mc1 = ali->components;
		if(mc1 == NULL) {
			printf("Wrong orthologs file\n");
			exit(1);
		}
		mc2 = mc1->next;
		if(mc2 == NULL) {
			printf("Wrong orthologs file\n");
			exit(1);
		}
		net[index][net_num[index]].block_index = ali->score;
		net[index][net_num[index]].start = mc1->start + 1;
		net[index][net_num[index]].end = mc1->start + mc1->size;
		strcpy(net[index][net_num[index]].chr1, mc1->name);
		strcpy(net[index][net_num[index]].contig1, mc1->contig);
		id = add_ctg_name(net_num[index], mc1->name, mc1->contig, mc1->srcSize, *contigs1, num1);
		net[index][net_num[index]].ctg_id1 = id;

		if( id == num1 ) {
			num1++;
		}

		if( num1 >= ((*num_blocks1) * ALLOC_UNIT) ) {
			*num_blocks1 = *num_blocks1 + 1;
			*contigs1 = (struct n_pair *) ckrealloc(*contigs1, sizeof(struct n_pair) * (*num_blocks1) * ALLOC_UNIT);
			init_pair(*contigs1, ((*num_blocks1)-1) * ALLOC_UNIT, (*num_blocks1) * ALLOC_UNIT);
		}

		if(mc2->strand == '+') {
			net[index][net_num[index]].start2 = mc2->start + 1;
			net[index][net_num[index]].end2 = mc2->start + mc2->size;
		}
		else {
			net[index][net_num[index]].start2 = mc2->srcSize - mc2->start - mc2->size + 1;
			net[index][net_num[index]].end2 = mc2->srcSize - mc2->start;
		}
		net[index][net_num[index]].len = mc2->srcSize;
		strcpy(net[index][net_num[index]].chr2, mc2->name);
		strcpy(net[index][net_num[index]].contig2, mc2->contig);
		id = add_ctg_name(net_num[index], mc2->name, mc2->contig, mc2->srcSize, *contigs2, num2);

		net[index][net_num[index]].ctg_id2 = id;
		if( id == num2 ) {
			num2++;
		}

		if( num2 >= ((*num_blocks2) * ALLOC_UNIT) ) {
			*num_blocks2 = *num_blocks2 + 1;
			*contigs2 = (struct n_pair *) ckrealloc(*contigs2, sizeof(struct n_pair) * (*num_blocks2) * ALLOC_UNIT);
			init_pair(*contigs2, ((*num_blocks2)-1) * ALLOC_UNIT, (*num_blocks2) * ALLOC_UNIT);
		}
 
		net[index][net_num[index]].orient = mc2->strand;
		net[index][net_num[index]].net1 = (char *)ckalloc(sizeof(char) * (strlen(mc1->text) + 1));
		strcpy(net[index][net_num[index]].net1, mc1->text);
		net[index][net_num[index]].net2 = (char *)ckalloc(sizeof(char) * (strlen(mc2->text) + 1));
		strcpy(net[index][net_num[index]].net2, mc2->text);
		net_num[index]++;
	}
	
	*num_contigs1 = num1;
	*num_contigs2 = num2;

	if( ( org_num1 == 0 ) && (num1 > 0 ) ) {
		*len_sum1 = (int *) ckalloc(sizeof(int) * (*num_contigs1) );	
		(*len_sum1)[0] = 0;
		org_num1++;
	}
	else if( num1 > org_num1 ) {
		*len_sum1 = (int *) ckrealloc(*len_sum1, sizeof(int) * (*num_contigs1) );	
	}

	for( i = org_num1; i < num1; i++ ) {
		(*len_sum1)[i] = (*len_sum1)[i-1] + (*contigs1)[i-1].len;
	}

	if( ( org_num2 == 0 ) && (num2 > 0 ) ) {
		*len_sum2 = (int *) ckalloc(sizeof(int) * (*num_contigs2) );	
		(*len_sum2)[0] = 0;
		org_num2++;
	}
	else if( num2 > org_num2 ) {
		*len_sum2 = (int *) ckrealloc(*len_sum2, sizeof(int) * (*num_contigs2) );	
	}

	for( i = org_num2; i < num2; i++ ) {
		(*len_sum2)[i] = (*len_sum2)[i-1] + (*contigs2)[i-1].len;
	}

	for( i = 0; i < net_num[index]; i++ ) {
		id = net[index][i].ctg_id1;
		net[index][i].start = net[index][i].start + (*len_sum1)[id];
		net[index][i].end = net[index][i].end + (*len_sum1)[id];
		id = net[index][i].ctg_id2;
		net[index][i].start2 = net[index][i].start2 + (*len_sum2)[id];
		net[index][i].end2 = net[index][i].end2 + (*len_sum2)[id];
	}

	qsort((void*)net[index], net_num[index], sizeof(net_t), (void *)sort_net);
	/*
	printf("net_num = %d\n", net_num[0]);
	int i;
	for(i=0; i<net_num[0]; i++) 
		printf("%d\t%d\n", net[0][i].start, net[0][i].end);
	*/
	return 0;
}

void find_IS(int ref) {
	int i;
	
	IS_num[ref-1] = 0;
	/*
	for(j=0; j<3; j++) {
		for(i=0; i<sites_num; i++)
			printf("%c", seqs[j][i]);
		printf("\n\n");
	}
	*/
	
	m[ref-1] = 0;
	n[ref-1] = 0;
	
	for(i=site_beg[ref-1]; i<site_end[ref-1]; i++) {
		//there are gaps in a site
		if((seqs[ref-1][0][i] == '-') || (seqs[ref-1][1][i] == '-') || (seqs[ref-1][2][i] == '-') || (seqs[ref-1][0][i] == 'N') || (seqs[ref-1][1][i] == 'N') || (seqs[ref-1][2][i] == 'N') || (seqs[ref-1][0][i] == 'n') || (seqs[ref-1][1][i] == 'n') || (seqs[ref-1][2][i] == 'n'))
			continue;
		//all are the same in a site
		if((toupper(seqs[ref-1][0][i]) == toupper(seqs[ref-1][1][i])) && (toupper(seqs[ref-1][0][i]) == toupper(seqs[ref-1][2][i])))
			continue;
		// all are different in a site
		else if((toupper(seqs[ref-1][0][i]) != toupper(seqs[ref-1][1][i])) && (toupper(seqs[ref-1][0][i]) != toupper(seqs[ref-1][2][i])) && (toupper(seqs[ref-1][1][i]) != toupper(seqs[ref-1][2][i])))
			continue;
		//type P
		else if(toupper(seqs[ref-1][0][i]) == toupper(seqs[ref-1][2][i])) {
			IS[ref-1][IS_num[ref-1]] = 1;
			m[ref-1]++;
		}
		//type Q	
		else if(toupper(seqs[ref-1][0][i]) == toupper(seqs[ref-1][1][i])) {
			IS[ref-1][IS_num[ref-1]] = 2;
			n[ref-1]++;
		}
		else
			continue;
		
		//printf("%d\n", i+1);
		IS_sites[ref-1][IS_num[ref-1]] = i;
		IS_num[ref-1]++;
	}
	
	/*
	int walk;
	walk = 0;
	for(i=IS_num[ref-1]-1; i>=0; i--) {
		if(IS[ref-1][i] == 1)
			walk++;
		else
			walk--;
		printf("%d\t%d\n", site_end[ref-1] - IS_sites[ref-1][i] - 1, walk);
	}
	
	printf("\n\n");
	
	for(i=0; i<IS_num[ref-1]; i++)
		printf("%d", IS[ref-1][i]);
	printf("\n\n");
	*/
}

void combine_IS(char orient) {
	int IS_index[2], i;

	if(orient == '+') {
		IS_index[0] = IS_index[1] = 0;
		IS_num[2] = 0;
		while(IS_index[0] < IS_num[0] || IS_index[1] < IS_num[1]) {
			//printf("%d\t%d\n", IS_index[0], IS_index[1]);
			if(IS_index[0] >= IS_num[0]) {
				for(i=IS_index[1]; i<IS_num[1]; i++) {
					IS[2][IS_num[2]] = IS[1][i];
					IS_ref[IS_num[2]] = 1;
					IS_sites[2][IS_num[2]] = seq2_pos[1][IS_sites[1][i]];
					IS_sites[3][IS_num[2]] = start[1] + IS_sites[1][i];
					IS_num[2]++;
				}
				break;
			}
			else if(IS_index[1] >= IS_num[1]) {
				for(i=IS_index[0]; i<IS_num[0]; i++) {
					IS[2][IS_num[2]] = IS[0][i];
					IS_ref[IS_num[2]] = 0;
					IS_sites[2][IS_num[2]] = start[0] + IS_sites[0][i];
					IS_sites[3][IS_num[2]] = seq2_pos[0][IS_sites[0][i]];
					IS_num[2]++;
				}
				break;
			}
			else {
				if((start[0] + IS_sites[0][IS_index[0]]) < seq2_pos[1][IS_sites[1][IS_index[1]]]) {
					IS[2][IS_num[2]] = IS[0][IS_index[0]];
					IS_ref[IS_num[2]] = 0;
					IS_sites[2][IS_num[2]] = start[0] + IS_sites[0][IS_index[0]];
					IS_sites[3][IS_num[2]] = seq2_pos[0][IS_sites[0][IS_index[0]]];
					IS_index[0]++;
				}
				else if ((start[0] + IS_sites[0][IS_index[0]]) > seq2_pos[1][IS_sites[1][IS_index[1]]]) {
					IS[2][IS_num[2]] = IS[1][IS_index[1]];
					IS_ref[IS_num[2]] = 1;
					IS_sites[2][IS_num[2]] = seq2_pos[1][IS_sites[1][IS_index[1]]];
					IS_sites[3][IS_num[2]] = start[1] + IS_sites[1][IS_index[1]];
					IS_index[1]++;
				}
				else {	
					IS[2][IS_num[2]] = IS[0][IS_index[0]];
					IS_ref[IS_num[2]] = 0;
					IS_sites[2][IS_num[2]] = start[0] + IS_sites[0][IS_index[0]];
					IS_sites[3][IS_num[2]] = seq2_pos[0][IS_sites[0][IS_index[0]]];
					IS_num[2]++;
					IS[2][IS_num[2]] = IS[1][IS_index[1]];
					IS_ref[IS_num[2]] = 1;
					IS_sites[2][IS_num[2]] = seq2_pos[1][IS_sites[1][IS_index[1]]];
					IS_sites[3][IS_num[2]] = start[1] + IS_sites[1][IS_index[1]];
					IS_index[0]++;
					IS_index[1]++;
				}
				IS_num[2]++;
			}   
		}
	} 
	else {
		IS_index[0] = 0;
		IS_index[1] = IS_num[1] - 1;
		IS_num[2] = 0;
		while(IS_index[0] < IS_num[0] || IS_index[1] >= 0) {
			//printf("%d\t%d\n", IS_index[0], IS_index[1]);
			if(IS_index[0] >= IS_num[0]) {
				for(i=IS_index[1]; i>=0; i--) {
					IS[2][IS_num[2]] = IS[1][i];
					IS_ref[IS_num[2]] = 1;
					IS_sites[2][IS_num[2]] = seq2_pos[1][IS_sites[1][i]];
					IS_sites[3][IS_num[2]] = start[1] + IS_sites[1][i];
					IS_num[2]++;
				}
				break;
			}
			else if(IS_index[1] < 0) {
				for(i=IS_index[0]; i<IS_num[0]; i++) {
					IS[2][IS_num[2]] = IS[0][i];
					IS_ref[IS_num[2]] = 0;
					IS_sites[2][IS_num[2]] = start[0] + IS_sites[0][i];
					IS_sites[3][IS_num[2]] = seq2_pos[0][IS_sites[0][i]];
					IS_num[2]++;
				}
				break;
			}
			else {
				if((start[0] + IS_sites[0][IS_index[0]]) < seq2_pos[1][IS_sites[1][IS_index[1]]]) {
					IS[2][IS_num[2]] = IS[0][IS_index[0]];
					IS_ref[IS_num[2]] = 0;
					IS_sites[2][IS_num[2]] = start[0] + IS_sites[0][IS_index[0]];
					IS_sites[3][IS_num[2]] = seq2_pos[0][IS_sites[0][IS_index[0]]];
					IS_index[0]++;
				}
				else if ((start[0] + IS_sites[0][IS_index[0]]) > seq2_pos[1][IS_sites[1][IS_index[1]]]) {
					IS[2][IS_num[2]] = IS[1][IS_index[1]];
					IS_ref[IS_num[2]] = 1;
					IS_sites[2][IS_num[2]] = seq2_pos[1][IS_sites[1][IS_index[1]]];
					IS_sites[3][IS_num[2]] = start[1] + IS_sites[1][IS_index[1]];
					IS_index[1]--;
				}
				else {
					IS[2][IS_num[2]] = IS[0][IS_index[0]];
					IS_ref[IS_num[2]] = 0;
					IS_sites[2][IS_num[2]] = start[0] + IS_sites[0][IS_index[0]];
					IS_sites[3][IS_num[2]] = seq2_pos[0][IS_sites[0][IS_index[0]]];
					IS_num[2]++;
					IS[2][IS_num[2]] = IS[1][IS_index[1]];
					IS_ref[IS_num[2]] = 1;
					IS_sites[2][IS_num[2]] = seq2_pos[1][IS_sites[1][IS_index[1]]];
					IS_sites[3][IS_num[2]] = start[1] + IS_sites[1][IS_index[1]];
					IS_index[0]++;
					IS_index[1]--;
				}
				IS_num[2]++;
			}   
		}
	}
	/*
	for(i=0; i<IS_num[2]; i++)
		printf("%d", IS[2][i]);
	printf("\n\n");
	*/
}

void find_max_descent(int ref) {
	int max_height = 0, height = 0, descent = 0, i, max_height_index = 0;
	
	max_descent[ref-1] = 0;
	beg_index = end_index = 0;
	for(i=0; i<IS_num[ref-1]; i++) {
		if(IS[ref-1][i] == 1) {
			height++;
			if( height > max_height ) { 
				max_height = height;
				max_height_index = i;
			}
			descent = max_height - height;
			if( descent >= max_descent[ref-1] ) { 
				max_descent[ref-1] = descent;
				beg_index = max_height_index;
				end_index = i;
			}
		}
		else if(IS[ref-1][i] == 2) {
			height--;
			if( height > max_height ) { 
				max_height = height;
				max_height_index = i;
			}
			descent = max_height - height;
			if( descent >= max_descent[ref-1] ) {
				max_descent[ref-1] = descent;
				beg_index = max_height_index;
				end_index = i;
			}
		}
	}
	
	if(long_seqs == 1 || ref < 3)
		return;
	
	min_len = max_len = 0;
	/*
	for(i=IS_sites[beg_index+1]; i<=IS_sites[end_index]; i++)
			if(seqs[0][i] != '-')
				min_len++;
	*/
	min_len = IS_sites[2][end_index] - IS_sites[2][beg_index+1] + 1;
	min_beg1 = IS_sites[2][beg_index+1];
	min_beg2 = IS_sites[3][beg_index+1];
	min_end1 = IS_sites[2][end_index];
	min_end2 = IS_sites[3][end_index];

	if(beg_index > 0 && end_index < IS_num[2] - 1) {
		/*
		for(i=IS_sites[beg_index]+1; i<IS_sites[end_index+1]; i++)
			if(seqs[0][i] != '-')
				max_len++;
		*/
		max_len = IS_sites[2][end_index+1] - IS_sites[2][beg_index];
						
		max_beg1 = IS_sites[2][beg_index]+1;
		max_beg2 = IS_sites[3][beg_index]+1;
		max_end1 = IS_sites[2][end_index+1]-1;
		max_end2 = IS_sites[3][end_index+1]-1;
	}
	else if(beg_index == 0 && end_index < IS_num[2] - 1) {
		/*
		for(i=0; i<IS_sites[end_index+1]; i++)
			if(seqs[0][i] != '-')
				max_len++;
		*/
		max_len = IS_sites[2][end_index+1] - start[0];
						
		max_beg1 = start[0];
		max_beg2 = start[1];
		max_end1 = IS_sites[2][end_index+1]-1;
		max_end2 = IS_sites[3][end_index+1]-1;
	}
	else if(beg_index > 0 && end_index == IS_num[2] - 1) {
		/*
		for(i=IS_sites[beg_index]+1; i<sites_num; i++)
			if(seqs[0][i] != '-')
				max_len++;
		*/
		max_len = start[0] + site_end[0] - IS_sites[2][beg_index]-1;
						
		max_beg1 = IS_sites[2][beg_index]+1;
		max_beg2 = IS_sites[3][beg_index]+1;
		max_end1 = start[0] + site_end[0] - 1;
		max_end2 = start[1] + site_end[1] - 1;
	}
	else {
		/*
		for(i=0; i<sites_num; i++)
			if(seqs[0][i] != '-')
				max_len++;
		*/
		max_len = site_end[0] + 1;
						
		max_beg1 = start[0];
		max_beg2 = start[1];
		max_end1 = start[0] + site_end[0] - 1;
		max_end2 = start[1] + site_end[1] - 1;
	}
}

char Complement(char c) {
	if(c == 'A')
		return 'T';
	else if(c == 'C')
		return 'G';
	else if(c == 'G')
		return 'C';
	else if(c == 'T')
		return 'A';
	else if(c == 'a')
		return 't';
	else if(c == 'c')
		return 'g';
	else if(c == 'g')
		return 'c';
	else if(c == 't')
		return 'a';
	else if(c == 'N')
		return 'N';
	else if(c == 'n')
		return 'n';
	else if(c == '-')
		return '-';
	else {
		printf("wrong dna sequence : %c\n", c);
		exit(1);
	}
}

bool check_orthologs(int start, int end) {
	int start_net = -1, end_net = -1, i, total_net_len, orthologs_start, orthologs_end;

	mismatch_net_num = 0;
	for(i=0; i<net_num[0]; i++) {
		if(net[0][i].end < start)
			continue;
		if(net[0][i].start > end)
			break;
	
		if(start_net == -1)
			start_net = i;
		end_net = i;
	}
	
	if(start_net == -1)
		return false;

	//check coverage of all nets
	if(net[0][start_net].start > start)
		orthologs_start = net[0][start_net].start;
	else
		orthologs_start = start;
	if(net[0][end_net].end < end)
		orthologs_end = net[0][end_net].end;
	else
		orthologs_end = end;

	if(start_net == end_net) // only one net
		total_net_len = orthologs_end - orthologs_start + 1;
	else {
		total_net_len = net[0][start_net].end - orthologs_start + 1;
		for(i=start_net+1; i<end_net; i++)
			total_net_len += net[0][i].end - net[0][i].start + 1;
		total_net_len += orthologs_end - net[0][end_net].start + 1;
	}
	if(((float)total_net_len / (float)(end - start + 1)) < MIN_COVERAGE) {
		return false;
	}
	else if(start_net == end_net) {
		return true;
	}

	return true;	
}

//int process_chain(chain_t c, int ref, int num_contigs, int *len_sum) {
int process_chain(chain_t c, int ref) {
	int i, j, k, l, old_aln_end1, old_aln_end2;
	char aln_text_temp[BIG], *aln_text_ptr1, *aln_text_ptr2;
	aln_t *aln;
//	int temp_id = -1;
//	int len_diff = 0;
	
	//check orthologs
	if(((ref == 1) && (check_orthologs(c.start1, c.end1) == false)) || ((ref == 2) && (check_orthologs(c.start2, c.end2) == false)))
		return 0;
	
	if(ref == 1)
		for(i=0; i<c.end1 - c.start1 + 1; i++) {
			seqs[ref-1][0][i] = '-';
			seqs[ref-1][1][i] = '-';
			seqs[ref-1][2][i] = '-';
		}
	else
		for(i=0; i<c.end2 - c.start2 + 1; i++) {
			seqs[ref-1][0][i] = '-';
			seqs[ref-1][1][i] = '-';
			seqs[ref-1][2][i] = '-';
		}

	old_aln_end1 = old_aln_end2 = 0; 
	aln = c.aln;
	while(aln != NULL) {
		// remove overlapped part
		aln_text_ptr1 = &(aln->text1[0]);
		aln_text_ptr2 = &(aln->text2[0]);
		if(aln->start1 <= old_aln_end1) {
			aln->len1 -= old_aln_end1 - aln->start1 + 1;
			i = aln->start1;
			while(i <= old_aln_end1) {
				if(*aln_text_ptr1 != '-')
					i++;
				if(*aln_text_ptr2 != '-') {
					aln->start2++;
					aln->len2--;
				}
				aln_text_ptr1++;
				aln_text_ptr2++;
			}
			aln->start1 = old_aln_end1 + 1;
		}
		if(aln->start2 <= old_aln_end2) {
			aln->len2 -= old_aln_end2 - aln->start2 + 1;
			i = aln->start2;
			while(i <= old_aln_end2) {
				if(*aln_text_ptr2 != '-')
					i++;
				if(*aln_text_ptr1 != '-') {
					aln->start1++;
					aln->len1--;
				}
				aln_text_ptr2++;
				aln_text_ptr1++;
			}
			aln->start2 = old_aln_end2 + 1;
		}
		strcpy(aln->text1, aln_text_ptr1);
		strcpy(aln->text2, aln_text_ptr2);
		old_aln_end1 = aln->start1 + aln->len1 - 1;
		old_aln_end2 = aln->start2 + aln->len2 - 1;

		if(ref == 1) {
			for(i=j=k=0; i<aln->len1; i++) {
				while(aln->text1[j] == '-') {
					j++;
					if(aln->text2[j] != '-')
						k++;
				}
				if(seqs[ref-1][0][i + aln->start1+1 - c.start1] == '-') {
					seqs[ref-1][0][i + aln->start1+1 - c.start1] = aln->text1[j];
					seqs[ref-1][1][i + aln->start1+1 - c.start1] = aln->text2[j];
					if(c.orient == '+')
						seq2_pos[ref-1][i + aln->start1+1 - c.start1] = aln->start2 + 1 + k;
					else
						seq2_pos[ref-1][i + aln->start1+1 - c.start1] = aln->src_len2 - aln->start2 - k;
				}
				j++;
				if(aln->text2[j] != '-')
					k++;
			}
		}
		else {
			//reverse and complement sequences if ref = 2
			if(c.orient == '-') {
				strcpy(aln_text_temp, aln->text1);
				for(i=0; i<(int)strlen(aln_text_temp); i++)
					aln->text1[i] = Complement(aln_text_temp[strlen(aln_text_temp) - i - 1]);
				aln->start1 = aln->src_len1 - (aln->start1 + aln->len1);
				strcpy(aln_text_temp, aln->text2);
				for(i=0; i<(int)strlen(aln_text_temp); i++)
					aln->text2[i] = Complement(aln_text_temp[strlen(aln_text_temp) - i - 1]);
				aln->start2 = aln->src_len2 - (aln->start2 + aln->len2);
			}
			for(i=j=k=0; i<aln->len2; i++) {
				while(aln->text2[j] == '-') {
					j++;
					if(aln->text1[j] != '-')
						k++;
				}
				if(seqs[ref-1][0][i + aln->start2+1 - c.start2] == '-') {
					seqs[ref-1][0][i + aln->start2+1 - c.start2] = aln->text2[j];
					seqs[ref-1][1][i + aln->start2+1 - c.start2] = aln->text1[j];
					if(c.orient == '+')
						seq2_pos[ref-1][i + aln->start2+1 - c.start2] = aln->start1 + 1 + k;
					else 					
						seq2_pos[ref-1][i + aln->start2+1 - c.start2] = aln->src_len1 - aln->start1 - k;
				}
				j++;
				if(aln->text1[j] != '-')
					k++;
			}
		}
		
		aln = aln->next;
	}

	//find orthologous sequence
	for(i=0; i<net_num[0]; i++) {
		if(net[0][i].end < start[ref-1])
			continue;
		if(net[0][i].start > end[ref-1])
			break;
				
		for(j=0; j<mismatch_net_num; j++)
			if(mismatch_net[j] == i)
				break;
		if(j<mismatch_net_num) {
			//printf("%d\n", i);
			continue;
		}
		
		orthologs_block[ref-1][orthologs_block_num[ref-1]] = net[0][i].block_index;
		orthologs_block_num[ref-1]++;
		
		if(if_rescoring && rescored[0][i] == 0) {
			rescoring(i);
			rescored[0][i] = 1;
		}

		strcpy(outgroup_name[ref-1], net[0][i].chr2);
/*		if( strcmp(net[0][i].contig2, net[0][i].chr2) != 0 ) {
			strcat(outgroup_name[ref-1], ".");
			strcat(outgroup_name[ref-1], net[0][i].contig2);
		}
		temp_id = net[0][i].ctg_id2;
		if( ( temp_id < 0 ) || ( temp_id >= num_contigs) ) {
			fatalf("error: %d illegal index of contig list\n", temp_id);
		}
		len_diff = len_sum[temp_id];
*/
		outgroup_orient[ref-1] = net[0][i].orient;
			
		k = l = 0;
		for(j=net[0][i].start; j<=net[0][i].end; j++) {			
			while(net[0][i].net1[k] == '-') {
				k++;
				if(net[0][i].net2[k] != '-')
					l++;
			}
			if(j >= start[ref-1] && j <= end[ref-1]) {
				if(net[0][i].net1[k] != 'X')
					seqs[ref-1][2][j - start[ref-1]] = net[0][i].net2[k];
				if(net[0][i].orient == '+') {
					if((outgroup_start[ref-1] == -1) || ((net[0][i].start2 + l) < outgroup_start[ref-1]))
						outgroup_start[ref-1] = net[0][i].start2 + l;
					if((outgroup_end[ref-1] == -1) || ((net[0][i].start2 + l) > outgroup_end[ref-1]))
						outgroup_end[ref-1] = net[0][i].start2 + l;
				}
				else {
/*
					if((outgroup_start[ref-1] == -1) || ((net[0][i].start2 - l) < outgroup_start[ref-1]))
						outgroup_start[ref-1] = net[0][i].start2 - l;
					if((outgroup_end[ref-1] == -1) || ((net[0][i].start2 - l) > outgroup_end[ref-1]))
						outgroup_end[ref-1] = net[0][i].start2 - l;
*/
					if((outgroup_start[ref-1] == -1) || ((net[0][i].end2 - l) < outgroup_start[ref-1]))
						outgroup_start[ref-1] = net[0][i].end2 - l;
					if((outgroup_end[ref-1] == -1) || ((net[0][i].end2 - l) > outgroup_end[ref-1]))
						outgroup_end[ref-1] = net[0][i].end2 - l;
				}

			}
			else if(j > end[ref-1]) {
				break;
			}
			k++;
			if(net[0][i].net2[k] != '-')
				l++;
		}
	}
	
/*
	if( outgroup_start[ref-1] != -1 ) {
		outgroup_start[ref-1] = outgroup_start[ref-1] - len_diff;
	}
	
	if( outgroup_end[ref-1] != -1 ) {
		outgroup_end[ref-1] = outgroup_end[ref-1] - len_diff;
	}
*/
	if( (outgroup_start[ref-1] < -1) || (outgroup_end[ref-1] < -1) ) {
		fatalf("unexpected negative values: %d %d\n", outgroup_start[ref-1], outgroup_end[ref-1]);
	}

	return 1;
}

void free_aln(aln_t *aln) {
	if(aln->next != NULL)
		free_aln(aln->next);
	
	free(aln);

	return;
}

int aln2_pos(chain_t c, int pos1) {
	aln_t *aln;
	int i, j, pos2;
	
	aln = c.aln;
	while(aln != NULL) {
		if(pos1 < aln->start1 + aln->len1) {
			j = aln->start1 - 1;
			if(c.orient == '+') {
				pos2 = aln->start2 - 1;
				for(i=0; j<pos1; i++) {
					if(aln->text1[i] != '-')
						j++;
					if(aln->text2[i] != '-')
						pos2++;
				}
			}
			else {
				pos2 = aln->src_len2 - aln->start2 + 1;
				for(i=0; j<pos1; i++) {
					if(aln->text1[i] != '-')
						j++;
					if(aln->text2[i] != '-')
						pos2--;
				}
			}
			return pos2;
		}
		else if(pos1 < aln->next->start1) {
			if(c.orient == '+')
				pos2 = aln->start2 + aln->len2 + (pos1 - (aln->start1 + aln->len1));
			else
				pos2 = aln->src_len2 - (aln->start2 + aln->len2) - (pos1 - (aln->start1 + aln->len1));
			return pos2;
		}
		
		aln = aln->next;
	}
	
	return -1;
}

void divided_index(chain_t c) {
	int len;
	
	len = c.end1 - c.start1 + 1;
	if(len <= 5000) {
		d_beg_index[0] = c.start2;
		d_end_index[0] = c.end2;
		d_index_num = 1;
	}
	else {
		d_index_num = 0;
		site_beg[0] = 0;
		site_end[0] = 4999;
		while(site_end[0] < len) {
			d_beg_index[d_index_num] = aln2_pos(c, site_beg[0] + c.start1 - 1);
			d_end_index[d_index_num] = aln2_pos(c, site_end[0] + c.start1 - 1);
			d_index_num++;
			if(site_end[0] == len - 1)
				break;
			else if(site_end[0] + 4000 >= len - 1) {
				site_beg[0] = len - 5000;
				site_end[0] = len - 1  ;
			}
			else {
				site_beg[0] = site_end[0] - 1000;
				site_end[0] = site_end[0] + 4000;
			}
		}
	}
}

float determine_directionality() {
	int md_m1, md_n1, md_m2, md_n2, i;
	float p, p1, p2, V, T;
	
	md_m1 = md_n1 = md_m2 = md_n2 = 0;
	for(i=beg_index+1; i<=end_index; i++) {
		if(IS_ref[i] == 0) {
			if(IS[2][i] == 1)
				md_m1++;
			else if(IS[2][i] == 2)
				md_n1++;
		}
		else if(IS_ref[i] == 1) {
			if(IS[2][i] == 1)
				md_m2++;
			else if(IS[2][i] == 2)
				md_n2++;
		}
	}

	p = (float)(md_n1 + md_n2) / (float)(md_m1 + md_m2 + md_n1 + md_n2);
	p1 = (float)md_n1 / (float)(md_m1 + md_n1);
	p2 = (float)md_n2 / (float)(md_m2 + md_n2);
	V = ((1.0/(float)(md_m1 + md_n1)) + (1.0/(float)(md_m2 + md_n2))) * p * (1.0 - p);
	if(V == 0.0 || (md_m1 + md_n1) == 0 || (md_m2 + md_n2) == 0)
		T = 0.0;
	else
		T = ((float)(p1 - p2)) / sqrtf(V);
	
	//printf("\n\nm1 = %d; n1 = %d; m2 = %d; n2 = %d; V = %f; T = %f\n\n", md_m1, md_n1, md_m2, md_n2, V, T);
	
	return T;
}

void print_block_index() {
	int i, j;

	for(i=0; i<2; i++) {
		if(orthologs_block_num[i] == 0)
			printf("\tNAN");
		else {
			printf("\t");
			for(j=0; j<orthologs_block_num[i]-1; j++) {
				printf("%d,", orthologs_block[i][j]);
			}
			printf("%d", orthologs_block[i][j]);
		}
	}
}

int check_orthologs_overlap(chain_t c) {
	int i, j, k, i_beg[10000], i_end[10000], i_num = 0, i_beg2, i_end2, o_beg, o_end;

	//printf("\n\n");
	//get all orthologs for paralogs 1
	for(i=0; i<net_num[1]; i++) {
		if(c.start1 <= net[1][i].end && c.end1 >=net[1][i].start) {			
			if(c.start1 > net[1][i].start) {
				i_beg[i_num] = net[1][i].start2 - 1;
				for(j=k=0; k<=c.start1 - net[1][i].start; j++) {
					if(net[1][i].net1[j] != '-')
						k++;
					if(net[1][i].net2[j] != '-')
						i_beg[i_num]++;
				}
			}
			else 
				i_beg[i_num] = net[1][i].start2;

			if(c.end1 < net[1][i].end) {
				i_end[i_num] = net[1][i].end2 + 1;
				j=strlen(net[1][i].net1) - 1;	
				for(k=0; k<=net[1][i].end - c.end1; j--) {
					if(net[1][i].net1[j] != '-')
						k++;
					if(net[1][i].net2[j] != '-')
						i_end[i_num]--;
				}
			}
			else 
				i_end[i_num] = net[1][i].end2;

			if(i_beg[i_num] > i_end[i_num]) {
				printf("Wrong many-to-many orthologs file: %d\t%d\n", i_beg[i_num], i_end[i_num]);
				exit(1);
			}
			//printf("[%d, %d]\t", i_beg[i_num], i_end[i_num]);

			i_num++;
		}
	}
	//printf("\n");
	
	//check if any orthologs for paralogs 2 overlaps with orthologs for paralogs 1
	for(i=0; i<net_num[1]; i++) {
		if(c.start2 <= net[1][i].end && c.end2 >=net[1][i].start) {			
			if(c.start2 > net[1][i].start) {
				i_beg2 = net[1][i].start2 - 1;
				for(j=k=0; k<=c.start2 - net[1][i].start; j++) {
					if(net[1][i].net1[j] != '-')
						k++;
					if(net[1][i].net2[j] != '-')
						i_beg2++;
				}
			}
			else 
				i_beg2 = net[1][i].start2;

			if(c.end2 < net[1][i].end) {
				i_end2 = net[1][i].end2 + 1;
				j=strlen(net[1][i].net1) - 1;	
				for(k=0; k<=net[1][i].end - c.end2; j--) {
					if(net[1][i].net1[j] != '-')
						k++;
					if(net[1][i].net2[j] != '-')
						i_end2--;
				}
			}
			else 
				i_end2 = net[1][i].end2;

			if(i_beg2 > i_end2) {
				printf("Wrong many-to-many orthologs file: %d\t%d\n", i_beg2, i_end2);
				exit(1);
			}

			//printf("[%d, %d]\n\n", i_beg2, i_end2);
			//check overlapping
			for(j=0; j<i_num; j++) {
				if(i_beg2 <= i_end[j] && i_end2 >= i_beg[j]) {
					if(i_beg2 >= i_beg[j])
						o_beg = i_beg2;
					else
						o_beg = i_beg[j];

					if(i_end2 <= i_end[j])
						o_end = i_end2;
					else
						o_end = i_end[j];

					if(((float)(o_end - o_beg + 1) / (float)(i_end2 - i_beg2 + 1) > ORTHOLOGS_OVERLAP_RATIO) && ((float)(o_end - o_beg + 1) / (float)(i_end[j] - i_beg[j] + 1) > ORTHOLOGS_OVERLAP_RATIO))
						return 0;
				}
			}
		}
	}

	return 1;
}

int read_chain(char *chain_fn, struct n_pair **contigs1, int *num_contigs1, int *num_blocks1, int **len_sum1) {
	FILE *fp = NULL;
	char *not_eof = NULL;
	chain_t c;
	int index = 0;
	aln_t **aln_ptr = NULL, *aln = NULL;
	int i = 0, match = 0, mismatch = 0, triplet_status = 0, aln_len = 0;
	float d_conf = (float)0;
	char name1[LEN_NAME] = "", name2[LEN_NAME] = "";
	int num1 = 0;
	int count = 0;
	int temp_id1 = 0, temp_id2 = 0;

	num1 = *num_contigs1;
	fp = fopen(chain_fn, "r");
	while((not_eof = fgets(buf, BIG, fp)) && buf[0] == '#')
		;
	while(not_eof) {
		sscanf(buf, "%d %s %d %d %s %d %d %c", &index, name1, &c.start1, &c.end1, name2, &c.start2, &c.end2, &c.orient);

		concat_ctg_name(name2, c.chr2, c.contig2);
		// test length
		if((c.end1 - c.start1 + 1) > BIG || (c.end2 - c.start2 + 1) > BIG) {
			printf("chain is too long : %d - %d\n", c.end1 - c.start1 + 1, c.end2 - c.start2 + 1);
			exit(1);
		}
		not_eof = fgets(buf, BIG, fp);
		aln_ptr = &c.aln;
		count = 0;
		while((not_eof) && (buf[0] == 'a')) {
			*aln_ptr = (aln_t *)ckalloc(sizeof(aln_t));
			not_eof = fgets(buf, BIG, fp);
			sscanf(buf, "s %s %d %d %c %d %s", name1, &(*aln_ptr)->start1, &(*aln_ptr)->len1, &c.orient, &(*aln_ptr)->src_len1, (*aln_ptr)->text1);
			
			not_eof = fgets(buf, BIG, fp);
			sscanf(buf, "s %s %d %d %c %d %s", name2, &(*aln_ptr)->start2, &(*aln_ptr)->len2, &c.orient, &(*aln_ptr)->src_len2, (*aln_ptr)->text2);

			if( count == 0 ) {
				concat_ctg_name(name1, c.chr1, c.contig1);
				temp_id1 = add_ctg_name(index, c.chr1, c.contig1, (*aln_ptr)->src_len1, *contigs1, num1);
				c.ctg_id1 = temp_id1;
				if( temp_id1 == num1 ) {
					num1++;

					if( num1 == 1 ) {
						*len_sum1 = (int *) ckalloc(sizeof(int));
						(*len_sum1)[0] = 0;
					}
					else {
						*len_sum1 = (int *) ckrealloc(*len_sum1, num1 * sizeof(int));
						(*len_sum1)[num1-1] = (*len_sum1)[num1-2] + (*contigs1)[num1-2].len;
					}
				}
				else if( temp_id1 > num1 ) {
					fatalf("warning: %d should not excceed %d\n", temp_id1, num1);
				}

				if( num1 >= ((*num_blocks1) * ALLOC_UNIT) ) {
					*num_blocks1 = *num_blocks1 + 1;
					*contigs1 = (struct n_pair *) ckrealloc(*contigs1, sizeof(struct n_pair) * (*num_blocks1) * ALLOC_UNIT);
				}

				c.start1 = c.start1 + (*len_sum1)[temp_id1];
				c.end1 = c.end1 + (*len_sum1)[temp_id1];

				concat_ctg_name(name2, c.chr2, c.contig2);
				temp_id2 = add_ctg_name(index, c.chr2, c.contig2, (*aln_ptr)->src_len2, *contigs1, num1);
				c.ctg_id2 = temp_id2;
				if( temp_id2 == num1 ) {
					num1++;

					if( num1 == 1 ) {
						*len_sum1 = (int *) ckalloc(sizeof(int));
						(*len_sum1)[0] = 0;
					}
					else { 
						*len_sum1 = (int *) ckrealloc(*len_sum1, num1 * sizeof(int));
						(*len_sum1)[num1-1] = (*len_sum1)[num1-2] + (*contigs1)[num1-2].len;
					}
				}
				else if( temp_id2 > num1 ) {
					fatalf("warning: %d should not excceed %d\n", temp_id1, num1);
				}

				if( num1 >= ((*num_blocks1) * ALLOC_UNIT) ) {
					*num_blocks1 = *num_blocks1 + 1;
					*contigs1 = (struct n_pair *) ckrealloc(*contigs1, sizeof(struct n_pair) * (*num_blocks1) * ALLOC_UNIT);
					init_pair(*contigs1, ((*num_blocks1)-1) * ALLOC_UNIT, (*num_blocks1) * ALLOC_UNIT);
				}
				c.start2 = c.start2 + (*len_sum1)[temp_id2];
				c.end2 = c.end2 + (*len_sum1)[temp_id2];
			}

			(*aln_ptr)->start1 = (*aln_ptr)->start1 + (*len_sum1)[temp_id1];
			(*aln_ptr)->start2 = (*aln_ptr)->start2 + (*len_sum1)[temp_id2];
			(*aln_ptr)->ctg_id1 = temp_id1;
			(*aln_ptr)->ctg_id2 = temp_id2;
			count++;

			(*aln_ptr)->next = NULL;
			aln_ptr = &((*aln_ptr)->next);
			not_eof = fgets(buf, BIG, fp);
			not_eof = fgets(buf, BIG, fp);
		}
		
		start[0] = c.start1;
		end[0] = c.end1;
		sites_num[0] = c.end1 - c.start1 + 1;
		start[1] = c.start2;
		end[1] = c.end2;
		sites_num[1] = c.end2 - c.start2 + 1;
		
		//calculate identity
		aln = c.aln;
		match = mismatch = 0;
		while(aln != NULL) {
			aln_len = strlen(aln->text1);
			for(i=0; i<aln_len; i++) {
				if(toupper(aln->text1[i]) == toupper(aln->text2[i]))
					match++;
				else
					mismatch++;
			}
			aln = aln->next;
		}
		identity = (double)match / (double)(match + mismatch);
		
		divided_index(c);

		site_beg[0] = 0;
		site_beg[1] = 0;
		long_seqs = 0;
		outgroup_start[0] = outgroup_start[1] = outgroup_end[0] = outgroup_end[1] = -1;
		strcpy(outgroup_name[0], "NAN");
		strcpy(outgroup_name[1], "NAN");
		outgroup_orient[0] = outgroup_orient[1] = 'N';
		orthologs_block_num[0] = orthologs_block_num[1] = 0;
		if(process_chain(c, 1) == 1) {	 
			if(process_chain(c, 2) == 1) 
				triplet_status = 3;
			else 
				triplet_status = 1;
		}		
		else if(process_chain(c, 2) == 1) 		
			triplet_status = 2;
		else
			triplet_status = 0;
		
		// test conversion ratio
		site_end[0] = sites_num[0];
		if(triplet_status == 1 || triplet_status == 3) {
			find_IS(1);
			find_max_descent(1);
		}
		else {
			IS_num[0] = 0;
			m[0] = 0;
			n[0] = 0;
			max_descent[0] = 0;
		}
		site_end[1] = sites_num[1];
		if(triplet_status == 2 || triplet_status == 3) {
			find_IS(2);
			find_max_descent(2);
		}
		else {
			IS_num[1] = 0;
			m[1] = 0;
			n[1] = 0;
			max_descent[1] = 0;
		}
		d_conf = 0.0;
		if(triplet_status != 0) {
			combine_IS(c.orient);
			find_max_descent(3);
			if(triplet_status == 3)
				d_conf = determine_directionality();
		}

/*
		if( strcmp(c.chr1, c.contig1) != 0 ) {
			strcat(c.chr1, ".");
			strcat(c.chr1, c.contig1);
		}
		c.start1 = c.start1 - (*len_sum1)[temp_id1];
		c.end1 = c.end1 - (*len_sum1)[temp_id1];

		if( strcmp(c.chr2, c.contig2) != 0 ) {
			strcat(c.chr2, ".");
			strcat(c.chr2, c.contig2);
		}
		c.start2 = c.start2 - (*len_sum1)[temp_id2];
		c.end2 = c.end2 - (*len_sum1)[temp_id2];
*/

		if(triplet_status == 0 || (float)(abs(max_end1 - max_beg1) + 1) / (float)sites_num[0] > ((float)max_gc_ratio)) {
			if((triplet_status == 3) && (n[0] + n[1] >= m[0] + m[1]) && check_orthologs_overlap(c) == 1) {
				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c", index, c.chr1, c.start1, c.end1, c.chr2, c.start2, c.end2, c.orient, sites_num[0], 0, 0, 0, 0, 0, 0, m[0] + m[1], n[0] + n[1], max_descent[2], identity, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1]);
				print_block_index();
				printf("\t4\n");
				for(i=5000; i<sites_num[0]; i+=4000) {
					printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c", index, c.chr1, c.start1, c.end1, c.chr2, c.start2, c.end2, c.orient, sites_num[0], 0, 0, 0, 0, 0, 0, m[0] + m[1], n[0] + n[1], max_descent[2], identity, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1]);
					print_block_index();
					printf("\t4\n");
				}
			}
			else {
				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c", index, c.chr1, c.start1, c.end1, c.chr2, c.start2, c.end2, c.orient, sites_num[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, identity, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1]);			// add block index
				print_block_index();
				printf("\t0\n");
				for(i=5000; i<sites_num[0]; i+=4000) {
					printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c", index, c.chr1, c.start1, c.end1, c.chr2, c.start2, c.end2, c.orient, sites_num[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, identity, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1]);
					print_block_index();
					printf("\t0\n");
				}
			}
		}
		else if(sites_num[0] <= 5000) {
			printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c", index, c.chr1, c.start1, c.end1, c.chr2, c.start2, c.end2, c.orient, sites_num[0], m[0], n[0], max_descent[0], m[1], n[1], max_descent[1], m[0] + m[1], n[0] + n[1], max_descent[2], identity, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1]);
			print_block_index();
			printf("\t%d\n", triplet_status);
		}
		else {
			long_seqs = 1;
			site_end[0] = 5000;
			i = 0;
			while(site_end[0] <= sites_num[0]) {
				if(triplet_status == 1 || triplet_status == 3) {
					find_IS(1);
					find_max_descent(1);
				}
				else {
					IS_num[0] = 0;
					m[0] = 0;
					n[0] = 0;
					max_descent[0] = 0;
				}
				if(triplet_status == 2 || triplet_status == 3) {
					if(c.orient == '+') {
						site_beg[1] = d_beg_index[i] - c.start2;
						site_end[1] = d_end_index[i] - c.start2;
					}
					else {
						site_beg[1] = d_end_index[i] - c.start2;
						site_end[1] = d_beg_index[i] - c.start2;
					}
					find_IS(2);
					find_max_descent(2);
				}
				else {
					IS_num[1] = 0;
					m[1] = 0;
					n[1] = 0;
					max_descent[1] = 0;
				}
				combine_IS(c.orient);
				find_max_descent(3);
				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c", index, c.chr1, c.start1, c.end1, c.chr2, c.start2, c.end2, c.orient, sites_num[0], m[0], n[0], max_descent[0], m[1], n[1], max_descent[1], m[0] + m[1], n[0] + n[1], max_descent[2], identity, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1]);
				print_block_index();
				printf("\t%d\n", triplet_status);
				
				if(site_end[0] == sites_num[0])
					break;
				else if(site_end[0] + 4000 >= sites_num[0]) {
					site_beg[0] = sites_num[0] - 5000;
					site_end[0] = sites_num[0];
				}
				else {
					site_beg[0] = site_end[0] - 1000;
					site_end[0] = site_end[0] + 4000;
				}
				i++;
			}
		}	
		
		//free aln
		free_aln(c.aln);
	}

	fclose(fp);
	
	return 0;
}

void init_contigs_list(FILE *ctg_f, struct n_pair **contigs, int *num_contigs, int *num_blocks, int **len_sum)
{
	int i = 0, j = 0;
	char name1[1000], name2[1000];
	int len1 = 0, len2 = 0;

	while (fgets(buf, 1000, ctg_f)) i++;
	if( i > 0 ) {
		*num_contigs = i;
		*num_blocks = (i / ALLOC_UNIT) + 1;
		*contigs = (struct n_pair *) ckalloc(sizeof(struct n_pair) * (*num_blocks) * ALLOC_UNIT);
		init_pair(*contigs, 0, (*num_blocks) * ALLOC_UNIT);
		*len_sum = (int *) ckalloc(sizeof(int) * (*num_contigs));
		for( j = 0; j < (*num_contigs); j++ ) (*len_sum)[j] = 0;
		fseek(ctg_f, 0, SEEK_SET);
		i = 0;
		while( fgets(buf, 1000, ctg_f) ) {
			sscanf(buf, "%s %s %d %d", name1, name2, &len1, &len2);

			if( i < (*num_contigs) ) {
				strcpy( (*contigs)[i].name1, name1);
				strcpy( (*contigs)[i].name2, name2);
				(*contigs)[i].len = len1;
				(*len_sum)[i] = len2;
			}
			else {
				fatalf("counting error: (md_quadruplet_entire_region.c) %d\n", i);
			}
			i++;
		}
	}
	else {
		*contigs = (struct n_pair *) ckalloc(sizeof(struct n_pair) * (*num_blocks) * ALLOC_UNIT);
		init_pair(*contigs, 0, (*num_blocks) * ALLOC_UNIT);
		*num_contigs = 0;
	}
}

void print_contig_list(FILE *ctg_f, struct n_pair *contigs, int *len_sum, int num_org, int num)
{
	int i = 0;

	for( i = num_org; i < num; i++ )
	{
		fprintf(ctg_f, "%s %s %d %d\n", contigs[i].name1, contigs[i].name2, contigs[i].len, len_sum[i]);
	}
}

int main(int argc, char *argv[]) {
	struct n_pair *contigs1, *contigs2;
	int *num_contigs1 = NULL, *num_contigs2 = NULL;
	int *num_blocks1 = NULL, *num_blocks2 = NULL;
	int *len_sum1, *len_sum2;
	FILE *ctg_f1, *ctg_f2;
	int num_org_contigs1 = 0, num_org_contigs2 = 0;

	if ( (argc != 6) && (argc != 7)) {
		printf("md_quadruplet_entire_region chain-file many-to-one-orthologs_file many-to-many-orthologs_file contigs_list <max_gc_ratio> \n");
		return 1;
	}

	num_contigs1 = (int *) ckalloc(sizeof(int));
	num_contigs2 = (int *) ckalloc(sizeof(int));
	num_blocks1 = (int *) ckalloc(sizeof(int));
	num_blocks2 = (int *) ckalloc(sizeof(int));
	*num_contigs1 = 0;
	*num_contigs2 = 0;
	*num_blocks1 = 1;
	*num_blocks2 = 1;

	ctg_f1 = ckopen(argv[4], "a+");
	ctg_f2 = ckopen(argv[5], "a+");
	
	init_contigs_list(ctg_f1, &contigs1, num_contigs1, num_blocks1, &len_sum1);
	num_org_contigs1 = *num_contigs1;
	init_contigs_list(ctg_f2, &contigs2, num_contigs2, num_blocks2, &len_sum2);
	num_org_contigs2 = *num_contigs2;

/*
	// read input file
	if(argc == 5) {
		read_score(argv[4]);
		if_rescoring = 1;
	}
*/

	max_gc_ratio = MAX_GC_RATIO;
	if(argc == 7) {
		max_gc_ratio = atof(argv[6]);
	}

	read_orthologs(argv[2], 0, &contigs1, num_contigs1, num_blocks1, &len_sum1, &contigs2, num_contigs2, num_blocks2, &len_sum2);
	read_orthologs(argv[3], 1, &contigs1, num_contigs1, num_blocks1, &len_sum1, &contigs2, num_contigs2, num_blocks2, &len_sum2);
	read_chain(argv[1], &contigs1, num_contigs1, num_blocks1, &len_sum1);
	//printf("d_0 = %d, d_1 = %d, d_2 = %d, d_3 = %d\n", d_0, d_1, d_2, d_3);

	print_contig_list(ctg_f1, contigs1, len_sum1, num_org_contigs1, *num_contigs1);
	print_contig_list(ctg_f2, contigs2, len_sum2, num_org_contigs2, *num_contigs2);
	
	fclose(ctg_f1);
	fclose(ctg_f2);

	free(contigs1);
	free(contigs2);
	if( (*num_contigs1) > 0 ) free(len_sum1);
	if( (*num_contigs2) > 0 ) free(len_sum2);
	free(num_contigs1);
	free(num_contigs2);
	free(num_blocks1);
	free(num_blocks2);
	return 0;
}
