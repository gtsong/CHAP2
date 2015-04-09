/* 3seq_2D_cluster_fast - calculate p-value for each triplet test by looking up p_value_table. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_DIM 1501

typedef struct P_Value_Table {
	int k;
	float p_value;
	struct P_Value_Table *next;
} p_value_table_t;

p_value_table_t *p_value_table[MAX_DIM][MAX_DIM];

float get_pvalue(int m, int n, int k) {
	p_value_table_t *p_value_table_ptr;
	
	if(p_value_table[m][n] == NULL)
		return 1.0;

	p_value_table_ptr = p_value_table[m][n];
	if(k < p_value_table_ptr->k)
		return 1.0;
	while(p_value_table_ptr != NULL && k > p_value_table_ptr->k)
		p_value_table_ptr = p_value_table_ptr->next;

	if(p_value_table_ptr == NULL)
		return 0.0;
	else
		return p_value_table_ptr->p_value;
}

float get_pvalue_entire_region(int m, int n) {
	int i, j, k;
	float p_value = 0.0, p;

	for(i=n; i<=n+m; i++) {
		p = 1.0;
		k = 1;
		//avoid overflow
		for(j=i+1; j<=n+m; j++) {
			p = p * ((float)j / (float)(j - i));
			while((p > 10000.0) && (k <= n + m)) {
				p = p * 0.5;
				k++;
			}
		}
		for(; k<=n+m; k++)
			p = p * 0.5;

		p_value += p;
	}

	return p_value;
}

int main(int argc, char *argv[]) {
	FILE *fp;
	char buf[5000], chr1[100], chr2[100], orient;
	int index, b1, e1, b2, e2, length, m1, n1, k1, m2, n2, k2, m3, n3, k3, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2;
	char identity[100], d_conf[100], outgroup_name[2][100], outgroup_orient[2], orthologs_block[2][100], pvalue[100];
	int i, j, outgroup_start[2], outgroup_end[2], m, n, k_num, k, triplet_status;
	p_value_table_t *p_value_table_ptr = NULL;
	
	if ( argc != 3) {
		printf("3seq_2D_cluster_fast md-file p_value_table\n");
		return 1;
	}

	for(i=0; i<MAX_DIM; i++)
		for(j=0; j<MAX_DIM; j++)
			p_value_table[i][j] = NULL;

	//read p_value_table file
	fp = fopen(argv[2], "r");
	while(fscanf(fp, "%d %d %d %d", &m, &n, &k_num, &k) != EOF) {
		for(i=0; i<k_num; i++) {
			if(fscanf(fp, "%s", pvalue) != 1) {
				printf("wrong p_value_table file\n");
				exit(1);
			}
			if(p_value_table[m][n] == NULL) {
				p_value_table_ptr = p_value_table[m][n] = (p_value_table_t *)malloc(sizeof(p_value_table_t));
				
			}
			else {
				p_value_table_ptr->next = (p_value_table_t *)malloc(sizeof(p_value_table_t));
				p_value_table_ptr = p_value_table_ptr->next;
			}
			p_value_table_ptr->k = k;
			p_value_table_ptr->p_value = atof(pvalue);
			p_value_table_ptr->next = NULL;
			k++;
		}
	}

	fp = fopen(argv[1], "r");
	while(fgets(buf, 5000, fp)) {
		if(buf[0] == '#') {
			printf("%s", buf);
			continue;
		}
		sscanf(buf, "%d %s %d %d %s %d %d %c %d %d %d %d %d %d %d %d %d %d %s %d %d %d %d %d %d %d %d %d %d %s %s %d %d %c %s %d %d %c %s %s %d", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, &m1, &n1, &k1, &m2, &n2, &k2, &m3, &n3, &k3, identity, &min_len, &max_len, &max_beg1, &min_beg1, &min_end1, &max_end1, &max_beg2, &min_beg2, &min_end2, &max_end2, d_conf, outgroup_name[0], &outgroup_start[0], &outgroup_end[0], &outgroup_orient[0], outgroup_name[1], &outgroup_start[1], &outgroup_end[1], &outgroup_orient[1], orthologs_block[0], orthologs_block[1], &triplet_status);
		printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t", index, chr1, b1, e1, chr2, b2, e2, orient, length);
		if(m1<MAX_DIM && n1<MAX_DIM && k1<MAX_DIM) {
			if(m1 > 0 && n1 > 0)
				printf("%d\t%d\t%d\t%e\t", m1, n1, k1, get_pvalue(m1, n1, k1));
			else
				printf("%d\t%d\t%d\t%e\t", m1, n1, k1, 1.0);
		}
		else
			printf("%d\t%d\t%d\tunknown\t", m1, n1, k1);
		if(m2<MAX_DIM && n2<MAX_DIM && k2<MAX_DIM) {
			if(m2 > 0 && n2 > 0)
				printf("%d\t%d\t%d\t%e\t", m2, n2, k2, get_pvalue(m2, n2, k2));
			else
				printf("%d\t%d\t%d\t%e\t", m2, n2, k2, 1.0);
		}
		else
			printf("%d\t%d\t%d\tunknown\t", m2, n2, k2);
		if(m3<MAX_DIM && n3<MAX_DIM && k3<MAX_DIM) {
			if(m3 > 0 && n3 > 0) {
				if(m1 == 0 && n1 == 0 && m2 == 0 && n2 == 0) { //for entire region
					printf("%d\t%d\t%d\t%e\t", m3, n3, k3, get_pvalue_entire_region(m3, n3));
				}
				else
					printf("%d\t%d\t%d\t%e\t", m3, n3, k3, get_pvalue(m3, n3, k3));
			}
			else
				printf("%d\t%d\t%d\t%e\t", m3, n3, k3, 1.0);
		}
		else
			printf("%d\t%d\t%d\tunknown\t", m3, n3, k3);
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%s\t%d\n", identity, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1], orthologs_block[0], orthologs_block[1], triplet_status);
	}
	
	fclose(fp);

	return 0;
}
