/* 3seq_2D_cluster - calculate p-value for each triplet test. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_DIM 2000

typedef struct Data_Info {
	int k;
	float pvalue;
	struct Data_Info *next;
} data_info_t;

data_info_t *data_info_table[MAX_DIM][MAX_DIM];
int MAX_M_DIM = 1, MAX_N_DIM = 1;

float A[MAX_DIM][MAX_DIM][2], Y[2][MAX_DIM][MAX_DIM], X[2][MAX_DIM][MAX_DIM];

/* Insert a triplet test into the 3-D linked list. */
void insert_data(int m, int n, int k) {
	data_info_t *(*di_ptr), *next_ptr;
	
	di_ptr = &data_info_table[m][n];
	if((*di_ptr) != NULL) {
		while((*di_ptr) != NULL && k > ((*di_ptr))->k)
			di_ptr = &((*di_ptr)->next);
		if((*di_ptr) != NULL) {
			if(k == (*di_ptr)->k) 
				return;
			else {
				next_ptr = (*di_ptr);
				(*di_ptr) = (data_info_t *)malloc(sizeof(data_info_t));
				(*di_ptr)->k = k;
				(*di_ptr)->pvalue = 0.0;
				(*di_ptr)->next = next_ptr;
			}
		}
		else { //insert in the end of linked list
			(*di_ptr) = (data_info_t *)malloc(sizeof(data_info_t));
			(*di_ptr)->k = k;
			(*di_ptr)->pvalue = 0.0;
			(*di_ptr)->next = NULL;
		}
	}
	else {
		(*di_ptr) = (data_info_t *)malloc(sizeof(data_info_t));
		(*di_ptr)->k = k;
		(*di_ptr)->pvalue = 0.0;
		(*di_ptr)->next = NULL;
	}
}

void free_data_info(data_info_t *data_info) {
	if(data_info->next != NULL)
		free_data_info(data_info->next);
	free(data_info);
}


/* Release the memory allocation of the 3-D linked list. */
void free_data_info_table() {
	int i, j;
	
	for(i=0; i<MAX_DIM; i++)
		for(j=0; j<MAX_DIM; j++)
			if(data_info_table[i][j] != NULL)
				free_data_info(data_info_table[i][j]);
}

/* Calculate p-values for all combination of m, n and k. */
void pvalue_table() {
	int m, n, k, j, pct = -1;
	data_info_t (*di_ptr);

	for(m=0; m<MAX_M_DIM; m++) {
		A[m][0][0] = 1.0;
		for(n=1; n<MAX_N_DIM; n++)
			A[m][n][0] = 0.0;
	}

	for(n=0; n<MAX_N_DIM; n++)
		for(k=0; k<MAX_N_DIM; k++)
			if(n==k)
				X[0][n][k] = 1.0;
			else
				X[0][n][k] = 0.0;

	for(k=1; k<MAX_N_DIM; k++) {
		if((k * 1000 / MAX_N_DIM) > pct) {
			pct = (k * 1000 / MAX_N_DIM);
			fprintf(stderr, "calculated p-values : %3.1f%%\r", (float)pct / 10.0);
		}
		//boundary conditions			
		for(n=0; n<MAX_N_DIM; n++)
			for(j=0; j<MAX_N_DIM; j++)
				if(k==n && j==n)
					Y[0][n][j] = 1.0;
				else
					Y[0][n][j] = 0.0;
				
		for(m=1; m<MAX_M_DIM; m++) {			
			for(j=0; j<MAX_N_DIM; j++)
				Y[1][0][j] = 0.0;

			X[1][0][0] = 1.0;
			X[1][0][k] = 0.0;
				
			for(n=1; n<MAX_N_DIM; n++) {
				for(j=0; j<MAX_N_DIM; j++)
					if(k>n || k<(n-m))
						Y[1][n][j] = 0.0;
					else if(j>k || j>n || j<(n-m))
						Y[1][n][j] = 0.0;
					else if (j == 0)
						Y[1][n][j] = ((float)m / (float)(m + n)) * (Y[0][n][1] + Y[0][n][0]);
					else
						Y[1][n][j] = ((float)m / (float)(m + n)) * Y[0][n][j+1] + ((float)n / (float)(m + n)) * Y[1][n-1][j-1];
										
				if(k>n || k<(n-m))
					A[m][n][1] = 0.0;
				else
					A[m][n][1] = ((float)n / (float)(m + n)) * (A[m][n-1][0] + Y[1][n-1][k-1]);
				Y[1][n][k] = A[m][n][1];
								
				X[1][n][k] = ((float)m / (float)(m + n)) * X[0][n][k] + ((float)n / (float)(m + n)) * (X[1][n-1][k] + A[m][n-1][0] - A[m][n-1][1]);
				if(k <= n) {
					di_ptr = data_info_table[m][n];
					while(di_ptr != NULL && di_ptr->k <= k) {
						di_ptr->pvalue += X[1][n][k];
						di_ptr = di_ptr->next;
					}
				}
			}
			
			for(n=0; n<MAX_N_DIM; n++) {
				A[m][n][0] = A[m][n][1];
				X[0][n][k] = X[1][n][k];
				for(j=0; j<MAX_N_DIM; j++)
					Y[0][n][j] = Y[1][n][j];
			}
			
		}
	}
	
	return;
}

/* Get p-value from the 3-D linked list for a particular triplet test. */
float get_pvalue(int m, int n, int k) {
	data_info_t (*di_ptr);
	
	di_ptr = data_info_table[m][n];
	while(di_ptr != NULL && di_ptr->k != k)
		di_ptr = di_ptr->next;
		
	if(di_ptr == NULL) {
		perror("Wrong pvalue!\n");
		exit(1);
	}
	else
		return di_ptr->pvalue;
}

int main(int argc, char *argv[]) {
	FILE *fp;
	char buf[5000], chr1[100], chr2[100], orient;
	int index, b1, e1, b2, e2, length, m1, n1, k1, m2, n2, k2, m3, n3, k3, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2;
	char identity[100], d_conf[100], outgroup_name[2][100], outgroup_orient[2], orthologs_block[2][100];
	int i, j, outgroup_start[2], outgroup_end[2];
	
	if ( argc != 2) {
		printf("3seq_2D_cluster MD-file\n");
		return 1;
	}

	//initialization of Hash table
	for(i=0; i<MAX_DIM; i++)
		for(j=0; j<MAX_DIM; j++)
			data_info_table[i][j] = NULL;

	fp = fopen(argv[1], "r");
	while(fgets(buf, 5000, fp)) {
		if(buf[0] == '#')
			continue;
		sscanf(buf, "%d %s %d %d %s %d %d %c %d %d %d %d %d %d %d %d %d %d %s %d %d %d %d %d %d %d %d %d %d %s %s %d %d %c %s %d %d %c %s %s", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, &m1, &n1, &k1, &m2, &n2, &k2, &m3, &n3, &k3, identity, &min_len, &max_len, &max_beg1, &min_beg1, &min_end1, &max_end1, &max_beg2, &min_beg2, &min_end2, &max_end2, d_conf, outgroup_name[0], &outgroup_start[0], &outgroup_end[0], &outgroup_orient[0], outgroup_name[1], &outgroup_start[1], &outgroup_end[1], &outgroup_orient[1], orthologs_block[0], orthologs_block[1]);
		if(m1 > 0 && n1 > 0) {
			if(m1 >= MAX_M_DIM)
				MAX_M_DIM = m1 + 1;
			if(n1 >= MAX_N_DIM)
				MAX_N_DIM = n1 + 1;
				
			if(m1 < MAX_DIM && n1 < MAX_DIM)
				insert_data(m1, n1, k1);
		}
		if(m2 > 0 && n2 > 0) {
			if(m2 >= MAX_M_DIM)
				MAX_M_DIM = m2 + 1;
			if(n2 >= MAX_N_DIM)
				MAX_N_DIM = n2 + 1;
				
			if(m2 < MAX_DIM && n2 < MAX_DIM)
				insert_data(m2, n2, k2);
		}
		if(m3 > 0 && n3 > 0) {
			if(m3 >= MAX_M_DIM)
				MAX_M_DIM = m3 + 1;
			if(n3 >= MAX_N_DIM)
				MAX_N_DIM = n3 + 1;
				
			if(m3 < MAX_DIM && n3 < MAX_DIM)
				insert_data(m3, n3, k3);
		}
	}

	fprintf(stderr, "MAX_M_DIM = %d; MAX_N_DIM = %d\n", MAX_M_DIM, MAX_N_DIM);
		
	if(MAX_M_DIM > MAX_DIM)
		MAX_M_DIM = MAX_DIM;
	if(MAX_N_DIM > MAX_DIM)
		MAX_N_DIM = MAX_DIM;
	
	pvalue_table();
	
	fseek(fp, 0, SEEK_SET);
	while(fgets(buf, 5000, fp)) {
		if(buf[0] == '#') {
			printf("%s", buf);
			continue;
		}
		sscanf(buf, "%d %s %d %d %s %d %d %c %d %d %d %d %d %d %d %d %d %d %s %d %d %d %d %d %d %d %d %d %d %s %s %d %d %c %s %d %d %c %s %s", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, &m1, &n1, &k1, &m2, &n2, &k2, &m3, &n3, &k3, identity, &min_len, &max_len, &max_beg1, &min_beg1, &min_end1, &max_end1, &max_beg2, &min_beg2, &min_end2, &max_end2, d_conf, outgroup_name[0], &outgroup_start[0], &outgroup_end[0], &outgroup_orient[0], outgroup_name[1], &outgroup_start[1], &outgroup_end[1], &outgroup_orient[1], orthologs_block[0], orthologs_block[1]);
		printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t", index, chr1, b1, e1, chr2, b2, e2, orient, length);
		if(m1<MAX_M_DIM && n1<MAX_N_DIM && k1<MAX_N_DIM) {
			if(m1 > 0 && n1 > 0)
				printf("%d\t%d\t%d\t%e\t", m1, n1, k1, get_pvalue(m1, n1, k1));
			else
				printf("%d\t%d\t%d\t%e\t", m1, n1, k1, 1.0);
		}
		else
			printf("%d\t%d\t%d\tunknown\t", m1, n1, k1);
		if(m2<MAX_M_DIM && n2<MAX_N_DIM && k2<MAX_N_DIM) {
			if(m2 > 0 && n2 > 0)
				printf("%d\t%d\t%d\t%e\t", m2, n2, k2, get_pvalue(m2, n2, k2));
			else
				printf("%d\t%d\t%d\t%e\t", m2, n2, k2, 1.0);
		}
		else
			printf("%d\t%d\t%d\tunknown\t", m2, n2, k2);
		if(m3<MAX_M_DIM && n3<MAX_N_DIM && k3<MAX_N_DIM) {
			if(m3 > 0 && n3 > 0)
				printf("%d\t%d\t%d\t%e\t", m3, n3, k3, get_pvalue(m3, n3, k3));
			else
				printf("%d\t%d\t%d\t%e\t", m3, n3, k3, 1.0);
		}
		else
			printf("%d\t%d\t%d\tunknown\t", m3, n3, k3);
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%s\n", identity, min_len, max_len, max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, d_conf, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1], orthologs_block[0], orthologs_block[1]);
	}
	
	fclose(fp);
	
	free_data_info_table();

	return 0;
}
