/* FDR_test - FDR test is applied. Cut off threshold is determined */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PVALUE 0.01

struct Pvalue {
	double pvalue;
};

struct Pvalue *pvalue[100000000];
int pvalue_num;

/* Sort all p-values in ascending order */
static int compar(const void *a, const void *b) {
  if((*(struct Pvalue **)a)->pvalue < (*(struct Pvalue **)b)->pvalue)
	  return -1;
  else 
	  return 1;
}


int main(int argc, char *argv[]) {
	FILE *fp;
	char buf[5000];
	char value[100];
	int i;

	if ( argc != 2) {
		printf("FDR_test FDR-file\n");
		return 1;
	}

	fp = fopen(argv[1], "r");
	pvalue_num = 0;
	while(fgets(buf, 5000, fp)) {
		sscanf(buf, "%s", value);
		pvalue[pvalue_num] = (struct Pvalue *)malloc(sizeof(struct Pvalue));
		pvalue[pvalue_num]->pvalue = atof(value);
		pvalue_num++;
	}

	qsort((void **)pvalue, pvalue_num, sizeof(struct Pvalue *), compar);
	for(i=0; i<pvalue_num; i++) {
		if(pvalue[i]->pvalue > ((float)((i + 1) * PVALUE) / (float)pvalue_num))
			break;
	}
	printf("%f\n", ((float)((i + 1) * PVALUE) / (float)pvalue_num));
		
	fclose(fp);

	return 0;
}

