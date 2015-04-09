#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
	FILE *fp;
	char buf[5000], chr1[100], chr2[100], orient;
	int index, b1, e1, b2, e2, length, status, status1; // 0 : unknown; 1 : no GC; 2 : GC
	char identity[100], pvalue[100], GC_len[100], outgroup_name[2][100], outgroup_orient[2], orthologs_block[2][100];
	int max_beg1, min_beg1, min_end1, max_end1, max_beg2, min_beg2, min_end2, max_end2, direction, outgroup_start[2], outgroup_end[2], triplet_status;
	
	if ( argc != 2) {
		printf("show_GC_cluster GC-file\n");
		return 1;
	}
	
	fp = fopen(argv[1], "r");
	while(fgets(buf, 5000, fp)) {
		sscanf(buf, "%d %s %d %d %s %d %d %c %d %s %d %d %s %s %d %d %d %d %d %d %d %d %d %s %d %d %c %s %d %d %c %s %s %d", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, identity, &status1, &status, GC_len, pvalue, &max_beg1, &min_beg1, &min_end1, &max_end1, &max_beg2, &min_beg2, &min_end2, &max_end2, &direction, outgroup_name[0], &outgroup_start[0], &outgroup_end[0], &outgroup_orient[0], outgroup_name[1], &outgroup_start[1], &outgroup_end[1], &outgroup_orient[1], orthologs_block[0], orthologs_block[1], &triplet_status);
		if(status == 2) {
			printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%s\t%d\t%1.2e\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%s\t%d\n", index, chr1, b1, e1, chr2, b2, e2, orient, length, identity, (int)atof(GC_len), atof(pvalue), min_beg1, min_end1, min_beg2, min_end2, direction, outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1], orthologs_block[0], orthologs_block[1], triplet_status);
		}
	}

	return 0;
}
