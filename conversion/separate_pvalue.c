/* separate_pvalue - separate results for different species */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
	FILE *fp, *fp1;
	char buf[10000], *str, *sp_name1, *sp_name2, filename[1000];
	int file_opend = 0;
	
	if ( argc != 2) {
		printf("separate_pvalue P-value-file\n");
		return 1;
	}
	
	fp = fopen(argv[1], "r");
	fp1 = NULL;

	while(fgets(buf, 10000, fp) != NULL) {
		if(buf[0] == '#') {
			str = buf + 1;
			sp_name1 = strtok(str, ".");
			sp_name2 = strtok(NULL, "\n");
			sprintf(filename, "%s.d/%s.d/all.pvalue", sp_name1, sp_name2);
			printf("%s\n", filename);
			if(file_opend)
				fclose(fp1);
			fp1 = fopen(filename, "w");
			file_opend = 1;
		}
		else
			fprintf(fp1, "%s", buf);
	}

	fclose(fp);
	if(file_opend)
		fclose(fp1);

	return 0;
}

