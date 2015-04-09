#include "main.h"
#include "util.h"
#include "seq.h"

#define LEN 10000000

char S[LEN], T[LEN];

int main(int argc, char **argv) {
	FILE *fp;
	SEQ *sf;
	char header[1000];
	int num_lower1 = 0, num_lower2 = 0;
	int num_len1 = 0, num_len2 = 0;
//	int col = 0;
	int i = 0;
	int N;
//	int len = 0;

	if (argc != 3)
		fatal("args: seq-file1 seq-file2");
	
	sf = seq_get(argv[1]);
	N = SEQ_LEN(sf);
	strcpy(S, (char *)SEQ_CHARS(sf));
	seq_close(sf);
	num_len1 = strlen(S);
	for( i = 0; i < num_len1; i++ ) {
		if(islower(S[i]) != 0) num_lower1++;
	}

	fp = fopen(argv[2], "r");
	fgets(header, 1000, fp);
	fclose(fp);

	sf = seq_get(argv[2]);
	N = SEQ_LEN(sf);
	strcpy(T, (char *)SEQ_CHARS(sf));
	seq_close(sf);
	num_len2 = strlen(T);
	for( i = 0; i < num_len2; i++ ) {
		if(islower(T[i]) != 0) num_lower2++;
	}

	if( num_lower1 > num_lower2 ) printf("first\n");
	else printf("second\n");
/*
	len = strlen(header);
	printf("%s", header);
	if( header[len-1] != '\n' ) printf("\n");

	if( num_lower1 > num_lower2 ) {
		for (col = i = 0; S[i] != '\0'; ++i) {
			putchar(S[i]);
			if (++col == 50) {
				putchar('\n');
				col = 0;
			}
		}
		if (col != 0)
			putchar('\n');
	}
	else {
		for (col = i = 0; T[i] != '\0'; ++i) {
			putchar(T[i]);
			if (++col == 50) {
				putchar('\n');
				col = 0;
			}
		}
		if (col != 0)
			putchar('\n');
	}
*/
	return EXIT_SUCCESS;
}
