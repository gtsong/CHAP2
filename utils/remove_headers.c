#include "main.h"
#include "util.h"

int main(int argc, char **argv) {
	FILE *fp;
	char buf[BIG];
//	int num_scafs = 0;
	int len = 0;

	if (argc != 2)
		fatal("args: seq-file(fasta)");
	
	fp = ckopen(argv[1], "r");

	while(fgets(buf, BIG, fp))
	{
		if( (buf[0] == '>') || (buf[0] == '<') ) {}
		else {
			len = strlen(buf);
			if( buf[len-1] != '\n' ) printf("%s\n", buf);
			else printf("%s", buf);
		}
	}
	
	fclose(fp);
	return EXIT_SUCCESS;
}
