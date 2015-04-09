#include "main.h"
#include "util.h"

int main(int argc, char **argv) {
	FILE *fp;
	char buf[BIG];
	char name[100];
	int num_scafs = 0;
	int len = 0;

	if (argc != 3)
		fatal("args: seq-file(fasta) seq_name");
	
	strcpy(name, argv[2]);

	fp = ckopen(argv[1], "r");

	while(fgets(buf, 1000, fp))
	{
		if( (buf[0] == '>') || (buf[0] == '<') ) {
			num_scafs++;	
		}
	}
	
	fseek(fp, 0, SEEK_SET);

	if( num_scafs <= 1 ) {
		while(fgets(buf, 1000, fp))
		{
			if( (buf[0] == '>') || (buf[0] == '<') ) {
				if( strstr(buf, name) != NULL ) {
					printf(">%s\n", name);
				}
			}
			else {
				len = strlen(buf);
				if( buf[len-1] != '\n' ) printf("%s\n", buf);
				else printf("%s", buf);
			}
		}
	}
	else {
		if( (buf[0] == '>') || (buf[0] == '<') ) {
			len = strlen(buf);
			if( strstr(buf, name) != NULL ) {
				if( buf[len-1] != '\n' ) printf(">%s\n", buf);
				else printf(">%s", buf);
			}
			else {
				if( buf[len-1] != '\n' ) printf(">%s.%s\n", name, buf);
				else printf(">%s.%s", name, buf);
			}
		}
		else {
			len = strlen(buf);
			if( buf[len-1] != '\n' ) printf("%s\n", buf);
			else printf("%s", buf);
		}	
	}

	return EXIT_SUCCESS;
}
