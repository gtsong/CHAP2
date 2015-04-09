#include "main.h"
#include "util.h"

int main(int argc, char **argv) {
	FILE *fp;
	char buf[BIG];
	char name[100], scaf_name[100];
//	int num_scafs = 0;
	int len = 0;

	if (argc != 3)
		fatal("args: seq-file(fasta) seq_name");
	
	strcpy(name, argv[2]);

	if( strchr(name, '.') != NULL ) {
		fatalf("%s: sequence (species) name does not allow to include \".\"", name);
	}

	fp = ckopen(argv[1], "r");

//	while(fgets(buf, BIG, fp))
//	{
//		if( (buf[0] == '>') || (buf[0] == '<') ) {
//			num_scafs++;	
//		}
//	}
	
	fseek(fp, 0, SEEK_SET);

//	if( num_scafs <= 1 ) {
//		while(fgets(buf, BIG, fp))
//		{
//			if( (buf[0] == '>') || (buf[0] == '<') ) {
//				sscanf(buf+1, "%s %*s", scaf_name);
//				if( strstr(scaf_name, name) != NULL ) {
//					printf(">%s\n", scaf_name);
//				}
//				else {
//					printf(">%s.%s\n", name, scaf_name);
//				}
//			}
//			else {
//				len = strlen(buf);
//				if( buf[len-1] != '\n' ) printf("%s\n", buf);
//				else printf("%s", buf);
//			}
//		}
//	}
//	else {
		while(fgets(buf, BIG, fp))
		{
			if( (buf[0] == '>') || (buf[0] == '<') ) {
				sscanf(buf+1, "%s %*s", scaf_name);
				if( strstr(scaf_name, name) != NULL ) {
					printf(">%s\n", scaf_name);
				}
				else {
					printf(">%s.%s\n", name, scaf_name);
				}
			}
			else {
				len = strlen(buf);
				if( buf[len-1] != '\n' ) printf("%s\n", buf);
				else printf("%s", buf);
			}
		}	
//	}

	fclose(fp);
	return EXIT_SUCCESS;
}
