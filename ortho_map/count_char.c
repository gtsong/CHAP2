#include "main.h"
#include "util.h"

int debug_mode;
int count_node; 
char S[BIG], T[BIG];

int main(int argc, char **argv)
{
	char buf[LEN_NAME];
	FILE *f;
	int len = 0;
	int max_len = 0;
	
	if( argc == 2 ) {}
	else {
		fatalf("%d arguments are too many or few\n", argc);
	}

	f = ckopen(argv[1], "r");
	while(fgets(buf, LEN_NAME, f)){
		len = strlen(buf);
		if( len > max_len ) max_len = len;
	}
	fclose(f);

	printf("%d",len);
	return EXIT_SUCCESS;
}

