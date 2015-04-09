#include "util.h"
#include "util_annot.h"

#define LEN_NAME 100

int main(int argc, char *argv[])
{
	char buf[5000], cur[5000];
	FILE *f;
	char chr_name[LEN_NAME], dir[LEN_NAME], chr[LEN_NAME], name[LEN_NAME], annot[LEN_NAME], gname[LEN_NAME];
	int b = 0, e = 0;
	
	strcpy(chr, "");
	strcpy(name, "");
	strcpy(annot, "");

	if( argc != 3 ) {
    fatal("gff2temp gff_file chr_name\n");
  }

	strcpy(chr_name, argv[2]);
	f = ckopen(argv[1], "r");

  while(fgets(buf, 5000, f))
  {
    if( (buf[0] == '#') || (buf[0] == '>' ) ) {}
    else if( sscanf(buf, "%s %*s %s %d %d %*s %s %*s %s", chr, annot, &b, &e, dir, cur) != 6 ) {
      fatalf("line in wrong gff format: %s\n", buf);
    }
    else {
			if( strcmp(chr, chr_name) == 0 ) {
 	     if( strcmp(annot, "gene") == 0 ) {
					get_gene_name(cur, gname);
					if( strcmp(dir, "+") == 0 ) {
						printf("> %d %d %s\n", b, e, gname);
					}
					else if( strcmp(dir, "-") == 0 ) {
						printf("< %d %d %s (complement)\n", b, e, gname);
					}
					else {
					}
 	     }
				else if( strcmp(annot, "CDS") == 0 ) {
					printf("%d %d\n", b, e);
				}
 	   }
		}
  }

  fclose(f);
	return EXIT_SUCCESS;
}
