/* remove_lower_tri - remove redundant paralogous pairs. */

#include <stdio.h>
#include "maf.h"
#include "mz_scores.h"

int main(int argc, char** argv) {
	struct mafFile* mf;
	struct mafAli* ali;
	struct mafComp* mc;
	
	if ( argc != 2) {
		printf("remove_self maf-file\n");
		return 1;
	}
	
	init_scores70();
	
	mafWriteStart(stdout, 0);
	
	mf = mafOpen(argv[1], 0);
	while((ali = mafNext(mf)) != NULL) {
		mc = ali->components;

		if(mc->next->strand == '+' && mc->start > mc->next->start)
			continue;
		else if(mc->next->strand == '-' && mc->start > (mc->next->srcSize - mc->next->start - mc->next->size))
			continue;
		else
			mafWrite(stdout, ali);
	}

	mafFileFree(&mf);
	
	mafWriteEnd(stdout);

	return 0;
}

