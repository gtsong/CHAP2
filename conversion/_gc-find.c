/* gc-find -- find specified lines in a gc.all file produce by the PSU
*  gene-conversion package, where each line of gc.all reports an oberservation
*  supporting a conversion event (i.e., two paralogs in the primary species and
*  orthologs of one or both of them in an "outgroup" species, meeting certain
*  conditions).
*
* Sample commands:
*   1. gc-find all.gc species human
*   2. gc-find all.gc addr1 12345
*   3. gc-find all.gc addr2 12345
*   4. gc-find all.gc len 400
*   5. gc-find all.gc gc-len 300
*   6. gc-find all.gc pval 3.7e-12
*   7. gc-find all.gc conv1 12345
*   8. gc-find all.gc conv2 12345
*   9. gc-find all.gc event 3
*  10. gc-find all.gc branch 14
*
* Specifying a DNA position (samples 2, 3, 7 and 8) asks for all intervals
* containing that position; addr1 refers to the lower-position paralog,
* addr2 to the higher-position paralog, conv1, conv2 to the converted
* subintervals. With len and gc-len (length of the paralogs or the converted
* subinterval), all observations with at least the specified length are
* reported. In sample 6, all observations with the specified p-value or less
* are reported. The first line of all.gc is a phylogenetic tree for the
* analyzed species, with a unique number assigned to each branch; sample 10
* selected observations where the putative conversion event is on that branch.
*
* Keywords can be truncated, e.g. "spe" in place of "species". The command
* "gc-find" (no arguments) produces a list of the permissible keywords.
*/

#include "lib.h"

char W[100][100];
int nwords, interval, n;

char buf[10000];

char *N[] = {
	"0",
	"species",
	"addr1",
	"3",
	"4",
	"addr2",
	"6",
	"7",
	"len",
	"9",
	"gc-len",
	"pval",
	"conv1",
	"13",
	"conv2",
	"15",
	"direction",
	"outgroup",
	"18",
	"19",
	"20",
	"21",
	"22",
	"23",
	"24",
	"event",
	"branch"
};

int test(char *val) {
	int n1, n2, i;

	if (nwords < n)
		return 0;
	if (interval) {
		i = atoi(val);
		n1 = atoi(W[n]);
		n2 = atoi(W[n+1]);
		return (n1 <= i && i <= n2);
	}
	if (n == 11) // pval
		return (atof(W[n]) <= atof(val));
	if (n == 8 || n == 10) // len or gc-len
		return atoi(W[n]) >= atoi(val);
	return same_string(W[n], val);
}

void get_words() {
	char temp[10000], *x;

	nwords = 0;
	strcpy(temp, buf);
	if ((x = strtok(temp, "\t\n")) == NULL)
		return;
	strcpy(W[0], x);

	for (nwords = 1; (x = strtok(NULL, "\t\n")) != NULL; ++nwords)
		strcpy(W[nwords], x);
}

int main(int argc, char **argv) {
	FILE *fp;
	int i;

	if (argc < 3) {
		fprintf(stderr, "options:");
		for (i = 0; i < 27; ++i)
			if (!isdigit(N[i][0]))
				fprintf(stderr, " %s", N[i]);
			fatal("");
	}
	fp = ckopen(argv[1], "r");
	for (n = 0; n <= 26 && strncmp(argv[2], N[n], strlen(argv[2])); ++n)
		;
	if (n > 26)
		fatalf("field-name %s not found", argv[2]);
	interval = (n == 2 || n == 5 || n == 12 || n == 14);
	//printf("n = %d\n", n);
	while (fgets(buf, 10000, fp)) {
		get_words();
		if (test(argv[3]))
			printf("%s", buf);
	}
	return 0;
}    // added newline to avoid compiler warning  [-CR, 2/2011]
