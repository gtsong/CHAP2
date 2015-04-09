/* exons.c -- read the file of gene/exon positions
 *
 *  Entry points:
 *
 *	int get_exons (char *filename, int **Beg, int **End, int **Dir,
 *	               char ***Gene_name, char ***Exon_name, int max_errs);
 *		For each exon specified in the file, get_exons creates a record
 *		(implemented as parallel arrays) with a begin point, an end
 *		point, and a transcription direction (1 = left-to-right,
 *		-1 = right-to-left, 0 = not applicable).  For each gene, there
 *		are two such records, one for the gene (distinguished from all
 *		other records by having gene_name != NULL) followed immediately
 *		by one for the CDS limits.  If the file contains no "+" line
 *		for that gene, then the limits are both 0.  Up to max_errs
 *		error messages will be printed before execution is terminated;
 *		however a value of 0 is treated as 1 for backward compatibility.
 *		A negative value for max_errs is equivalent to infinity, i.e.
 *		all messages will be printed, without premature termination.
 *
 *	void free_exons (int *beg, int *end, int *dir, char **gene_name,
 *	                 char **exon_name, int n);
 *
 *	char *exon_header (char *filename);
 *
 * $Revision: 1.25 $$Date: 2003/11/04 21:17:04 $
 */

#include "exons.h"
#include "util.h"

enum { BUFSIZE = 1000 };

static int nerrs = 0;

static char *next_nonempty_line (FILE *fp, char *buf, int bufsize, int *lnbr);
static int test_isheader (char *str);    /* "is*" names are reserved by ansii rules :-( */
static int test_isnumber (char *str);
static char *get2numbers (char *p, int *i, int *j, char *buf, int line_nbr,
                          int max_errs);
static void set_exon (int n, int *beg, int *end, int *dir, char **gene_name,
                      char **exon_name, int b, int e, int d, char *g_name,
                      char *e_name);
static void error (int line_nbr, char *buf, const char *msg, int max_errs);


int get_exons (char *filename, int **Beg, int **End, int **Dir,
               char ***Gene_name, char ***Exon_name, int max_errs) {
	FILE *fp;
	int n, Nexons, nexon, ngene, i, j, line_nbr, first_exon, direc;
	int *beg, *end, *dir;
	char buf [BUFSIZE], **gene_name, **exon_name, *p;

	fp = ckopen (filename, "r");

	line_nbr = 0;
	if (next_nonempty_line (fp, buf, BUFSIZE, &line_nbr) && !test_isheader (buf))
		rewind (fp);
	nexon = ngene = 0;
	while (next_nonempty_line (fp, buf, BUFSIZE, &line_nbr))
		if (strchr ("<|>", buf[0]))
			ngene++;
		else if (buf[0] != '+')
			nexon++;
	Nexons = 2*ngene + nexon;

	*Beg = beg = (int *) ckalloc (Nexons * sizeof (int));
	*End = end = (int *) ckalloc (Nexons * sizeof (int));
	*Dir = dir = (int *) ckalloc (Nexons * sizeof (int));
	*Gene_name = gene_name = (char **) ckalloc (Nexons * sizeof (char *));
	*Exon_name = exon_name = (char **) ckalloc (Nexons * sizeof (char *));
	for (n = 0; n < Nexons; n++)
		exon_name[n] = NULL;

	rewind (fp);
	line_nbr = 0;
	if (next_nonempty_line (fp, buf, BUFSIZE, &line_nbr) && !test_isheader (buf)) {
		rewind (fp);
		line_nbr = 0;	/* this time it counts */
	}
	first_exon = 0;
	direc = 0;    /* pacify lint */
	for (n = 0; next_nonempty_line (fp, buf, BUFSIZE, &line_nbr); n++)
		if (strchr ("<|>", buf[0])) {
			if (n >= Nexons-1)
				fatal ("Internal error #1.");
			direc = (buf[0] == '>' ? 1 : (buf[0] == '<' ? -1 : 0));
			p = get2numbers (buf+1, &i, &j, buf, line_nbr, max_errs);
			set_exon (n, beg, end, dir, gene_name, exon_name,
			          i, j, direc, p, NULL);
			n++;  /* CDS information */
			set_exon (n, beg, end, dir, gene_name, exon_name,
			          0, 0, direc, NULL, NULL);
			first_exon = n+1;
		} else if (buf[0] == '+') {
			if (*get2numbers (buf+1, &i, &j, buf, line_nbr, max_errs))
				error (line_nbr, buf, "Expecting just two numbers",
				       max_errs);
			if (!first_exon)
				error (line_nbr, buf, "No gene named for CDS positions.",
				       max_errs);
			if (i < beg[first_exon-2] || j > end[first_exon-2])
				error (line_nbr, buf, "CDS positions outside of gene.",
				       max_errs);
			if (beg[first_exon-1] != 0)
				error (line_nbr, buf, "Redundant CDS positions.",
				       max_errs);
			set_exon (first_exon-1, beg, end, dir, gene_name,
			          exon_name, i, j, direc, NULL, NULL);
			n--;  /* information recorded in an earlier record */
		} else {
			if (n >= Nexons)
				fatal ("Internal error #2.");
			p = get2numbers (buf, &i, &j, buf, line_nbr, max_errs);
			if (n > first_exon && i <= end[n-1])
				error (line_nbr, buf, "Exons out of order.", max_errs);
			set_exon (n, beg, end, dir, gene_name, exon_name,
			          i, j, direc, NULL, (!*p ? NULL : p));
		}
	fclose (fp);
	if (nerrs)
		exit (nerrs);
	return Nexons;
}

void free_exons (int *beg, int *end, int *dir, char **gene_name,
                 char **exon_name, int n) {
	free (beg);
	free (end);
	free (dir);
	while (--n >= 0) {
		if (gene_name[n])
			free (gene_name[n]);
		if (exon_name[n])
			free (exon_name[n]);
	}
	free (gene_name);
	free (exon_name);
}

char *exon_header (char *filename) {
	FILE *fp;
	char buf [BUFSIZE], *p;
	int line_nbr = 0;

	fp = ckopen (filename, "r");

	p = next_nonempty_line (fp, buf, BUFSIZE, &line_nbr);
	fclose (fp);
	return ((p && test_isheader (buf)) ? copy_string (buf) : NULL);
}

static char *next_nonempty_line (FILE *fp, char *buf, int bufsize, int *lnbr) {
	char *status, *p;
	int c;

	(*lnbr)++;
	while ((status = fgets (buf, bufsize, fp))) {
		p = skip_ws (buf);
		/* if comment line is longer than buf, skip the rest of it */
		if (*p == '#' && !strchr (p, '\n')) {
			while ((c = getc (fp)) != '\n' && c != EOF)
				;
		}
		if (*p && *p != '#')
			break;
		(*lnbr)++;
	}

	if (status && (p = strchr (buf, '\n')))
		*p = '\0';

	return status;
}

static int test_isheader (char *str) {
	return (!test_isnumber (str) && !strchr (">|<+", str[0]));
}

static int test_isnumber (char *str) {
	long num;
	char *endptr;

	errno = 0;
	num = strtol (str, &endptr, 10);
	if (errno == EINVAL)
		return 0;
	return (!*endptr || isspace (*endptr));
}

static char *get2numbers (char *p, int *i, int *j, char *buf, int line_nbr,
                          int max_errs) {
	while (*p == ' ' || *p == '\t')
		p++;
	*i = atoi (p);
	while (isdigit (*p))
		p++;
	while (*p == ' ' || *p == '\t' || *p == '-' || *p == '.')
		p++;
	*j = atoi (p);
	if (*i <= 0 || *i > *j)
		error (line_nbr, buf, "Improper endpoints.", max_errs);
	while (isdigit (*p))
		p++;
	while (*p == ' ' || *p == '\t')
		p++;

	return p;
}

static void set_exon (int n, int *beg, int *end, int *dir, char **gene_name,
                      char **exon_name, int b, int e, int d, char *g_name,
                      char *e_name) {
	beg[n] = b;
	end[n] = e;
	dir[n] = d;
	gene_name[n] = (g_name ? copy_string (g_name) : NULL);
	exon_name[n] = (e_name ? copy_string (e_name) : NULL);
}

static void error (int line_nbr, char *buf, const char *msg, int max_errs) {
	nerrs++;
	fprintf (stderr, "\nError in exons file at line %d.\n%s\n",
	         line_nbr, buf);
	if (nerrs < max_errs || max_errs < 0)
		fprintf (stderr, "%s\n", msg);
	else
		fatal (msg);
}

