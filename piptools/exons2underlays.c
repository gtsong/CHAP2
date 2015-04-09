/* exons2underlays.c -- convert an exons file to underlay format
 *
 * This program was suggested by Mike Wilson.
 *
 * Based on $Revision: 1.24 $$Date: 2001/06/26 23:51:03 $ from piptools.2003-11-04
 */

#include "exons.h"
#include "util.h"

#ifndef PIPTOOLS_EXONS
	fatal ("Missing or incorrect header file \"exons.h\"." \
	       "  (Are there several with the same name?)");
#endif
#ifndef PIPTOOLS_UTIL
	fatal ("Missing or incorrect header file \"util.h\"." \
	       "  (Are there several with the same name?)");
#endif

#define SYNTAX  "exons_file [-laj] [-split] " \
                  "-fexon fwd_exon_color [-futr fwd_utr_color] " \
                  "[-rexon rev_exon_color] [-rutr rev_utr_color] " \
                  "[-pexon pseudogene_color] [-intron intron_color] > output_file"

#define PSEUDO_SUFFIX         "_ps"
#define DEFAULT_PSEUDO_COLOR  "PaleGray"
#define DEFAULT_INTRON_COLOR  "White"

static char *Types[6] = {
	"fwd_exon", "fwd_UTR", "rev_exon", "rev_UTR", "pseudo", "gene"
};

static void get_args (int argc, char **argv, char **filename, bool *laj,
                      bool *split, char **colors);
static void print_underlay (int begin, int end, char *gname, char *item_name,
                            int exon_num, int dir, bool laj, bool split);
static bool pseudo (char *gname);
static bool ends_with (char *s, char *t);
static char *strtrim (char *s);
static char *one_word (char *s);

int main (int argc, char **argv) {
	char *filename = NULL;                  /* name of exons file           */
	bool laj = 0;                           /* make labeled underlays?      */
	bool split = 0;                         /* split strands?               */
	char *colors[6] = {NULL, NULL, NULL,    /* color names, in SYNTAX order */
	                   NULL, NULL, NULL};
	char *header;                           /* header line from exons file  */
	char **gene_name, **exon_name;          /* arrays of gene/exon info     */
	int *begin, *end, *dir, *exon_num;      /* arrays of gene/exon info     */
	int length;                             /* length of arrays             */
	int i, j;                               /* array indices                */
	char *gname;                            /* name of current gene         */
	int cds_begin, cds_end;                 /* translated region            */
	int num;                                /* exon number within gene      */

	/* get arguments */
	set_argv0 (argv[0]);
	get_args (argc, argv, &filename, &laj, &split, colors);
	if (!filename)
		fatal (usage (SYNTAX));
	if (!colors[0] && !colors[2])     /* exons */
		fatal (usage (SYNTAX));
	else if (!colors[0])
		colors[0] = colors[2];
	else if (!colors[2])
		colors[2] = colors[0];
	if (!colors[1] && !colors[3]) {   /* utrs */
		colors[1] = colors[0];
		colors[3] = colors[2];
	}
	else if (!colors[1])
		colors[1] = colors[3];
	else if (!colors[3])
		colors[3] = colors[1];
	if (!colors[4])                   /* pseudo-exons */
		colors[4] = DEFAULT_PSEUDO_COLOR;
	if (!colors[5])                   /* genes, pseudo-genes */
		colors[5] = DEFAULT_INTRON_COLOR;

	/* print header */
	header = exon_header (filename);
	printf ("# Underlays converted from \"%s\" by %s.\n", filename, argv0);
	printf ("# %s\n", (header ? header : "(no header specified)"));

	/* print color legend */
	putchar ('\n');
	for (i = 0; i < 6; i++) {
		assert (colors[i]);
		printf ("%s %s\n", colors[i], Types[i]);
	}

	/* read gene and exon info */
	length = get_exons (filename, &begin, &end, &dir, &gene_name,
	                    &exon_name, -1);    /* print all errors */

	/* compute exon numbers */
	exon_num = (int *) ckalloc (length * sizeof (int));
	assert (!length || gene_name[0]);
	for (i = 0; i < length; ) {
		num = 0;
		if (dir[i] < 0) {
			for (j = i+2; j < length && !gene_name[j]; j++)
				num++;
		}
		for (j = i+2; j < length && !gene_name[j]; j++)
			exon_num[j] = (dir[j] < 0 ? num-- : ++num);
		i = j;
	}

	/* print underlays */
	gname = NULL; cds_begin = cds_end = 0;    /* pacify lint */
	for (i = 0; i < length; i++) {
		if (gene_name[i]) {     /* gene */
			gname = one_word (strtrim (gene_name[i]));
			printf ("\n");
			print_underlay (begin[i], end[i], gname, gname,
			                0, dir[i], laj, split);
			assert (++i < length);     /* cds */
			cds_begin = begin[i];
			cds_end = end[i];
		}
		else {     /* exon */
			print_underlay (begin[i], end[i], gname, strtrim (exon_name[i]),
			                exon_num[i], dir[i], laj, split);
			if (begin[i] < cds_begin)
				print_underlay (begin[i], MIN (cds_begin-1, end[i]),
				                gname, (dir[i] < 0 ? "3'UTR" : "5'UTR"),
				                0, dir[i], laj, split);
			if (cds_end > 0 && cds_end < end[i])
				print_underlay (MAX (cds_end+1, begin[i]), end[i],
				                gname, (dir[i] < 0 ? "5'UTR" : "3'UTR"),
				                0, dir[i], laj, split);
		}
	}

	/* clean up */
	free_exons (begin, end, dir, gene_name, exon_name, length);
	free (exon_num);
	return 0;
}

static void get_args (int argc, char **argv, char **filename, bool *laj,
                      bool *split, char **colors) {
	int i;

	for (i = 1; i < argc; i++) {
		if (strsame (argv[i], "-laj"))
			*laj = 1;
		else if (strsame (argv[i], "-split"))
			*split = 1;
		else if (strsame (argv[i], "-fexon"))
			colors[0] = optval (argc, argv, &i, SYNTAX);
		else if (strsame (argv[i], "-futr"))
			colors[1] = optval (argc, argv, &i, SYNTAX);
		else if (strsame (argv[i], "-rexon"))
			colors[2] = optval (argc, argv, &i, SYNTAX);
		else if (strsame (argv[i], "-rutr"))
			colors[3] = optval (argc, argv, &i, SYNTAX);
		else if (strsame (argv[i], "-pexon"))
			colors[4] = optval (argc, argv, &i, SYNTAX);
		else if (strsame (argv[i], "-intron"))
			colors[5] = optval (argc, argv, &i, SYNTAX);
		else if (!*filename)
			*filename = argv[i];
		else
			fatal (usage (SYNTAX));
	}
}

static void print_underlay (int begin, int end, char *gname, char *item_name,
                            int exon_num, int dir, bool laj, bool split) {
	int type;
	char *half;

	if (exon_num > 0)    /* exon */
		type = (pseudo (gname) ? 4 : (dir < 0 ? 2 : 0));
	else if (ends_with (item_name, "'UTR"))
		type = (pseudo (gname) ? 4 : (dir < 0 ? 3 : 1));
	else
		type = 5;    /* gene */

	half = (split && dir ? (dir < 0 ? "-" : "+") : "");

	if (laj) {
		if (item_name)
			printf ("%d %d (%s) %s %s\n",
			        begin, end, item_name, Types[type], half);
		else
			printf ("%d %d (%s.%d) %s %s\n",
			        begin, end, gname, exon_num, Types[type], half);
	}
	else
		printf ("%d %d %s %s\n", begin, end, Types[type], half);
}

static bool pseudo (char *gname) {
	return ends_with (gname, PSEUDO_SUFFIX);
}

static bool ends_with (char *s, char *t) {
	return (strsame (s + strlen (s) - strlen (t), t));
}

static char *strtrim (char *s) {
	char *p;

	if (!s) return s;
	p = s + strlen (s) - 1;
	while (p >= s && isspace (*p))
		p--;
	*++p = '\0';
	for (p = s; isspace (*p); p++)
		;
	return p;
}

static char *one_word (char *s) {
	char *p;

	if (!s) return s;
	p = strchr (s, ' ');
	if (p)
		*p = '\0';
	return s;
}
