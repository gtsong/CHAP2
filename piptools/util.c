/* util.c -- general utility routines to augment the standard C library
 *
 * $Revision: 1.22 $$Date: 2001/10/24 22:38:25 $
 */

#include "util.h"

char *argv0;

/* set_argv0 -------------------------------------------- set name of program */
void set_argv0 (char *name) {
	char *p = strrchr (name, '/');
	argv0 = (p ? p+1 : name);
}

/* print_argv0 ---------------------------------------- print name of program */
void print_argv0 (void) {
	if (argv0) {
		char *p = strrchr (argv0, '/');
		(void) fprintf (stderr, "%s: ", p ? p+1 : argv0);
	}
}

/* optval -------------------------------- get next argument; increment index */
char *optval (int argc, char **argv, int *i, char *syntax) {
	if (++(*i) < argc)
		return argv[*i];
	else
		fatal (usage (syntax));
}

/* usage --------------------------------------------- format a usage message */
char *usage (char *syntax) {
	int bufsize = strlen (argv0) + strlen (syntax) + 20;
	char *buf = ckalloc (bufsize);
	sprintf (buf, "Usage:\n\n%s %s", argv0, syntax);
	return buf;
}

/* fatal ---------------------------------------------- print message and die */
void fatal (const char *msg) {
	fatalf ("%s", msg);
}

/* fatalf --------------------------------- format message, print it, and die */
void fatalf (const char *fmt, ...) {
	va_list ap;
	va_start (ap, fmt);
	fflush (stdout);
	print_argv0 ();
	(void) vfprintf (stderr, fmt, ap);
	(void) fputc ('\n', stderr);
	va_end (ap);
	exit (1);
}

/* fatalfr ----------------- format message + error string, print it, and die */
void fatalfr (const char *fmt, ...) {
	va_list ap;
	va_start (ap, fmt);
	fflush (stdout);
	print_argv0 ();
	(void) vfprintf (stderr, fmt, ap);
	(void) fprintf (stderr, ": %s\n", strerror (errno));
	va_end (ap);
	exit (1);
}

/* ckopen -------------------------------------- open file; check for success */
FILE *ckopen (const char *name, const char *mode) {
	FILE *fp;
	if (!(fp = fopen (name, mode)))
		fatalf ("Cannot open %s.", name);
	return fp;
}

/* ckalloc -------------------------------- allocate space; check for success */
void *ckalloc (size_t amount) {
	void *p;
	if ((long) amount < 0)                                 /* was "<= 0" -CR */
		fatal ("Request for negative space in ckalloc().");
	if (amount == 0)
		amount = 1;   /* ANSI portability hack */
	if (!(p = malloc (amount)))
		fatalf ("Ran out of memory trying to allocate %lu.",
		        (unsigned long) amount);
	return p;
}

/* ckallocz -------------------- allocate space; zero fill; check for success */
void *ckallocz (size_t amount) {
	void *p = ckalloc (amount);
	memset (p, 0, amount);
	return p;
}

/* ckrealloc --------------------- re-size allocated space; check for success */
void *ckrealloc (void *p, size_t amount) {
	if ((long) amount < 0)
		fatal ("Request for negative space in ckrealloc().");
	p = (p ? realloc (p, amount) : malloc (amount));
	if (!p)
		fatalf ("Ran out of memory trying to re-allocate %lu.",
		        (unsigned long) amount);
	return p;
}

/* same_string ---------------------- determine whether two strings are equal */
bool same_string (const char *s, const char *t) {
	return (!strcmp (s, t));
}

/* strsame ------------ determine whether two strings are equal, even if null */
bool strsame (const char *s, const char *t) {
	return ((!s || !t) ? (s == t) : !strcmp (s, t));
}

/* starts ------------------------------ determine whether t is a prefix of s */
bool starts (const char *s, const char *t) {
	return (!strncmp (s, t, strlen (t)));
}

/* skip_ws ------------------- find the first non-whitespace char in a string */
char *skip_ws (const char *s) {
	while (isspace (*s))
		s++;
	return (char *) s;
}

/* copy_string ---------------------- save string s somewhere; return address */
char *copy_string (const char *s) {
	char *p = ckalloc (strlen (s) + 1);    /* +1 to hold '\0' */
	return strcpy (p, s);
}

/* copy_substring ------------ save first n chars of string s; return address */
char *copy_substring (const char *s, int n) {
	char *p = ckalloc ((size_t) n + 1);    /* +1 to hold '\0' */
	memcpy (p, s, (size_t) n);
	p[n] = 0;
	return p;
}

/* swap --------------------------------------------------- swap two integers */
void swap (int *a, int *b) {
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

/* overlaps ------------------------- determine whether two intervals overlap */
bool overlaps (int a0, int a1, int b0, int b1) {
	if (a1 < a0)
		swap (&a0, &a1);
	if (b1 < b0)
		swap (&b0, &b1);
	return (b0 <= a1 && b1 >= a0);
}

