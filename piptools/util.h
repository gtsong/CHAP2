/* $Revision: 1.21 $$Date: 2001/10/24 22:38:25 $ */

#ifndef PIPTOOLS_UTIL
#define PIPTOOLS_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>     /* HUGE_VAL */
#include <limits.h>   /* INT_MAX, INT_MIN, LONG_MAX, LONG_MIN, etc. */
#include <stdarg.h>
#include <assert.h>
#include <errno.h>

#undef isspace        /* use functions instead of macros, to handle char */
#undef isdigit

typedef int bool;

extern char *argv0;

void set_argv0 (char *name);
void print_argv0 (void);
char *optval (int argc, char **argv, int *i, char *syntax);
char *usage  (char *syntax);
#ifdef __GNUC__     /* avoid some "foo might be used uninitialized" warnings */
	void fatal (const char *msg) __attribute__ ((noreturn));
	void fatalf (const char *fmt, ...) __attribute__ ((noreturn));
	void fatalfr (const char *fmt, ...) __attribute__ ((noreturn));
#else
	void fatal (const char *msg);
	void fatalf (const char *fmt, ...);
	void fatalfr (const char *fmt, ...);
#endif
FILE *ckopen (const char *name, const char *mode);
void *ckalloc (size_t amount);
void *ckallocz (size_t amount);
void *ckrealloc (void *p, size_t amount);
bool same_string (const char *s, const char *t);
bool strsame (const char *s, const char *t);
bool starts (const char *s, const char *t);
char *skip_ws (const char *s);
char *copy_string (const char *s);
char *copy_substring (const char *s, int n);
void swap (int *a, int *b);
bool overlaps (int a0, int a1, int b0, int b1);

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#undef XOR
#define XOR(x,y) (((x) && !(y)) || ((y) && !(x)))

#endif
