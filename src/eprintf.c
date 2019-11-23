
/*
 * The functions handle errors and then terminate execution.  They are copied
 * from pp 109-114 of:
 *  Brian Kernighan and Rob Pike. 1997. The Practice of Programming. 
 *  Addison-Wesley.
 * I modified them so that progname is passed as an argument rather than being
 * set as a static external.  The static external seemed to invite errors.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include "eprintf.h"

/* eprintf: print error message and exit */
void         eprintf(const char *progname, const char *fmt,...)
{
  va_list      args;

  fflush(stdout);
  if (progname != NULL)
    fprintf(stdout, "%s: ", progname);

  va_start(args, fmt);
  vfprintf(stdout, fmt, args);
  va_end(args);

  if (fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
    fprintf(stdout, " %s", strerror(errno));
  fprintf(stdout, "\n");
  exit(2);			/* conventional value for failed execution */
}

/* estrdup: duplicate a string, report if error */
char        *estrdup(const char *progname, const char *s)
{
  char        *t;

  t = (char *) malloc(strlen(s) + 1);
  if (t == NULL)
    eprintf(progname, "estrdup(\"%.20s\") failed:", s);
  strcpy(t, s);
  return (t);
}

/* emalloc: malloc and report if error */
void        *emalloc(const char *progname, size_t n)
{
  void        *p;

  p = malloc(n);
  if (p == NULL)
    eprintf(progname, "emalloc of %u bytes failed:", n);
  return (p);
}
