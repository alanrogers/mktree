
#include <stdio.h>
#include <stdlib.h>
#include "mytypes.h"
#include "bye.h"

/****************************************************************
This file contains short functions of that print error messages
and abort.  They are in a separate file so that they can be easily
shared among different directories.
****************************************************************/

/** open a file, abort if unsuccessful **/
FILE        *mustopen(char *name, char *mode)
{
  FILE        *fp;

  fp = fopen(name, mode);
  if (fp == NULL) {
    fprintf(stderr, "\nCan't open file \"%s\" with mode \"%s\"\n",
	    name, mode);
    exit(1);
  }
  return (fp);
}

/** print an error message and quit **/
void         error(char *s)
{
  fflush(stdout);
  fprintf(stderr, "\nERROR: %s\n", s);
  exit(1);
}

/** allocate memory using malloc; abort on failure **/
char        *mustalloc(unsigned bytes)
{
  char        *p;

  p = (char *) malloc(bytes);
  if (p == NULL) {
    fflush(stdout);
    fprintf(stderr, "\nmustalloc: Can't allocate %d bytes of memory.\n",
	    bytes);
    exit(1);
  }
  return (p);
}

/****************************************************************
Allocate a vector of floats whose index runs from low to high.  On
successful completion, return pointer to new structure of type
FLOATVEC.  On return, abort.
****************************************************************/
FLOATVEC    *newfloatvec(int lo, int hi)
{
  FLOATVEC    *v;
  int          len;

  v = (FLOATVEC *) malloc(sizeof(FLOATVEC));
  if (v == NULL)
    error("newfloatvec: out of memory");
  len = hi - lo + 1;
  v->f = (float *) malloc(len * sizeof(float));
  if (v->f == NULL)
    error("newfloatvec: out of memory");
  v->f -= lo;			/* now f[lo] is 1st entry in vector */
  v->lo = lo;
  v->hi = hi;
  return (v);
}

/****************************************************************
Allocate a vector of doubles whose index runs from low to high.  On
successful completion, return pointer to new structure of type
DBLVEC.  On return, abort.
****************************************************************/
DBLVEC      *newdblvec(int lo, int hi)
{
  DBLVEC      *v;
  int          len;

  v = (DBLVEC *) malloc(sizeof(DBLVEC));
  if (v == NULL)
    error("newdblvec: out of memory");
  len = hi - lo + 1;
  v->f = (double *) malloc(len * sizeof(double));
  if (v->f == NULL)
    error("newdblvec: out of memory");
  v->f -= lo;			/* now f[lo] is 1st entry in vector */
  v->lo = lo;
  v->hi = hi;
  return (v);
}
