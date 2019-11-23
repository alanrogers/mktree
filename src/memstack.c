/****************************************************************
memstack.c: This file implements a very simple system for allocating
and deallocating storage.  stackalloc(10) returns a pointer to an
array 10 bytes long.  It works just like malloc().  stackfree()
frees ALL memory allocated by previous calls to stackalloc().  These
routines are much faster than malloc and free.

                                  Alan Rogers 8/17/92
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memstack.h"

MEMSTACK    *sstack[MAX_STACKS];	/* stack of pointers to memory stacks */
int          nstacks = 0;	/* # of allocated stacks */
int          sndx = -1;		/* index of current stack */

/*** Allocate a stack for subsequent calls to stackalloc() **/
MEMSTACK    *setmemstack(int size)
{
  if (nstacks == MAX_STACKS) {
    fprintf(stderr, "\nError: Can only allocate %d memory stacks\n",
	    MAX_STACKS);
    exit(1);
  }
  sstack[nstacks] = (MEMSTACK *) malloc(sizeof(MEMSTACK));
  if (sstack[nstacks] == NULL)
    return (NULL);
  sstack[nstacks]->stack = (char *) malloc(size);
  if (sstack[nstacks]->stack == NULL) {
    free(sstack[nstacks]);
    return (NULL);
  }
  sndx = nstacks++;
  sstack[sndx]->size = size;
  sstack[sndx]->ndx = 0;
  fprintf(stderr, "\n%d bytes allocated to stack %d\n", size, nstacks);
  return (sstack[sndx]);
}

/** allocate memory **/
char        *stackalloc(int size)
{
  int          i;
  MEMSTACK    *p;

  if (sstack[sndx]->ndx + size > sstack[sndx]->size) {	/* need a new stack */
    if (sndx < nstacks - 1) {	/* next stack already allocated */
      ++sndx;
      return (stackalloc(size));
    } else {
      if (size > sstack[sndx]->size)
	p = setmemstack(size);
      else
	p = setmemstack(sstack[sndx]->size);
      if (p == NULL)
	return (NULL);		/* can't allocate that much memory */
    }
  }
  i = sstack[sndx]->ndx;	/* index of next free memory */
  sstack[sndx]->ndx += size;	/* reset stack index */
  return (sstack[sndx]->stack + i);
}
/** free all memory in each stack **/
void         stackfree(void)
{
  int          i;

  for (i = 0; i < nstacks; i++)
    sstack[i]->ndx = 0;
  sndx = 0;
}
/*** deallocate everything **/
void         deallocstack(void)
{
  int          i;

  for (i = 0; i < nstacks; i++) {
    free(sstack[i]->stack);
    free(sstack[i]);
  }
  nstacks = 0;
  sndx = -1;
}
/** return number of available bytes of memory **/
int          checkstack(void)
{
  return (sstack[sndx]->size - sstack[sndx]->ndx);
}
/**** allocate memory, exit on failure ****/
char        *stackmustalloc(int size)
{
  char        *x;

  x = (char *) stackalloc(size);
  if (x == NULL) {
    fprintf(stderr, "\nstackmustalloc(): Can't allocate %d bytes\n", size);
    exit(1);
  }
  return (x);
}
/*** duplicate a string, allocating memory ****/
char        *stackstrdup(char *s)
{
  char        *s2, *s3;

  s3 = s2 = stackmustalloc((strlen(s) + 1) * sizeof(char));
  while (*s != '\0')
    *s2++ = *s++;
  *s2 = '\0';
  return (s3);
}

/*** report state of stack ****/
void         stackstatus(FILE * fp)
{
  int          i;

  if (nstacks == 0)
    fprintf(fp, "\nNo memstacks allocated: nstacks=%d sndx=%d",
	    nstacks, sndx);
  for (i = 0; i < nstacks; i++)
    fprintf(fp, "\nMemStack %d: %d/%d used; %d avail",
	    i, sstack[i]->ndx, sstack[i]->size,
	    sstack[i]->size - sstack[i]->ndx);
  fflush(fp);
}
