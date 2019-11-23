#include <stdio.h>
#include <stdlib.h>
#include "alloc2d.h"
#ifdef S_ALLOC
#define alloc(nelem,elsize) S_alloc((unsigned)(nelem),(unsigned)(elsize))
#define free(x)        /* NOOP */;
#else
#define alloc(nelem,elsize) malloc((nelem)*(elsize))
#endif

/*
 * Matrix allocation routines.  To allocate a 10 by 10 matrix of floats
 * use:
 * 
 * void **alloc2d(); float **x;
 * 
 * x = (float **) alloc2d(10, 10, sizeof(float));
 * 
 * To free this matrix use:
 * 
 * free2d(x);
 */
/*--Allocate & free 2 arrays--*/
void       **alloc2d(int dim1, int dim2, int size)
{
  int          i;
  unsigned     nelem;
  char        *p, **pp;

  nelem = dim1 * dim2;
  p = (void *) alloc((unsigned) nelem, size);
  if (p == NULL)
    return (NULL);
  pp = (char **) alloc((unsigned) dim1, (unsigned) sizeof(char *));
  if (pp == NULL) {
    free(p);
    return (NULL);
  }
  for (i = 0; i < dim1; i++)
    pp[i] = p + i * dim2 * size;

  return ((void **) pp);
}
int          free2d(void **mat)
{
  if (mat != NULL && *mat != NULL)
    free((void *) *mat);
  if (mat != NULL)
    free((void *) mat);
  return (0);
}
/*--Allocate 3 dimensional array----------*/
void      ***alloc3d(int dim1, int dim2, int dim3, int size)
{
  int          i;
  void       **pp, ***ppp;

  pp = (void **) alloc2d(dim1 * dim2, dim3, size);
  if (pp == NULL)
    return (NULL);
  ppp = (void ***) alloc((unsigned) dim1, (unsigned) sizeof(void **));
  if (ppp == NULL) {
    free2d(pp);
    return (NULL);
  }
  for (i = 0; i < dim1; i++)
    ppp[i] = pp + i * dim2;
  return (ppp);
}
int          free3d(void ***mat)
{
  free2d(*mat);
  if (mat != NULL)
    free((void *) mat);
  return (0);
}
/*---Allocate 4 dimensional array--------*/
void     ****alloc4d(int dim1, int dim2, int dim3, int dim4, int size)
{
  int          i;
  void      ***ppp, ****pppp;

  ppp = (void ***) alloc3d(dim1 * dim2, dim3, dim4, size);
  if (ppp == NULL)
    return (NULL);
  pppp = (void ****) alloc((unsigned) dim1, sizeof(void ****));
  if (pppp == NULL) {
    free3d(ppp);
    return (NULL);
  }
  for (i = 0; i < dim1; i++)
    pppp[i] = ppp + i * dim2;
  return (pppp);
}
int          free4d(void ****mat)
{
  free3d(*mat);
  if (mat != NULL)
    free((void *) mat);
  return (0);
}
/*--allocate lower triangule of square matrix------------*/
void       **alloclt(int dim, int size)
{
  int          i;
  char        *p, **pp;
  unsigned     nelem;

  nelem = (dim * (dim + 1)) / 2;
  p = (char *) alloc((unsigned) nelem, size);
  if (p == NULL)
    return (NULL);
  pp = (char **) alloc((unsigned) dim, (unsigned) sizeof(char *));
  if (pp == NULL) {
    free(p);
    return (NULL);
  }
  pp[0] = p;
  for (i = 1; i < dim; i++)
    pp[i] = pp[i - 1] + i * size;
  return ((void **) pp);
}

/*--allocate upper triangule of square matrix------------*/
void       **allocut(int dim, int size)
{
  int          i;
  char        *p, **pp;
  unsigned     nelem;

  nelem = (dim * (dim + 1)) / 2;
  p = (char *) alloc(nelem, (unsigned) size);
  if (p == NULL)
    return (NULL);
  pp = (char **) alloc((unsigned) dim, (unsigned) sizeof(char *));
  if (pp == NULL) {
    free(p);
    return (NULL);
  }
  pp[0] = p;
  for (i = 1; i < dim; i++)
    pp[i] = pp[i - 1] + (dim - i) * size;
  return ((void **) pp);
}
