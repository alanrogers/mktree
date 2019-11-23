#if !defined( FLT_MANT_DIG )
#include <float.h>
#endif
/*Must somewhere define a type called "real" by e.g. "typedef float real" */
/****************************************************************
randint(n) returns a random integer between 0 and n-1.  It is now
implemented as a macro, which should improve speed.  The float-int
conversion rounds down, thus giving an int uniformly distributed on
0,1,...,(n-1).  We will never get n because uni() is uniformly
distributed on [0,1), not [0,1].
****************************************************************/
#define randint(n)  (uni() * (n))
/*********Defined in unirand.c************/
int          getseed(void);
int          initrand(int seed);
int         *randperm(int *vec, int n);
int          multinomial(real * cum, int n);
real         uni(void);
int          ustart(int iseed);
