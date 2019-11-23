#if 1				/* single precision */
typedef float real;
#define RFMT "%f"
#else
typedef double real;
#define RFMT "%lf"
#endif

typedef struct {
  int          lo;		/* f[lo] is lowest permissible index */
  int          hi;		/* f[hi] is highest permissible index */
  float       *f;
} FLOATVEC;

typedef struct {
  int          lo;		/* f[lo] is lowest permissible index */
  int          hi;		/* f[hi] is highest permissible index */
  double      *f;
} DBLVEC;
