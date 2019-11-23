#define MAX 1000
#define INFTY -1.0

typedef struct list {
  int          n;		/* # of descendants */
  int         *ndx;		/* array of indices of descendants */
} LIST;
typedef char SITE;

/* state variable's interpretation depends on mutational model */
/* under infinite sites model, state variable = # of mutations */
typedef union {
  SITE        *seq;		/* sequence sites : finite sites models */
  int          steps;		/* # of steps: stepwise model */
} STATE;


typedef struct node {
  real         branch;
  int          pop;		/* which pop is this node currently in? */
  int          pop0;		/* which pop is did node start in? */
  STATE        state;		/* state of chromosome */
  int          mutations;
/*****************************************************************/
  struct list *d;		/* descendants of this node */
  struct node *left, *right, *ancestor;
} NODE;

typedef struct pophist {
  real         mn;		/* m*n, migrants per generation */
  real         theta;		/* 2*N*u where N is subpop size */
  real         tau;		/* 2*u*time in generations */
  /* for infinity, set tau=-1 */
  int          K;		/* # of subpopulations */
  struct pophist *next;		/* next pophist parameters */
} POPHIST;


/* prototypes defined in iscoales.h */
NODE        *newnode(NODE * n1, NODE * n2);
#ifndef NDEBUG
void         check_s(char *msg, int *s, int K, NODE ** node, int S, int SS);
#endif
real         getr(real S, real SS, real mn, real theta);
NODE        *iscoales(int init_nsubs, int *init_subsize, POPHIST * history,
		      int domut);
POPHIST     *gethistory(FILE * fp);
POPHIST     *newhistory(real theta, real mn, real tau, int K);
void         freehistory(POPHIST * ph);
POPHIST     *duphist(POPHIST * new, POPHIST * old);
void         writehistory(FILE * fp, POPHIST * history, char comment);
int          icompar(int *x, int *y);
int          get_collapse_vector(int K1, int K2, int *v);
void         prtree(NODE * node, int indent);
real         treedepth(NODE * node);
real         treelength(NODE * node);
void         clear_externals(void);
LIST        *newlist(int n, int *ndx);
void         prlist(LIST * l);
LIST        *mergelist(LIST * l1, LIST * l2);
