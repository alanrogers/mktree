#define Pi 3.14159265358979323846
typedef enum {
  infinite, finite_gamma, finite_flat, stepwise
} MUTATION_MODEL;

/* defined in mutation.c */
int          getmatch(int init_msize, NODE * tree);
void         crosstab(int depth, NODE * node);
real         init_mutation(MUTATION_MODEL init_mut_model, int init_n_sites,
			   real init_shape, int init_reset);
void         set_gamma_rates(void);
void         clear_mutation(void);
void         mutate(NODE * node, STATE * inherited);
int          getsite(void);
SITE        *copyseq(SITE * s1);
int          ndiffs(NODE * s1, NODE * s2);
int          getsegregating(void);
real         lnfact(int n);
double       poisson(int x, double mean);
double       choose(int n, int k);
real         gammln(real xx);
real         poidev(real xm);
real         gamma_dev(double a);
real         igamma_dev(int ia);
int          is_even(int i);
int          treesumsteps(int *sum, NODE * node);
int          treemoments(DBLVEC * m, double origin, NODE * node);
void         pairwise_moments(DBLVEC * G);
void         stepwise_mom(real * theta0, real * tau, DBLVEC * G);
int          getsteps(int mutations);
real         prob_step(int s, int x);
void         prstepwise(FILE * fp);
int          treemutations(NODE * node);
void         init_node_mutations(NODE * new);
void         prtree(NODE * node, int indent);
