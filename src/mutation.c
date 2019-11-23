#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "mytypes.h"
#include "bye.h"
#include "unirand.h"
#include "alloc2d.h"
#include "getcic.h"
#include "memstack.h"
#include "iscoales.h"
#include "mutation.h"
#define INIT_SITEVAL 0          /* initial value of each site */
#define MUTATE(x) ((x)==0)      /* converts 0 to 1 and vice versa */
#define INBUFF 100
#define KEEPSIZE 100

/*** externals ***/
extern NODE *nodevec;
extern int  sampsize;           /* size of sample, all subdivisions included */
int        *subsize;            /* subsize[i] = size of i'th subdivision */
extern int  nsubs;              /* number of subdivisions in sample */

/** for mutations **/
int         n_sites = 0;
int       **xtab = NULL, xtab_dim = 0;
real       *cumprob = NULL;     /* cumprob[i] = prob that site is <= i */
char       *hit = NULL;         /* hit[i] = 1 (0) if site i has (has not) mutated */
MUTATION_MODEL mut_model = infinite;    /* model of mutation */
int         reset_mut = 0;      /* should mutation rates be reset each time? */
real        gamma_shape;        /* gamma shape parameter */
real        sim_mpd;            /* mean pairwise difference. Set by getmatch() */

/* externals initialized by getmatch() */
int         msize = 0;
int        *mismatch = NULL, ***intermatch = NULL;

/**** Get mismatch and intermatch distributions ******************/
int getmatch(int init_msize, NODE * tree) {
    register int dif, sum;
    int         i, j, ipop, jpop;
    int         size;

    assert(sampsize > 0);
    assert(mut_model == infinite || mut_model == finite_flat
           || mut_model == finite_gamma || mut_model == stepwise);

    /* allocate arrays on first call */
    if(init_msize > msize) {
        if(mismatch != NULL)
            free(mismatch);
        if(intermatch != NULL)
            free(intermatch);
        msize = init_msize;
        mismatch = (int *) mustalloc(msize * sizeof(int));
        intermatch = (int ***) alloc3d(nsubs, nsubs, msize, sizeof(int));
        if(intermatch == NULL)
            error("iscoales: memory");
    }
    memset((void *) mismatch, 0, sizeof(int) * msize);  /* initialize */
    memset((void *) (intermatch[0][0]), 0,
           sizeof(int) * nsubs * nsubs * msize);

    sim_mpd = 0.0;
    switch (mut_model) {
    case finite_flat:
    case finite_gamma:
        for(i = 1; i < sampsize; i++) {
            assert(nodevec[i].pop0 >= 0);
            for(j = 0; j < i; j++) {
                sim_mpd += dif = ndiffs(nodevec + i, nodevec + j);
                if(dif >= msize)
                    dif = msize - 1;    /* dump extreme vals into last entry */
                ++mismatch[dif];    /* accumulate mismatch counts */
                assert(nodevec[j].pop0 >= 0);
                if(nodevec[i].pop0 > nodevec[j].pop0) {
                    ipop = nodevec[i].pop0;
                    jpop = nodevec[j].pop0;
                } else {
                    jpop = nodevec[i].pop0;
                    ipop = nodevec[j].pop0;
                }
                ++intermatch[ipop][jpop][dif];
            }
        }
        break;
    case stepwise:
    case infinite:
        /* put counts in xtab matrix: */
        crosstab(0, tree);

    /****************************************************************
    use xtab matrix to calculate mismatch distribution and mean
    pairwise difference
    ****************************************************************/
        for(i = 1; i < sampsize; i++) {
            for(j = 0; j < i; j++) {
                fflush(stdout);
                sim_mpd += dif = xtab[i][j];
                if(dif >= msize)
                    dif = msize - 1;
                ++mismatch[dif];
                if(nodevec[i].pop0 > nodevec[j].pop0) {
                    ipop = nodevec[i].pop0;
                    jpop = nodevec[j].pop0;
                } else {
                    jpop = nodevec[i].pop0;
                    ipop = nodevec[j].pop0;
                }
                ++intermatch[ipop][jpop][dif];
            }
        }
        break;
#ifndef NDEBUG
    default:
        fprintf(stderr, "\nIllegal switch value in getmatch\n");
        exit(1);
#endif
    }
    sim_mpd /= (sampsize * (sampsize - 1)) / 2;
    /* discard trailing zeroes */
    for(size = msize; size > 0 && mismatch[size - 1] == 0; --size) ;

#if 1
    /* check mismatch distribution */
    for(sum = i = 0; i < size; i++)
        sum += mismatch[i];
    if(sum != (sampsize * (sampsize - 1)) / 2)
        printf("\nmismatch:sum(m)=%d, should = %d",
               sum, (sampsize * (sampsize - 1)) / 2);
    for(ipop = 0; ipop < nsubs; ipop++)
        for(jpop = 0; jpop <= ipop; jpop++) {
            for(sum = i = 0; i < size; i++)
                sum += intermatch[ipop][jpop][i];
            if(ipop == jpop) {
                if(sum != (subsize[ipop] * (subsize[ipop] - 1) / 2)) {
                    fflush(stdout);
                    fprintf(stderr,
                            "\nw/i group mismatch sum=%d; should be %d\n",
                            sum, (subsize[ipop] * (subsize[ipop] - 1) / 2));
                    exit(1);
                }
            } else {
                if(sum != subsize[ipop] * subsize[jpop]) {
                    fflush(stdout);
                    fprintf(stderr,
                            "\nw/i group mismatch sum=%d; should be %d\n",
                            sum, (subsize[ipop] * subsize[jpop]));
                    exit(1);
                }
            }
        }
#endif
    return (size);
}

/****************************************************************
 Create a matrix whose ij'th entry is the difference between individuals
 i and j.
 ****************************************************************/
void crosstab(int depth, NODE * node) {
    int         nd, *descendant, nc, complement[MAX];
    int         increment = 0;
    register int i, j;

    if(node == NULL)
        return;

    if(depth == 0) {            /* initialize xtab array */
        if(sampsize != xtab_dim || xtab == NULL) {  /* allocate */
            if(xtab != NULL)
                free2d((void **) xtab);
            xtab = (int **) alloclt(sampsize, sizeof(int));
            if(xtab == NULL)
                error("crosstab: memory");
            xtab_dim = sampsize;
        }

    /****************************************************************
    fill xtab w/ zeroes: Note that we have
    sampsize*(sampsize*(sampsize+1))/2 entries rather than
    sampsize*(sampsize*(sampsize-1))/2.  This is because alloclt() allocates
    a full lower triangular matrix (including the diagonal) even though
    the diagonal entries are not used.  If I were a better person I would
    eliminate this inefficiency.
    ****************************************************************/
        memset((void *) xtab[0], 0,
               sizeof(int) * (sampsize * (sampsize + 1)) / 2);
    }
    switch (mut_model) {
    case stepwise:
        increment = node->state.steps;
        break;
    case infinite:
        increment = node->mutations;
        break;
    default:
        error("Wrong mutation model in crosstab");
    }

    if(increment > 0) {
        nd = node->d->n;
        descendant = node->d->ndx;
        /* Get complement of list of descendants */
        for(nc = i = j = 0; i < nd; i++) {
            while(j < descendant[i])
                complement[nc++] = j++;
            j++;
        }
        while(j < sampsize)
            complement[nc++] = j++;
        assert(nc + nd == sampsize);

        for(i = 0; i < nd; i++)
            for(j = 0; j < nc; j++) {
                if(descendant[i] > complement[j])
                    xtab[descendant[i]][complement[j]] += increment;
                else
                    xtab[complement[j]][descendant[i]] += increment;
            }
    }
    crosstab(depth + 1, node->left);
    crosstab(depth + 1, node->right);
}

/****************************************************************
init_mutation: call this to initialize mutation model
Parameters
 init_mut_model: specifies mutation model
 init_n_sites  : the number of sites.
 init_shape    : shape parameter of gamma distribution
 init_reset    : if nonzero, gamma rates are reset by each call to iscoales 
init_n_sites and init_shape are both ignored when mut_model==infinite.
****************************************************************/
real init_mutation(MUTATION_MODEL init_mut_model, int init_n_sites,
                   real init_shape, int init_reset) {

/** if mutation vectors are already set then unset them **/
    if(n_sites > 0 && n_sites != init_n_sites) {
        if(cumprob != NULL) {
            free(cumprob);
            cumprob = NULL;
        }
        if(hit != NULL) {
            free(hit);
            hit = NULL;
        }
    }
    mut_model = init_mut_model;
    switch (mut_model) {
    case finite_gamma:
        reset_mut = init_reset;
        gamma_shape = init_shape;
        cumprob = (real *) mustalloc(init_n_sites * sizeof(real));
        n_sites = init_n_sites;
        hit = (char *) mustalloc(n_sites * sizeof(char));

    /****************************************************************
    If reset_mut is nonzero, gamma rates will be set at top of iscoales.
    Otherwise, set them here.
    ****************************************************************/
        if(!reset_mut)
            set_gamma_rates();
        break;
    case finite_flat:
        n_sites = init_n_sites;
        hit = (char *) mustalloc(n_sites * sizeof(char));
        break;
    case stepwise:
    case infinite:
        n_sites = 0;
        break;
    default:
        fflush(stdout);
        fprintf(stderr, "\ninit_mutation: illegal value of mut_model");
        exit(1);
    }
    return (0.0);               /* meaningless return value */
}

/****************************************************************
SET GAMMA MUTATION RATES:
The mutational time scale requires that the sum of mutation rates
across sites equal unity so that there will on average be one mutation
per unit of mutational time.  To accomplish this goal, we first define
$y\equiv x/\beta$, where $x$ is gamma-distributed with scale parameter
$\beta$ and shape parameter $\alpha$.  The variable $y$ is also
gamma-distributed, with density 
\begin{displaymath} 
      y^{\alpha-1} e^{-y} / \Gamma(\alpha) 
\end{displaymath} 
We first use this density to generate a vector
$(y_1,y_2,\ldots{},y_K)$.  These variates are each equal to $x/\beta$,
where $x$ is the gamma-distributed variate we seek and $\beta$ is
arbitrary.  To satisfy the constraint imposed by the mutational time
scale, we set $\beta = 1/\sum_i y_i$.  In other words, we set $x_i =
y_i / \sum_i y_i$.  This provides a vector of gamma-distributed
mutation rates that sums to unity, as required.
****************************************************************/
void set_gamma_rates(void) {
    int         i;
    real        b;

    b = 0.0;
    for(i = 0; i < n_sites; i++)
        b += cumprob[i] = gamma_dev((double) gamma_shape);
    cumprob[0] /= b;
    for(i = 1; i < n_sites; i++)
        cumprob[i] = cumprob[i - 1] + cumprob[i] / b;
}

/**** call before building each tree ****/
void clear_mutation(void) {

    sim_mpd = -1.0;             /* indicates that sim_mpd has not been calculated */
    switch (mut_model) {
    case finite_gamma:
        if(cumprob == NULL) {
            fprintf(stderr,
                    "\nerror: call init_mutation before clear_mutation\n");
            exit(1);
        }
        if(reset_mut)           /* set gamma mutation rates on each call */
            set_gamma_rates();

/** fall through to finite_flat **/
    case finite_flat:
        if(hit == NULL) {
            fprintf(stderr,
                    "\nerror: call init_mutation before clear_mutation\n");
            exit(1);
        }
        memset(hit, 0, n_sites);    /* set to zero */
        break;
    case stepwise:
    case infinite:
        n_sites = 0;
        break;
#ifndef NDEBUG
    default:
        fflush(stdout);
        fprintf(stderr, "\nclear_mutation: illegal value of mut_model");
        exit(1);
#endif
    }
}

/* put mutations onto tree */
void mutate(NODE * node, STATE * inherited) {
    int         site, nmut;

    if(node == NULL)
        return;
    /* # of mutations */
    node->mutations = nmut = poidev(0.5 * node->branch);

    switch (mut_model) {
    case stepwise:
        if(inherited == NULL)
            node->state.steps = getsteps(node->mutations);
        else
            node->state.steps = inherited->steps + getsteps(node->mutations);
        break;
    case infinite:
        n_sites += nmut;
        break;
    case finite_gamma:
    case finite_flat:
        if(node->mutations == 0 && inherited == NULL)
            node->state.seq = NULL;
        else
            node->state.seq = copyseq(inherited->seq);
        while(nmut-- > 0) {
            /* choose a site to mutate */
            site = getsite();
            hit[site] = 1;
            node->state.seq[site] = MUTATE(node->state.seq[site]);
        }
        break;
    default:
        error("mutate: Illegal mutation model");
    }
    mutate(node->left, &(node->state));
    mutate(node->right, &(node->state));
}

/* choose a site to mutate under one of the finite sites models */
int getsite(void) {
    real        r;
    int         site, lo, mid, hi;

    switch (mut_model) {
    case finite_flat:
        site = randint(n_sites);
        break;
    case finite_gamma:
        /* choose a site to mutate by binary search */
        r = uni();              /* uniform random variable */
        lo = 0;
        hi = n_sites - 1;
        while(lo < hi) {
            mid = lo + (hi - lo) / 2;
            if(mid == lo) {
                if(cumprob[lo] < r)
                    lo = hi;
                break;
            }
            if(cumprob[mid] <= r)
                lo = mid;
            else
                hi = mid;
        }
        return (lo);
        break;
    default:
        fflush(stdout);
        fprintf(stderr, "\nIllegal mutation model in getsite()\n");
        exit(1);
    }
    return (site);
}

/****************************************************************
Allocate and initialize a new DNA (or restriction site) sequence.  
If s1==NULL then on return each site of the new sequence will equal 
INIT_SITEVAL.  Otherwise, sequence s1 will be copied into the new
sequence.
****************************************************************/
SITE       *copyseq(SITE * s1) {
    int         i;
    SITE       *s;

    s = (SITE *) stackalloc(n_sites * sizeof(SITE));
    if(s == NULL) {
        fflush(stdout);
        fprintf(stderr, "\ncopyseq: stackalloc failed\n");
        exit(1);
    }
    if(s1 == NULL)
        for(i = 0; i < n_sites; i++)
            s[i] = INIT_SITEVAL;
    else
        for(i = 0; i < n_sites; i++)
            s[i] = s1[i];
    return (s);
}

/* count diffs between two nodes */
int ndiffs(NODE * s1, NODE * s2) {
    register int i, diffs = 0.0;

    if(s1 != NULL && s2 != NULL) {  /* neither == NULL */
        switch (mut_model) {
        case infinite:
            i = s1->mutations - s2->mutations;
            return (i >= 0 ? i : -i);
        case finite_gamma:
        case finite_flat:
            for(i = 0; i < n_sites; i++)
                if(s1->state.seq[i] != s2->state.seq[i])
                    diffs++;
            return (diffs);
        case stepwise:
            i = s1->state.steps - s2->state.steps;
            return (i >= 0 ? i : -i);
        }
    }
    /* if we get this far then at least one is NULL */
    if(s2 != NULL) {            /* s1 == NULL; s2 != NULL */
        switch (mut_model) {
        case infinite:
            return (s2->mutations);
        case finite_gamma:
        case finite_flat:
            for(i = 0; i < n_sites; i++)
                if(s2->state.seq[i] != INIT_SITEVAL)
                    diffs++;
            return (diffs);
        case stepwise:
            return (s2->state.steps);
        }
    }
    if(s1 != NULL) {            /* s1 != NULL; s2==NULL */
        switch (mut_model) {
        case infinite:
            return (s1->mutations);
        case finite_gamma:
        case finite_flat:
            for(i = 0; i < n_sites; i++)
                if(s1->state.seq[i] != INIT_SITEVAL)
                    diffs++;
            return (diffs);
        case stepwise:
            return (s1->state.steps);
        }
    }
    /* both are NULL */
    return (0);
}

/* count segregating sites */
int getsegregating(void) {
    int         i, segregating = 0;

    switch (mut_model) {
    case stepwise:
        break;
    case infinite:
        segregating = n_sites;
        break;
    case finite_gamma:
    case finite_flat:
        for(i = 0; i < n_sites; i++)
            if(hit[i] > 0)
                ++segregating;
        break;
    }
    return (segregating);
}

/* log of factorial */
real lnfact(int n) {
    static int  first_call = 1;
    static real a[KEEPSIZE];
    int         j;

    if(first_call) {
        first_call = 0;
        for(j = 0; j < KEEPSIZE; j++)
            a[j] = 0.0;
    }
    assert(n >= 0);

    if(n <= 1)
        return (0.0);
    if(n >= KEEPSIZE)
        return (gammln(n + 1.0));
    if(a[n] == 0.0)
        a[n] = gammln(n + 1.0);
    return a[n];
}

/* poisson probability function */
double poisson(int x, double mean) {
    double      lnprob;

    assert(x >= 0);
    assert(mean >= 0.0);
    switch (x) {                /* for small x, do it the simple way */
    case 0:
        return (exp(-mean));
    case 1:
        return (mean * exp(-mean));
    case 2:
        return (mean * mean * exp(-mean) / 2.0);
    case 3:
        return (mean * mean * mean * exp(-mean) / 6.0);
    }

    /* for larger x, work with logs to avoid large numbers */
    lnprob = x * log(mean) - mean - lnfact(x);
    return (exp(lnprob));
}

/* binomial coefficient */
double choose(int n, int k) {
    assert(n >= 0);
    assert(k >= 0);
    assert(k <= n);
    return (floor(0.5 + exp(lnfact(n) - lnfact(k) - lnfact(n - k))));
}

/* log gamma, this came from Numerical Recipes */
real gammln(real xx) {
    double      x, tmp, ser;

    static double cof[6] = { 76.18009173, -86.50532033, 24.01409822,
        -1.231739516, 0.120858003e-2, -0.536382e-5
    };
    int         j;

    x = xx - 1.0;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.0;
    for(j = 0; j <= 5; j++) {
        x += 1.0;
        ser += cof[j] / x;
    }
    return -tmp + log(2.50662827465 * ser);
}

/* Poisson random deviate, from Numerical Recipes */
real poidev(real xm) {
    static real sq, alxm, g, oldm = (-1.0);
    real        em, t, y;

    if(xm < 12.0) {
        if(xm != oldm) {
            oldm = xm;
            g = exp(-xm);
        }
        em = -1;
        t = 1.0;
        do {
            em += 1.0;
            t *= uni();
        } while(t > g);
    } else {
        if(xm != oldm) {
            oldm = xm;
            sq = sqrt(2.0 * xm);
            alxm = log(xm);
            g = xm * alxm - gammln(xm + 1.0);
        }
        do {
            do {
                y = tan(Pi * uni());
                em = sq * y + xm;
            } while(em < 0.0);
            em = floor(em);
            t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - g);
        } while(uni() > t);
    }
    return (em);
}

/****************************************************************
Random deviates from standard gamma distribution with density
         a-1
        x    exp[ -x ]
f(x) = ----------------
         Gamma[a]

where a is the shape parameter.  The algorithm for integer a comes
from numerical recipes, 2nd edition, pp 292-293.  The algorithm for
a<1 uses code from p 213 of Statistical Computing, by Kennedy and
Gentle, 1980 edition.  This algorithm was originally published in:

Ahrens, J.H. and U. Dieter (1974), "Computer methods for sampling from
Gamma, Beta, Poisson, and Binomial Distributions".  COMPUTING
12:223-246.

The mean and variance of these values are both supposed to equal a.
My tests indicate that they do.

This algorithm has problems when a is small.  In single precision, the
problem  arises when a<0.1, roughly.  That is why I have declared
everything as double below.  Trouble is, I still don't know how small 
a can be without causing trouble.  f(x) doesn't seem to have a series
expansion around x=0.
****************************************************************/
real gamma_dev(double a) {
    int         ia;
    double      u, b, p, x, y = 0.0, recip_a;

    if(a <= 0) {
        fprintf(stderr, "\ngamma_dev: parameter must be positive\n");
        exit(1);
    }
    ia = floor(a);              /* integer part */
    a -= ia;                    /* fractional part */
    if(ia > 0) {
        y = igamma_dev(ia);     /* gamma deviate w/ integer argument ia */
        if(a == 0.0)
            return (y);
    }
    /* get gamma deviate with fractional argument "a" */
    b = (M_E + a) / M_E;
    recip_a = 1.0 / a;
    for(;;) {
        u = uni();
        p = b * u;
        if(p > 1) {
            x = -log((b - p) / a);
            if(uni() > pow(x, a - 1))
                continue;
            break;
        } else {
            x = pow(p, recip_a);
            if(uni() > exp(-x))
                continue;
            break;
        }
    }
    return (x + y);
}

/****************************************************************
gamma deviate for integer shape argument.  Code modified from pp
292-293 of Numerical Recipes in C, 2nd edition.
****************************************************************/
real igamma_dev(int ia) {
    int         j;
    real        am, e, s, v1, v2, x, y;

    if(ia < 1) {
        fprintf(stderr, "\nerror: arg of igamma_dev was <1\n");
        exit(1);
    }
    if(ia < 6) {
        x = 1.0;
        for(j = 0; j < ia; j++)
            x *= uni();
        x = -log(x);
    } else {
        do {
            do {
                do {            /* next 4 lines are equivalent */
                    v1 = 2.0 * uni() - 1.0; /* to y = tan(Pi * uni()).     */
                    v2 = 2.0 * uni() - 1.0;
                } while(v1 * v1 + v2 * v2 > 1.0);
                y = v2 / v1;
                am = ia - 1;
                s = sqrt(2.0 * am + 1.0);
                x = s * y + am;
            } while(x <= 0.0);
            e = (1.0 + y * y) * exp(am * log(x / am) - s * y);
        } while(uni() > e);
    }
    return (x);
}

/* is the integer even? */
int is_even(int i) {
    return (2 * (i / 2) == i);
}

/**************************************************************** 
Sum step variable across all nodes in tree.  On return, m.f[i] = sum
of steps and function returns number of terms in sum.
****************************************************************/
int treesumsteps(int *sum, NODE * node) {
    int         n;

    if(node == NULL)
        return (0);

    n = treesumsteps(sum, node->left) + treesumsteps(sum, node->right);

    if(n > 0)                   /* We're not at tip. */
        return (n);             /* Only tip nodes count. */

    *sum += node->state.steps;

    return (1);
}

/**************************************************************** 
Calculate moments about the value "origin" by traversing the tree.
On return, m.f[i] = sum of (steps-origin)^i and function returns number
of terms in sum.  Moments are calculated by:

  DBLVEC *m;
  double origin;

  m = newdblvec(1,4);
  for(i=m->lo; i<=m->hi; i++)  
    m->f[i] = 0.0;
  origin = 0.0       
  n = treemoments(m, origin, tree);
  for(i=m->lo; i<=m->hi; i++)
    m->f[i] /= n;

****************************************************************/
int treemoments(DBLVEC * m, double origin, NODE * node) {
    int         i, lo = 1, n = 0;
    double      v = 1.0, diff;

    if(node == NULL)
        return (0);

    if(node->left != NULL)
        n += treemoments(m, origin, node->left);

    if(node->right != NULL)
        n += treemoments(m, origin, node->right);

    if(n > 0)                   /* We're not at tip. */
        return (n);             /* Only tip nodes count. */

    if(m->lo > 1)
        lo = m->lo;
    assert(lo <= m->hi);
    diff = node->state.steps - origin;

    for(i = 1; i < lo; i++)     /* set v = (steps^(lo-1) */
        v *= diff;

    for(i = lo; i <= m->hi; i++) {

    /******************************************
    on successive iterations of loop, v equals
    diff, diff^(lo), diff^(lo+1), and so on
    *******************************************/
        v *= diff;
        m->f[i] += (double) v;  /* sum of steps^(i+1) */
    }
    return (1);
}

/****************************************************************
Convert moments (about origin) of the distribution across chromosomes
of repeat counts into central moments of the distribution of 
pairwise differences in repeat counts.  On entry, G->f[i] is the
expectation of steps^i.  On return, 

  G->f[1] = 0
  G->f[2] = E[ (x-y)^2 ]
  G->f[3] = 0
  G->f[4] = E[ (x-y)^4 ]
  G->f[5] = 0
  G->f[6] = E[ (x-y)^6 ]

where x and y are the repeat counts of random individuals in the
populations.
****************************************************************/
void pairwise_moments(DBLVEC * G) {
    assert(G->lo <= 1);
    assert(G->hi >= 2);

    if(G->hi >= 6) {            /* 6th pairwise moment */
        G->f[6] = 2.0 * G->f[6] - 12.0 * G->f[1] * G->f[5]
            + 30.0 * G->f[2] * G->f[4]
            - 20.0 * (G->f[3]) * (G->f[3]);
    }
    if(G->hi >= 4) {            /* 4th pairwise moment: */
        G->f[4] = 6.0 * G->f[2] * G->f[2]
            - 8.0 * G->f[1] * G->f[3] + 2.0 * G->f[4];
    }
    /* 2 * 2nd moment: */
    G->f[2] = 2.0 * G->f[2] - 2.0 * G->f[1] * G->f[1];

    /* odd moments */
    G->f[5] = G->f[3] = G->f[1] = 0.0;

}

/**************************************************************** 
On input: G->f[2] and G->f[4] contain the 2nd and 4th moments of
the distribution of (x-y) the pairwise differences in step counts.
On return, F[0] and F[1] contain method-of-moments estimates 
of theta0 and tau, assuming the one-step mutational model.
****************************************************************/
void stepwise_mom(float *theta0, float *tau, DBLVEC * G) {
    assert(G->lo <= 2);
    assert(G->hi >= 4);

    /* theta0 = sqrt( (G4 - G2)/3 - G2^2 ) */
    *theta0 = ((G->f[4]) - (G->f[2])) / 3.0 - (G->f[2]) * (G->f[2]);
    if(*theta0 < 0.0)
        *theta0 = 0.0;          /* avoid imaginary result */
    else
        *theta0 = sqrt(*theta0);

    /* tau = theta0 - F->f[0] */
    *tau = G->f[2] - *theta0;
}

/* convert mutations into steps */
int getsteps(int mutations) {
    int         steps;
    real        u, sum, p;

    if(mutations == 0)          /* shortcut */
        return (0);

    u = uni();                  /* uniform r.v. */
    sum = 0.0;
    steps = -(mutations + 1);
    do {
        sum += p = prob_step(++steps, mutations);
    } while(u > sum && steps < mutations);
    return (steps);
}

/****************************************************************
Return the probability of observing a step value of s given that
x mutations have occurred and each mutational change is either
+1 or -1, the two possibilities have equal probability 1/2.

Under this symetrical 1-step model, the change after one mutation
is 2*b - 1, where b is a Bernoulli random variable.  The change
after x mutations is

  s(x) = 2*B(x) - x

where B(x) is binomial with parameters x and 1/2.  Thus, the
probability that s(x) equals, say, s is equal to the probability that
B(x) equals (x + s)/2.

Note that s(x) + x = 2*B and must therefore be even.  The probability
of s(x) is 0 if this sum is odd.
****************************************************************/
real prob_step(int s, int x) {
    static int  length = 1;
    static double pow2[KEEPSIZE] = { 1.0 }; /* pow2[x] = 1/2^x */
    float       p2;

    assert(x >= 0);
    /* return Pr[s|x] = 0 unless s+x is even */
    if(!is_even(s + x))
        return (0.0);

    if(s < -x || s > x)         /* out of range */
        return (0.0);

    /* store 1/2^x values to minimize calls to pow function */
    while(length <= x && x < KEEPSIZE) {
        pow2[length] = 0.5 * pow2[length - 1];
        length += 1;
    }
    if(x < length)
        p2 = pow2[x];
    else
        p2 = pow(0.5, (double) x);

    return (choose(x, (x + s) / 2) * p2);
}

/* describe the model of mutation */
void prstepwise(FILE * fp) {
    fprintf(fp, "\n#Stepwise mutation model with mutational distribution:");
    fprintf(fp, "\n Pr[-1] = Pr[+1] = 1/2");
}

/**** count mutations in tree ****/
int treemutations(NODE * node) {
    if(node == NULL)
        return (0);
    return (node->mutations + treemutations(node->left)
            + treemutations(node->right));
}

/*** initialize mutations of a new node ***/
void init_node_mutations(NODE * new) {
    new->mutations = 0;
    switch (mut_model) {
    case infinite:             /* do nothing */
        break;
    case finite_flat:
    case finite_gamma:
        new->state.seq = NULL;
        break;
    case stepwise:
        new->state.steps = 0;
        break;
    default:
        error("bad switch value in init_node_mutations");
    }
}

void prtree(NODE * node, int indent) {
    int         i;

    if(node == NULL)
        return;
    putchar('\n');
    for(i = 0; i < indent; i++)
        putchar('-');
    printf("lngth=%f pop=%d mut=%d", node->branch, node->pop,
           node->mutations);
    switch (mut_model) {
    case stepwise:
        printf(" steps=%d", node->state.steps);
        break;
    case finite_gamma:
    case finite_flat:
        putchar(':');
        for(i = 0; i < n_sites; i++)
            printf("%c", node->state.seq[i]);
        break;
    default:
        break;
    }
    prtree(node->left, indent + 2);
    prtree(node->right, indent + 2);
}
