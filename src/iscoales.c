#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "mytypes.h"
#include "bye.h"
#include "unirand.h"
#include "memstack.h"
#include "alloc2d.h"
#include "getcic.h"
#include "iscoales.h"
#include "mutation.h"
#define INBUFF 100
#define KEEPSIZE 100

/*** externals ***/
NODE       *nodevec = NULL;
int         nnodes = 0, nextnode = 0;
int         sampsize = 0;       /* size of sample, all subdivisions included */
int         nsubs = 0;          /* number of subdivisions in sample */
int        *subsize = NULL;     /* subsize[i] = size of i'th subdivision */
NODE     ***smpl = NULL;        /* smpl[i][j] -> j'th sample from i'th subdivision */

/* Join two nodes into a single node */
NODE       *newnode(NODE * n1, NODE * n2) {
    NODE       *new;

#ifndef NDEBUG
    if(n1 != NULL && n2 != NULL && n1->pop != n2->pop) {
        fprintf(stderr,
                "\nError in newnode: %d = n1->pop != n2->pop=%d\n",
                n1->pop, n2->pop);
        exit(1);
    }
#endif
    new = (NODE *) nodevec + nextnode;
    ++nextnode;
    new->branch = 0.0;
    new->left = n1;
    new->right = n2;
    new->d = NULL;
    new->ancestor = NULL;
    new->mutations = 0;
    init_node_mutations(new);
    if(n1 != NULL)
        n1->ancestor = new;
    if(n2 != NULL)
        n2->ancestor = new;
    if(n1 != NULL)
        new->pop = n1->pop;
    else if(n2 != NULL)
        new->pop = n2->pop;
    else
        new->pop = -1;
    if(n1 != NULL)
        new->d = n1->d;
    if(n2 != NULL)
        new->d = mergelist(new->d, n2->d);
    return (new);
}

#ifndef NDEBUG
void check_s(char *msg, int *s, int K, NODE ** node, int S, int SS) {
    int         sum, sum2, count[MAX], i, bad = 0;

    for(i = 0; i < K; i++)
        count[i] = 0;
    for(i = 0; i < S; i++) {
        if((node[i]->pop) < 0 || (node[i]->pop) >= K) {
            fprintf(stderr, "\ncheck_s:node[%d]->pop=%d. Legal=0..%d",
                    i, node[i]->pop, K - 1);
            bad = 1;
        }
        count[node[i]->pop] += 1;
    }
    for(sum = sum2 = i = 0; i < K; i++) {
        sum += count[i];
        sum2 += count[i] * count[i];
    }
    if(sum != S) {
        fprintf(stderr, "\ncheck_s: sum=%d != %d=S", sum, S);
        bad = 1;
    }
    if(sum2 != SS) {
        fprintf(stderr, "\ncheck_s: sum2=%d != %d=SS", sum2, SS);
        bad = 1;
    }
    for(i = 0; i < K; i++) {
        if(count[i] != s[i]) {
            fprintf(stderr, "\ncheck_s: count[%d]=%d!=%d=s[%d]",
                    i, count[i], s[i], i);
            bad = 1;
        }
    }
    if(bad) {
        fprintf(stderr, "\ncheck_s: %s, S=%d\ns=", msg, S);
        for(i = 0; i < K; i++)
            fprintf(stderr, " %d", s[i]);
        exit(1);
    }
}
#endif

real getr(real S, real SS, real mn, real theta) {
    return ((S * mn + 0.5 * (SS - S)) / theta);
}

/****************************************************************
Build a random tree using coalescent w/ geographic structure and
island model migration.

Let s[i] = # of nodes in pop i
    S = sum(s[i]) = total nodes
    SS = sum(s[i]*s[i]), the sum of squares
    n = size of a subpopulation
    K = # of subpopulations
    m = island model migration rate: Pr a haploid individual is an
        immigrant from some other subdivision.
    a = the hazard of an event at time t (i.e. the conditional 
        probability density that an event will happen at t given 
        that it did not happen before that) 

       K-1
       ---
       \
    a = >   (  s[i]*m + (s[i]*(s[i]-1)/2) / n )
       /__
       i=0
    
      = S*m + ( SS - S)/(2*n)

I will only need a in calculating a*t, where t is time in generations.
This is equal to

      2*u*t
a*t = ----- * [ S*(mn) + (SS - S)/2 ]
      2*u*n

Thus, if I use mutational scale, defined by tau = 2*u*t, then the
parameters are:  theta=2*u*n and mn = m*n.  In mutational time, the
hazard "a" becomes 

r = [ S*(mn) + (SS - S)/2 ]/theta

The time of the next event is an exponential r.v. with rate r.  The
cdf is 1 - exp(-r*x).  Since any cdf is uniformly distributed on
(0,1), one can generate the exponential variate as follows.  First
generate a uniform variate u, and then solve u = 1 - exp(-r*x) to get 
x = -ln(1-u)/r.

Given that an event has occurred, it is a migration with probability 

 pmig = S*(mn)/[ S*(mn) + (SS - S)/2]

When a coalescent event occurs in group i, the following adjustments
are made:
1. 2 nodes are joined to form a single new node
2. S -= 1;
3. SS -= 2*s[i] - 1;
4. s[i] -= 1;
Here, step 3 is equivalent to

  SS' = SS - s[i]^2 + (s[i] - 1)^2

When a migration from i to j occurs, the following adjustments are
made:
1. SS += 2*(s[j] - s[i] + 1);
2. s[i] -= 1;
3. s[j] += 1;
Here, step 1 is equivalent to

  SS' = SS - s[i]^2 + (s[i]-1)^2 - s[j]^2 + (s[j]+1)^2

ON ENTRY:
  sampsize: # of nodes per subpopulation
  domut   : If domut is nonzero, mutations will be placed on tree.
****************************************************************/
NODE       *iscoales(int init_nsubs, int *init_subsize, POPHIST * history,
                     int domut) {
    POPHIST    *ph;
    NODE      **node;
    int         s[MAX];         /*s[i]=#of nodes in pop i */
    int         S, SS;
    int         i, j, k, ii, jj, ipop, jpop, oldK, kk, ndx;
    int         v[MAX];
    real        t_g;            /*time since last migration or coalescence */
    real        t_p;            /*time since last change in population parameters */
    real        c[MAX];         /* cumulative distribution function */
    real        pmig;
    real        x, r;

    clear_mutation();

  /****************************************************************
  Set permanent arrays, which do *not* get cleared each time.  If
  nsubs hasn't changed, I assume w/o checking that
  subdivision sizes haven't changed either.
  ****************************************************************/
    if(nsubs != init_nsubs) {
        clear_externals();
        nsubs = init_nsubs;
        subsize = (int *) mustalloc(nsubs * sizeof(int));
        sampsize = 0;
        for(i = 0; i < nsubs; i++)
            sampsize += subsize[i] = init_subsize[i];
        nnodes = 2 * sampsize - 1;  /* total nodes in tree */
        /* nodevec: pointers that WILL get reshuffled */
        nodevec = (NODE *) mustalloc(nnodes * sizeof(NODE));

        /* smpl: pointers that will not be reshuffled by coalescent events */
        smpl = (NODE ***) mustalloc(nsubs * sizeof(NODE **));
        for(i = 0; i < nsubs; i++)
            smpl[i] = (NODE **) mustalloc(subsize[i] * sizeof(NODE *));
    }

  /****************************************************************
  The number of subdivisions in the sample may not match the number
  in the population.  If not, I add another epoch to the population
  history with length 0.  The values of mn and theta don't matter
  since this epoch has 0 length anyway.
  ****************************************************************/
    if(history->K != nsubs) {
        ph = newhistory(1.0, 1.0, 0.0, nsubs);  /*new epoch */
        ph->next = history;     /*attach old history */
        history = ph;
    }

  /****************************************************************
  The newnode function doesn't allocate memory; it just grabs the
  next available node from the nodevec vector.  Setting nextnode=0
  tells newnode to start at the beginning again, effectively freeing
  all the nodes that were used last time.
  ****************************************************************/
    nextnode = 0;
    node = (NODE **) stackmustalloc(sampsize * sizeof(NODE *));
    for(k = j = 0; j < nsubs; j++) {
        for(i = 0; i < subsize[j]; i++) {
            smpl[j][i] = node[k] = newnode(NULL, NULL);
            node[k]->d = newlist(1, &k);
            node[k]->pop = j;
            node[k]->pop0 = j;
            k += 1;
        }
        s[j] = subsize[j];
    }

  /****************************************************************
  At this point we can access the nodes in the sample either through
  smpl[j][i] or through node[k].  After the building the tree, the
  entries of node[] will will be scrambled but those of smpl[][] will
  still point to the sample.
  ****************************************************************/

    S = sampsize;               /* number of nodes */
    for(SS = 0.0, i = 0; i < nsubs; i++)
        SS += subsize[i] * subsize[i];  /*sum of squared sample sizes */
    if(S > MAX) {
        fprintf(stderr, "\nError: trees can only have %d tips\n", MAX);
        exit(1);
    }
    assert(S == sampsize);
#ifndef NDEBUG
    check_s("top of iscoales", s, history->K, node, S, SS);
#endif
    if(history->tau == 0.0)
        history->theta = 1.0;

    /* initial value of rate r of exponential random variable */
    r = getr(S, SS, history->mn, history->theta);

    /* Iterate until sample size (S) equals 1 */
    t_g = t_p = 0.0;
    while(S > 1) {
        x = -log(1.0 - uni());  /* An exponential random variate. */
        /* translate x into units of 1/(2u) generations */
        while(history->next != NULL && (r == 0 || t_p + x / r > history->tau)) {

      /****************************************************************
      1st line below subtracts off portion of x that was "used up" by the
      epoch we are about to leave.
      ****************************************************************/
            x -= r * (history->tau - t_p);

      /****************************************************************
      The rest of this loop resets parameters to reflect the new epoch
      and recalculates r.
      ****************************************************************/
            t_g += history->tau - t_p;
            t_p = 0.0;
            oldK = history->K;
            history = history->next;
            if(history->tau == 0.0)
                history->theta = 1.0;
            if(history->K < oldK) { /*  reduce # of subpopulations */
                /* get vector to map oldK old pops into history->K new ones */
                get_collapse_vector(oldK, history->K, v);
                for(i = 0; i < oldK; i++)
                    s[i] = 0;
                for(i = 0; i < S; i++) {
                    node[i]->pop = v[node[i]->pop]; /* reset population values */
                    s[node[i]->pop]++;  /* recalc pop sample sizes */
                }
                SS = 0;         /* reset sum of squares, SS */
                for(i = 0; i < history->K; i++)
                    SS += s[i] * s[i];
            } else if(history->K > oldK) {  /* increase # of supopopulations */

    /***************************************************
	Initially, the new groups will be empty: s[i]=0.
	These groups gain members only by migration.   
	This assumes that, in forwards time, the reduction in group
	numbers was caused by group extinction rather than
	coalescence. 
        ****************************************************/
                for(i = oldK; i < history->K; i++)
                    s[i] = 0;
            }
            /* reset r for new history parameters */
            r = getr(S, SS, history->mn, history->theta);
        }
        if(r == 0)
            error("iscoales: NO COALESCENCE\n");
        t_p += x / r;           /* time since last change in history parameters */
        t_g += x / r;           /* time since last coalescence */

#ifndef NDEBUG
        check_s("after setting t_p & t_g", s, history->K, node, S, SS);
#endif

        /* Classify current event */
        pmig = (real) S *(history->mn);

        assert(pmig + 0.5 * (SS - S) > 0.0);
        pmig = pmig / (pmig + 0.5 * (SS - S));  /* Pr[event is a migration] */
        if(history->mn > 0 && uni() < pmig) {
            /* event is a migration: choose migrant */
            i = randint(S);
            ipop = node[i]->pop;
            assert(s[ipop] > 0);
            /* choose group of origin */
            jpop = randint((history->K) - 1);
            if(jpop >= ipop)
                jpop++;
            assert(ipop < history->K);
            assert(jpop < history->K);
            /* move node i */
            node[i]->pop = jpop;
            /* reset parameters */
            SS += 2 * (s[jpop] - s[ipop] + 1);
            s[ipop] -= 1;
            s[jpop] += 1;
#if 0
            check_s("after migration", s, history->K, node, S, SS);
#endif
            continue;
        } else {                /* Event is a coalescence. */
            /* record time */
            for(i = 0; i < S; i++)
                node[i]->branch += t_g;
            t_g = 0.0;
            assert(SS > S);
            /* Choose group */
            x = SS - S;
            c[0] = s[0] * (s[0] - 1) / x;
            for(i = 1; i < history->K; i++)
                c[i] = c[i - 1] + s[i] * (s[i] - 1) / x;
            c[history->K - 1] = 1.0;
            x = uni();
            for(ipop = 0; ipop < history->K && c[ipop] <= x; ipop++) ;
            /* Choose two nodes w/i group ipop */
            assert(s[ipop] >= 2);
            i = randint(s[ipop]);
            j = randint(s[ipop] - 1);
            if(j >= i)
                j += 1;
            /* find nodes i & j w/i group ipop */
            ndx = ii = jj = -1;
            for(kk = 0; kk < S; kk++) {
                if(node[kk]->pop == ipop) {
                    ndx++;
                    if(ndx == i)
                        ii = kk;
                    else if(ndx == j)
                        jj = kk;
                    if(ii >= 0 && jj >= 0)
                        break;
                }
            }
#ifndef NDEBUG
            if(ii < 0 || jj < 0) {
                fprintf(stderr, "\nbad ii or jj: ii=%d jj=%d ipop=%d", ii, jj,
                        ipop);
                fprintf(stderr, "  s[%d]=%d", ipop, s[ipop]);
                fprintf(stderr, "\n   node[*]->pop=");
                for(i = 0; i < S; i++) {
                    fprintf(stderr, " %d", node[i]->pop);
                    if((i + 1) % 20 == 0)
                        putc('\n', stderr);
                }
                exit(1);
            }
            if(node[ii]->pop != node[jj]->pop) {
                fprintf(stderr,
                        "\niscoales: node[%d]->pop=%d != %d=node[%d]->pop",
                        ii, node[ii]->pop, node[jj]->pop, jj);
                fprintf(stderr, " S=%d\n", S);
                exit(1);
            }
#endif
            /* join the two nodes and shorten vector */
            node[ii] = newnode(node[ii], node[jj]);
            node[jj] = node[S - 1];

            /* reset parameters */
            S -= 1;
            SS -= 2 * s[ipop] - 1;
            s[ipop] -= 1;
        }
        /* reset r because SS and maybe S have changed */
        r = getr(S, SS, history->mn, history->theta);
    }
    if(domut)
        mutate(node[0], NULL);    /* add mutations to tree */
    return (node[0]);
}

/****************************************************************
Read file fp, which should contain one line of data for each time
period, with the most recent time period first.  Each line should
contain four items of data: theta, mn, tau, and K.  Here theta=4*u*n*K
where n is the size of the population as a whole.

The function returns a pointer to the base of a linked list of type
POPHIST. 
****************************************************************/
POPHIST    *gethistory(FILE * fp) {
    real        mn, theta, tau;
    int         i, K, got_infinity = 0;
    POPHIST    *root = NULL, *ph, *last = NULL;
    char        buff[INBUFF];

    while(getrealic(&theta, buff, INBUFF, fp) != EOF) {
        if(got_infinity)
            error("pophist.ini continues after an infinite epoch");
        if(getrealic(&mn, buff, INBUFF, fp) == EOF)
            error("Unexpected EOF reading mn in pophist.ini");
        if(getwordic(buff, INBUFF, fp) == NULL)
            error("Unexpected EOF reading tau in pophist.ini");
        /* check for word "inf" */
        for(i = 0; i < INBUFF && buff[i] != '\0'; i++)
            buff[i] = tolower(buff[i]);
        if(strncmp(buff, "inf", (sizeof buff)-1) == 0) {
            /* doesn't matter what number we put here */
            tau = INFTY;        /* because initial value of tau is ignored */
            got_infinity = 1;
        } else if(!isdigit(*buff)) {
            fprintf(stderr, "\nIllegal value of tau: \"%s\"\n", buff);
            exit(1);
        } else
            sscanf(buff, "%f", &tau);
        if(getintic(&K, buff, INBUFF, fp) == EOF)
            error("Unexpected EOF reading K in pophist.ini");

        if(K < 1) {
            fprintf(stderr, "\nerror in population history data: K=%d\n", K);
            exit(1);
        }
        /* On input theta is the aggregate size of whole population.  But
         * the value stored is the size of a single subdivision. */
        ph = newhistory(theta / K, mn, tau, K);
        if(last == NULL)
            root = ph;          /* set root on first pass through loop */
        else
            last->next = ph;    /* append to list on later passes */
        last = ph;              /* always set last to end of list */
    }
    if(!got_infinity) {
        fprintf(stderr,
                "\nWarning: tau was not infinite in initial epoch of pophist.int.");
        fprintf(stderr, "\n         I'm treating it as infinite anyway.");
    }
    return (root);
}

/*** create a new structure of type POPHIST ****/
POPHIST    *newhistory(real theta, real mn, real tau, int K) {
    POPHIST    *ph;

    ph = (POPHIST *) mustalloc(sizeof(POPHIST));
    ph->next = NULL;
    ph->mn = mn;
    ph->theta = theta;
    ph->tau = tau;
    ph->K = K;
    if(K == 1)
        ph->mn = 0.0;           /* can't have migration w/ only 1 group */
    return (ph);
}

/*** free a list of type POPHIST ***/
void freehistory(POPHIST * ph) {
    if(ph == NULL)
        return;
    freehistory(ph->next);
    free(ph);
}

/*** duplicate a list of type POPHIST ****/
POPHIST    *duphist(POPHIST * new, POPHIST * old) {
    if(old == NULL) {
        freehistory(new);
        return (NULL);
    }
    if(new == NULL) {
        new = newhistory(old->theta, old->mn, old->tau, old->K);
        new->next = NULL;
    } else {
        new->theta = old->theta;
        new->mn = old->mn;
        new->tau = old->tau;
        new->K = old->K;
    }
    new->next = duphist(new->next, old->next);
    return (new);
}
void writehistory(FILE * fp, POPHIST * history, char comment) {
    fprintf(fp, "\n%c%10s %10s %10s %10s",
            comment, "theta", "mn", "tau", "K");
    if(history == NULL)
        return;
    while(history->next != NULL) {
        fprintf(fp, "\n%c%10.4f %10.4f %10.4f %10d",
                comment,
                (history->theta) * (history->K),
                history->mn, history->tau, history->K);
        history = history->next;
    }
    fprintf(fp, "\n%c%10.4f %10.4f %10s %10d",
            comment,
            (history->theta) * (history->K), history->mn, "Inf", history->K);
}

/********integer comparison function used by qsort()************/
int icompar(int *x, int *y) {
    if(*x > *y)
        return (1);
    if(*x < *y)
        return (-1);
    return (0);
}

/****************************************************************
The collapse_vector, v, is used to reduce the number of groups from K1
to K2.  Each individual in old pop i is assigned to new pop v[i].  To
generate v, I calculate a random partition of the old groups.
****************************************************************/
int get_collapse_vector(int K1, int K2, int *v) {
    int         i, j, r[MAX], p[MAX];

    (void) randperm(r, K1);     /* randomly reorders groups */
    (void) randperm(p, K1 - 1); /* vector of random partition points */
    if(K2 > 2)
        qsort(p, K2 - 1, sizeof(int),   /* sort 1st K2-1 partition points */
              (int (*)(const void *, const void *)) icompar);
    for(i = j = 0; i < K2 - 1; i++) {   /* j indexes partitions */
        while(j <= p[i])        /* we are in the i'th partition */
            v[r[j++]] = i;      /* assign i to all members of this partition */
    }                           /* r[] randomizes partition membership */
    while(j < K1)               /* the last partition goes to the end */
        v[r[j++]] = i;
    return (0);
}

/***** measure depth of tree *****/
real treedepth(NODE * node) {
    if(node == NULL)
        return (0.0);
    return (node->branch + treedepth(node->left));
}

/**** measure length of tree = sum of all branch lengths ****/
real treelength(NODE * node) {
    if(node == NULL)
        return (0.0);
    return (node->branch + treelength(node->left) + treelength(node->right));
}

/* set external defs back to original state */
void clear_externals(void) {
    int         i;

    if(nodevec != NULL) {
        free(nodevec);
        nodevec = NULL;
    }
    if(smpl != NULL) {
        for(i = 0; i < nsubs; i++)
            free(smpl[i]);
        free(smpl);
    }
    nodevec = (NODE *) NULL;
    smpl = (NODE ***) NULL;
    sampsize = nsubs = nnodes = nextnode = 0;
}

LIST       *newlist(int n, int *ndx) {
    int         i;
    LIST       *l;

    l = (LIST *) stackmustalloc(sizeof(LIST));
    l->n = n;
    l->ndx = (int *) stackmustalloc(n * sizeof(int));
    for(i = 0; i < n; i++)
        l->ndx[i] = ndx[i];
    return (l);
}

/* print list */
void prlist(LIST * l) {
    int         i;

    for(i = 0; i < l->n; i++)
        printf(" %d", l->ndx[i]);
}

/* merge two lists */
LIST       *mergelist(LIST * l1, LIST * l2) {
    LIST       *new;
    int         i1, i2, inew;

    if(l2 == NULL)
        return (l1);
    if(l1 == NULL)
        return (l2);

    new = (LIST *) stackmustalloc(sizeof(LIST));
    new->n = l1->n + l2->n;
    new->ndx = (int *) stackmustalloc(new->n * sizeof(int));
    i1 = i2 = inew = 0;
    while((i1 < l1->n) && (i2 < l2->n)) {
        if(l1->ndx[i1] < l2->ndx[i2])
            new->ndx[inew++] = l1->ndx[i1++];
        else
            new->ndx[inew++] = l2->ndx[i2++];
    }
    while(i1 < l1->n) {
        new->ndx[inew++] = l1->ndx[i1++];
    }
    while(i2 < l2->n) {
        new->ndx[inew++] = l2->ndx[i2++];
    }
    assert(i1 == l1->n && i2 == l2->n && inew == new->n);

    return (new);
}
