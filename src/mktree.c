#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "mytypes.h"
#include "memstack.h"
#include "iscoales.h"
#include "textree.h"
#include "unirand.h"
#include "mutation.h"
#include "bye.h"
#define YES(x) ((x) ? "Yes" : "No")

/******* prototypes ************/
void        usage(void);
int         fold_spectrum(int *spec, int n);
int         getspec(NODE * node, int *spec, int len);

extern int *mismatch, ***intermatch;

int         domutation = 1;
int         sampsize;           // all subdivisions included
int         dohistory = 1;
int         dotree = 1;
int         domismatch = 1;
int         dospectrum = 1;
int         folded_spectrum = 0;

const char *msg = "Edit the file pophist.ini to make it reflect the history of the\n\
hypothetical population in which you are interested.  The pophist.ini file\n\
has a row for each epoch of the population's history.  There can\n\
be as many rows (corresponding to as many epochs) as you want.  Within each\n\
epoch, all parameters are constant.  Within each row, there are four\n\
columns:\n\
  column 1: the variable theta = 2*N*u, where N is haploid population\n\
            size and u is the mutation rate per generation.\n\
            In other words, theta measures population size in units\n\
            of 1/(2*u) individuals.\n\
  column 2: mn = m*n, where \"m\" is the immigration rate and \"n\" is\n\
            the haploid size of each subdivision. This parameter equals\n\
            the number of immigrants per generation into each population.\n\
  column 3: The length of the epoch in units of 1/(2*u) generations.\n\
  column 4: The number of sub-populations, all of equal size.\n\
The input routine ignores comments, which begin with '%' and end with\n\
end-of-line.\n";

void usage(void) {
    fflush(stdout);
    fprintf(stderr, "\nusage: mktree [options]");
    fprintf(stderr, "\n  where options may include");
    fprintf(stderr, "\n  -n<x>  set aggregate sample size. Def=%d",
            sampsize);
    fprintf(stderr, "\n  -u     put mutations on tree? Def=%s", YES(domutation));
    fprintf(stderr, "\n  -f     folded spectrum. Def=%s", YES(folded_spectrum));
    fprintf(stderr, "\n  -m     print mismatch? Def=%s", YES(domismatch));
    putc('\n', stderr);

    fputs(msg,stderr);
    exit(1);
}

int main(int argc, char **argv) {
    int         i, fold;
    int         msize = 500;    /* max size of a mismatch distribution */
    int         mmlength = 0;   /* realized size of mismatch */
    int        *spec;
    NODE       *tree;
    POPHIST    *history;
    real        hinches, vinches, spinches;
    real        maxtree, maxtau, maxmm, maxx;
    int         imaxx, imaxtic, ibyx;
    FILE       *fp;
    int         nsubs;          /* number of subdivisions in sample */
    int        *subsize;        /* subsize[i] = size of i'th subdiv */
    real        mpd;            /* mean pairwise difference */
    char        message[200];

    hinches = 4.0;              /* horizontal size of plots in inches */
    vinches = 1.4;              /* vertical size */
    spinches = 0.5;             /* spacing btw plots */
    sampsize = 50;

/****** Command line arguments *********/
    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            switch (argv[i][1]) {
            case 'f':
                folded_spectrum = 1;
                break;
            case 'n':
                sampsize = atoi(argv[i] + 2);
                break;
            case 'u':
                domutation = 0;
                break;
            case 'm':
                domismatch = !domismatch;
                break;
            default:
                usage();
                break;
            }
        } else
            usage();
    }

  /**** get history *****/
    fp = fopen("pophist.ini", "r");

    if(fp == NULL) {
        fprintf(stderr, "\nWarning: Can't open pophist.ini\n");
        exit(1);
    }
    history = gethistory(fp);

    (void) setmemstack(4000 + 108.5447 * sampsize
                       + 2.03351 * sampsize * sampsize);

    nsubs = history->K;
    subsize = (int *) mustalloc(nsubs * sizeof(int));

    /* set subdivision sizes */
    for(i = 0; i < nsubs; i++)
        subsize[i] = sampsize / nsubs;

    /* put the remainder into subdivision 0 */
    subsize[0] += sampsize % nsubs;

    /* initialize random number generator */
    initrand(0);

    /* echo parameters */
    fputs("% -*-latex-*-", stdout);
    bold_comment("Tree in PicTeX Format", stdout);
    printf("\n%%Sample size   = %d", sampsize);
    printf("\n%%Do mutation = %d", domutation);
    writehistory(stdout, history, '%');

    if(sampsize > MAX) {
        fflush(stdout);
        fprintf(stderr, "\nERROR: The sample size (%d) is too big.  Max=%d\n",
                sampsize, MAX);
        exit(1);
    }
    tree = iscoales(nsubs, subsize, history, domutation);

    /* calculate histogram */
    crosstab(0, tree);
    mmlength = getmatch(msize, tree);

    /* x axis upper limits: */
    maxmm = mmlength - 1;       /* required by mismatch */
    maxtau = 1.1 * total_tau(history);  /* required by history */
    maxtree = 1.1 * treedepth(tree);    /* required by tree */

    /* maxx is max(maxmm, maxtau, maxtree) */
    maxx = maxmm;
    if(maxtau > maxx)
        maxx = maxtau;
    if(maxtree > maxx)
        maxx = maxtree;
    printf("\n%%maxmm=%f maxtau=%f maxtree=%f maxx=%f",
           maxmm, maxtau, maxtree, maxx);

    int_axis(maxx, &imaxx, &imaxtic, &ibyx);

    setup_pictex(hinches, vinches, spinches, stdout);

    if(dohistory) {
        print_history(history, imaxx, imaxtic, ibyx, stdout);
    }

    if(dotree) {
        print_textree(tree, (real) imaxx, hinches, vinches, stdout);
    }

    if(domismatch) {
        mpd =
            print_mismatch(mismatch, mmlength, imaxx, imaxtic, ibyx, stdout);
        sprintf(message, "mean pairwise difference: %f", mpd);
        fflush(stdout);
    }
    if(dospectrum) {
        /* allocate spectrum */
        spec = (int *) mustalloc(sampsize * sizeof(int));
        memset((void *) spec, 0, sampsize * sizeof(int));
        i = getspec(tree, spec, sampsize);
        assert(i == sampsize);
        if(folded_spectrum) {
            fold = 1 + fold_spectrum(spec, sampsize);
        } else {
            fold = sampsize;
        }
        print_spectrum(spec, fold, sampsize, folded_spectrum, stdout);
    }
    if(domismatch)
        wrapup_pictex(message, stdout);
    else
        wrapup_pictex(NULL, stdout);
    putchar('\n');
    exit(0);
}

/* get unfolded site frequency spectrum */
int getspec(NODE * node, int *spec, int len) {
    int         descendants;

    if(node == NULL) {
        return (0);
    }
    /* how many descendants of this node? */
    descendants = getspec(node->left, spec, len)
        + getspec(node->right, spec, len);
    /* is this node a leaf? */
    if(descendants == 0)
        descendants = 1;

    assert(node->mutations >= 0);

    /* add mutated sites to spectrum */
    if(node->mutations > 0) {
        if(descendants >= len) {
            printf("\n%% Warning: truncating descendants from %d to %d",
                   descendants, len - 1);
            spec[len - 1] += node->mutations;
        } else
            spec[descendants] += node->mutations;
    }
    /* returned value is number of descendants of current node */
    return descendants;
}

/*
 * fold spectrum about middle
 * spectrum should have length n
 * 0'th entry of spectrum is unused
 * returns length of folded spectrum (including unused 0'th position)
 */
int fold_spectrum(int *spec, int n) {
    int         i;

    for(i = 1; i <= n / 2; i++)
        spec[i] += spec[n - i];

    /* find index of maximum entry in folded spectrum */
    if(n % 2 == 0)
        i = n / 2;
    else
        i = n / 2 + 1;

    /* rtn length including unused 0'th position */
    return i;
}
