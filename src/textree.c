#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mytypes.h"
#include "iscoales.h"
#include "unirand.h"
#include "textree.h"

/*
 * Print initial part of PicTeX output, specifying:
 * hinches: horizontal dimension of each plot in inches
 * vinches: vertical dimension of each plot in inches
 * spinches: spacing between plots.  ofp is output file
 */
void setup_pictex(real hinches, real vinches, real spinches, FILE * ofp) {
    fputs("\n\\newdimen\\offsety", ofp);
    fputs("\n\\newdimen\\yunit", ofp);
    fputs("\n\\newdimen\\xunit", ofp);
    fputs("\n\\newdimen\\thusfar", ofp);
    fputs("\n\\newdimen\\plotht", ofp);
    fputs("\n\\newdimen\\plotwd", ofp);
    fputs("\n\\newdimen\\plotsp", ofp);
    fprintf(ofp, "\n\\thusfar=%fin  %s",
            0.0, "% Keeps track of what's above");
    fprintf(ofp, "\n\\plotht=%fin   %s", vinches, "% Height of each plot");
    fprintf(ofp, "\n\\plotwd=%fin   %s", hinches, "% Width of each plot");
    fprintf(ofp, "\n\\plotsp=%fin   %s", spinches, "% Spacing between plots");
    fputs("\n\\begin{figure}", ofp);
    fputs("\n\\begin{center}", ofp);
    fputs("\n\\mbox{\\beginpicture", ofp);
    fputs("\n\\def\\mutation{\\tiny$\\bullet$}", ofp);
    fputs("\n\\small", ofp);
    fputs("\n\\valuestolabelleading=.4\\baselineskip", ofp);
    fputs("\n\\headingtoplotskip=0.4\\baselineskip", ofp);
}

void wrapup_pictex(char *message, FILE * ofp) {
    fputs("\n\\endpicture}", ofp);
    fputs("\n\\end{center}", ofp);
    if(message != NULL)
        fprintf(ofp, "\n%s", message);
    fputs("\n\\end{figure}", ofp);
}

void do_coordinates(real lo_x, real hi_x, real lo_y, real hi_y, FILE * ofp) {
    real        rx, ry;

    rx = hi_x - lo_x;
    ry = hi_y - lo_y;

    fputs("\n\\yunit=\\plotht", ofp);
    fputs("\n\\xunit=\\plotwd", ofp);
    fprintf(ofp, "\n\\Divide <\\xunit> by <%fpt> forming <\\xunit>", rx);
    fprintf(ofp, "\n\\Divide <\\yunit> by <%fpt> forming <\\yunit>", ry);
    fputs("\n\\advance\\thusfar by \\plotht", ofp);
    fputs("\n\\advance\\thusfar by \\plotsp", ofp);
    fputs("\n\\Divide <\\thusfar> by <\\yunit> forming <\\offsety>", ofp);
    fputs("\n\\placevalueinpts of <\\offsety> in {\\mktree_tmp}", ofp);
    fprintf(ofp, "\n\\setcoordinatesystem units <\\xunit, \\yunit>");
    fputs(" point at 0 {\\mktree_tmp}", ofp);
    fprintf(ofp, "\n\\setplotarea x from %f to %f, y from %f to %f",
            lo_x, hi_x, lo_y, hi_y);
}

void print_textree(NODE * tree, real maxx, real hinches,
                   real vinches, FILE * ofp) {
    real        vpos, mid, hpos;

    hpos = treedepth(tree);
    vpos = countleaves(tree);

    bold_comment("Gene Genealogy", ofp);
    do_coordinates(0.0, maxx, 0.0, vpos - 1.0, ofp);
    fputs("\n\\axis left invisible label {\\lines{Gene\\cr genealogy}} /",
          ofp);
    fputs("\n\\axis bottom invisible", ofp);
    fputs("\n     label {Mutational time before present} /", ofp);
    vpos = mid = 0.0;
    hpos = textree(tree, &vpos, &mid, ofp);
    putrule(hpos, mid, maxx, mid, ofp);
}

void bold_comment(char *str, FILE * ofp) {
    const char *b = "%%%%%%%%%%%%%%";

    fprintf(ofp, "\n%%\n%s %s %s\n%%", b, str, b);
}

/*
 * Recursive algorithm for drawing a tree.
 *
 * On entry:
 *  node points to current node of tree
 *  *vpos should equal 0.0.
 *  *mid's value doesn't matter.
 *  ofp is the output file
 *
 * On return:
 *  *node is unchanged
 *  *vpos  equals the number of leaves in the tree.
 *  *mid is the vertical coordinate of the middle of the segment separating
 *     the two top-level clades.
 * The function returns the horizontal coordinate of the tree's root.
 */

real textree(NODE * node, real * vpos, real * mid, FILE * ofp) {
    real        top, bottom, hpos;
    int         i;

    if(node == NULL) {
        *mid = *vpos;
        return (0);
    }
    hpos = textree(node->left, vpos, mid, ofp);
    top = *mid;
    hpos = textree(node->right, vpos, mid, ofp);
    bottom = *mid;
    if(hpos > 0 && top != bottom) {
        putrule(hpos, bottom, hpos, top, ofp);
    }
    if(top == bottom)
        *mid = top;
    else
        *mid = bottom + 0.5 * (top - bottom);
    if(node->branch > 0)
        putrule(hpos, *mid, hpos + node->branch, *mid, ofp);
    for(i = 0; i < node->mutations; i++)
        putmutation(hpos + uni() * (node->branch), *mid, ofp);
    if(hpos == 0.0)
        *vpos += 1.0;
    return (hpos + node->branch);
}

void putrule(real x1, real y1, real x2, real y2, FILE * ofp) {
    fprintf(ofp, "\n\\putrule from %f %f to %f %f", x1, y1, x2, y2);
}

void putmutation(real x, real y, FILE * ofp) {
    fprintf(ofp, "\n\\put {\\mutation} at %f %f", x, y);
}
void texscatterplot(int *h, int max, FILE * ofp) {
    int         i;
    real        sum = 0.0;

    for(i = 0; i <= max; i++)
        sum += h[i];
    fprintf(ofp, "\n\\multiput {$\\circ$} at");
    if(sum > 0.0)
        for(i = 0; i <= max; i++) {
            if(i % 5 == 0)
                putc('\n', ofp);
            fprintf(ofp, " %d %f", i, h[i] / sum);
        }
    fprintf(ofp, "\n/");
}
void texlineplot(int *h, int max, FILE * ofp) {
    int         i;
    real        sum = 0.0;

    for(i = 0; i <= max; i++)
        sum += h[i];
    fprintf(ofp, "\n\\plot");
    if(sum > 0.0) {
        for(i = 0; i <= max; i++) {
            if(i % 5 == 0)
                putc('\n', ofp);
            fprintf(ofp, " %d %f", i, h[i] / sum);
        }
        fprintf(ofp, "\n/");
    } else {
        fflush(stdout);
        fprintf(stderr, "\nWarning: texlineplot was given an empty vector");
    }
}

void texhist(int *h, int lo, int hi, FILE * ofp) {
    int         i;
    real        sum = 0.0;

    for(i = lo; i <= hi; i++)
        sum += h[i];
    fprintf(ofp, "\n\\sethistograms");
    fprintf(ofp, "\n\\plot 0 0");
    if(sum > 0.0)
        for(i = lo; i <= hi; i++) {
            if(i % 5 == 0)
                putc('\n', ofp);
            fprintf(ofp, " %d %f", i, h[i] / sum);
        }
    fprintf(ofp, "\n/");
    fprintf(ofp, "\n\\setlinear");
}

int countleaves(NODE * node) {
    if(node == NULL)
        return (0);

    if(node->left == NULL && node->right == NULL)
        return (1);

    return (countleaves(node->left)
            + countleaves(node->right));
}

real
print_mismatch(int *m, int length, int imaxx, int imaxtic, int ibyx,
               FILE * ofp) {
    int         i, isum;
    real        hi_y, mpd = 0;

    bold_comment("Mismatch Distribution", stdout);
    isum = 0;
    for(i = 0; i < length; i++) {
        isum += m[i];
        mpd += i * m[i];
        if(m[i] > hi_y)
            hi_y = m[i];
    }
    hi_y /= isum;
    mpd /= isum;
    do_coordinates(0.0, (real) imaxx, 0.0, hi_y, ofp);
    fputs("\n\\axis left label {\\lines{Mismatch\\cr distribution}} /", ofp);
    fprintf(ofp, "\n\\axis bottom");
    fprintf(ofp, "\n   label {Pairwise differences}");
    fprintf(ofp, "\n   ticks numbered from 0 to %d by %d /", imaxtic, ibyx);
    texscatterplot(m, length - 1, ofp);
    texlineplot(m, length - 1, ofp);
    return (mpd);
}

void
print_spectrum(int *s, int length, int sampsize, int folded_spectrum,
               FILE * ofp) {
    int         i, sumspec, lim;
    real        x, a, maxy, maxe, maxo;
    real        stretch = 1.1;

    /* make room for bold line */
    fprintf(ofp, "\n\\plotht=%f\\plotht", stretch);

    /* find y dimension */
    maxe = a = 0.0;
    for(i = 1; i < sampsize; i++) {
        x = 1.0 / i;
        if(x > maxe)
            maxe = x;
        a += x;
    }
    maxe /= a;                  /* maximal expected value */

    maxo = sumspec = 0.0;
    for(i = 1; i < length; i++) {
        if(s[i] > maxo)
            maxo = s[i];
        sumspec += s[i];
    }
    maxo /= sumspec;            /* maximal observed value */

    maxy = (maxo > maxe ? maxo : maxe); /* maximal y value */
    maxy *= stretch;

    bold_comment("Simulated site frequency spectrum", ofp);
    do_coordinates(0.0, length - 1.0, 0.0, maxy, ofp);
    fputs("\n\\axis left invisible", ofp);
    fputs(" label {\\lines{Site\\cr frequency\\cr spectrum}} /", ofp);
    fprintf(ofp, "\n\\axis bottom");
    if(folded_spectrum) {
        fprintf(ofp, "\n   label {Frequency of minor allele}");
        fprintf(ofp, "\n   ticks withvalues 0 {1/2} / at 0 %d / /",
                length - 1);
    } else {
        fprintf(ofp, "\n   label {Frequency of mutant allele}");
        fprintf(ofp, "\n   ticks withvalues 0 1 / at 0 %d / /", length - 1);
    }
    fputs("\n\\linethickness 1.2pt", ofp);
    fputs("\n\\axis top /", ofp);
    fputs("\n\\linethickness 0.4pt", ofp);
    texhist(s, 1, length - 1, stdout);

    bold_comment("Neutral expectation of site freq spectrum", ofp);
    fputs("\n\\multiput {$\\bullet$} at ", ofp);
    lim = (folded_spectrum ? sampsize / 2 : sampsize - 1);
    for(i = 1; i <= lim; i++) {
        if((i - 1) % 5 == 0)
            fputs("\n   ", stdout);
        if(folded_spectrum)
            printf("  %f %g", i - 0.5,
                   1.0 / (i * a) + 1.0 / ((sampsize - i) * a));
        else
            printf("  %f %g", i - 0.5, 1.0 / (i * a));
    }
    if(folded_spectrum && sampsize % 2 != 0) {
        i = sampsize / 2 + 1;
        printf("\n   %f %f", i - 0.5, 1.0 / (i * a));
    }
    fputs("\n/", stdout);
    fflush(ofp);
}
void
print_history(POPHIST * history, int imaxx, int imaxtic, int ibyx,
              FILE * ofp) {
    real        tau, theta;
    const real  inflate = 1.2;  /* y dim inflation factor */
    const real  downshift = 0.03;   /* x axis downshift factor */

    tau = total_tau(history);
    assert(tau <= imaxx);
    theta = max_theta(history);
    bold_comment("Population size", ofp);
    do_coordinates(0.0, (real) imaxx, 0.0, inflate * theta, ofp);
    fprintf(ofp, "\n\\advance\\thusfar by %f\\plotht", downshift);
    fputs("\n\\axis left label {\\lines{Population\\cr size}} /", ofp);
    fprintf(ofp, "\n\\axis bottom shiftedto y=-%f",
            downshift * inflate * theta);
    fprintf(ofp, "\n    label {Mutational time before present}");
    fprintf(ofp, "\n    ticks numbered from %d to %d by %d /",
            0, imaxtic, ibyx);
    print_ph_node(0.0, history, (real) imaxx, ofp);
}

void print_ph_node(real from, POPHIST * ph, real maxtau, FILE * ofp) {
    real        tau;
    real        totalTheta, nextTotalTheta; /* aggregate across subpopulations */

    if(ph == NULL)
        return;

    if(ph->tau == INFTY)
        tau = maxtau;
    else
        tau = ph->tau;

    totalTheta = ph->theta * ph->K;

    fprintf(ofp, "\n\\putrule from %f %f to %f %f",
            from, totalTheta, from + tau, totalTheta);
    if(ph->next != NULL) {
        nextTotalTheta = ph->next->theta * ph->next->K;
        if(totalTheta != nextTotalTheta)
            fprintf(ofp, "\n\\putrule from %f %f to %f %f",
                    from + tau, totalTheta, from + tau, nextTotalTheta);
        print_ph_node(from + ph->tau, ph->next, maxtau - ph->tau, ofp);
    }
}

/*
 * Find the total amount of mutational time in a history list
 * excluding the final infinite epoch.
 */
real total_tau(POPHIST * ph) {
    if(ph == NULL || ph->tau == INFTY)
        return (0.0);

    return (ph->tau + total_tau(ph->next));
}

/*
 * Find the maximal theta in a history list
 */
real max_theta(POPHIST * ph) {
    real        max;
    real        totalTheta;

    if(ph == NULL)
        return (0.0);

    max = max_theta(ph->next);  /* max in rest of list */

    totalTheta = ph->theta * ph->K; /* aggregate theta across subpopulations */

    if(max > totalTheta)
        return (max);

    return (totalTheta);
}

/* set axis limit and increment */
enum { NIV = 3 };
int         iv[] = { 2, 3, 5 };
void int_axis(real max, int *imax, int *imaxtic, int *iby) {
    int         i;

    *imax = (int) ceil(max);

    /* find maximal tick value */
    *imaxtic = find_pretty_max(*imax, 1);

    for(i = 0; i < NIV; i++) {
        if(*imaxtic % iv[i] == 0)
            break;
    }
    if(i == NIV) {
        fprintf(stderr, "\nerror in int_axis");
        exit(1);
    }
    *iby = *imaxtic / iv[i];
}

/* find largest pretty number <= need */
int find_pretty_max(int need, int got) {
    int         i, try, max;

    max = got;
    for(i = 0; i < NIV; i++) {
        if(got * iv[i] == need)
            return (got * iv[i]);
        if(got * iv[i] < need) {
            try = find_pretty_max(need, got * iv[i]);
            if(try > max)
                max = try;
        }
    }

    return (max);
}
