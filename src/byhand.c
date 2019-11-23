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
void         usage(void);
int          sampsize;		/* all subdivisions included */
int          dotree = 1;

void         usage(void)
{
  fflush(stdout);
  fprintf(stderr, "\nusage: byhand");
  putc('\n', stderr);
  exit(1);
}

int            main(int argc, char **argv)
{
  extern int   nnodes;
  extern NODE  *nodevec;
  int          ntips;
  real         root10;
  real         hinches, vinches, spinches;
  real         maxx;
  int          imaxx, imaxtic, ibyx;
  NODE        *tree, *n1, *n2, *n3, *n4, *n5, *clade1, *clade2;


  root10 = sqrt(10.0);
  hinches = 4.0;		/* horizontal size of plots in inches */
  vinches = 1.4;		/* vertical size */
  spinches = 0.5;		/* spacing btw plots */

  /* echo parameters */
  fputs("% -*-latex-*-", stdout);
  bold_comment("Tree in PicTeX Format", stdout);

  /* allocate nodevec */
  ntips = 9;
  nnodes = 2*ntips - 1;
  nodevec = (NODE *) mustalloc(nnodes * sizeof(NODE));

  /* construct tree */
  n1 = newnode(NULL,NULL); /* Hawaii */
  n2 = newnode(NULL,NULL); /* Maui */
  n1->branch = n2->branch = 1.0;
  n3 = newnode(NULL,NULL); /* Oahu */
  n3->branch = 2.0;
  n4 = newnode(n1,n2); /* (Hawaii,Maui) */
  n4->branch = 1.0;
  clade1 = newnode(n4, n3); /*((Hawaii,Maui),Oahu) */
  clade1->branch = 3.0;

  n1 = newnode(NULL,NULL); /* Hawaii */
  n2 = newnode(NULL,NULL); /* Maui */
  n1->branch = n2->branch = 1.0;
  n3 = newnode(NULL,NULL); /* Maui again */
  n3->branch = 2.0;
  n4 = newnode(NULL,NULL); /* Molokai */
  n4->branch = 3.0;
  n5 = newnode(NULL,NULL); /* Oahu */
  n5->branch = 4.0;

  clade2 = newnode(n1,n2); /* (Hawaii,Maui)*/
  clade2->branch = 1.0;
  clade2 = newnode(clade2, n3); /*((Ha,Ma),Ma)*/
  clade2->branch = 1.0;
  clade2 = newnode(clade2,n4); /*(((Ha,Ma),Ma),Mo)*/
  clade2->branch=1.0;
  clade2 = newnode(clade2,n5); /*((((Ha,Ma),Ma),Mo),Oa)*/
  clade2->branch=1.0;
  
  clade1 = newnode(clade1, clade2);
  clade1->branch = 1.0;
  n1 = newnode(NULL,NULL); /* Kaui */
  n1->branch = 6.0;
  tree = newnode(clade1, n1);

  maxx = 1.1*treedepth(tree);   /* required by tree */
  printf("\n%%maxx=%f", maxx);
  int_axis(maxx, &imaxx, &imaxtic, &ibyx);
  setup_pictex(hinches, vinches, spinches, stdout);
  print_textree(tree, (real) imaxx, hinches, vinches, stdout);
  wrapup_pictex(NULL, stdout);
  putchar('\n');
  exit(0);
}

