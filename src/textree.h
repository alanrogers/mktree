void         bold_comment(char *str, FILE * ofp);
void         print_textree(NODE * tree, real maxx, real hinches,
			   real vinches, FILE * ofp);
void         setup_pictex(real hinches, real vinches, real spinches,
			  FILE * ofp);
void         wrapup_pictex(char *message, FILE * ofp);
void         do_coordinates(real lo_x, real hi_x, real lo_y, real hi_y,
			    FILE * ofp);
real         textree(NODE * node, real * vpos, real * mid, FILE * ofp);
void         putrule(real x1, real y1, real x2, real y2, FILE * ofp);
void         putmutation(real x, real y, FILE * ofp);
void         texscatterplot(int *h, int max, FILE * ofp);
void         texlineplot(int *h, int max, FILE * ofp);
void         texhist(int *h, int lo, int hi, FILE * ofp);
int          countleaves(NODE * node);
real         print_mismatch(int *m, int length, int imaxx, int imaxtic,
			    int ibyx, FILE *ofp);
void         print_spectrum(int *s, int length, int sampsize, int folded_spectrum, FILE *ofp);
real total_tau(POPHIST *ph);
real max_theta(POPHIST *ph);
void print_ph_node(real from, POPHIST *ph, real maxtau, FILE *ofp);
void print_history(POPHIST *history, int imaxx, int imaxtic, int ibyx,
		   FILE *ofp);
void int_axis(real max, int *imax, int *imaxtic, int *iby);
int find_pretty_max(int need, int got);
