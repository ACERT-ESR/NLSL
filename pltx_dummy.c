/* dummy function to stand in for pltx functions and allow compilation on windows
                                          */

/* Structure definitions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fortrancall.h"

#define MAX_WINDOWS 5

typedef char word[20];		/* Defined sting type word */

/* Global variable declarations */
int numwndws = 0;		/* Initialize number of windows to 0 */
int done = 0;			/* Initialize "done" to "not done" */
int displayOK = 0;		/* Variable to keep wpoll from attempting to
				   check events for nonexistant windows */

/********************************/
void FORTRAN (getwndws) (int *n, word * title)
{				/* Begin getwndws */
  /* declarations */
  //printf("Received a call to getwndws(%d,%s)",n,title);
}				/* End getwndws */

/********************************/
void FORTRAN (fstplt) (double *y1, double *y2, double *xmin1,
		       double *xstep1, int *indx1, int *wnum)
/* double *y1, *y2;               y1 is the experimental, y2 is the
   double *xmin1, *xstep1;        difference (experimental-simulated)
   int *indx1, *wnum; */
{				/* Begin fstplt */
  /* declarations */
    //printf("received a call to fstplt(%g,%g,%g,%g,%d,%d)",y1,y2,xmin1,xstep1,indx1,wnum);
}				/* End fstplt */

/********************************/
void FORTRAN (wpoll) ()
{				/* Begin wpoll */
    //printf("received a call to wpoll");
}				/* End wpoll */

/********************************/
void FORTRAN (shtwndws) ()
{				/* Begin shtwndws */
    //printf("Received a call to shtwndws");
}				/* End shtwndws */

/********************************/
