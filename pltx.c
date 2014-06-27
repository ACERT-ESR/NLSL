/* Pltx.c is an XWindows application capable of handling several 
   independant windows and their contents simultaneously.  Each
   window contains a plot of two spectra (experimental and simulated).
   A line representing the difference of these two spectra can also
   be plotted, plus an axes box (without labels) drawn enclosing
   the plot.  

   While running, the function will perform certain operations on the 
   plotted data in each window as follows:

   KEYSTROKE OPERATIONS:

   "a": Turns the axes box on or off, switch-fashion.

   "o": Displays the graphs of the two arrays sent to the function,
        the original plot.  The simulation curve is shown in green,
        the experimental in orange.

   "d": Displays the curve that represents the difference between
        the simulation and the difference.  Shown in magenta.

   "s": Displays the simulation curve only.

   "e": Displays the experimental curve only.

   "b": Displays all three curves.

   "q": Exits the function and returns to the calling program.

   MOUSE-BUTTON PRESSES: 

   Right button:  When pressed, the rightheand mouse button will highlight
        a pixel point and display the data coordinates (x, y).

   Left button:  When pressed, moved, and released, the lefthand mouse 
        button will highlight the pixel points at the spots where the 
        button was pressed and released, and will display the difference 
        between the two points (x, y).

   The function will also automatically rescale the data if the
   window size is changed.

   Written by Diane Crepeau             
   MODIFIED by David Budil 2/8/93                                       */

/* Structure definitions */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fortrancall.h"

#define MAX_WINDOWS 5

struct wstruct			/* Structure to contain window information */
{
  Display *dspy;		/* Connection between workstation and window */
  Drawable wndw;		/* Window ID */
  GC gc, gc1, gc2, gc3;		/* Graphics context IDs */
  Colormap map;			/* Window colormap */
  XPoint xye[4096], xys[4096], ydif[4096];	/* Pixel arrays */
  float ye[4096], ys[4096];	/* Data arrays associated with the window */
  float xmin, xmax, xstep, ymin, ymax;	/* Minima, maxima, and X step size */
  int indx;			/* Number of data points in arrays */
  int scrn;			/* Screen number */
  int mode;			/* Plotting mode, origin is corner of screen */
  unsigned int fgnd, bgnd;	/* Default foregraound and background */
  int flag, flag2;		/* Flags determining which spectra to plot */
  int yzer;			/* Pixel y-coordinate of zero data value */
  int nx, ny;			/* Dimensions of window in pixels */
  float er, el, et, eb;		/* Fractions of screen for edges (blank) */
  int led, red, bed, ted;	/* Pixel coordinates of window inner edges */
};

struct xarry			/* Structure for linked list of x-values */
{
  float x;			/* Data value */
  struct xarry *next;		/* Pointer to next 'node' in list */
};

typedef char word[20];		/* Defined sting type word */

/* Global variable declarations */
struct wstruct wds[MAX_WINDOWS + 1];	/* Maximum 5 windows */
int numwndws = 0;		/* Initialize number of windows to 0 */
int done = 0;			/* Initialize "done" to "not done" */
int displayOK = 0;		/* Variable to keep wpoll from attempting to
				   check events for nonexistant windows */

/********************************/
/* Function findzero */
/* Calculates pixel value for '0' data value */
int
findzero (wn, ymax, ymin)
     float ymax, ymin;		/* Max and min for calculating value */
     int wn;			/* Window number */
{				/* Begin findzero */
  int yzer1;			/* yzer contains the pixel y-coordinate of the
				   data value closest to zero, i is the loop
				   control index */

  yzer1 = (wds[wn].ny * (1 - wds[wn].eb) - wds[wn].ny * wds[wn].et) /
    (ymin - ymax) * (-ymax) + (wds[wn].ny * wds[wn].et);
  /* Calculate zero pixel value */

  return yzer1;			/* Return pixel coordinate of zero */

}				/* End findzero */

/********************************/
/* Function drawln */
/* Accesses external variables to draw spectra in window */
void
drawln (int wn)			/* wn=Window number */
{				/* Begin drawln */
  /* Access variables needed */
  XClearWindow (wds[wn].dspy, wds[wn].wndw);	/*  Clear window */

  if (wds[wn].flag2 == 1)	/* Check flag2 for axes on or off */
    {				/* Begin if */
      XDrawLine (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, wds[wn].led,
		 wds[wn].ted, wds[wn].red, wds[wn].ted);
      XDrawLine (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, wds[wn].led,
		 wds[wn].bed, wds[wn].red, wds[wn].bed);
      XDrawLine (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, wds[wn].led,
		 wds[wn].ted, wds[wn].led, wds[wn].bed);
      XDrawLine (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, wds[wn].red,
		 wds[wn].ted, wds[wn].red, wds[wn].bed);
      XDrawLine (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, wds[wn].led,
		 wds[wn].yzer, wds[wn].red, wds[wn].yzer);
    }				/* End if */

  switch (wds[wn].flag)		/* Check which spectrum to draw */
    {				/* Begin switch */
    case 0:			/* Draw experimental and simulated spectra */
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc1, wds[wn].xye,
		  wds[wn].indx, wds[wn].mode);
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc2, wds[wn].xys,
		  wds[wn].indx, wds[wn].mode);
      break;
    case 1:			/* Draw difference line */
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc3, wds[wn].ydif,
		  wds[wn].indx, wds[wn].mode);
      break;
    case 2:			/* Draw all three spectra */
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc1, wds[wn].xye,
		  wds[wn].indx, wds[wn].mode);
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc2, wds[wn].xys,
		  wds[wn].indx, wds[wn].mode);
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc3, wds[wn].ydif,
		  wds[wn].indx, wds[wn].mode);
      break;
    case 3:			/* Draw experimental spectrum */
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc1, wds[wn].xye,
		  wds[wn].indx, wds[wn].mode);
      break;
    case 4:			/* Draw simulated spectrum */
      XDrawLines (wds[wn].dspy, wds[wn].wndw, wds[wn].gc2, wds[wn].xys,
		  wds[wn].indx, wds[wn].mode);
      break;
    }				/* End switch */
  XFlush (wds[wn].dspy);	/* Clear event queue */

}				/* End drawln */

/********************************/
/* Function xnew */
/* Dynamically allocates space for the next node of the linked 
   list, initializes the values to zero, and returns a pointer
   to the new node.                                            */
struct xarry *
xnew (struct xarry *list)
{				/* Begin xnew */
  list = (struct xarry *) malloc (sizeof (struct xarry));
  /* Allocate new node */
  list->x = 0.0f;		/* Set internal values to zero */
  list->next = NULL;		/* Set next pointer to end of list */
  return list;			/* Return pointer to calling function */

}				/* End xnew */

/********************************/
void
xdump (struct xarry *list)
{				/* Begin xdump */
  struct xarry *ptr;

  while (list != NULL)
    {
      ptr = list;
      list = list->next;
      free (ptr);
    }
}				/* End xdump */

/********************************/
/* Function findmin */
/* Searches array for smallest value */
float
findmin (float *ary, int indx)
{				/* Begin findmin */
  float min;			/* Minimum */
  int i;			/* Loop index */

  min = ary[0];			/* Initialize min */

  for (i = 1; i < indx; i++)	/* Loop control */
    {				/* Begin for */
      if (min > ary[i])		/* Compare min with next array value */
	min = ary[i];		/* Assign smallest value to min */
    }				/* End for */

  return min;			/* Return smallest value in array */
}				/* End findmin */

/********************************/
/* Function findmax */
/* Searches array for largest value */
float
findmax (float *ary, int indx)
{				/* Begin findmax */
  float max;			/* Maximum */
  int i;			/* Loop index */

  max = ary[0];			/* Initialize max */

  for (i = 1; i < indx; i++)	/* Loop control */
    {				/* Begin for */
      if (max < ary[i])		/* Compare max with next array value */
	max = ary[i];		/* Assign largest to max */
    }				/* End for */

  return max;			/* Return largest value in array */

}				/* End findmax */

/********************************/
/* Function differnce */

/*** Calculates difference spectrum pixel values
     given pixel values for the data and experimental curves ***/

void
differnce (XPoint * xy1, XPoint * xy2, int indx, XPoint * ydif, int yzer)
/* XPoint *xy1, *xy2, *ydif;      Pixel arrays to calculate with
   int indx, yzer;                Array index and zero coordinate */
{				/* Begin differnce */
  int i;			/* Loop index */

  for (i = 0; i < indx; ++i)	/* Loop control */
    {				/* Begin for */
      ydif[i].x = xy1[i].x;	/* Set x-coordinates the same */
      ydif[i].y = (xy1[i].y - xy2[i].y) + yzer;
      /* Set y-coordinates by taking difference of
         xy1 and xy2 and adding zero coordinate */
    }				/* End for */

}				/* End differnce */

/********************************/
/* Function pixary */
/* Translates data array into pixel array */
void
pixary (float *y, float min, float step, int n, XPoint * xy, float ymax,
	float ymin, int wn)
/* XPoint *xy;                 Pixel coordinate array to fill
  float *y, min, step, ymax, ymin;
                               Y-array,x-values,and y-maximum and minimum
  int n, wn;                   Index to arrays, window number */
{				/* Begin pixary */
  int i = 0, xer, xel, yet, yeb;	/* "For" loop index */
  float xmin, temp;		/* X minimum,and temp holding variable */
  struct xarry *xlist, *list;	/* 2 linked list pointers to a list for x */

  /* Assign values to xer, xel, yet, and yeb */
  xer = wds[wn].nx * (1 - wds[wn].er);	/* Right edge coordinate */
  xel = wds[wn].nx * wds[wn].el;	/* Left edge coordinate */
  yet = wds[wn].ny * wds[wn].et;	/* Top edge coordinate */
  yeb = wds[wn].ny * (1 - wds[wn].eb);	/* Bottom edge coordinate */

  /* translate min, step, and n into array of X values */
  xlist = xnew (xlist);		/* Initialize xlist with new node */
  list = xlist;			/* Set list to indicate same node */

  for (i = 0; i <= n; i++)	/* Loop control */
    {				/* Begin for */
      (*list).x = min + (i * step);	/* Calculate next x and assign to list */
      temp = (*list).x;		/* Set temp to the newest x value */
      (*list).next = xnew ((*list).next);	/* Get new node in list */
      list = (*list).next;	/* Increment list to point to new node */
    }				/* End for */

  wds[wn].xmax = temp;		/* Assign highest x value to xmax */
  xmin = (*xlist).x;		/* Assign 1st value in linked list to xmin */
  i = 0;			/* Re-initialize i */
  list = xlist;			/* Reset list to beginning of linked list */

  /* loop translates 2 data arrays into 1 */
  while (xlist != NULL)		/* Loop control */
    {				/* Begin while */
      xy[i].x =
	((xer - xel) / (wds[wn].xmax - xmin)) * ((*xlist).x - xmin) + xel;
      /* Calculate x-coordinate */
      xy[i].y = (yeb - yet) / (ymin - ymax) * (y[i] - ymax) + yet;
      /* Calculate y-coordinate */
      xlist = (*xlist).next;	/* Increment list pointer to next node */
      ++i;			/* Increment loop index */
    }				/* End while */
  xdump (xlist);

}				/* End pixary */

/********************************/
/* Function calcul */
/* Calls other functions to calculate the pixel values from the
   data arrays, adjust for vertical offset, and calculate the 
   difference between them                                     */
void
calcul (int wn)			/* wn=window number */
{				/* Begin calcul */
  float y1max, y1min, y2max, y2min;
  /* Y array maxima and minima */

  y1max = findmax (wds[wn].ye, wds[wn].indx);	/* Get maximum value of y1 array */
  y2max = findmax (wds[wn].ys, wds[wn].indx);	/* Get maximum value of y2 array */
  y1min = findmin (wds[wn].ye, wds[wn].indx);	/* Get minimum value of y1 array */
  y2min = findmin (wds[wn].ys, wds[wn].indx);	/* Get minimum value of y2 array */

  wds[wn].ymax = (y1max > y2max) ? y1max : y2max;
  wds[wn].ymin = (y1min < y2min) ? y1min : y2min;

  pixary (wds[wn].ye, wds[wn].xmin, wds[wn].xstep, wds[wn].indx,
	  wds[wn].xye, wds[wn].ymax, wds[wn].ymin, wn);
  /* Calculate pixel coordinates for y1 */
  pixary (wds[wn].ys, wds[wn].xmin, wds[wn].xstep, wds[wn].indx,
	  wds[wn].xys, wds[wn].ymax, wds[wn].ymin, wn);
  /* Calculate pixel coordinates for y1 */

  wds[wn].yzer = findzero (wn, wds[wn].ymax, wds[wn].ymin);
  /* Get zero coordinate for y1 */

  differnce (wds[wn].xye, wds[wn].xys, wds[wn].indx, wds[wn].ydif,
	     wds[wn].yzer);
  /* Calculate difference curve */
}				/* End calcul */

/********************************/
void
calcpt (int xpt, int ypt, float *xpt1, float *ypt1, int wn)
{				/* Begin calcpt */
  int xer, xel, yet, yeb;	/* Edge-of-plot coordinate variables */

  /* Assign values to xer, xel, yet, and yeb */
  xer = wds[wn].nx * (1 - wds[wn].er);	/* Right edge coordinate */
  xel = wds[wn].nx * wds[wn].el;	/* Left edge coordinate */
  yet = wds[wn].ny * wds[wn].et;	/* Top edge coordinate */
  yeb = wds[wn].ny * (1 - wds[wn].eb);	/* Bottom edge coordinate */

  *xpt1 =
    ((wds[wn].xmax - wds[wn].xmin) / (xer - xel)) * (xpt - xel) +
    wds[wn].xmin;
  /* Calculate X data value of pixel coord. */
  *ypt1 =
    ((wds[wn].ymax - wds[wn].ymin) / (yet - yeb) * (ypt - yet) +
     wds[wn].ymax);
  /* Calculate Y data value of pixel coord. */
}				/* End calcpt */

/********************************/
void
setxpts (XPoint * xypts, int xpt, int ypt)
{				/* Begin setxpts */
  xypts[0].x = xpt;
  xypts[0].y = ypt;
  xypts[1].x = xpt;
  xypts[1].y = ypt + 1;
  xypts[2].x = xpt;
  xypts[2].y = ypt - 1;
  xypts[3].x = xpt + 1;
  xypts[3].y = ypt;
  xypts[4].x = xpt - 1;
  xypts[4].y = ypt;
}				/* End setxpts */

/********************************/
void
getstr (char *pltstr, float xpt, float ypt)
{				/* Begin getstr */
  char cxstr[10], cystr[10];

  sprintf (cxstr, "%0.2f", xpt);
  sprintf (cystr, "%0.2e", ypt);

  strcpy (pltstr, "( ");
  strcat (pltstr, cxstr);
  strcat (pltstr, ", ");
  strcat (pltstr, cystr);
  strcat (pltstr, " ) ");
}				/* End getstr */

/********************************/
void
doevnt (int wn)
{				/* Begin doevnt */
  /* declarations */
  XEvent myevent;
  KeySym mykey;
  XWindowAttributes myattribs;
  char txt[10];
  int i;
  float xpt1, xpt2, xptd, ypt1, ypt2, yptd;
  int xpt, xptr, ypt, yptr, xlow, ylow;
  XPoint pts1[5], pts2[5];
  char cpltstr[34];

  XNextEvent (wds[wn].dspy, &myevent);	/* read the next event */

  switch (myevent.type)		/* Determine type of event */
    {				/* Begin switch */

      /* Repaint window on expose events */
    case Expose:		/* moved, resized, uncovered etc. window */
      if (myevent.xexpose.count == 0)
	{			/* Begin if */
	  XGetWindowAttributes (wds[wn].dspy, wds[wn].wndw, &myattribs);
	  /* Get information about window configuration */
	  wds[wn].ny = myattribs.height;	/* Set ny to window height in pixels */
	  wds[wn].nx = myattribs.width;	/* Set nx to window width in pixels */

	  calcul (wn);
	  /* Call calcul to fill pixel arrays */
	  wds[wn].led = wds[wn].nx * wds[wn].el;
	  /* Set left edge */
	  wds[wn].red = wds[wn].nx - (wds[wn].nx * wds[wn].er);
	  /* Set right edge */
	  wds[wn].ted = wds[wn].ny * wds[wn].et;
	  /* Set top edge */
	  wds[wn].bed = wds[wn].ny - (wds[wn].ny * wds[wn].eb);
	  /* Set bottom edge */
	  drawln (wn);		/* Call drawln to plot spectra */
	}			/* End if */
      break;			/* End case Expose */

      /* Process keyboard mapping changes */
    case MappingNotify:	/* changes in layout of windows on screen */
      XRefreshKeyboardMapping ((XMappingEvent *) & myevent);	/* make necessary changes */
      break;			/* End case MappingNotify */

      /* Process mouse-button presses */
    case ButtonPress:		/* any mouse-button pressed */
      drawln (wn);
      if (myevent.xbutton.button == 1)
	{
	  xpt = myevent.xbutton.x;
	  ypt = myevent.xbutton.y;

	  calcpt (xpt, ypt, &xpt1, &ypt1, wn);
	  setxpts (pts1, xpt, ypt);

	  XDrawPoints (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, pts1, 5,
		       wds[wn].mode);
	}
      if (myevent.xbutton.button == 3)
	{
	  xpt = myevent.xbutton.x;
	  ypt = myevent.xbutton.y;

	  calcpt (xpt, ypt, &xpt1, &ypt1, wn);
	  setxpts (pts1, xpt, ypt);

	  getstr (cpltstr, xpt1, ypt1);

	  XDrawPoints (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, pts1, 5,
		       wds[wn].mode);
	  if ((xpt + (strlen (cpltstr) * 6)) < wds[wn].nx)
	    /* Doesn't overlap right edge */
	    XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xpt + 6,
			 ypt + 3, cpltstr, strlen (cpltstr));
	  else if ((xpt - (strlen (cpltstr) * 6)) > 0)
	    {
	      XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xpt -
			   (strlen (cpltstr) * 6) + 6, ypt + 3, cpltstr,
			   strlen (cpltstr));
	    }
	  else if ((ypt + 10) > wds[wn].ny)
	    {
	      XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xpt -
			   (strlen (cpltstr) * 6) / 2 + 6,
			   ypt + 10, cpltstr, strlen (cpltstr));
	    }
	  else
	    {
	      XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xpt -
			   (strlen (cpltstr) * 6) / 2 + 6,
			   ypt - 10, cpltstr, strlen (cpltstr));
	    }
	}
      break;			/* End case ButtonPress */

    case ButtonRelease:	/* any mouse-button released */
      if (myevent.xbutton.button == 1)
	{
	  xptr = myevent.xbutton.x;
	  yptr = myevent.xbutton.y;

	  xlow = (xpt - xptr) / 2 + xptr;
	  ylow = (ypt - yptr) / 2 + yptr;

	  calcpt (xptr, yptr, &xpt2, &ypt2, wn);
	  setxpts (pts2, xptr, yptr);

	  xptd = xpt1 - xpt2;
	  yptd = ypt1 - ypt2;

	  getstr (cpltstr, xptd, yptd);

	  XDrawPoints (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, pts2, 5,
		       wds[wn].mode);

	  if ((xlow + 250) < wds[wn].nx)
	    XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xlow + 6,
			 ylow + 3, cpltstr, strlen (cpltstr));
	  else if ((xlow - 250) > 0)
	    XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xlow -
			 (strlen (cpltstr) * 6) + 6,
			 ylow + 3, cpltstr, strlen (cpltstr));
	  else if ((ylow + 10) > wds[wn].ny)
	    XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xlow -
			 (strlen (cpltstr) * 6) / 2 + 6,
			 ylow + 10, cpltstr, strlen (cpltstr));
	  else
	    XDrawString (wds[wn].dspy, wds[wn].wndw, wds[wn].gc, xlow -
			 (strlen (cpltstr) * 6) / 2 + 6,
			 ylow - 10, cpltstr, strlen (cpltstr));
	}
      break;

      /* Process keyboard input */
    case KeyPress:		/* any key pressed */
      i = XLookupString ((XKeyEvent *) & myevent, txt, 10, &mykey, 0);	/* get i and txt */
      if (i == 1)		/* number of characters returned equals 1 */
	{			/* Begin if */

	  switch (txt[0])	/* Determine key pressed */
	    {			/* Begin switch */

	    case 'q':		/* Lowercase 'q' pressed */
	      done = 1;		/* End of loop */
	      break;		/* End case 'q' */

	    case 'o':		/* Lowercase 'o' pressed */
	      wds[wn].flag = 0;	/* Flag set to 'original 2 spectra' */
	      break;		/* End case 'o' */

	    case 'd':		/* Lowercase 'd' pressed */
	      wds[wn].flag = 1;	/* Flag set to 'difference spectrum' */
	      break;		/* End case 'd' */

	    case 'b':		/* Lowercase 'b' pressed */
	      wds[wn].flag = 2;	/* Flag set to 'all three spectra' */
	      break;		/* End case 'b' */

	    case 'e':		/* Lowercase 'e' pressed */
	      wds[wn].flag = 3;	/* Flag set to 'exp. spectrum only' */
	      break;		/* End case 'e' */

	    case 's':		/* Lowercase 's' pressed */
	      wds[wn].flag = 4;	/* Flag set to 'sim. spectrum only' */
	      break;		/* End case 's' */

	    case 'a':		/* Lowercase 'a' pressed */
	      wds[wn].flag2 = wds[wn].flag2 * -1;	/* Switch flag2 off or on */
	      break;		/* End case 'a' */

	    }			/* End switch (txt[0]) */

	  drawln (wn);		/* Redraw spectra */
	}			/* End if */

      break;			/* End case KeyPress */

    }				/* End switch (myevent.type) */
  XFlush (wds[wn].dspy);
}				/* End doevnt */

/********************************/
void FORTRAN (getwndws) (int *n, word * title)
{				/* Begin getwndws */
  /* declarations */
  XSizeHints myhint;
  XColor orng, grn, mgnt, blk, wht;
  Font myfont;
  Status resultb, resultba, resultw, resultwa;
  Status result1a, result1, result2a;
  Status result2, result3a, result3;
  char *name1, *name2, *name3, *nameb, *namew;
  int wn;

  if (getenv ("DISPLAY"))
    displayOK = 1;
  if (!displayOK)
    return;

  if (*n != 0)
    {				/* Begin if */

      /* Clear windows from screen */
      for (wn = 1; wn <= numwndws; ++wn)
	{			/* Begin for */
	  XFreeGC (wds[wn].dspy, wds[wn].gc);	/* Explicitly deallocates mygc */
	  XFreeGC (wds[wn].dspy, wds[wn].gc1);	/* Explicitly deallocates mygc1 */
	  XFreeGC (wds[wn].dspy, wds[wn].gc2);	/* Explicitly deallocates mygc2 */
	  XFreeGC (wds[wn].dspy, wds[wn].gc3);	/* Explicitly deallocates mygc3 */
	  XDestroyWindow (wds[wn].dspy, wds[wn].wndw);
	  XCloseDisplay (wds[wn].dspy);
	}			/* End for */

      numwndws = *n;

#ifdef DEBUG
      printf ("Numwndws = %d\n", numwndws);
#endif

      for (wn = 1; wn <= numwndws; ++wn)
	{			/* Begin for */
	  /* Initialization */
	  wds[wn].dspy = XOpenDisplay ("");	/*  Open display connection */
	  wds[wn].scrn = DefaultScreen (wds[wn].dspy);	/* Get default screen */
	  wds[wn].map = XDefaultColormap (wds[wn].dspy, wds[wn].scrn);
	  /* Get colormap ID */
	  wds[wn].mode = CoordModeOrigin;

#ifdef DEBUG
	  printf ("Getwndws: Wds[%d].dspy = %d\n", wn, wds[wn].dspy);
#endif

	  /* Set default program-specified window position and size */
#if 0
	  myhint.x = 650;
	  myhint.y = (wn - 1) * 400;
	  myhint.width = 500;
	  myhint.height = 400;
#endif
	  myhint.x = (wn - 1) * 400;
	  myhint.y = 0;
	  myhint.width = 400;
	  myhint.height = 300;
	  myhint.flags = PPosition | PSize;

	  /* Create window with properties set in myhint  */
	  wds[wn].wndw = XCreateSimpleWindow (wds[wn].dspy,
					      DefaultRootWindow (wds[wn].
								 dspy),
					      myhint.x, myhint.y,
					      myhint.width, myhint.height, 5,
					      wds[wn].fgnd, wds[wn].bgnd);
	  XSetStandardProperties (wds[wn].dspy, wds[wn].wndw, title[wn - 1],
				  title[wn - 1], None, NULL, 0, &myhint);

	  /* Set name variables to color names */
	  name1 = "orange";
	  name2 = "green";
	  name3 = "magenta";
	  nameb = "black";
	  namew = "white";

	  /* Find color values */
	  result1 = XParseColor (wds[wn].dspy, wds[wn].map, name1, &orng);
	  result2 = XParseColor (wds[wn].dspy, wds[wn].map, name2, &grn);
	  result3 = XParseColor (wds[wn].dspy, wds[wn].map, name3, &mgnt);
	  resultb = XParseColor (wds[wn].dspy, wds[wn].map, nameb, &blk);
	  resultw = XParseColor (wds[wn].dspy, wds[wn].map, namew, &wht);

	  /* Allocate color cells */
	  result1a = XAllocColor (wds[wn].dspy, wds[wn].map, &orng);
	  result2a = XAllocColor (wds[wn].dspy, wds[wn].map, &grn);
	  result3a = XAllocColor (wds[wn].dspy, wds[wn].map, &mgnt);
	  resultba = XAllocColor (wds[wn].dspy, wds[wn].map, &blk);
	  resultwa = XAllocColor (wds[wn].dspy, wds[wn].map, &wht);

	  /* Access apropriate font for strings */
	  myfont = XLoadFont (wds[wn].dspy, "6x13");

	  /* GC creation and initialization */
	  /* Set GC for axes */
	  wds[wn].gc = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	  /* Create GC mygc */
	  XSetBackground (wds[wn].dspy, wds[wn].gc, blk.pixel);
	  /* Set background to black */
	  XSetForeground (wds[wn].dspy, wds[wn].gc, wht.pixel);
	  /* Set foreground to white */
	  XSetFont (wds[wn].dspy, wds[wn].gc, myfont);
	  /* Set font */

	  /* Set GC for experimental spectrum */
	  if ((result1 != 0) && (result1a != 0))	/* Check if colors were 
							   successfully allocated */
	    {			/* Begin if */
	      wds[wn].gc1 = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	      /* Create GC mygc1 */
	      XSetBackground (wds[wn].dspy, wds[wn].gc1, blk.pixel);
	      /* Set background to black */
	      XSetForeground (wds[wn].dspy, wds[wn].gc1, orng.pixel);
	      /* Set foreground to orange */
	    }			/* End if */
	  else
	    {			/* Begin else */
	      wds[wn].gc1 = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	      /* Create GC mygc1 */
	      XSetBackground (wds[wn].dspy, wds[wn].gc1, blk.pixel);
	      /* Set background to black */
	      XSetForeground (wds[wn].dspy, wds[wn].gc1, wht.pixel);
	      /* Set foreground to white */
	    }			/* End else */

	  /* Set GC for simulated spectrum */
	  if ((result2 != 0) && (result2a != 0))
	    {			/* Begin if */
	      wds[wn].gc2 = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	      /* Create GC mygc2 */
	      XSetBackground (wds[wn].dspy, wds[wn].gc2, blk.pixel);
	      /* Set background to black */
	      XSetForeground (wds[wn].dspy, wds[wn].gc2, grn.pixel);
	      /* Set foreground to green */
	    }			/* End if */
	  else
	    {			/* Begin else */
	      wds[wn].gc2 = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	      /* Create GC mygc2 */
	      XSetBackground (wds[wn].dspy, wds[wn].gc2, blk.pixel);
	      /* Set background to black */
	      XSetForeground (wds[wn].dspy, wds[wn].gc2, wht.pixel);
	      /* Set foreground to white */
	    }			/* End else */

	  /* Set GC for difference graph */
	  if ((result3 != 0) && (result3a != 0))
	    {			/* Begin if */
	      wds[wn].gc3 = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	      /* Create GC mygc3 */
	      XSetBackground (wds[wn].dspy, wds[wn].gc3, blk.pixel);
	      /* Set background to black */
	      XSetForeground (wds[wn].dspy, wds[wn].gc3, mgnt.pixel);
	      /* Set foreground to magenta */
	    }			/* End if */
	  else
	    {			/* Begin else */
	      wds[wn].gc3 = XCreateGC (wds[wn].dspy, wds[wn].wndw, 0, 0);
	      /* Create GC mygc3 */
	      XSetBackground (wds[wn].dspy, wds[wn].gc3, blk.pixel);
	      /* Set background to black */
	      XSetForeground (wds[wn].dspy, wds[wn].gc3, wht.pixel);
	      /* Set foreground to white */
	    }			/* End else */

	  /* input event selection */
	  XSelectInput (wds[wn].dspy, wds[wn].wndw,
			/* Accept key or mouse button presses */
			/* or window exposures as event types */
			ButtonPressMask | ButtonReleaseMask |
			OwnerGrabButtonMask | KeyPressMask | ExposureMask);

	  /* window mapping */
	  XMapRaised (wds[wn].dspy, wds[wn].wndw);
	  /* Map window "above" existing windows */

	  /* Initialize external variables */
	  wds[wn].er = 0.1;	/* Set right edge to be .1 of window */
	  wds[wn].el = 0.1;	/* Set left edge to be .1 of window */
	  wds[wn].et = 0.1;	/* Set top edge to be .1 of window */
	  wds[wn].eb = 0.1;	/* Set bottom edge to be .1 of window */

	  /* Initialize indx and flags */
	  wds[wn].flag = 0;	/* Set spectra-to-plot to "original" */
	  wds[wn].flag2 = 1;	/* Set axes to "on" */

#ifdef DEBUG
	  printf ("Getwndws: End of for, wds[%d].dspy = %d\n", wn,
		  wds[wn].dspy);
#endif
	}			/* End for */
    }				/* End if */
}				/* End getwndws */

/********************************/

/* NOTE: Modified by DEB to accept experimental spectrum in y1
   and difference spectrum (simulated - experimental) in y2
*/
void
putary (double *y1, double *y2, int wn)
{				/* Begin putary */
  int i;

  for (i = 0; i < wds[wn].indx; ++i)
    {
      wds[wn].ye[i] = y1[i];
      wds[wn].ys[i] = (y1[i] - y2[i]);
    }
}				/* End putary */

/********************************/
void FORTRAN (fstplt) (double *y1, double *y2, double *xmin1,
		       double *xstep1, int *indx1, int *wnum)
/* double *y1, *y2;               y1 is the experimental, y2 is the
   double *xmin1, *xstep1;        difference (experimental-simulated)
   int *indx1, *wnum; */
{				/* Begin fstplt */
  /* declarations */
  XWindowAttributes myattribs;
  int wn;

#ifdef DEBUG
  printf ("Wnum = %d\n", *wnum);
  printf ("Fstplt: Before setting variables to parameters\n");
#endif
  if (!displayOK)
    return;
  wn = *wnum;

#ifdef DEBUG
  printf ("After wn, wn = %d, before indx\n", wn);
#endif
  wds[wn].indx = *indx1;	/* Set global indx to equal parameter indx1 */

#ifdef DEBUG
  printf ("After indx, before xmin\n");
#endif
  wds[wn].xmin = /**/ *xmin1;

#ifdef DEBUG
  printf ("After xmin, before xstep\n");
#endif
  wds[wn].xstep = /**/ *xstep1;

#ifdef DEBUG
  printf ("Before XGetWindowAttributes, dspy = %d, wndw = %d\n",
	  wds[wn].dspy, wds[wn].wndw);
#endif

  XGetWindowAttributes (wds[wn].dspy, wds[wn].wndw, &myattribs);
  /* Get information about window configuration */

  wds[wn].ny = myattribs.height;	/* Set ny to window height in pixels */
  wds[wn].nx = myattribs.width;	/* Set nx to window width in pixels */

  putary (y1, y2, wn);
  calcul (wn);			/* Call calcul to fill pixel arrays */

  wds[wn].led = wds[wn].nx * wds[wn].el;	/* Set left edge */
  wds[wn].red = wds[wn].nx - (wds[wn].nx * wds[wn].er);	/* Set right edge */
  wds[wn].ted = wds[wn].ny * wds[wn].et;	/* Set top edge */
  wds[wn].bed = wds[wn].ny - (wds[wn].ny * wds[wn].eb);	/* Set bottom edge */
  drawln (wn);			/* Call drawln to plot spectra */
}				/* End fstplt */

/********************************/
void FORTRAN (wpoll) ()
{				/* Begin wpoll */
  int wn;

  if (displayOK == 1)
    {
      for (wn = 1; wn <= numwndws; ++wn)
	{
	  while (XEventsQueued (wds[wn].dspy, QueuedAfterReading) > 0)
	    doevnt (wn);
	}
    }
}				/* End wpoll */

/********************************/
void FORTRAN (shtwndws) ()
{				/* Begin shtwndws */
  /* declarations */
  int wn;
  if (displayOK)
    {
      for (wn = 1; wn <= numwndws; ++wn)
	{
	  /* Terminate window */
	  XFreeGC (wds[wn].dspy, wds[wn].gc);	/* Explicitly deallocates mygc */
	  XFreeGC (wds[wn].dspy, wds[wn].gc1);	/* Explicitly deallocates mygc1 */
	  XFreeGC (wds[wn].dspy, wds[wn].gc2);	/* Explicitly deallocates mygc2 */
	  XFreeGC (wds[wn].dspy, wds[wn].gc3);	/* Explicitly deallocates mygc3 */
	  XDestroyWindow (wds[wn].dspy, wds[wn].wndw);
	  XCloseDisplay (wds[wn].dspy);
	}
      numwndws = 0;
    }
}				/* End shtwndws */

/********************************/
