/*
 *	SIGINT handler
 *	rewritten to use portable standard C routines.
 *	Note that handlers cannot nest!
 *	- WH 20020210
 */

#include <stddef.h>
#include <signal.h>
#include "fortrancall.h"

static int *flagptr=NULL;
void (*prev_handler)(int);


/***********************************************************************
 *
 *   ctrlc_handler: simple handler for receiving SIGINT signal
 *
 ***********************************************************************/
void sigint_handler(int sig)
{
	*flagptr=1;
}

/* **********************************************************************
 *    subroutine catchc
 *
 *       Installs a SIGINT handler that sets flag non-zero
 *
 ***********************************************************************/

void FORTRAN(catchc)(int *flag)
{
	if (flagptr==NULL) {
		flagptr=flag;
		*flagptr=0;
		prev_handler=signal(SIGINT, &sigint_handler);
	}
}

/* **********************************************************************
 *
 * subroutine uncatchc()
 *
 *    Restores default action for handling SIGINT signal
 *
 ********************************************************************** */
void FORTRAN(uncatchc)(int *flag)
{
	if (flagptr==flag) {
		if (flagptr) *flagptr=0; /* test */
		signal(SIGINT, prev_handler);
		flagptr=NULL;
	}
}

