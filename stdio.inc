c   Version 1.5  5/2/94
c*********************************************************************
c
c	declarations of	standard I/O parameters	
c
c	Notes:
c		1) The parameters define here are used as
c
c			LUTTYI : logical unit for terminal input
c			LUTTYO : logical unit for terminal output
c			LUDISK : logical unit for disk files
c			LULOG  : logical unit for disk log files
c
c		2) The appropriate values for these parameters are
c		   highly operating system and compiler	dependent.
c
c	written	by DJS 11-Nov-86
c
c*********************************************************************
c
	integer	LUTTYI,LUTTYO,LUDISK,LULOG
	parameter (LUTTYI=5,LUTTYO=6,LUDISK=3,LULOG=4)
c
