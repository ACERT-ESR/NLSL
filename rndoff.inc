c      Version 1.5  5/2/94
c*********************************************************************
c
c	declarations of	machine	epsilon	parameter
c
c	Notes:
c		1) The machine epsilon,	eps, is	defined	to be the 
c		   largest number such that the	result of the 
c		   expression
c
c			delta=(1.0D0+eps)-1.0D0
c
c		   returns the value delta=0.0D0.  In other words,
c		   1.0D0+eps is	smallest double	precision floating
c		   point number	greater	than 1.0D0.
c
c		2) A FORTRAN 77	program	to determine the machine 
c		   epsilon is included below for completeness.
c		   This	program	should be compiled with	all 
c		   compiler optimization disabled.  The	parameter
c		   rndoff declared in this file	should assigned	to 
c		   be a	small multiple of the larger of	the two	
c		   numbers printed by this program.  This program 
c		   is an adaptation of the function subroutine epslon 
c		   supplied with EISPACK.
c
c		3) On most DEC machines, the two methods of 
c		   calculating the machine epsilon give	the nearly
c		   same	the result, eps=2.77555756D-17.	 An
c		   appropriate value of	rndoff is then 1.0D-16.	 An
c		   example of a	machine	which gives grossly 
c		   different results is	the Prime 9955.	
c
c	written	by DJS 11-SEP-87
c
c*********************************************************************
c	
c	program	meps
c
c	double precision a,b,c,eps
c
c#######################################################################
c
c	a=1.0D0
c10	b=a+1.0D0
c	eps=c
c	c=b-1.0D0
c	if (c.eq.0.0D0)	go to 20
c	a=a/2.0D0
c	go to 10
c
c20	write (*,*) 'first epsilon = ',eps
c
c	a=4.0D0/3.0D0
c30	b=a-1.0D0
c	c=b+b+b
c	eps=dabs(c-1.0D0)
c
c	write (*,*) 'second epsilon = ',eps
c
c	stop
c	end
c
c***********************************************************************
c
	double precision RNDOFF
	parameter (RNDOFF=3.0D-15)
c
