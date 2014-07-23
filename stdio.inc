c  NLSPMC VERSION 1.0 2/5/99
c*********************************************************************
c
c	declarations of	standard I/O parameters	
c
c	Notes:
c		1) The parameters define here are used as
c
c			luttyi : logical unit for terminal input
c			luttyo : logical unit for terminal output
c                       lulog  : logical unit for log file
c                       lutrc  : logical unit for trace file
c			ludisk : logical unit for other disk files
c			ludiskb : logical unit for unformatted disk files
c                       
c                       NOTE: The units that will be used for additional
c                       disk input depend on the MXSPEC parameter defined
c                       in "nlsdim.inc". If MXSPEC is greater than 1, 
c                       units <ludisk+1>, <ludisk+2>, ... <ludisk+MXSPEC>
c                       are assumed to be available.
c
c		2) The appropriate values for these parameters are
c		   highly operating system and compiler	dependent.
c
c*********************************************************************
c
      integer luttyi,luttyo,ludisk,ludiskb,lulog,lutrc,ludeb
      parameter (ludeb=3,lulog=4,luttyi=5,luttyo=6,lutrc=7,ludisk=9,
     #		ludiskb=8)
c
      integer lucmd,luecho,luout
      common /stdio/ lucmd,luecho,luout

      external mclock
      integer time1, time2, mclock
      real seconds
      common /atimes/seconds,time1,time2
