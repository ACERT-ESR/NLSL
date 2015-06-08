c      Version 1.5  5/2/94
c----------------------------------------------------------------------
c
c        subroutine for doing the complex vector scaling and
c        swapping at the end of every Lanczos step.
c
c        written by DJS 11-Nov-86
c
c       Includes:
c               stddim.inc
c               rnfoff.inc
c
c       Uses:
c
c----------------------------------------------------------------------
c
      subroutine zscsw(x,y,scale,ndim)
c
      include 'stddim.inc'
c
      integer ndim
      complex*16 x,y,scale
      dimension x(MXDIM),y(MXDIM)
c
      integer iel
      complex*16 scalex,scaley,tempc
c
c######################################################################
c
      scaley=1.0D0/scale
      scalex=-scale
c
      do 10 iel=1,ndim
        tempc=y(iel)
        y(iel)=scaley*x(iel)
        x(iel)=scalex*tempc
 10   continue
c
      return
      end
