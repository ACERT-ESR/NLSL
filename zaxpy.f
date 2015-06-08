c     Version 1.5   5/2/94 
c**********************************************************************
c
c             complex double precision vector scale and add
c             ---------------------------------------------
c
c        This subroutine will update a complex double precision
c        vector Y by the formula,
c
c                          Y=a*X+Y ,
c
c        where a is a complex double precision constant and X is a
c        complex double presion vector.
c
c       Includes:
c               stddim.inc
c               rndoff.inc
c               
c       Uses:
c
c**********************************************************************
c
      subroutine zaxpy(x,y,scale,ndim)
c
      include 'stddim.inc'
c
      integer ndim,iel
      complex*16 x,y,scale
      dimension x(MXDIM),y(MXDIM)
c
c*djs      include 'rndoff.inc'
c*djs      complex*16 CZERO
c*djs      parameter (CZERO=(0.0D0,0.0D0))
c
c######################################################################
c
      do iel=1,ndim
        y(iel)=scale*x(iel)+y(iel)
c*djs        if (abs(y(iel)).lt.RNDOFF) y(iel)=CZERO
      end do
      return
      end
