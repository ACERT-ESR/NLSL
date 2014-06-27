c VC        NLS VERSION 20 May 1992
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
c               nlsdim.inc
c               rndoff.inc
c               
c       Uses:
c
c**********************************************************************
c
      subroutine zaxpy(x,y,scale,ndim)
c
      include 'nlsdim.inc'
c*djs      include 'rndoff.inc'
c
      integer ndim
      double complex scale
c
      double complex x,y
      dimension x(MXDIM),y(MXDIM)
c
      integer iel
c
c*djs      double complex CZERO
c*djs      parameter (CZERO=(0.0D0,0.0D0))
c
c######################################################################
c
      do iel=1,ndim
        y(iel)=scale*x(iel)+y(iel)
c*djs        if (abs(y(iel)).lt.rndoff) y(iel)=CZERO
      end do
c
      return
      end
