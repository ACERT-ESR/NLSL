! VC        NLS VERSION 20 May 1992
!**********************************************************************
!
!             complex double precision vector scale and add
!             ---------------------------------------------
!
!        This subroutine will update a complex double precision
!        vector Y by the formula,
!
!                          Y=a*X+Y ,
!
!        where a is a complex double precision constant and X is a
!        complex double presion vector.
!
!       Includes:
!               nlsdim.inc
!               rndoff.inc
!               
!       Uses:
!
!**********************************************************************
!
      subroutine zaxpy(x,y,scale,ndim)
!
      use nlsdim
!*djs      use rnddbl
!
      integer ndim
      double complex scale
!
      double complex x,y
      dimension x(MXDIM),y(MXDIM)
!
      integer iel
!
!*djs      double complex CZERO
!*djs      parameter (CZERO=(0.0D0,0.0D0))
!
!######################################################################
!
      do iel=1,ndim
        y(iel)=scale*x(iel)+y(iel)
!*djs        if (abs(y(iel)).lt.rndoff) y(iel)=CZERO
      end do
!
      return
      end
