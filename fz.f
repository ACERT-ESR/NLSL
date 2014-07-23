c NLSMC VERSION 1.0   2/5/99
c**********************************************************************
c               
c       This double precision function evaluates the integrand of 
c       the orientational integral in the definition of the starting 
c       vector.  The integral over the gamma Euler angle is done 
c       analytically.
c
c       Notes:
c               1) The L and K quantum numbers are passed via a 
c                  common block.
c               2) For more information on the evaluation of the 
c                  modified Bessel and associated Legendre functions
c                  see the respective routines.
c               3) The often used special case of lptmx=2, kptmx=0
c                  is handled separately for faster execution.
c
c       written by DJS 10-SEP-87
c
c       Includes:
c               nlsdim.inc
c               eprprm.inc
c
c       Uses:
c               bessi.f
c               plgndr.f
c
c**********************************************************************
c
      double precision function fz(z)
c
      double precision z
c
      include 'limits.inc'
      include 'simparm.inc'
c
      integer lr,kr
      common/ifzdat/lr,kr
c
      double precision a,b
      integer k
c
      double precision dsq24,dsq360
      parameter (dsq24=4.89897948556635619639,
     #           dsq360=18.97366596101027599199)
c
      double precision bessi,plgndr
      external bessi,plgndr
c
c######################################################################
c
      if((lptmx.eq.2).and.(kptmx.eq.0)) then
        if(kr.eq.0) then
          fz=dexp(0.5D0*cpot(2,1)*plgndr(2,0,z))*
     #         plgndr(lr,kr,z)
        else
          fz=0.0D0
        end if
      else 
        a=0.5D0*(cpot(2,1)*plgndr(2,0,z)
     #           +cpot(3,1)*plgndr(4,0,z))
        if (kptmx.ne.0) then
          b=cpot(2,2)*plgndr(2,2,z)/dsq24
     #     +cpot(3,2)*plgndr(4,2,z)/dsq360
        else
          b=0.0D0
        end if
        k=kr/2
        fz=bessi(k,b)*dexp(a)*plgndr(lr,kr,z)
      end if
c
c----------------------------------------------------------------------
c     return to caller
c----------------------------------------------------------------------
c
      return
      end
