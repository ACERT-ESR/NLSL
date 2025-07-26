! NLS VERSION
!**********************************************************************
!               
!       This double precision function evaluates the integrand of 
!       the orientational integral in the definition of the starting 
!       vector.  The integral over the gamma Euler angle is done 
!       analytically.
!
!       Notes:
!               1) The L and K quantum numbers are passed via a 
!                  common block.
!               2) For more information on the evaluation of the 
!                  modified Bessel and associated Legendre functions
!                  see the respective routines.
!               3) The often used special case of lptmx=2, kptmx=0
!                  is handled separately for faster execution.
!
!       written by DJS 10-SEP-87
!
!       Includes:
!               eprprm.inc
!
!       Uses:
!               bessi.f
!               plgndr.f
!
!**********************************************************************
!
      double precision function fz(z)
!
      use eprprm
!
      implicit none
      double precision z
!
      integer lr,kr,iberr
      common/ifzdat/lr,kr,iberr
!
      double precision a,b
      integer k,ierr
!
      double precision DSQ24,DSQ360
      parameter (DSQ24=4.89897948556635619639d0,
     #           DSQ360=18.97366596101027599199d0)
!
      double precision bessi,plgndr
      external bessi,plgndr
!
!######################################################################
!
      iberr=0
!
      if((lptmx.eq.2).and.(kptmx.eq.0)) then
         if(kr.eq.0) then
            fz=dexp(0.5D0*cpot(2,1)*plgndr(2,0,z))*
     #           plgndr(lr,kr,z)
         else
            fz=0.0D0
         end if
      else 
         a=0.5D0*(cpot(2,1)*plgndr(2,0,z)
     #        +cpot(3,1)*plgndr(4,0,z))
         if (kptmx.ne.0) then
            b=cpot(2,2)*plgndr(2,2,z)/DSQ24
     #       +cpot(3,2)*plgndr(4,2,z)/DSQ360
         else
            b=0.0D0
         end if
         k=kr/2
         fz=bessi(k,b,ierr)*dexp(a)*plgndr(lr,kr,z)
         if (ierr.ne.0) iberr=-1
      end if
!
      return
      end
