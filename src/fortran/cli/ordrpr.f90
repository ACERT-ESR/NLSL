! NLSL Version 1.5 7/2/95
!---------------------------------------------------------------------
!                      =========================
!                          subroutine ORDRPR
!                      =========================
!
! This routine calculates the molecular order parameters corresponding 
! to a specified axial orienting potential. The potential is 
! expressed as
!
!                                L   L
!    U(phi,theta,0) = Sum   Sum c   D   (phi,theta,0)
!                      L     K   K   0K
! 
! where D^L_{0K} is a Wigner rotation matrix element (equivalent to
! the spherical harmonic Y^L_K), c^L_K are user-specified coefficients,
! and the sums are restricted to L=0 to 4 and K=0 to L (for even L,K only)
!
! The coefficients are specified in the cf array in the order 
! c20, c22, c40, c42, c44. The routine calculates an order parameter for 
! each nonzero coefficient specififed. The order parameter is the average
! of <D^L_{0K}(phi,theta)>  weighted by the value of exp(-U(phi,theta)) 
! and integrated over angles theta and phi.
!
!
! Includes:
!   rndoff.inc
! 
! Uses:
!    ccrints.f  (Integration routines ccrint ccrin1)
!    ftht       (Function returning theta-dependent part of U)
!    fphi       (Function returning phi-dependent part of U)
!
!    The last two functions are included in this file.
!    Written by David E. Budil 
!    Cornell Department of Chemistry, March 1992
! 
!----------------------------------------------------------------------

      subroutine ordrpr(cf,d)
!
      use rnddbl
!
      implicit none
      double precision cf(5),d(5)
!
      double precision ZERO,ONE,SMALL
      parameter(ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-16)
!
      integer npt,nup,id
      double precision fz,pfn
!
      integer kount
      double precision acc,c,d20,d40,bd22,bd42,bd44
      common/potprm/acc,c(5),kount
      common/func/d20,d40,bd22,bd42,bd44
!
      double precision ftht
      external ftht
!......................................................................
      acc=0.001
      npt=0
!
      do kount=1,5
         c(kount)=cf(kount)
         if (dabs(c(kount)).gt.rndoff) npt=kount
      end do
!
      if (npt.gt.0) then        
!
!--------------------------------------------------
! Integrate unweighted potential function (kount=0)
!--------------------------------------------------
         kount=0
         call ccrint(ZERO,ONE,acc,SMALL,pfn,nup,ftht,id)
!
!-------------------------------------------------------------------
! Integrate potential weighted by D20,D22,D40,D42,D44 for kount=1..5
! Only calculate up to order of highest-order potential coefficient (npt)
!-------------------------------------------------------------------
         do kount=1,npt
            call ccrint(ZERO,ONE,acc,SMALL,fz,nup,ftht,id)
            d(kount)=fz/pfn
         end do
      endif
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                          function FTHT
!                    =========================
!
!     Contains theta-dependence of orienting pseudopotential
!----------------------------------------------------------------------
!
      function ftht(ctht)
!
      use pidef
!
      implicit none
      double precision ftht,ctht
!
      integer nup,id
      double precision ctht2,stht2,result
!
      integer kount
      double precision acc,c,d20,d40,bd22,bd42,bd44
      common/potprm/acc,c(5),kount
      common/func/d20,d40,bd22,bd42,bd44
!
      double precision fphi
      external fphi
!
!---------------------------------------------------------
!  definition of some constants
!    A22 = sqrt(3/2), A42 = sqrt(5/8),  A44 = sqrt(35/8)
!---------------------------------------------------------
      double precision A22,A42,A44,ONE,ZERO,SMALL
      parameter (A22=1.22474487139159d0,
     2           A42=0.790569415042095d0,
     3           A44=1.04582503316759d0 )
      parameter(ONE=1.0D0,ZERO=0.0D0,SMALL=1.0D-16 )
!
!......................................................................
!
      ctht2=ctht*ctht
      stht2=ONE-ctht2
      d20 =1.5d0*ctht2-0.5d0
      d40 =( (4.375d0*ctht2)-3.75d0)*ctht2+0.375d0
      bd22=A22*stht2
      bd42=A42*stht2*(7.0d0*ctht2-ONE)
      bd44=A44*stht2*stht2
!
       call ccrin1(ZERO,PI,acc,SMALL,result,nup,fphi,id)
       ftht=result
       return
       end
                               
!----------------------------------------------------------------------
!                    =========================
!                          function FPHI
!                    =========================
!
!     Contains phi-dependence of orienting pseudopotential
!----------------------------------------------------------------------
      function fphi(phi)
      implicit none
      double precision fphi,phi,c2phi,c4phi
!
      integer kount
      double precision acc,c,d20,d40,bd22,bd42,bd44
      common/potprm/acc,c(5),kount
      common/func/d20,d40,bd22,bd42,bd44
!
      double precision ONE,TWO
      parameter ( ONE=1.0D0,
     #            TWO=2.0D0 )
!
      c2phi=dcos(phi+phi)
      c4phi=TWO*c2phi*c2phi - ONE
!
       fphi=dexp(  c(1)*d20
     2           + c(2)*bd22*c2phi
     3           + c(3)*d40
     4           + c(4)*bd42*c2phi
     5           + c(5)*bd44*c4phi )
!
      if(kount.eq.0) return
      if(kount.eq.1) fphi=d20*fphi
      if(kount.eq.2) fphi=bd22*fphi*c2phi
      if(kount.eq.3) fphi=d40*fphi
      if(kount.eq.4) fphi=bd42*fphi*c2phi
      if(kount.eq.5) fphi=bd44*fphi*c4phi
      return
      end

!----------------------------------------------------------------------
!                    =========================
!                         function FU20
!                    =========================
!  Function defining Boltzmann population at a given orientation
!  resulting from an orienting pseudopotential using only the c20
!  term. Used for calculating weighting as a function of angle
!  between the average director and the spectrometer field, for
!  the "partial MOMD" model. The potential is intended to describe
!  partial macroscopic orientation of the directors of local domains.
!
!  This is required to be a single-argument function so that it 
!  may be called by ccrint to calculate the normalization factor
!  for orientational weighting coefficients. The second parameter, the 
!  director c20 coefficient, is passed in common /dfunc/
!----------------------------------------------------------------------
      function fu20(cost)
!
      use dfunc
!
      implicit none
      double precision fu20,cost
!
      fu20 = exp( dlam*(1.5d0*cost*cost-0.5d0 ) )
      return
      end
!
!
!  phi is expected in radians
!
      function fu20phi(phi)
!
      use dfunc
      use pidef
!
      implicit none
      double precision costht,fu20phi,phi
!
      costht = ss*cos(phi)+cc
      fu20phi = exp( dlam*(1.5d0*costht*costht-0.5d0 ) )
      return
      end
