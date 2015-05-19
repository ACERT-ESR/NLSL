c NLSL Version 1.5 7/2/95
c---------------------------------------------------------------------
c                      =========================
c                          subroutine ORDRPR
c                      =========================
c
c This routine calculates the molecular order parameters corresponding 
c to a specified axial orienting potential. The potential is 
c expressed as
c
c                                L   L
c    U(phi,theta,0) = Sum   Sum c   D   (phi,theta,0)
c                      L     K   K   0K
c 
c where D^L_{0K} is a Wigner rotation matrix element (equivalent to
c the spherical harmonic Y^L_K), c^L_K are user-specified coefficients,
c and the sums are restricted to L=0 to 4 and K=0 to L (for even L,K only)
c
c The coefficients are specified in the cf array in the order 
c c20, c22, c40, c42, c44. The routine calculates an order parameter for 
c each nonzero coefficient specififed. The order parameter is the average
c of <D^L_{0K}(phi,theta)>  weighted by the value of exp(-U(phi,theta)) 
c and integrated over angles theta and phi.
c
c
c Includes:
c   rndoff.inc
c 
c Uses:
c    ccrints.f  (Integration routines ccrint ccrin1)
c    ftht       (Function returning theta-dependent part of U)
c    fphi       (Function returning phi-dependent part of U)
c
c    The last two functions are included in this file.
c    Written by David E. Budil 
c    Cornell Department of Chemistry, March 1992
c 
c----------------------------------------------------------------------

      subroutine ordrpr(cf,d)
c
      use rnddbl
c
      implicit none
      double precision cf(5),d(5)
c
      double precision ZERO,ONE,SMALL
      parameter(ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-16)
c
      integer i,npt,nup,id
      double precision fz,pfn
c
      integer kount
      double precision acc,c,d20,d40,bd22,bd42,bd44
      common/potprm/acc,c(5),kount
      common/func/d20,d40,bd22,bd42,bd44
c
      double precision ftht
      external ftht
c......................................................................
      acc=0.001
      npt=0
c
      do kount=1,5
         c(kount)=cf(kount)
         if (dabs(c(kount)).gt.rndoff) npt=kount
      end do
c
      if (npt.gt.0) then        
c
c--------------------------------------------------
c Integrate unweighted potential function (kount=0)
c--------------------------------------------------
         kount=0
         call ccrint(ZERO,ONE,acc,SMALL,pfn,nup,ftht,id)
c
c-------------------------------------------------------------------
c Integrate potential weighted by D20,D22,D40,D42,D44 for kount=1..5
c Only calculate up to order of highest-order potential coefficient (npt)
c-------------------------------------------------------------------
         do kount=1,npt
            call ccrint(ZERO,ONE,acc,SMALL,fz,nup,ftht,id)
            d(kount)=fz/pfn
         end do
      endif
      return
      end

c----------------------------------------------------------------------
c                    =========================
c                          function FTHT
c                    =========================
c
c     Contains theta-dependence of orienting pseudopotential
c----------------------------------------------------------------------
c
      function ftht(ctht)
c
      use pidef
c
      implicit none
      double precision ftht,ctht
c
      integer nup,id
      double precision ctht2,stht2,result
c
      integer kount
      double precision acc,c,d20,d40,bd22,bd42,bd44
      common/potprm/acc,c(5),kount
      common/func/d20,d40,bd22,bd42,bd44
c
      double precision fphi
      external fphi
c
c---------------------------------------------------------
c  definition of some constants
c    A22 = sqrt(3/2), A42 = sqrt(5/8),  A44 = sqrt(35/8)
c---------------------------------------------------------
      double precision A22,A42,A44,ONE,ZERO,SMALL
      parameter (A22=1.22474487139159d0,
     2           A42=0.790569415042095d0,
     3           A44=1.04582503316759d0 )
      parameter(ONE=1.0D0,ZERO=0.0D0,SMALL=1.0D-16 )
c
c......................................................................
c
      ctht2=ctht*ctht
      stht2=ONE-ctht2
      d20 =1.5d0*ctht2-0.5d0
      d40 =( (4.375d0*ctht2)-3.75d0)*ctht2+0.375d0
      bd22=A22*stht2
      bd42=A42*stht2*(7.0d0*ctht2-ONE)
      bd44=A44*stht2*stht2
c
       call ccrin1(ZERO,PI,acc,SMALL,result,nup,fphi,id)
       ftht=result
       return
       end
                               
c----------------------------------------------------------------------
c                    =========================
c                          function FPHI
c                    =========================
c
c     Contains phi-dependence of orienting pseudopotential
c----------------------------------------------------------------------
      function fphi(phi)
      implicit none
      double precision fphi,phi,c2phi,c4phi
c
      integer kount
      double precision acc,c,d20,d40,bd22,bd42,bd44
      common/potprm/acc,c(5),kount
      common/func/d20,d40,bd22,bd42,bd44
c
      double precision ONE,TWO
      parameter ( ONE=1.0D0,
     #            TWO=2.0D0 )
c
      c2phi=dcos(phi+phi)
      c4phi=TWO*c2phi*c2phi - ONE
c
       fphi=dexp(  c(1)*d20
     2           + c(2)*bd22*c2phi
     3           + c(3)*d40
     4           + c(4)*bd42*c2phi
     5           + c(5)*bd44*c4phi )
c
      if(kount.eq.0) return
      if(kount.eq.1) fphi=d20*fphi
      if(kount.eq.2) fphi=bd22*fphi*c2phi
      if(kount.eq.3) fphi=d40*fphi
      if(kount.eq.4) fphi=bd42*fphi*c2phi
      if(kount.eq.5) fphi=bd44*fphi*c4phi
      return
      end

c----------------------------------------------------------------------
c                    =========================
c                         function FU20
c                    =========================
c  Function defining Boltzmann population at a given orientation
c  resulting from an orienting pseudopotential using only the c20
c  term. Used for calculating weighting as a function of angle
c  between the average director and the spectrometer field, for
c  the "partial MOMD" model. The potential is intended to describe
c  partial macroscopic orientation of the directors of local domains.
c
c  This is required to be a single-argument function so that it 
c  may be called by ccrint to calculate the normalization factor
c  for orientational weighting coefficients. The second parameter, the 
c  director c20 coefficient, is passed in common /dfunc/
c----------------------------------------------------------------------
      function fu20(cost)
c
      use dfunc
c
      implicit none
      double precision fu20,cost
c
      fu20 = exp( dlam*(1.5d0*cost*cost-0.5d0 ) )
      return
      end
c
c
c  phi is expected in radians
c
      function fu20phi(phi)
c
      use dfunc
      use pidef
c
      implicit none
      double precision costht,fu20phi,phi
c
      costht = ss*cos(phi)+cc
      fu20phi = exp( dlam*(1.5d0*costht*costht-0.5d0 ) )
      return
      end
