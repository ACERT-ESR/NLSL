! NLSL Version 1.9.0 beta 2/13/15
!***********************************************************************
!                   =========================
!                          module PIDEF
!                   =========================
!
!         Define the constant pi and conversion factor
!         for degrees-to-radians
!
!***********************************************************************
!
      module pidef
      implicit none
!
      double precision , parameter ::
     #    PI=3.1415926535897932384D0, RADIAN=PI/180.0D0
!
      end module pidef
