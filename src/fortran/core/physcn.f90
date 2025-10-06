! NLSL Version 1.9.0 beta 2/13/15
!----------------------------------------------------------------------
!                   =========================
!                         module PHYSCN
!                   =========================
!
! Physical constants (in cgs/electromagnetic units)
!    ge:     free-electron g-value
!    betae:  Bohr magneton
!    hbar:   Planck's constant
!    gammae: electronic gyromagnetic ratio
!    kb:     Boltzmann's constant
!
! --- updated to CODATA 1998 values 20020213 WH
!----------------------------------------------------------------------
!
      module physcn
      implicit none
! 
      double precision, parameter ::
     #             GE=2.002 319 304 3737 D0,
     #             BETAE=927.400 899 D-23,
     #             HBAR=1.054 571 596 D-27,
     #             GAMMAE=1.760 859 794 D7,
     #             KB=1.380 6503 D-16
!
      end module physcn
