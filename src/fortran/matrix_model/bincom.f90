! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                    =========================
!                          module BINCOM
!                    =========================
!
!     Storage for pre-calculated binomial coefficients used by
!     subroutine wig3j (included in file 'w3j.f') for calculation
!     of Wigner 3J symbols. Used exclusively by routines in that file.
!
!----------------------------------------------------------------------
!
      module bincom
      use maxl
      implicit none
!
      integer, parameter :: nb=2*(mxlval+8)+2
      integer, parameter :: nbncf=nb*(nb+1)+1
!
      integer, dimension(nb), save :: bncfx
      double precision, dimension(nbncf), save :: bncf
!
      end module bincom
