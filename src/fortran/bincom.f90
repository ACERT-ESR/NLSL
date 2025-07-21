c NLSL Version 1.9.0 beta 2/12/15
c----------------------------------------------------------------------
c                    =========================
c                          module BINCOM
c                    =========================
c
c     Storage for pre-calculated binomial coefficients used by
c     subroutine wig3j (included in file 'w3j.f') for calculation
c     of Wigner 3J symbols. Used exclusively by routines in that file.
c
c----------------------------------------------------------------------
c
      module bincom
      use maxl
      implicit none
c
      integer, parameter :: nb=2*(mxlval+8)+2
      integer, parameter :: nbncf=nb*(nb+1)+1
c
      integer, dimension(nb), save :: bncfx
      double precision, dimension(nbncf), save :: bncf
c
      end module bincom
