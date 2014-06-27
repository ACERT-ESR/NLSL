c Version 1.1  12/19/91
c----------------------------------------------------------------------
c                    =========================
c                          common BINCOM
c                    =========================
c
c     Storage for pre-calculated binomial coefficients used by
c     subroutine wig3j (included in file 'w3j.f') for calculation
c     of Wigner 3J symbols. Used exclusively by routines in that file.
c
c----------------------------------------------------------------------
c
      integer nb,nbncf
      parameter (nb=2*(mxlval+8)+2)
      parameter (nbncf=nb*(nb+1)+1)
c
      integer bncfx
      double precision bncf
      common /bincom/bncf(nbncf),bncfx(nb)
c
