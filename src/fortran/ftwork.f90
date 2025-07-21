c NLSL Version 1.9.0 beta 2/12/15
c----------------------------------------------------------------------
c                    =========================
c                          module FTWORK
c                    =========================
c
c Defines and saves working arrays used for correlation and
c convolution function calculations by Fourier transform methods
c
c    tmpclc     Temporary storage array for calculated spectrum 
c               used by function sshift (calls subroutine correl)
c               and by subroutine eprls (and subroutine convlv)
c
c    tmpdat     Temporary storage array for data
c               used by function sshift (calls subroutine correl)
c               and by subroutine eprls (and subroutine convlv)
c
c    work,work2 Temporary work array used by correl (and convlv)
c
c----------------------------------------------------------------------
c
      module ftwork
      use nlsdim
      implicit none
c
      double precision, dimension (MXPT), save :: 
     #   work, work2, tmpclc, tmpdat
c
      end module ftwork

