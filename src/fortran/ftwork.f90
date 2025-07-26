! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                    =========================
!                          module FTWORK
!                    =========================
!
! Defines and saves working arrays used for correlation and
! convolution function calculations by Fourier transform methods
!
!    tmpclc     Temporary storage array for calculated spectrum 
!               used by function sshift (calls subroutine correl)
!               and by subroutine eprls (and subroutine convlv)
!
!    tmpdat     Temporary storage array for data
!               used by function sshift (calls subroutine correl)
!               and by subroutine eprls (and subroutine convlv)
!
!    work,work2 Temporary work array used by correl (and convlv)
!
!----------------------------------------------------------------------
!
      module ftwork
      use nlsdim
      implicit none
!
      double precision, dimension (MXPT), save :: 
     #   work, work2, tmpclc, tmpdat
!
      end module ftwork

