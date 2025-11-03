! NLSL Version 1.9.0 beta 2/3/15
!----------------------------------------------------------------------
!                    =========================
!                          module ERRMSG
!                    =========================
!
! Contains text for error messages and definitions of symbols used
! to define error message numbers. Errors after ZEROB0=16 are defined
! as fatal, and will cause the calculation to halt.
!
! This module replaces errmsg.inc and the second half of nlstxt.f.
!
! For eprerr, the definitions of the symbols are as follows:
!
! ZEROIN2     Zero nuclear spin assumed
! POSB0       Positive B0 assumed
! BADHIFLD    High-field approx. may not apply
! BADJMP      Nonzero potential with jump/free: ipdf set to 0
! BADAV       Zero potential with aniso. viscos.: ipdf set to 0
! BADDJ       Discrete jumps with aniso. viscos.: ipdf set 0
! LEMXHI      lemx must be 48 or less with potenti
! BSSADJ      Basis set adjusted to (iii,iii,iii,iii,iii,iii,ii)
! NSTEPHI     NSTEP too large: adjusted to iiiii
! CGTOLHI     CGTOL too small: adjusted to xxxxxxxxx
! SHIFTHI     SHIFTR too small: adjusted to xxxxxxxxx
! NOCONVRG    CG did not converge after iiiii steps
! NOALPHAM    Angle alm not needed with axial g-matrix
! NOALPHAD    Angle ald not needed with axial diffusion tensor
! ZEROB0      Zero B0
! ZEROG       Zero g-tensor element
! TDGERR      Error in tridiagonal matrix storage
! MTXBIG      Too many matrix elements
! DIMBIG      Matrix dimension too large
! MTXHLT      User halt during matrix calculation
! CGHLT       User halt during CG tridiagonalization
! TDGBIG      Too many tridiagonal matrix elements
! NOMIN       Bracketing routine (search command) did not find a minimum
! ZERONFLD    NFIELD is zero (no points in spectrum)
! ZRORANGE    RANGE is zero 
! 
! For minerr, the meanings of the INFO index values are found in the
! comments at the beginning of lmnls.f.
!
!----------------------------------------------------------------------

      module errmsg
      implicit none

      integer, parameter :: FATAL=16,NEPERR=26,NMNERR=12
!
      integer, parameter :: ZEROIN2=1,POSB0=2,
     #            BADHIFLD=3,BADJMP=4,BADAV=5,BADDJ=6,
     #            LEMXHI=7,BSSADJ=8,NSTEPHI=9,CGTOLHI=10,SHIFTHI=11,
     #            NOCONVRG=12,NOALPHAM=13,NOALPHAD=14,BADBESS=15,
     #            ZEROB0=16,ZEROG=17,TDGERR=18,MTXBIG=19,
     #            DIMBIG=20,MTXHLT=21,CGHLT=22,TDGBIG=23,NOMIN=24,
     #            ZERONFLD=25,ZRORANGE=26

      character*50, dimension(NEPERR), save ::
     #  eprerr =
     1  (/ 'Zero nuclear spin assumed                         ' ,
     2     'Positive B0 assumed                               ' ,
     3     'High-field approx. may not apply                  ' ,
     4     'Nonzero potential with jump/free: ipdf set to 0   ' ,
     5     'Zero potential with aniso. viscos.: ipdf set to 0 ' ,
     6     'Discrete jumps with aniso. viscos.: ipdf set 0    ' ,
     7     'lemx must be 48 or less with potential            ' ,
     8     'Basis set adjusted to (iii,iii,iii,iii,iii,iii,ii)' ,
     9     'NSTEP too large: adjusted to iiiii                ' ,
     A     'CGTOL too small: adjusted to xxxxxxxxx            ' ,
     1     'SHIFTR too small: adjusted to xxxxxxxxx           ' ,
     2     'CG did not converge after iiiii steps             ' ,
     3     'Angle alpham not needed w/ axial g-matrix         ' ,
     4     'Angle alphad not needed w/ axial diffusion tensor ' ,
     5     'Taylor series for Bessel function did not converge' ,
     6     'Zero B0                                           ' ,
     7     'Zero g-tensor element                             ' ,
     8     'Error in tridiagonal matrix storage               ' ,
     9     'Too many matrix elements                          ' ,
     A     'Matrix dimension too large                        ' ,
     1     'User halt during matrix calculation               ' ,
     2     'User halt during CG tridiagonalization            ' ,
     3     'Too many tridiagonal matrix elements              ' ,
     4     'Bracketing did not find minimum within step bound ' ,
     5     'Number of points in spectrum is zero              ' ,
     6     'Field range of spectrum is zero                   ' /)
!
      character*32, dimension(0:NMNERR-1), save ::
     #  minerr =  (/ 'Bad input parameters            ' ,
     1               'Chi-squared convergence         ' ,
     2               'X vector convergence            ' ,
     3               'Chi-squared/X vector convergence' ,
     4               'Gradient convergence            ' ,
     5               'Max function evaluations reached' ,
     6               'Max iterations reached          ' ,
     7               'FTOL too small                  ' ,
     8               'XTOL too small                  ' ,
     9               'Zero Jacobian or GTOL too small ' ,
     A               'Terminated internally           ' ,
     1               'Single calculation completed    ' /)

      end module errmsg
