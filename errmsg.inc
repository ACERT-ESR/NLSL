c Version 1.5.1 beta 03/07/96
c----------------------------------------------------------------------
c                    =========================
c                         common ERRMSG
c                    =========================
c
c Contains text for error messages and definitions of symbols used
c to define error message numbers. Errors after ZEROB0 are defined
c as fatal, and will cause the calculation to halt.
c
c Error message texts are set in a block data section in nlstxt.f
c Definitions of the symbols are as follows:
c
c ZEROIN2     Zero nuclear spin assumed
c POSB0       Positive B0 assumed
c BADHIFLD    High-field approx. may not apply
c BADJMP      Nonzero potential with jump/free: ipdf set to 0
c BADAV       Zero potential with aniso. viscos.: ipdf set to 0
c BADDJ       Discrete jumps with aniso. viscos.: ipdf set 0
c LEMXHI      lemx must be 48 or less with potenti
c BSSADJ      Basis set adjusted to (iii,iii,iii,iii,iii,iii,ii)
c NSTEPHI     NSTEP too large: adjusted to iiiii
c CGTOLHI     CGTOL too small: adjusted to xxxxxxxxx
c SHIFTHI     SHIFTR too small: adjusted to xxxxxxxxx
c NOCONVRG    CG did not converge after iiiii steps
c NOALPHAM    Angle alm not needed with axial g-matrix
c NOALPHAD    Angle ald not needed with axial diffusion tensor
c ZEROB0      Zero B0
c ZEROG       Zero g-tensor element
c TDGERR      Error in tridiagonal matrix storage
c MTXBIG      Too many matrix elements
c DIMBIG      Matrix dimension too large
c MTXHLT      User halt during matrix calculation
c CGHLT       User halt during CG tridiagonalization
c TDGBIG      Too many tridiagonal matrix elements
c NOMIN       Bracketing routine (search command) did not find a minimum
c ZERONFLD    NFIELD is zero (no points in spectrum)
c ZRORANGE    RANGE is zero 
c----------------------------------------------------------------------

      integer FATAL,NEPERR
      parameter (FATAL=16,NEPERR=26)
c
       integer ZEROIN2,POSB0,BADHIFLD,BADJMP,BADAV,BADDJ,LEMXHI,BSSADJ,
     #         NSTEPHI,CGTOLHI,SHIFTHI,NOCONVRG,NOALPHAM,NOALPHAD,
     #         BADBESS,ZEROB0,ZEROG,TDGERR,MTXBIG,DIMBIG,MTXHLT,CGHLT,
     #         TDGBIG,NOMIN,ZERONFLD,ZRORANGE
c
       parameter (ZEROIN2=1,POSB0=2,BADHIFLD=3,BADJMP=4,BADAV=5,BADDJ=6,
     #            LEMXHI=7,BSSADJ=8,NSTEPHI=9,CGTOLHI=10,SHIFTHI=11,
     #            NOCONVRG=12,NOALPHAM=13,NOALPHAD=14,BADBESS=15,
     #            ZEROB0=16,ZEROG=17,TDGERR=18,MTXBIG=19,
     #            DIMBIG=20,MTXHLT=21,CGHLT=22,TDGBIG=23,NOMIN=24,
     #            ZERONFLD=25,ZRORANGE=26)

      character*50 eprerr
      character*32 minerr
      common /errmsg/ eprerr(NEPERR),minerr(0:11)



