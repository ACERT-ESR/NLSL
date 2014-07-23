c----------------------------------------------------------------------
c                    =========================
c                       include file LMTXT
c                    =========================
c  Contains definitions of strings identifying return codes from
c  the subroutine LMDER in MINPACK
c----------------------------------------------------------------------
      character*32 minerr(0:10)
c
       data minerr / 'Bad input parameters',
     1               'Chi-squared convergence',
     2               'X vector convergence',
     3               'Chi-squared/X vector convergence',
     4               'Gradient convergence',
     5               'Max function evaluations reached',
     5               'Max iterations reached',
     6               'FTOL too small',
     7               'XTOL too small' ,
     8               'Zero Jacobian or GTOL too small',
     9               'Terminated internally'/
