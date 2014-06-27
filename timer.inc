c Version 1.4 10/10/94
c----------------------------------------------------------------------
c                    =========================
c                         common TIMER
c                    =========================
c
c     Defines array used by UNIX system function dtime
c----------------------------------------------------------------------
      real tarray(8)
      common /timer/ tarray
