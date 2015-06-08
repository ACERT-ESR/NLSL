c    Version 1.5  5/2/94
c----------------------------------------------------------------------
c                    =========================
c                         common SPECTR
c                    =========================
c
c   Storage for calculated spectrum and field array
c
c   Note: stddim.inc must be included before this file
c----------------------------------------------------------------------
      double precision spect
      common /spectr/ spect(2,MXPT)
