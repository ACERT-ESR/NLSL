c----------------------------------------------------------------------
c                    =========================
c                      include file PARMEQU
c                    =========================
c   Defines equivalence of arrays fepr and iepr with common block 
c   /eprprm/ contained in sims.inc
c
c   fepr    is a vector corresponding to all the user-settable floating
c           point EPRLL parameters
c   iepr    is a vector corresponding to all the user-settable integer
c           parameters
c
c   NOTE: both nlsdim.inc and eprprm.inc must be included before
c   this file
c----------------------------------------------------------------------
      double precision fepr(NFPRM)
      integer iepr(NIPRM)
      equivalence (fepr(1),gib),(iepr(1),in2)
