c Version 1.6  8/12/94
c----------------------------------------------------------------------
c Physical constants (in cgs/electromagnetic units)
c    ge:     free-electron g-value
c    betae:  Bohr magneton
c    hbar:   Planck's constant
c    gammae: electronic gyromagnetic ratio
c    kb:     Boltzmann's constant
c----------------------------------------------------------------------
        double precision GE,BETAE,GAMMAE,HBAR,KB
        parameter( GE=2.002322D0, 
     #             BETAE=9.2731D-21, 
     #             HBAR=1.05443D-27,
     #             GAMMAE=GE*BETAE/HBAR,
     #             KB=1.380622D-16 )
