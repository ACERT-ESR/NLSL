c----------------------------------------------------------------------
c                    ====================
c                    include file MTSDEF
c                    ====================
c----------------------------------------------------------------------
c
c    Defines positions of elements of mts array used in lbasix and  setmts
c    routines
c
c          mts(1)  lemx    Maximum even L
c          mts(2)  lomx    Maximum odd L
c          mts(3)  kmn     Minimum K (K<0 => jK = -1)
c          mts(4   kmx     Maximum K
c          mts(5)  mmn     Minimum M (M<0 => jM = -1)
c          mts(6)  mmx     Maximum M
c          mts(7)  ipnmx   Maximum pI
c          mts(8)  ldelta  L increment
c          mts(9)  kdelta  K increment
c          mts(10) ipsi0   Tilt angle flag
c          mts(11) in2     Nuclear spin
c          mts(12) jkmn    Minimum jK
c          mts(13) jmmn    Minimum jM

      integer NLEMX,NLOMX,NMMX,NKMX,NIPNMX,NLDEL,NKDEL,NIPSI0,NIN2,
     #        NKMN,NMMN,NJKMN,NJMMN,NTRC   
      parameter (NLEMX=1,NLOMX=2,NKMN=3,NKMX=4,NMMN=5,NMMX=6,NIPNMX=7,
     #           NLDEL=8,NKDEL=9,NIPSI0=10,NIN2=11,NJKMN=12,NJMMN=13,
     #           NTRC=7)
