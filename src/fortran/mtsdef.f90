! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                    ====================
!                       module MTSDEF
!                    ====================
!
!    Defines positions of elements of mts array used in lbasix and setmts
!    routines
!
!          mts(1)  lemx    Maximum even L
!          mts(2)  lomx    Maximum odd L
!          mts(3)  kmn     Minimum K (K<0 => jK = -1)
!          mts(4   kmx     Maximum K
!          mts(5)  mmn     Minimum M (M<0 => jM = -1)
!          mts(6)  mmx     Maximum M
!          mts(7)  ipnmx   Maximum pI
!          mts(8)  ldelta  L increment
!          mts(9)  kdelta  K increment
!          mts(10) ipsi0   Tilt angle flag
!          mts(11) in2     Nuclear spin
!          mts(12) jkmn    Minimum jK
!          mts(13) jmmn    Minimum jM
!
!----------------------------------------------------------------------
!
      module mtsdef
      implicit none
!
      integer, parameter ::
     #           NLEMX=1,NLOMX=2,NKMN=3,NKMX=4,NMMN=5,NMMX=6,NIPNMX=7,
     #           NLDEL=8,NKDEL=9,NIPSI0=10,NIN2=11,NJKMN=12,NJMMN=13,
     #           NTRC=7
!
      end module mtsdef
