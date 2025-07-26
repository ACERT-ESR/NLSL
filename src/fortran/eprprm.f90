! NLSL Version 1.9.0 beta 2/6/15
!----------------------------------------------------------------------
!                    =========================
!                          module EPRPRM
!                    =========================
!
!  Definitions of variables and pointer-aliases for NLSL fitting of
!  experimental spectra using LBLL-family slow-motional calculations.
!
!  This module supersedes various F77 common blocks and equivalences,
!  in particular those that were in the files eprprm.inc and prmeqv.inc.
!
!
! IMPORTANT NOTES (mostly taken from version 1.5.1 beta, 2/3/96): 
!
!   The order in which the parameters appear in fparm and iparm is critical 
!   to the proper functioning of the NLSL program. If the parameter order is
!   to be changed, the following rules must be observed:
!
!   (1) To permit alias names to be used for the g, A, and diffusion
!       tensors, the parameters gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz,wxx,wyy,wzz
!       should appear contiguously in that order in fparm.
!
!   (2) The parameters igsph,iasph,irsph,iwsph should appear contiguously
!       in that order in iparm. (See function tcheck in file tensym.f90.)
!
!   (3) The order in which the parameters appear in module eprprm must
!       be consistent with the order of names defined in module lpnam.
!       This enables function ipfind to work properly, among other things.
!
!   (4) The residence times tl,tkxy,tkzz have been replaced with the
!       parameters pml, pmxy, pmzz. This is to facilitate fitting using
!       non-Brownian models where it is desirable to vary the residence
!       time and the diffusion rate constant together (i.e. keeping the
!       product, which is related to the rms jump angle, constant).
!
!   (5) The non-Brownian model flags have been changed to ml, mzz, and
!       mxy for perpendicular, parallel, and internal motions, respectively.
!       The flags may be set to 0, 1, and 2 for Brownian, free, and
!       jump diffusion, respectively.
!
!   (6) The old "wint" parameters have been changed to "gib", and
!       now the additional width refers to Gaussian inhomogeneous width
!       rather than an additional Lorentzian width
!
!   (7) A new spectral linewidth broadening parameter, "lb", has been
!       added
!
!***********************************************************************
!
      module eprprm
      use parcom
      use nlsdim
      implicit none
!
      double precision, target, save :: fepr(NFPRM)
      character*10, dimension(NFPRM), save:: fepr_name
      integer, target, save :: iepr(NIPRM)
      character*10, dimension(NIPRM), save:: iepr_name

      double precision, pointer, save :: phase,gib0,gib2,
     #         wxx,wyy,wzz,
     #         gxx,gyy,gzz,
     #         axx,ayy,azz,
     #         dx,dy,dz,
     #         pml,pmxy,pmzz,
     #         djf,djfprp,oss,
     #         psi,ald,bed,gad,alm,bem,gam,
     #         c20,c22,c40,c42,c44,lb,dc20,b0,gamman,
     #         cgtol,shiftr,shifti,range,fldi,dfld

      integer, pointer, save :: in2,ipdf,ist,
     #         ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #         nort,nstep,nfld,ideriv,iwflg,igflg,iaflg,irflg,jkmn,jmmn,
     #         ndim
!
!     The following are not declared with the pointer attribute,
!     as they are not special names that are used to point into
!     another array, like the variables above
!
      double precision, save ::
     #         a0,g0,w0,expl,expkxy,expkzz,faa(5),fgm(5),fwm(5),
     #         fam(2,5),fgd(2,5),fad(2,5),fwd(2,5),cpot(5,5),xlk(5,5)

      integer, save ::
     #         itype,ipt,itm,itd,ipsi0,lband,kband,ldelta,
     #         kdelta,lptmx,kptmx,neltot,nelv,nelre,nelim,ncgstp

!
! *** The following constants identify the position of many of the
!     parameters within the fepr/iepr (and fparm/iparm) arrays.
!     THESE CONSTANTS *MUST* BE REDEFINED IF THE PARAMETER ORDER 
!     IN EPRPRM IS CHANGED!!!
!
      integer, parameter :: IPHASE=1,IGIB0=2,IGIB2=3,
     #         IWXX=4,IWZZ=6,IGXX=7,IGZZ=9,
     #         IAXX=10,IAZZ=12,IDX=13,
     #         IDZ=15,IPML=16,IPMXY=17,IPMZZ=18,
     #         IDJF=19,IDJFPRP=20,IOSS=21,IPSI=22,IALD=23,IBED=24,
     #         IGAD=25,IALM=26,IGAM=28,IC20=29,IC44=33,ILB=34,IDC20=35,
     #         IB0=36,IGAMAN=37,ICGTOL=38,ISHIFT=39,IRANGE=41,IFLDI=42,
     #         IDFLD=43
!
      integer, parameter :: IIN2=1,IIPDF=2,IIST=3,IML=4,
     #         ILEMX=7,INORT=14,INSTEP=15,
     #         INFLD=16,IIDERV=17,IIWFLG=18,IIGFLG=19,IIAFLG=20,
     #         IIRFLG=21,INDIM=24

!     The following constants were absent from the original lists.
!     They are now included for consistency.  
      integer, parameter :: IWYY=5,IGYY=8,IAYY=11,IDY=14,
     #         IBEM=27,IC22=30,IC40=31,IC42=32,ISHIFTI=40,
     #         IMXY=5,IMZZ=6,
     #         ILOMX=8,IKMN=9,IKMX=10,IMMN=11,IMMX=12,IIPNMX=13,
     #         IJKMN=22,IJMMN=23

!     In the original lists, IGAMAN, ISHIFT, IIDERV had odd spellings.
!     The following extra constants conform to the typical pattern.
      integer, parameter :: IGAMMAN=IGAMAN,ISHIFTR=ISHIFT,
     #         IIDERIV=IIDERV
!
!      integer, save :: current_site

      contains

      subroutine prm_ptr_init
      implicit none

      phase => fepr(IPHASE)
      fepr_name(IPHASE) = "phase"
      gib0 => fepr(IGIB0)
      fepr_name(IGIB0) = "gib0"
      gib2 => fepr(IGIB2)
      fepr_name(IGIB2) = "gib2"
      wxx => fepr(IWXX)
      fepr_name(IWXX) = "wxx"
      wyy => fepr(IWYY)
      fepr_name(IWYY) = "wyy"
      wzz => fepr(IWZZ)
      fepr_name(IWZZ) = "wzz"
      gxx => fepr(IGXX)
      fepr_name(IGXX) = "gxx"
      gyy => fepr(IGYY)
      fepr_name(IGYY) = "gyy"
      gzz => fepr(IGZZ)
      fepr_name(IGZZ) = "gzz"
      axx => fepr(IAXX)
      fepr_name(IAXX) = "axx"
      ayy => fepr(IAYY)
      fepr_name(IAYY) = "ayy"
      azz => fepr(IAZZ)
      fepr_name(IAZZ) = "azz"
      dx => fepr(IDX)
      fepr_name(IDX) = "rx"
      dy => fepr(IDY)
      fepr_name(IDY) = "ry"
      dz => fepr(IDZ)
      fepr_name(IDZ) = "rz"
      pml => fepr(IPML)
      fepr_name(IPML) = "pml"
      pmxy => fepr(IPMXY)
      fepr_name(IPMXY) = "pmxy"
      pmzz => fepr(IPMZZ)
      fepr_name(IPMZZ) = "pmzz"
      djf => fepr(IDJF)
      fepr_name(IDJF) = "djf"
      djfprp => fepr(IDJFPRP)
      fepr_name(IDJFPRP) = "djfprp"
      oss => fepr(IOSS)
      fepr_name(IOSS) = "oss"
      psi => fepr(IPSI)
      fepr_name(IPSI) = "psi"
      ald => fepr(IALD)
      fepr_name(IALD) = "alphad"
      bed => fepr(IBED)
      fepr_name(IBED) = "betad"
      gad => fepr(IGAD)
      fepr_name(IGAD) = "gammad"
      alm => fepr(IALM)
      fepr_name(IALM) = "alpham"
      bem => fepr(IBEM)
      fepr_name(IBEM) = "betam"
      gam => fepr(IGAM)
      fepr_name(IGAM) = "gammam"
      c20 => fepr(IC20)
      fepr_name(IC20) = "c20"
      c22 => fepr(IC22)
      fepr_name(IC22) = "c22"
      c40 => fepr(IC40)
      fepr_name(IC40) = "c40"
      c42 => fepr(IC42)
      fepr_name(IC42) = "c42"
      c44 => fepr(IC44)
      fepr_name(IC44) = "c44"
      lb => fepr(ILB)
      fepr_name(ILB) = "lb"
      dc20 => fepr(IDC20)
      fepr_name(IDC20) = "dc20"
      b0 => fepr(IB0)
      fepr_name(IB0) = "b0"
      gamman => fepr(IGAMMAN)
      fepr_name(IGAMMAN) = "gamman"
      cgtol => fepr(ICGTOL)
      fepr_name(ICGTOL) = "cgtol"
      shiftr => fepr(ISHIFTR)
      fepr_name(ISHIFTR) = "shiftr"
      shifti => fepr(ISHIFTI)
      fepr_name(ISHIFTI) = "shifti"
      range => fepr(IRANGE)
      fepr_name(IRANGE) = "range"
      fldi  => fepr(IFLDI)
      fepr_name(IFLDI) = ""
      dfld => fepr(IDFLD)
      fepr_name(IDFLD) = ""

      in2 => iepr(IIN2)
      iepr_name(IIN2) = "in2"
      ipdf => iepr(IIPDF)
      iepr_name(IIPDF) = "ipdf"
      ist => iepr(IIST)
      iepr_name(IIST) = "ist"
      ml => iepr(IML)
      iepr_name(IML) = "ml"
      mxy => iepr(IMXY)
      iepr_name(IMXY) = "mxy"
      mzz => iepr(IMZZ)
      iepr_name(IMZZ) = "mzz"
      lemx => iepr(ILEMX)
      iepr_name(ILEMX) = "lemx"
      lomx => iepr(ILOMX)
      iepr_name(ILOMX) = "lomx"
      kmn => iepr(IKMN)
      iepr_name(IKMN) = "kmn"
      kmx => iepr(IKMX)
      iepr_name(IKMX) = "kmx"
      mmn => iepr(IMMN)
      iepr_name(IMMN) = "mmn"
      mmx => iepr(IMMX)
      iepr_name(IMMX) = "mmx"
      ipnmx => iepr(IIPNMX)
      iepr_name(IIPNMX) = "ipnmx"
      nort => iepr(INORT)
      iepr_name(INORT) = "nort"
      nstep => iepr(INSTEP)
      iepr_name(INSTEP) = "nstep"
      nfld => iepr(INFLD)
      iepr_name(INFLD) = "nfield"
      ideriv => iepr(IIDERV)
      iepr_name(IIDERV) = "ideriv"
      iwflg => iepr(IIWFLG)
      iepr_name(IIWFLG) = ""
      igflg => iepr(IIGFLG)
      iepr_name(IIGFLG) = ""
      iaflg => iepr(IIAFLG)
      iepr_name(IIAFLG) = ""
      irflg => iepr(IIRFLG)
      iepr_name(IIRFLG) = ""
      jkmn => iepr(IJKMN)
      iepr_name(IJKMN) = ""
      jmmn => iepr(IJMMN)
      iepr_name(IJMMN) = ""
      ndim => iepr(INDIM)
      iepr_name(INDIM) = ""

      end subroutine prm_ptr_init
	  
      end module eprprm
