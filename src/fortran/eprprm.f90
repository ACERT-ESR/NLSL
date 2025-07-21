c NLSL Version 1.9.0 beta 2/6/15
c----------------------------------------------------------------------
c                    =========================
c                          module EPRPRM
c                    =========================
c
c  Definitions of variables and pointer-aliases for NLSL fitting of
c  experimental spectra using LBLL-family slow-motional calculations.
c
c  This module supersedes various F77 common blocks and equivalences,
c  in particular those that were in the files eprprm.inc and prmeqv.inc.
c
c
c IMPORTANT NOTES (mostly taken from version 1.5.1 beta, 2/3/96): 
c
c   The order in which the parameters appear in fparm and iparm is critical 
c   to the proper functioning of the NLSL program. If the parameter order is
c   to be changed, the following rules must be observed:
c
c   (1) To permit alias names to be used for the g, A, and diffusion
c       tensors, the parameters gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz,wxx,wyy,wzz
c       should appear contiguously in that order in fparm.
c
c   (2) The parameters igsph,iasph,irsph,iwsph should appear contiguously
c       in that order in iparm. (See function tcheck in file tensym.f90.)
c
c   (3) The order in which the parameters appear in module eprprm must
c       be consistent with the order of names defined in module lpnam.
c       This enables function ipfind to work properly, among other things.
c
c   (4) The residence times tl,tkxy,tkzz have been replaced with the
c       parameters pml, pmxy, pmzz. This is to facilitate fitting using
c       non-Brownian models where it is desirable to vary the residence
c       time and the diffusion rate constant together (i.e. keeping the
c       product, which is related to the rms jump angle, constant).
c
c   (5) The non-Brownian model flags have been changed to ml, mzz, and
c       mxy for perpendicular, parallel, and internal motions, respectively.
c       The flags may be set to 0, 1, and 2 for Brownian, free, and
c       jump diffusion, respectively.
c
c   (6) The old "wint" parameters have been changed to "gib", and
c       now the additional width refers to Gaussian inhomogeneous width
c       rather than an additional Lorentzian width
c
c   (7) A new spectral linewidth broadening parameter, "lb", has been
c       added
c
c***********************************************************************
c
      module eprprm
      use parcom
      use nlsdim
      implicit none
c
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
c
c     The following are not declared with the pointer attribute,
c     as they are not special names that are used to point into
c     another array, like the variables above
c
      double precision, save ::
     #         a0,g0,w0,expl,expkxy,expkzz,faa(5),fgm(5),fwm(5),
     #         fam(2,5),fgd(2,5),fad(2,5),fwd(2,5),cpot(5,5),xlk(5,5)

      integer, save ::
     #         itype,ipt,itm,itd,ipsi0,lband,kband,ldelta,
     #         kdelta,lptmx,kptmx,neltot,nelv,nelre,nelim,ncgstp

c
c *** The following constants identify the position of many of the
c     parameters within the fepr/iepr (and fparm/iparm) arrays.
c     THESE CONSTANTS *MUST* BE REDEFINED IF THE PARAMETER ORDER 
c     IN EPRPRM IS CHANGED!!!
c
      integer, parameter :: IPHASE=1,IGIB0=2,IGIB2=3,
     #         IWXX=4,IWZZ=6,IGXX=7,IGZZ=9,
     #         IAXX=10,IAZZ=12,IDX=13,
     #         IDZ=15,IPML=16,IPMXY=17,IPMZZ=18,
     #         IDJF=19,IDJFPRP=20,IOSS=21,IPSI=22,IALD=23,IBED=24,
     #         IGAD=25,IALM=26,IGAM=28,IC20=29,IC44=33,ILB=34,IDC20=35,
     #         IB0=36,IGAMAN=37,ICGTOL=38,ISHIFT=39,IRANGE=41,IFLDI=42,
     #         IDFLD=43
c
      integer, parameter :: IIN2=1,IIPDF=2,IIST=3,IML=4,
     #         ILEMX=7,INORT=14,INSTEP=15,
     #         INFLD=16,IIDERV=17,IIWFLG=18,IIGFLG=19,IIAFLG=20,
     #         IIRFLG=21,INDIM=24

c     The following constants were absent from the original lists.
c     They are now included for consistency.  
      integer, parameter :: IWYY=5,IGYY=8,IAYY=11,IDY=14,
     #         IBEM=27,IC22=30,IC40=31,IC42=32,ISHIFTI=40,
     #         IMXY=5,IMZZ=6,
     #         ILOMX=8,IKMN=9,IKMX=10,IMMN=11,IMMX=12,IIPNMX=13,
     #         IJKMN=22,IJMMN=23

c     In the original lists, IGAMAN, ISHIFT, IIDERV had odd spellings.
c     The following extra constants conform to the typical pattern.
      integer, parameter :: IGAMMAN=IGAMAN,ISHIFTR=ISHIFT,
     #         IIDERIV=IIDERV
c
c      integer, save :: current_site

      contains

c      subroutine select_site(isite)
c      implicit none
c
c      Possible future subroutine...
c      Instead of copying a column of fparm into fepr, or copying a
c      column of iparm into iepr, it may be possible to make fepr
c      and iepr into pointers to the heads of the relevant columns.
c      Currently, though, fepr and iepr are treated as temp arrays
c      whose values can change, so making them into pointers would
c      have undesirable effects on global arrays fparm and iparm. 
c
c      integer :: isite
c      current_site = isite
c              ------------------------------------------------
c              After call momdls, all eprprm parameters point
c              to the current site via call select_site(isi)
c              ------------------------------------------------
c      fepr => fparm(:,isite)
c      iepr => iparm(:,isite)
c      end subroutine select_site

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
      fepr_name(IDX) = "dx"
      dy => fepr(IDY)
      fepr_name(IDY) = "dy"
      dz => fepr(IDZ)
      fepr_name(IDZ) = "dz"
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
      fepr_name(IALD) = "ald"
      bed => fepr(IBED)
      fepr_name(IBED) = "bed"
      gad => fepr(IGAD)
      fepr_name(IGAD) = "gad"
      alm => fepr(IALM)
      fepr_name(IALM) = "alm"
      bem => fepr(IBEM)
      fepr_name(IBEM) = "bem"
      gam => fepr(IGAM)
      fepr_name(IGAM) = "gam"
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
      fepr_name(IFLDI) = "fldi"
      dfld => fepr(IDFLD)
      fepr_name(IDFLD) = "dfld"

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
      iepr_name(INFLD) = "nfld"
      ideriv => iepr(IIDERV)
      iepr_name(IIDERV) = "ideriv"
      iwflg => iepr(IIWFLG)
      iepr_name(IIWFLG) = "iwflg"
      igflg => iepr(IIGFLG)
      iepr_name(IIGFLG) = "igflg"
      iaflg => iepr(IIAFLG)
      iepr_name(IIAFLG) = "iaflg"
      irflg => iepr(IIRFLG)
      iepr_name(IIRFLG) = "irflg"
      jkmn => iepr(IJKMN)
      iepr_name(IJKMN) = "jkmn"
      jmmn => iepr(IJMMN)
      iepr_name(IJMMN) = "jmmn"
      ndim => iepr(INDIM)
      iepr_name(INDIM) = "ndim"

      end subroutine prm_ptr_init
	  
      end module eprprm
