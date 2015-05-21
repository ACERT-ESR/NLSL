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
      implicit none
c
      double precision, target, save :: fepr(NFPRM)
      integer, target, save :: iepr(NIPRM)

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
      gib0 => fepr(IGIB0)
      gib2 => fepr(IGIB2)
      wxx => fepr(IWXX)
      wyy => fepr(IWYY)
      wzz => fepr(IWZZ)
      gxx => fepr(IGXX)
      gyy => fepr(IGYY)
      gzz => fepr(IGZZ)
      axx => fepr(IAXX)
      ayy => fepr(IAYY)
      azz => fepr(IAZZ)
      dx => fepr(IDX)
      dy => fepr(IDY)
      dz => fepr(IDZ)
      pml => fepr(IPML)
      pmxy => fepr(IPMXY)
      pmzz => fepr(IPMZZ)
      djf => fepr(IDJF)
      djfprp => fepr(IDJFPRP)
      oss => fepr(IOSS)
      psi => fepr(IPSI)
      ald => fepr(IALD)
      bed => fepr(IBED)
      gad => fepr(IGAD)
      alm => fepr(IALM)
      bem => fepr(IBEM)
      gam => fepr(IGAM)
      c20 => fepr(IC20)
      c22 => fepr(IC22)
      c40 => fepr(IC40)
      c42 => fepr(IC42)
      c44 => fepr(IC44)
      lb => fepr(ILB)
      dc20 => fepr(IDC20)
      b0 => fepr(IB0)
      gamman => fepr(IGAMMAN)
      cgtol => fepr(ICGTOL)
      shiftr => fepr(ISHIFTR)
      shifti => fepr(ISHIFTI)
      range => fepr(IRANGE)
      fldi  => fepr(IFLDI)
      dfld => fepr(IDFLD)

      in2 => iepr(IIN2)
      ipdf => iepr(IIPDF)
      ist => iepr(IIST)
      ml => iepr(IML)
      mxy => iepr(IMXY)
      mzz => iepr(IMZZ)
      lemx => iepr(ILEMX)
      lomx => iepr(ILOMX)
      kmn => iepr(IKMN)
      kmx => iepr(IKMX)
      mmn => iepr(IMMN)
      mmx => iepr(IMMX)
      ipnmx => iepr(IIPNMX)
      nort => iepr(INORT)
      nstep => iepr(INSTEP)
      nfld => iepr(INFLD)
      ideriv => iepr(IIDERV)
      iwflg => iepr(IIWFLG)
      igflg => iepr(IIGFLG)
      iaflg => iepr(IIAFLG)
      irflg => iepr(IIRFLG)
      jkmn => iepr(IJKMN)
      jmmn => iepr(IJMMN)
      ndim => iepr(INDIM)

      end subroutine prm_ptr_init
	  
      end module eprprm
