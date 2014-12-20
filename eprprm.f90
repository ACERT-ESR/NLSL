c Version 1.5.1 beta 2/3/96
c----------------------------------------------------------------------
c                    =========================
c                       common block EPRPRM
c                    =========================
c
c       Variable and common block definition for NLSL fitting of
c  experimental spectra using LBLL-family slow-motional calculations.
c
c
c IMPORTANT NOTES: 
c
c   The order in which the parameters appear in /eprprm/ is critical to 
c   the proper functioning of the NLSL program. If this common block is
c   to be changed, the following rules must be observed:
c
c   (1) To permit alias names to be used for the g,A, and diffusion
c       tensors, the parameters gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz,wxx,wyy,wzz
c       should appear contiguously in that order in the common block.
c
c   (2) The parameters igsph,iasph,irsph,iwsph should appear contiguously
c       in that order in the common block. (See function tcheck in file
c       sphten)
c
c   (3) The order in which the parameters appear in /eprprm/ must
c       be consistent with the names defined in the block data section
c       of file "ipfind.f"
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

c        double precision gxx,gyy,gzz,axx,ayy,azz,wxx,wyy,wzz,dx,dy,dz,
c     #         pml,pmxy,pmzz,djf,djfprp,oss,psi,ald,bed,gad,alm,bem,gam,
c     #         c20,c22,c40,c42,c44,gib0,gib2,lb,b0,gamman,cgtol,shiftr,
c     #         shifti,range,fldi,dfld,a0,g0,w0,expl,expkxy,expkzz,phase,
c     #         dc20
c
c        double precision faa,fgm,fam,fwm,fgd,fad,fwd,cpot,xlk
c
c        integer in2,ipdf,ist,ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
c     #          nort,igflg,iaflg,iwflg,irflg,nstep,ideriv,itype,ndim,
c     #          nfld,jkmn,jmmn,ipt,itm,itd,ipsi0,lband,kband,ldelta,
c     #          kdelta,lptmx,kptmx,neltot,nelv,nelre,nelim,ncgstp
c
        double precision, pointer, save :: phase,gib0,gib2,
     #         wxx,wyy,wzz,
     #         gxx,gyy,gzz,
     #         axx,ayy,azz,
     #         dx,dy,dz,
     #         pml,pmxy,pmzz,
     #         djf,djfprp,oss,
     #         psi,ald,bed,gad,alm,bem,gam,
     #         c20,c22,c40,c42,c44,lb,dc20,b0,gamman,
     #         cgtol,shiftr,shifti,range,fldi,dfld,
     #         a0,g0,w0,expl,expkxy,expkzz,faa(:),fgm(:),fwm(:),
     #         fam(:,:),fgd(:,:),fad(:,:),fwd(:,:),cpot(:,:),xlk(:,:)

c     #         a0,g0,w0,expl,expkxy,expkzz,faa(5),fgm(5),fwm(5),
c     #         fam(2,5),fgd(2,5),fad(2,5),fwd(2,5),cpot(5,5),xlk(5,5),

        integer, pointer, save :: in2,ipdf,ist,
     #         ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #         nort,nstep,nfld,ideriv,iwflg,igflg,iaflg,irflg,jkmn,jmmn,
     #         ndim,itype,ipt,itm,itd,ipsi0,lband,kband,ldelta,
     #         kdelta,lptmx,kptmx,neltot,nelv,nelre,nelim,ncgstp
c
c *** The following constants identify the position of several
c     important parameters within the fepr (and fparm) arrays.
c     THESE CONSTANTS *MUST* BE REDEFINED IF THE PARAMETER ORDER 
c     IN EPRPRM IS CHANGED!!!
c
c      integer IPHASE,IGIB0,IGIB2,IGXX,IGZZ,IAXX,IAZZ,IWXX,IWZZ,
c     #        IDX,IDZ,IPML,IPMXY,IPMZZ,IDJF,IDJFPRP,IOSS,IPSI,IALD,IBED,
c     #        IGAD,IALM,IGAM,IC20,IC44,ILB,IB0,IGAMAN,IML,ILEMX,IIPDF,
c     #        IIST,IIN2,IIGFLG,IIAFLG,IIWFLG,IIRFLG,INORT,INSTEP,ISHIFT,
c     #        ICGTOL,IFLDI,IDFLD,INFLD,IIDERV,INDIM,IDC20,IRANGE
c
      integer, parameter :: IPHASE=1,IGIB0=2,IGIB2=3,
     #         IWXX=4,IWZZ=6,IGXX=7,IGZZ=9,
     #         IAXX=10,IAZZ=12,IDX=13,
     #         IDZ=15,IPML=16,IPMXY=17,IPMZZ=18,
     #         IDJF=19,IDJFPRP=20,IOSS=21,IPSI=22,IALD=23,IBED=24,
     #         IGAD=25,IALM=26,IGAM=28,IC20=29,IC44=33,ILB=34,IDC20=35,
     #         IB0=36,IGAMAN=37,ICGTOL=38,ISHIFT=39,IRANGE=41,IFLDI=42,
     #         IDFLD=43,
     #         IIN2=1,IIPDF=2,IIST=3,IML=4,ILEMX=7,INORT=14,INSTEP=15,
     #         INFLD=16,IIDERV=17,IIWFLG=18,IIGFLG=19,IIAFLG=20,
     #         IIRFLG=21,INDIM=24

c     The following constants were absent from the original list.
c     They are now included for consistency.  
      integer, parameter :: IWYY=5,IGYY=8,IAYY=11,IDY=14,
     #         IBEM=27,IC22=30,IC40=31,IC42=32,ISHIFTI=40,
     #         IMXY=5,IMZZ=6,
     #         ILOMX=8,IKMN=9,IKMX=10,IMMN=11,IMMX=12,IIPNMX=13,
     #         IJKMN=22,IJMMN=23

c     In the original list, IGAMAN, ISHIFT, IIDERV had odd spellings.
c     The following extra constants conform to the typical pattern.
      integer, parameter :: IGAMMAN=IGAMAN,ISHIFTR=ISHIFT,
     #         IIDERIV=IIDERV
c
      integer, save :: current_site

      contains

      subroutine select_site(isite)
      implicit none
      integer :: isite
      current_site = isite
      phase => fparm(IPHASE,isite)
      gib0 => fparm(IGIB0,isite)
      gib2 => fparm(IGIB2,isite)
      wxx => fparm(IWXX,isite)
      wyy => fparm(IWYY,isite)
      wzz => fparm(IWZZ,isite)
      gxx => fparm(IGXX,isite)
      gyy => fparm(IGYY,isite)
      gzz => fparm(IGZZ,isite)
      axx => fparm(IAXX,isite)
      ayy => fparm(IAYY,isite)
      azz => fparm(IAZZ,isite)
      dx => fparm(IDX,isite)
      dy => fparm(IDY,isite)
      dz => fparm(IDZ,isite)
      pml => fparm(IPML,isite)
      pmxy => fparm(IPMXY,isite)
      pmzz => fparm(IPMZZ,isite)
      djf => fparm(IDJF,isite)
      djfprp => fparm(IDJFPRP,isite)
      oss => fparm(IOSS,isite)
      psi => fparm(IPSI,isite)
      ald => fparm(IALD,isite)
      bed => fparm(IBED,isite)
      gad => fparm(IGAD,isite)
      alm => fparm(IALM,isite)
      bem => fparm(IBEM,isite)
      gam => fparm(IGAM,isite)
      c20 => fparm(IC20,isite)
      c22 => fparm(IC22,isite)
      c40 => fparm(IC40,isite)
      c42 => fparm(IC42,isite)
      c44 => fparm(IC44,isite)
      lb => fparm(ILB,isite)
      dc20 => fparm(IDC20,isite)
      b0 => fparm(IB0,isite)
      gamman => fparm(IGAMMAN,isite)
      cgtol => fparm(ICGTOL,isite)
      shiftr => fparm(ISHIFTR,isite)
      shifti => fparm(ISHIFTI,isite)
      range => fparm(IRANGE,isite)
      fldi  => fparm(IFLDI,isite)
      dfld => fparm(IDFLD,isite)

      in2 => iparm(IIN2,isite)
      ipdf => iparm(IIPDF,isite)
      ist => iparm(IIST,isite)
      ml => iparm(IML,isite)
      mxy => iparm(IMXY,isite)
      mzz => iparm(IMZZ,isite)
      lemx => iparm(ILEMX,isite)
      lomx => iparm(ILOMX,isite)
      kmn => iparm(IKMN,isite)
      kmx => iparm(IKMX,isite)
      mmn => iparm(IMMN,isite)
      mmx => iparm(IMMX,isite)
      ipnmx => iparm(IIPNMX,isite)
      nort => iparm(INORT,isite)
      nstep => iparm(INSTEP,isite)
      nfld => iparm(INFLD,isite)
      ideriv => iparm(IIDERV,isite)
      iwflg => iparm(IIWFLG,isite)
      igflg => iparm(IIGFLG,isite)
      iaflg => iparm(IIAFLG,isite)
      irflg => iparm(IIRFLG,isite)
      jkmn => iparm(IJKMN,isite)
      jmmn => iparm(IJMMN,isite)
      ndim => iparm(INDIM,isite)

      end subroutine select_site

      end module eprprm
