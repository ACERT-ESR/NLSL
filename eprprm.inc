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
        double precision gxx,gyy,gzz,axx,ayy,azz,wxx,wyy,wzz,dx,dy,dz,
     #         pml,pmxy,pmzz,djf,djfprp,oss,psi,ald,bed,gad,alm,bem,gam,
     #         c20,c22,c40,c42,c44,gib0,gib2,lb,b0,gamman,cgtol,shiftr,
     #         shifti,range,fldi,dfld,a0,g0,w0,expl,expkxy,expkzz,phase,
     #         dc20
c
        double precision faa,fgm,fam,fwm,fgd,fad,fwd,cpot,xlk
c
        integer in2,ipdf,ist,ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #          nort,igflg,iaflg,iwflg,irflg,nstep,ideriv,itype,ndim,
     #          nfld,jkmn,jmmn,ipt,itm,itd,ipsi0,lband,kband,ldelta,
     #          kdelta,lptmx,kptmx,neltot,nelv,nelre,nelim,ncgstp
c
        common /eprprm/ phase,gib0,gib2,
     #         wxx,wyy,wzz,
     #         gxx,gyy,gzz,
     #         axx,ayy,azz,
     #         dx,dy,dz,
     #         pml,pmxy,pmzz,
     #         djf,djfprp,oss,
     #         psi,ald,bed,gad,alm,bem,gam,
     #         c20,c22,c40,c42,c44,lb,dc20,b0,gamman,
     #         cgtol,shiftr,shifti,range,fldi,dfld,
     #         a0,g0,w0,expl,expkxy,expkzz,faa(5),fgm(5),fwm(5),
     #         fam(2,5),fgd(2,5),fad(2,5),fwd(2,5),cpot(5,5),xlk(5,5),
     #         in2,ipdf,ist,ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #         nort,nstep,nfld,ideriv,iwflg,igflg,iaflg,irflg,jkmn,jmmn,
     #         ndim,itype,ipt,itm,itd,ipsi0,lband,kband,ldelta,
     #         kdelta,lptmx,kptmx,neltot,nelv,nelre,nelim,ncgstp
c
c *** The following constants identify the position of several
c     important parameters within the fepr (and fparm) arrays.
c     THESE CONSTANTS *MUST* BE REDEFINED IF THE PARAMETER ORDER 
c     IN EPRPRM IS CHANGED!!!
c
      integer IPHASE,IGIB0,IGIB2,IGXX,IGZZ,IAXX,IAZZ,IWXX,IWZZ,
     #        IDX,IDZ,IPML,IPMXY,IPMZZ,IDJF,IDJFPRP,IOSS,IPSI,IALD,IBED,
     #        IGAD,IALM,IGAM,IC20,IC44,ILB,IB0,IGAMAN,IML,ILEMX,IIPDF,
     #        IIST,IIN2,IIGFLG,IIAFLG,IIWFLG,IIRFLG,INORT,INSTEP,ISHIFT,
     #        ICGTOL,IFLDI,IDFLD,INFLD,IIDERV,INDIM,IDC20,IRANGE
c
      parameter(IPHASE=1,IGIB0=2,IGIB2=3,IWXX=4,IWZZ=6,IGXX=7,IGZZ=9,
     #         IAXX=10,IAZZ=12,IDX=13,IDZ=15,IPML=16,IPMXY=17,IPMZZ=18,
     #         IDJF=19,IDJFPRP=20,IOSS=21,IPSI=22,IALD=23,IBED=24,
     #         IGAD=25,IALM=26,IGAM=28,IC20=29,IC44=33,ILB=34,IDC20=35,
     #         IB0=36,IGAMAN=37,ICGTOL=38,ISHIFT=39,IRANGE=41,IFLDI=42,
     #         IDFLD=43,
     #         IIN2=1,IIPDF=2,IIST=3,IML=4,ILEMX=7,INORT=14,INSTEP=15,
     #         INFLD=16,IIDERV=17,IIWFLG=18,IIGFLG=19,IIAFLG=20,
     #         IIRFLG=21,INDIM=24)
c





