c Version 1.6  8/12/94
**********************************************************************
c       
c                       ==================
c                       SUBROUTINE : WRDAT
c                       ==================
c
c       This writes the parameters contained in the common block
c       defined in eprdat.inc into a file in a form that can be read 
c       by rddat.f. An identification string for the current program
c       version is also written with the parameter file.
c
c       Version 1.6 ald,gad (diffusion tilt) angles added
c                   alm,bem,gam (magnetic tilt) angles added
c                   gamman added
c
c       Version 1.5 has wint0,wint2 revised to gib0, gib2
c                   wxx,wyy,wzz have been added
c 
c       Version 1.4 included parameters needed for the MOMD calculation.
c
c       Notes:
c               Care should be taken to ensure that the common 
c               block definition in the include file eprdat.inc
c               is consistent with the information stored in the
c               disk files by this routine to avoid corrupting 
c               the parameter file.
c
c       written by DJS 15-NOV-87
c       modified by DEB 13-MAR-92 to use longer filenames
c                       22-OCT-92 to write identification string
c
c
c       Includes:
c               stdio.inc
c               eprdat.inc
c		version.inc
c
c       Uses:
c
c**********************************************************************
c
      subroutine wrdat(prname,namlth)
      implicit none
c
      include 'stdio.inc'
      include 'eprdat.inc'
c
      character prname*30,eprID*5
      integer i,j,namlth
c
      include 'version.inc'
c
      data eprID / 'EPRLL' /
c
c######################################################################
c
c----------------------------------------------------------------------
c     open parameter file
c----------------------------------------------------------------------
c
      open (unit=ludisk,file=prname(:namlth),status='unknown',
     #     access='direct',form='unformatted',recl=1024)
c
c     -----------------------------------------------------------------
c     New floating parameters in Version 1.6 are:
c       alm,bem,gam,ald,gad,gamman
c     Also note that fam,fgd,fad,fwd are now complex-valued and
c     the order of parameters is somewhat rearranged from version 1.5
c     and the information has been redistributed between the two
c     records in the file.
c     -----------------------------------------------------------------
c
      write (ludisk,rec=1) eprID,version,
     #     gxx,gyy,gzz,axx,ayy,azz,gamman,pl,pkxy,pkzz,b0,g0,a0,
     #     dx,dy,dz,tl,tkxy,tkzz,djf,oss,psi,ald,bed,gad,alm,bem,gam,
     #     cgtol,shiftr,shifti,fieldi,fieldf,djfprp,gib0,gib2,btol,
     #     (faa(i),i=1,5),(fgm(i),i=1,5),(fam(1,i),fam(2,i),i=1,5),
     #     (fgd(1,i),fgd(2,i),i=1,5),(fad(1,i),fad(2,i),i=1,5)
c
c     New parameters in version 1.6 are kmn and mmn
c
      write (ludisk,rec=2)
     #     (cxp(i),i=1,6),((cpot(i,j),j=1,5),i=1,5),
     #     ((xlk(i,j),j=1,5),i=1,5),wxx,wyy,wzz,w0,
     #     (fwm(i),i=1,5),(fwd(1,i),fwd(2,i),i=1,5),
     #     mpl,mpkxy,mpkzz,in2,ipdf,ist,ipt,lptmx,kptmx,
     #     lband,kband,itd,ipsi0,lemx,lomx,kmx,kmn,mmx,mmn,ipnmx,
     #     ldelta,kdelta,jkmn,jmmn,nstep,ndim,neltot,itype,nfield,
     #     ((ixp(i,j),j=1,2),i=1,6),nelre,nelim,nort,ideriv
c
c----------------------------------------------------------------------
c     close parameter file
c----------------------------------------------------------------------
c
      close (unit=ludisk)
c
c----------------------------------------------------------------------
c     return to caller
c----------------------------------------------------------------------
c
      return
      end
