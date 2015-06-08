c Version 1.6  8/12/94
c***********************************************************************
c
c	Variable and common block definition for EPRLL/EPRCG/EPRBL series
c
c	Notes:
c		Any changes made to this file should also be made 
c		to the subroutines which read and write	this 
c		information to disk files.  See	the subroutines	
c		wrdat.f	and rddat.f.
c
c	written	by DJS 29-AUG-87
c	modified by DJS	24-OCT-87 to allow use with eprcgl (itype,cgtol,
c				  shiftr and shifti added)
c	modified by DJS	15-NOV-87 to allow use with eprbl (nfield,
c				  fieldi and fieldf added)
c
c       modified by DEB The following parameters have been added:
c                                 gib0, gib2     (Gaussian line broadening)
c                                 nelv           (number of vector elements)
c                                 ideriv         (derivative mode flag)
c                                 igsph,iasph    (spherical g,A tensor flags)
c
c                        23-Oct-92 Dx and Dy replace Dxy
c                        2-May-93 btol added: pruning tolerance
c                        2-May-94 wxx,wyy,wzz (Lorentzian widths) added
c                        12-Aug-94 alm,bem,gam (magnetic tilt) added
c                                  ald,gad (diffusion tilt) added
c                                  gamman added
c                                  truncation indices kmn, mmn 
c                                  and symmetry flags jkmn, jmmn added
c
c***********************************************************************
c
	double precision gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz,wxx,wyy,wzz,
     #			 tl,tkxy,tkzz,djf,djfprp,oss,psi,ald,bed,gad,
     #                   alm,bem,gam,c20,c22,c40,c42,c44,bshift,gib0,
     #			 gib2,b0,gamman,fieldi,fieldf,cgtol,shiftr,
     #                   shifti,a0,g0,w0,pl,pkxy,pkzz,btol

	double precision faa,fgm,fam,fgd,fad,fwm,fwd,cxp,cpot,xlk
c
	integer	in2,ipdf,ist,itype,ideriv,igsph,iasph,idrat,
     #          lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,mpl,mpkxy,mpkzz,
     #		nfield,nstep,nort,ipt,itm,itd,ipsi0,lband,kband,
     #          ldelta,kdelta,lptmx,kptmx,jkmn,jmmn,ndim,
     #          neltot,nelv,nelre,nelim
	integer	ixp
c
	common /eprdat/	 gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz,wxx,wyy,wzz,
     #			 tl,tkxy,tkzz,djf,djfprp,oss,psi,ald,bed,gad,
     #                   alm,bem,gam,c20,c22,c40,c42,c44,bshift,
     #			 gib0,gib2,b0,gamman,fieldi,fieldf,cgtol,
     #			 shiftr,shifti,a0,g0,w0,pl,pkxy,pkzz,btol,
     #                   faa(5),fgm(5),fam(2,5),fgd(2,5),fad(2,5),
     #                   fwm(5),fwd(2,5),cxp(6),cpot(5,5),xlk(5,5),
     #			 in2,ipdf,ist,itype,ideriv,igsph,iasph,
     #                   idrat,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #                   mpl,mpkxy,mpkzz,nfield,nstep,nort,
     #			 ipt,itm,itd,ipsi0,lband,kband,ldelta,kdelta,
     #                   lptmx,kptmx,jkmn,jmmn,ndim,neltot,nelv,
     #                   nelre,nelim,ixp(6,2)
c

