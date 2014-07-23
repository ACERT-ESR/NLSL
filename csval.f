c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       ================
c                       SUBROUTINE CSVAL
c                       ================
c
c     This subroutine calculates the distinct eigenvalues of a given
c     tridiagonal matrix.  The identification test for good versus
c     spurious eigenvalues is done according to the procedure in
c     Cullum & Willoughby's book.
c
c       Arguments:
c           alpha    diagonal T matrix elements
c           beta     extradiagonal T matrix elements
c           v1,v2    work spaces
c           gr,gc    work spaces
c           vs       work space, on exit = eigenvalues of T-matrix sorted
c                    in the order of decreasing weighting factor
c           w        weighting factor
c           ngood    # of good eigenvalues (passed via eprdat.inc common)
c           ndiag    dimension of T-matrix used in the calculation of the
c                    eigenvalues (passed via eprdat.inc common)
c           ndone    dimension of T-matrix generated in Lanczos routine
c           ierr     error flag for failure in cmtqli routine
c
c     modified from CSLEVAL.F by SHL 17-AUG-92
c
c     Includes :
c               nlsdim.inc
c               stdio.inc
c               eprprm.inc
c               rndoff.inc
c     Uses :
c               cmtqli.f
c               compev.f
c               lump.f
c               comgap.f
c               isoev.f
c               inverr.f
c
c*********************************************************************
c
      subroutine csval(alpha,beta,v1,v2,gr,gc,vs,w,ndone,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'simparm.inc'
      include 'rndoff.inc'
c
      integer  mxinit,mev,ndis,loop,ndone
      integer  i,j,niso,it,iwrite,l,m,mmb,ierr
      integer  mp(mxstep),mp2(mxstep)
c
      real*8  g(mxstep),gg(mxstep),gtp,ggtp
c
      double precision  gr(mxstep),gc(mxstep)
      double precision  gaptol,ttol,epsm,reltol,evmax,t0,t1
      double precision  scale1,scale2,sputol,contol,multol
      double precision  one,zero
      parameter (one=1.0D0,zero=0.0D0)
c
      complex*16  alpha(mxstep),beta(mxstep),vs(mxstep),w(mxstep)
      complex*16  v1(mxstep),v2(mxstep)
      complex*16  betam,wt,egv
      complex*16  czero,ci
      parameter (czero=(0.0D0,0.0D0),ci=(0.0D0,1.0D0))
c
c#####################################################################
c
      if (sptst.eq.1) go to 99
c
c---------------------------------------------------------------------
c     Eigenvalue calculation without test for spurious eigenvalues
c---------------------------------------------------------------------
c
	if (ndone.gt.mxstep) then
	end if
c      stop
      do 82 i=1,ndone
        vs(i)=alpha(i)
        v1(i)=beta(i)
 82   continue
c
      call cmtqli(ndone,vs,v1,w,1,ierr)
c
      if (ierr.ne.0) then
        write(lulog,1000)
        if (lulog.ne.luttyo) write(luttyo,1000)
        ierr=-abs(ierr)
        return
      end if
c
      ngood=ndone
      ndiag=ndone
c
      return
c
c---------------------------------------------------------------------
c     Eigenvalue calculation with test for spurious eigenvalues
c---------------------------------------------------------------------
c
 99   epsm=100.0D0*rndoff
      scale1=5.0D2
      scale2=5.0D0
      gaptol=1.0D-7
      reltol=1.0D-8
      mxinit=10      
      mev=ndone
c
      betam=beta(mev)
      multol=500.D0*dfloat(mev+1000)*epsm
      sputol=multol
c
c---------------------------------------------------------------------
c     calculate the eigenvalues and test for the spurious values
c---------------------------------------------------------------------
c
      call compev(alpha,beta,v1,v2,vs,w,gr,multol,sputol,mp,mp2,
     #            mev,ndis,ierr)
c
      if (ierr.ne.0) then
        write(lulog,1000)
        if (lulog.ne.luttyo) write(luttyo,1000)
        ierr=-abs(ierr)
        return
      end if
c
      if (ndis.eq.0) then
         write (lulog,1010)
         if (lulog.ne.luttyo) write (luttyo,1010)
         ierr=-abs(ierr)
         return
      end if
c
      evmax=gr(ndis)
      loop=ndis
c
c---------------------------------------------------------------------
c
      call lump(vs,v1,gr,w,reltol,sputol,scale2,mp,mp2,loop)
c
      ndis=loop
c
c---------------------------------------------------------------------
c     calculate mingaps for distinct t(1,mev) eigenvalues.
c---------------------------------------------------------------------
c
      call comgap(vs,gr,gg,mp,mp2,ndis)
c
c     set convergence critierion
c
      ttol=epsm*evmax
      contol=cdabs(betam)*1.D-10
c
      beta(mev)=betam
c
c---------------------------------------------------------------------
c
      call isoev(vs,gr,gg,gaptol,sputol,scale1,mp,ndis,niso)
c      stop
c
      if (niso.eq.0) go to 370
c
c---------------------------------------------------------------------
c     error estimates for isolated good T-eigenvalues
c---------------------------------------------------------------------
c
      it=mxinit
c
      do 330 i=1,ndone
         gr(i)=zero
         gc(i)=zero
 330  continue
c
      gr(1)=one
      gc(1)=one
c
      call inverr(alpha,beta,v1,v2,vs,epsm,gr,gc,g,gg,
     #            mp,mp2,mev,ndis,niso,it)
c
c      write(luttyo,340) contol
c
      do 350 i=1,niso
         if (abs(g(i)).gt.contol) then
            write (lulog,1020)
            if (lulog.ne.luttyo) write (luttyo,1020)
            go to 370
         end if
 350  continue
c
c      write(luttyo,360) contol
c
 370  ngood=0
      do 390 i=1,ndis
      if (mp(i).eq.0) go to 390
      ngood=ngood+1
      mp(ngood)=mp(i)
      vs(ngood)=vs(i)
      w(ngood)=w(i)
      g(ngood)=g(i)
      gg(ngood)=gg(i)
 390  continue
c
c---------------------------------------------------------------------
c     sort goodev in the order of decreasing weighting factor
c---------------------------------------------------------------------
c
c      write (luttyo,440) ngood,niso,ndis
c
      do 450 l=1,ngood-1
         t0=cdabs(w(l))
         do 450 m=l+1,ngood
            t1=cdabs(w(m))
            if (t1.gt.t0) then
               t0=t1
               wt=w(m)
               w(m)=w(l)
               w(l)=wt
               egv=vs(m)
               vs(m)=vs(l)
               vs(l)=egv
               gtp=g(m)
               g(m)=g(l)
               g(l)=gtp
               ggtp=gg(m)
               gg(m)=gg(l)
               gg(l)=ggtp
            end if
 450  continue
c
      ndiag=mev
c
      beta(mev)=betam
c
c     end of loop on different size T-matrices allowed.
c
c      write(luttyo,750) ngood,ndiag
c      write(lulog,750) ngood,ndiag
c
      return
c
c=====================================================================
c     format statements
c=====================================================================
c
 340  format(' convergence is tested using the convergence tolerance',
     #        e13.4)
 360  format(' all computed error estimates were less than',e15.4/
     #       ' therefore procedure terminates')
 440  format(/i6,' good T-eigenvalues have been computed'/
     #       i6,' of these are isolated'/
     #       i6,'=number of distinct T-eigenvalues computed'/)
 750  format(/,2x,'# of good eigenvalues : ',i4,/,2x,'Actual ',
     #     'Lanczos steps used in the eigenvalue estimation : ',i4)
 1000 format(2x,'on return from CMTQLI, error flag was not zero'/)
 1010 format(2x,'on return from COMPEV no distinct eigenvalues ',
     #     'were found'/)
 1020 format(2x,'inverse iteration did not converge'/)
c
      end
