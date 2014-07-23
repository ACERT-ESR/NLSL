c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       =================
c                       SUBROUTINE COMPEV
c                       =================
c
c     Subroutine COMPEV calls subroutine CMTQLI to compute the
c     T-eigenvalues.  COMPEV then applies the T-multiplicity and
c     spurious tests to the computed T-eigenvalues.
c
c     on return from COMPEV:
c        ndis=number of distinct eigenvalues of T(1,mev)
c        vs=distinct T-eigenvalues in increasing order of magnitude
c        gr(k)=|vs(k)|, k=1,ndis,  gr(k).le.gr(k+1)
c        mp=T-multiplicities of the T-eigenvalues in vs
c        mp(i)=(0,1,mi), mi>1, i=1,ndis  means:
c           (0)  vs(i) is spurious
c           (1)  vs(i) is simple and good
c           (mi) vs(i) is T-multiple and is therefore not only good but
c                also a converged good T-eigenvalue.
c
c     includes:
c           nlsdim.inc
c           eprprm.inc
c     uses:
c           cmtqli.f
c
c*********************************************************************
c
      subroutine compev(alpha,beta,v1,v2,vs,w,evmag,multol,sputol,
     #                  mp,t2flag,mev,ndis,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
c
      complex*16 alpha,beta,vs,v1,v2,w,w2,eval,ctemp,wtemp,czero
      dimension alpha(mxstep),beta(mxstep),vs(mxstep),v1(mxstep),
     #            v2(mxstep),w(mxstep),w2(mxstep)
      parameter (czero=(0.0D0,0.0D0))
      double precision evmag(mxstep)
      double precision temp,dgap,tol,delmin
      double precision multol,sputol,evalr,evalc
      integer i,j,mev,mev1,imin,index,jp1,isort,ndis,ierr,mp,t2flag
      dimension mp(mxstep),t2flag(mxstep)
c
c---------------------------------------------------------------------
c
      mev1=mev - 1
c
      do 50 j=1,mev
        vs(j)=alpha(j)
 50     v1(j)=beta(j)
c
      v1(mev)=czero
c
      isort=0
      call cmtqli(mev,vs,v1,w,isort,ierr)
      if (ierr.ne.0) return
c
c     T-eigenvalues are in vs in increasing order of magnitude
c
      do 100 j=1,mev
 100  evmag(j)=cdabs(vs(j))
c
      multol=multol*evmag(mev)
      sputol=sputol*evmag(mev)
      tol=1000.0d0*sputol
c
c----------------------------------------------------------------------
c     T-multiplicity determination
c----------------------------------------------------------------------
c
      j=0
      ndis=0
      do 150 i=1,mev
 150  t2flag(i)=0
c
 160  j=j+1
c
      if (j.le.mev) then
        if (t2flag(j).eq.1) then
          go to 160
        else
          wtemp=w(j)
          ctemp=vs(j)
          eval=ctemp
          temp=evmag(j)
          ndis=ndis+1
          index=1
          t2flag(j)=1
          i=j
 170      i=i+1
          if (i.gt.mev) go to 180
          if (t2flag(i).eq.1) go to 170
          dgap=evmag(i)-temp
          if (dgap.gt.multol) go to 180
          dgap=cdabs(eval-vs(i))
          if (dgap.gt.multol) go to 170
c
c     T-multiplicity increases
c
          index=index+1
          ctemp=ctemp+vs(i)
          wtemp=wtemp+w(i)
          t2flag(i)=1
          go to 170
c
c     T-multiplicity for vs(ndis) has been determined
c
 180      vs(ndis)=ctemp/dble(index)
          w(ndis)=wtemp/dble(index)
          mp(ndis)=index
          go to 160
        end if
      end if
c
c----------------------------------------------------------------------
c     T(2,mev) eigenvalue calculation
c----------------------------------------------------------------------
c
      do 210 j=1,mev1
        jp1=j+1
        v2(j)=alpha(jp1)
 210    v1(j)=beta(jp1)
c
      v1(mev1)=czero
c
      isort=0
      call cmtqli(mev1,v2,v1,w2,isort,ierr)
      if (ierr.ne.0) return
c
      do 250 j=1,mev1
 250  evmag(j)=cdabs(v2(j))
c
c----------------------------------------------------------------------
c     test for the spurious eigenvalues
c----------------------------------------------------------------------
c
      do 280 i=1,mev1
 280  t2flag(i)=0
c
c     go through the eigenvalues of t2-hat.  find the closest eigenvalue
c     of T(1,mev).  if it is T-multiple go on.  if it is simple declare
c     spurious whenever delmin < sputol by setting mp(i)=0
c
      j=0
 290  j=j+1
      if (j.gt.mev1) go to 390
c
      temp=evmag(j)
      eval=v2(j)
      evalr=temp+sputol
      evalc=temp-sputol
      delmin=2.D0*cdabs(vs(mev))
      imin=0
c
c     backward search
c
      i=j+1
 310  i=i-1
      if(i.lt.1) go to 320
      if(i.gt.ndis) i=ndis
c
      temp=cdabs(vs(i))
      if (temp.lt.evalc) go to 320
      if(mp(i).eq.0) go to 310
      dgap=cdabs(vs(i)-eval)
      if (dgap.ge.delmin) go to 310
      delmin=dgap
      imin=i

      go to 310
c
c     forward search
c
 320  i=j
 330  i=i+1
      if(i.gt.ndis) go to 340
c
      temp=cdabs(vs(i))
      if (temp.gt.evalr) go to 340
      if(mp(i).eq.0) go to 330
      dgap=cdabs(vs(i)-eval)
      if (dgap.ge.delmin) go to 330
      delmin=dgap
      imin=i
c
      go to 330
c
 340  continue
c
      if(imin.eq.0) go to 370
c
      if(delmin.gt.sputol) go to 290
      if(mp(imin).gt.1)  go to 290
      mp(imin)=0
c
      go to 290
c
 370  continue
c
      go to 290
c
 390  continue
c
      do 400 j=1,ndis
 400  evmag(j)=cdabs(vs(j))
c
      return
c
      end
