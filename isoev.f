c NLSPMC Version 1.0 2/5/99
c*********************************************************************
c
c                       ================
c                       SUBROUTINE ISOEV
c                       ================
c
c     Subroutine ISOEV uses tmingaps to label the isolated good
c     T-eigenvalues that are very close to spurious ones.
c     Error estimates will not be computed for these T-eigenvalues.
c
c     on entry and exit:
c        vs contains the computed distinct eigenvalues of t(1,mev)
c        gr(k) = |vs(k)|, k = 1,ndis, gr(k).le.gr(k+1)
c        gg(k) = min(j.ne.k,|vs(k)-vs(j)|)   mingap
c        mp contains the corresponding T-multiplicities
c        ndis = number of distinct T-eigenvalues
c        gaptol = relative gap tolerance set in main
c
c     on exit:
c        mp(j) is not changed except that  mp(j)=-1, if mp(j)=1,
c        and a spurious T-eigenvalue is too close.
c
c        if mp(i)=-1 that simple good T-eigenvalue will be skipped
c        in the subsequent error estimate computations in inverr
c        that is, we compute error estimates only for those good
c        T-eigenvalues with mp(j)=1.
c
c     includes :
c
c        nlsdim.inc
c        eprprm.inc
c
c*********************************************************************
c
      subroutine isoev(vs,gr,gg,gaptol,sputol,scale1,mp,ndis,niso)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
c
      complex*16 vs(mxstep),t0
      double precision  sputol,gaptol,scale1,temp,tol,tj,dgap,one
      double precision  gr(mxstep)
      real gg(mxstep)
      integer  i,j,ndis,niso,mp(mxstep)
c
c---------------------------------------------------------------------
c
      one=1.0D0
      dgap=scale1*sputol
      niso=0
      do 40 j=1,ndis
      if (mp(j).ne.1) go to 40
      tj=gr(j)
      t0=vs(j)
      tol=dmax1(dgap,gaptol*tj)
c     tol=dmax1(one,tj)*gaptol
c
c     vs(j) is next simple good T-eigenvalue
c
      niso=niso+1
      if (abs(gg(j)).gt.tol) go to 40
      i=j
 10   i=i-1
      if (i.lt.1) go to 20
      if (tj-gr(i).gt.tol) go to 20
      if (mp(i).ne.0) go to 10
      temp=cdabs(t0-vs(i))
      if (temp.gt.tol) go to 10
      mp(j)=-mp(j)
      niso=niso-1
      go to 40
 20   i=j
 30   i=i+1
      if (i.gt.ndis) go to 40
      if (gr(i)-tj.gt.tol) go to 40
      if (mp(i).ne.0) go to 30
      temp=cdabs(t0-vs(i))
      if (temp.gt.tol) go to 30
      mp(j)=-mp(j)
      niso=niso-1
 40   continue
c
      return
c
      end
