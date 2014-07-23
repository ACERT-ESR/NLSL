c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       =================
c                       SUBROUTINE COMGAP
c                       =================
c
c     Subroutine COMGAP computes the minimal gaps between
c     distinct T-eigenvalues supplied.
c
c     on return :
c        |gg(k)| = min(j.ne.k,|vs(k)-vs(j)|)
c        mp2(k)=j index for minimum.
c
c        gg(k) < 0 means nearest neighbor is spurious.
c
c     includes :
c
c        nlsdim.inc
c
c*********************************************************************
c
      subroutine comgap(vc,va,gg,mp,ind,m)
c
      implicit none
c
      include 'limits.inc'
c
      complex*16 vc(mxstep),z
      double precision  va(mxstep),t0,t1,tu,tk
      real  gg(mxstep),gtemp
      integer  index,itemp,j,k,kk,km1,kp1,l,m
      integer  mp(mxstep),ind(mxstep)
c
c---------------------------------------------------------------------
c     va(k) = |vc(k)|  va(k) <= va(k+1)
c     gg(k) = min |vc(k)-vc(j)|  j .ne. k.
c---------------------------------------------------------------------
c
      tu=va(m)+va(m)
      k=0
 10   k=k+1
      if (k.gt.m) return
      index=0
      t1=tu
      tk=va(k)
      z=vc(k)
      j=k
c     backwards
 20   j=j-1
      if (j.lt.1) go to 30
      t0=tk-va(j)
      if (t0.gt.t1) go to 30
      t0=cdabs(z-vc(j))
      if (t1.le.t0) go to 20
      t1=t0
      index=j
      go to 20
c     forwards
 30   j=k
 40   j=j+1
      if (j.gt.m) go to 50
      t0=va(j)-tk
      if (t0.gt.t1) go to 50
      t0=cdabs(z-vc(j))
      if (t1.le.t0) go to 40
      t1=t0
      index=j
      go to 40
 50   ind(k)=index
      gg(k)=t1
      if(mp(index).eq.0)  gg(k)=-gg(k)
      go to 10
c
      end
