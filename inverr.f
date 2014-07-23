c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c
c                       =================
c                       SUBROUTINE INVERR
c                       =================
c
c     Subroutine INVERR computes error estimates for computed isolated
C     good T-eigenvalues in vs and writes these eigenvalues and
c     estimates to file.  By definition a good T-eigenvalue is
c     isolated if its closest neighbor is also good, or if its closest
c     neighbor is spurious but that neighbor is far enough away.  So
c     in particular, we will compute estimates for any good
c     T-eigenvalue that is in a cluster of good T-eigenvalues.
c
c     uses inverse iteration on T(1,mev) solving the equation
c     (T - x*I)v2 = right-hand side for each such good T-eigenvalue x.
c
c     program refactors T-x*I on each iteration of inverse iteration.
c     typically only one iteration is needed per T-eigenvalue x.
c
c     on entry and exit:
c
c        mev = order of T
c        n = order of original matrix A
c        vs = computed distinct eigenvalues of T(1,mev)
c        mp = T-multiplicity of each T-eigenvalue in vs.
c             mp(i) = -1 means vs(i) is a good T-eigenvalue but that
c             it is sitting close to a spurious T-eigenvalue.
c             mp(i) = 0 means vs(i) is spurious.  Estimates are
c             computed only for those T-eigenvalues with mp(i) = 1.
c             Flagging was done in subroutine isoev prior to entering
c             INVERR.
c        niso = number of isolated good T-eigenvalues contained in vs
c        ndis =  number of distinct T-eigenvalues in vs
c
c     In program:
c
c        iter = maximum number of inverse iteration steps allowed for
c             each x.  iter = it on entry.
c        gr,gc = arrays of dimension at least mev+niso.  used to
c             store randomly-generated right-hand side.  This is not
c             regenerated for each x. g is also used to store error
c             estimates as they are computed for later printout.
c        v1,v2 = work spaces used in the factorization of T(1,mev).
c             At the end of the inverse iteration computation for x,
c             v2 contains the unit eigenvector of T(1,mev) corresponding
c             to x.  v1 and v2 must be of dimension at least mev.
c
c     on exit:
c
c        gg(j) = minimum gap in T(1,mev) for each vs(j), j=1,ndis
c        g(i) = |betam|*|v2(mev)| = error estimate for isolated good
c             T-eigenvalues, where i = 1,niso  and  betam = beta(mev)
c             T(1,mev) corresponding to ith isolated good T-eigenvalue.
c             If for some x it.gt.iter then the error estimate in g
c             is marked with a - sign.
c        v2 = isolated good T-eigenvalues
c        v1 = minimal T-gaps for the T-eigenvalues in v2.
c             These are constructed for write-out purposes only and not
c             needed elsewhere in the program.
c
c     includes :
c
c        nlsdim.inc
c        eprprm.inc
c
c*********************************************************************
c
      subroutine inverr(a,b,v1,v2,vs,eps,gr,gc,g,gg,mp,
     #                  intc,mev,ndis,niso,it)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
c
      integer  i,ii,iso,igood,ispur,it,iter
      integer  j,jev,mev,mp1,mm1,ndis,ng,niso
      integer  mp(mxstep), intc(mxstep)
c
      real g(mxstep),gg(mxstep)
c
      double precision est,estr,estc,sum,xu,norm,tsum,gsum,gap
      double precision  eps,eps3,eps4,zero,one,gr(mxstep),gc(mxstep)
c
      complex*16  a(mxstep),b(mxstep),v1(mxstep),v2(mxstep)
      complex*16  u,z,x,ratio,betam,temp,czero,csum,znormu,vs(mxstep)
c
      parameter (one=1.0D0,zero=0.0D0,czero=(0.0D0,0.0D0))
c
      external znormu
c
c#######################################################################
c
c     initialization and parameter specification
c
      ng=0
      niso=0
      iter=it
      mp1=mev+1
      mm1=mev-1
      betam=b(mev)
      b(mev)=czero
c
c     calculate scale and tolerances
      tsum=zero
      do 30 i=1,mev
        tsum=tsum+cdabs(a(i))+cdabs(b(i))
 30   continue
c
      eps3=eps*tsum
      eps4=dble(mev)*eps3
c
c     generate scaled right-hand side
      gsum=zero
      do 60 i=1,mev
        gsum=gsum+dabs(gr(i))+dabs(gc(i))
 60   continue
      gsum=eps4/gsum
c
      do 70 i=1,mev
      gr(i)=gsum*gr(i)
 70   gc(i)=gsum*gc(i)
c
c     loop on isolated good T-eigenvalues in vs (mp(i)=1) to
c     calculate corresponding unit eigenvector of T(1,mev)
c
      do 200 jev=1,ndis
      if (mp(jev).eq.0) go to 200
      ng=ng+1
      if (mp(jev).ne.1) go to 200
      it=1
      niso=niso+1
      x=vs(jev)
c
c     initialize right hand side for inverse iteration
c     and the flag on which rows are interchanged
      do 80 i=1,mev
        intc(i)=0
        v2(i)=dcmplx(gr(i),gc(i))
 80   continue
c
c-----------------------------------------------------------------------
c      triangular factorization
c-----------------------------------------------------------------------
c
 90   continue
      u=a(1)-x
      z=b(1)
c
      do 110 i=2,mev
        if (cdabs(b(i-1)).le.cdabs(u)) then
          v1(i-1)=z/u
          v2(i-1)=v2(i-1)/u
          v2(i)=v2(i)-b(i-1)*v2(i-1)
          ratio=b(i-1)/u
          u=a(i)-x-z*ratio
          z=b(i)
        else
          ratio=u/b(i-1)
          intc(i)=1
          v1(i-1)=a(i)-x
          u=z-ratio*v1(i-1)
          z=-ratio*b(i)
          temp=v2(i-1)
          v2(i-1)=v2(i)
          v2(i)=temp-ratio*v2(i)
        end if
 110  continue
c
      if (cdabs(u).eq.zero) u=dcmplx(eps3,eps3)
c
c     back substitution
      v2(mev)=v2(mev)/u
      do 130 ii=1,mm1
      i=mev-ii
      if (intc(i+1).ne.1) then
        v2(i)=v2(i)-v1(i)*v2(i+1)
      else
        v2(i)=(v2(i)-v1(i)*v2(i+1)-b(i+1)*v2(i+2))/b(i)
      end if
 130  continue
c
c-----------------------------------------------------------------------
c      tests for convergence of inverse iteration
c-----------------------------------------------------------------------
c
      norm=cdabs(v2(mev))
      do 140 ii=1,mm1
        i=mev-ii
        norm=norm+cdabs(v2(i))
 140  continue
c
      if (norm.ge.one) go to 160
c
      it=it+1
      if (it.gt.iter) go to 160
c
      xu=eps4/norm
c
      do 150 i=1,mev
        intc(i)=0
        v2(i)=v2(i)*xu
 150  continue
c
      go to 90
c     another inverse iteration step
c
c     inverse iteration finished
c     normalize computed T-eigenvector : v2=v2/||v2||
 160  continue
c
       csum=znormu(v2,mev)
c
       do 170 ii=1,mev
          v2(ii)=v2(ii)/csum
 170   continue
c
c     save error estimate for later output
      est=cdabs(betam)*cdabs(v2(mev))
      estr=dabs(dreal(v2(mev)))
      estc=dabs(dimag(v2(mev)))
      gsum=cdabs(betam)
      if (it.gt.iter) est=-est
      g(niso)=est
c
 200  continue
c
c     end error estimate loop on isolated good T-eigenvalues.
c     Generate distinct mingaps for T(1,mev).  This is useful as an
c     indicator of the goodness of the inverse iteration estimates.
c     Transfer isolated good T-eigenvalues and corresponding tmingaps
c     to v2 and v1 for output purposes only.
c
      iso=0
      do 210 j=1,ndis
      if (mp(j).ne.1) go to 210
      iso=iso+1
      gr(iso)=gg(j)
      v2(iso)=vs(j)
 210  continue
      if(niso.eq.0) go to 270
c
      ispur=0
      i=0
      do 260 j=1,ndis
      if(mp(j).ne.0) go to 240
      ispur=ispur+1
      go to 260
  240 if(mp(j).ne.1) go to 260
      i=i+1
      igood=j-ispur
  260 continue
c
c     restore b(mev+1)=betam
 270  b(mev)=betam
c
      return
c
      end
