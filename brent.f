c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                         =======================
c                            function BRENT
c                         =======================
c
c   Modified from the routine in NUMERICAL RECIPES for the non-linear
c   least squares iteration.  Notice that it uses the absolute tolerance
c   instead of the fractional tolerence.
c
c   Error flag of negative value means the function "func" failed.
c
c   Given a function func, and given a bracketing triplet of abscissas ax,
c   bx, cx (such that bx is between ax and cx, and func(bx) is less than
c   both func(ax) and func(cx), this routine isolates the minimum to
c   a ABSOLUTE precision of about tol using Brent's method.
c
c   Includes :
c       stdio.inc, iterat.f
c----------------------------------------------------------------------
c
      function brent(ax,bx,cx,func,tol,xmin,prmID,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'parms.inc'
c
      character prmID*6
      double precision brent,ax,bx,cx,func,tol,xmin,cgold
      double precision a,b,d,e,etemp,fu,fv,fw,fx,r,q,p,tol1,tol2,u,v,w,
     #     x,xm
      integer iterr,itmax,ierr
      parameter (itmax=100,cgold=.3819660d0)
c
      external func
c
      tol1=tol
      tol2=2.*tol1
c
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=func(x,ierr)
      write(*,*)'brent1 call func, x= ',x,fx
      if (ierr.lt.0 .or .ihltcmd.ne.0) return
c      seteval=.true.	! got eigenvalues - not used anymore
      fv=fx
      fw=fx
      do 11 iterr=1,itmax
        xm=0.5*(a+b)
c
        if (iterr.eq.1) then
          write (luout,1000) prmID
          if (luout.ne.luttyo) write (luttyo,1000) prmID
        end if
c
        write (luout,1010) iterr,fx,x
        if (luout.ne.luttyo) write (luttyo,1010) iterr,fx,x
c
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or. 
     *        p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=cgold*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=func(u,ierr)
      write(*,*)'brent2 call func, x= ',x,fu
        if (ierr.lt.0 .or. ihltcmd.ne.0) return
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'Brent exceed maximum iterations.'
3     xmin=x
      brent=fx
      write (luout,1020)
      if (luout.ne.luttyo) write (luttyo,1020)
      return
c
c ### format statements ##############################################
c
 1000 format(/,28('-'),/'Iter',5x,'RmsDv',6x,a,/28('-'))
 1010 format(i4,4x,g10.4,2x,g24.18)
 1020 format(28('-')/)
c
      end
