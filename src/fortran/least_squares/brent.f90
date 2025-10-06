      function brent(ax,bx,cx,fb0,f,tol,ftol,xmin,iflag)
      implicit none
!
      integer iflag,iter
      double precision brent,ax,bx,cx,f,tol,ftol,xmin
      double precision a,b,d,e,etemp,fu,fv,fw,fx,r,q,p,tol1,tol2,u,v,w,
     #     x,xm,fb0
!
      external f
!
      integer ITMAX
      double precision CGOLD,ZEPS
      parameter (ITMAX=100,CGOLD=.3819660d0,ZEPS=1.0d-10)
!
!     a is lower bound
!     b is upper bound 
!     x is function minimum so far
!     w is second-least function value
!     v is previous w
!     u is most recent function eval
!
      iflag=1
!
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
!
!     Use initial function value from mnbrak
!
      fx = fb0
!
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
         tol1=tol
         tol2=2.*tol1
         xm=0.5*(a+b)
!
!        Check for parameter convergence
!
         if(abs(x-xm).le.(tol2-.5d0*(b-a)))goto 3
!
         if(abs(e).gt.tol1) then
!
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            if(q.gt.0.) p=-p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or. 
     *           p.ge.q*(b-x)) goto 1
!
!           Take parabolic step
!
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
         endif
!
!        Take golden section step into larger of [a,x], [x,b]
!
 1       if(x.ge.xm) then
            e=a-x
         else
            e=b-x
         endif
         d=CGOLD*e
!
!        Evaluate function at new point if it is at least
!        tol away from existing minimum
! 
 2       if(abs(d).ge.tol1) then
            u=x+d
         else
            u=x+sign(tol1,d)
         endif
!
         fu=f(u,iflag)
         if (iflag.lt.0) go to 3
!
!        New step is the new minimum: old minimum is now a bracket
!
         if(fu.le.fx) then
            if(u.ge.x) then
               a=x
            else
               b=x
            endif
!
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
!
!        New step increased function: found a new bracket
! 
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
!
!     Check for function convergence
!
         if (abs(fx-fw).lt.ftol) goto 3
!
 11   continue
!
      print *, 'Brent exceed maximum iterations.'
      print *, '[execution paused, press enter to continue]'
      read (*,*)
!
 3    xmin=x
      brent=fx
      return
      end
