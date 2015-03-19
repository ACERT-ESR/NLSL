      function brent(ax,bx,cx,fb0,f,tol,ftol,xmin,iflag)
      implicit none
c
      integer iflag,iter
      double precision brent,ax,bx,cx,f,tol,ftol,xmin
      double precision a,b,d,e,etemp,fu,fv,fw,fx,r,q,p,tol1,tol2,u,v,w,
     #     x,xm,fdmy,fb0
c
      external f
c
      integer ITMAX
      double precision CGOLD,ZEPS
      parameter (ITMAX=100,CGOLD=.3819660d0,ZEPS=1.0d-10)
c
c     a is lower bound
c     b is upper bound 
c     x is function minimum so far
c     w is second-least function value
c     v is previous w
c     u is most recent function eval
c
      iflag=1
c
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
c
c     Use initial function value from mnbrak
c
      fx = fb0
c
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
         tol1=tol
         tol2=2.*tol1
         xm=0.5*(a+b)
c
c        Check for parameter convergence
c
         if(abs(x-xm).le.(tol2-.5d0*(b-a)))goto 3
c
         if(abs(e).gt.tol1) then
c
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
c
c           Take parabolic step
c
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
         endif
c
c        Take golden section step into larger of [a,x], [x,b]
c
 1       if(x.ge.xm) then
            e=a-x
         else
            e=b-x
         endif
         d=CGOLD*e
c
c        Evaluate function at new point if it is at least
c        tol away from existing minimum
c 
 2       if(abs(d).ge.tol1) then
            u=x+d
         else
            u=x+sign(tol1,d)
         endif
c
         fu=f(u,iflag)
         if (iflag.lt.0) go to 3
c
c        New step is the new minimum: old minimum is now a bracket
c
         if(fu.le.fx) then
            if(u.ge.x) then
               a=x
            else
               b=x
            endif
c
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
c
c        New step increased function: found a new bracket
c 
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
c
c     Check for function convergence
c
         if (abs(fx-fw).lt.ftol) goto 3
c
 11   continue
c
      pause 'Brent exceed maximum iterations.'
c
 3    xmin=x
      brent=fx
      return
      end
