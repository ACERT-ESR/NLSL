      subroutine mnbrak(ax,bx,cx,fa,fb,fc,func,iflag,bound)
      implicit none
      double precision ax,bx,cx,fa,fb,fc,func,bound
      integer iflag
      external func
c
      double precision dum,fu,plim,r,q,u,ulim
c
      double precision GLIMIT,GOLD,TINY
      parameter (GOLD=1.618034D0,GLIMIT=100.0D0,TINY=1.0D-20)
c
      include 'errmsg.inc'
c
c  -- Evaluate function at points a and b
c
      fa=func(ax,iflag)
      if (iflag.lt.0) return
c
      fb=func(bx,iflag)
      if (iflag.lt.0) return
c
c  -- Rearrange the two points so that a to b is always "downhill"
c     Set parameter search boundary along this direction
c
      if(fb.gt.fa)then
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
         plim=bx+sign(bound,bx-ax)
      else
         plim=ax+sign(bound,bx-ax)
      end if
c
c     Choose a third point using golden ratio
c
      cx=bx+GOLD*(bx-ax)
c
      fc=func(cx,iflag)
      if (iflag.lt.0) return
c
c==================================================
c     Main loop: do while f(B) > f(C)
c                (if f(B) < f(C), we're done!)
c==================================================
c
 1    if (fb.ge.fc) then
c
c        ----------------------------------------------------
c        Check that the point with the lowest function value
c        does not go beyond the initial step bound
c        ----------------------------------------------------
         if ((cx-ax)*(plim-cx).le.0.) then
            iflag=NOMIN
            return
         end if
c
c        -----------------------------------------------
c        Take parabolic step based on three known points
c        -----------------------------------------------
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
c
c        ---------------------------------------------------------------
c        Choose a limit that is the smaller of the parabolic step limit
c        and the step bound from the initial point (ax)
c        ----------------------------------------------------------------
         ulim=bx+GLIMIT*(cx-bx)
         if ((ulim-ax)*(plim-ulim).lt.0.) ulim=plim
c
c        --------------------------------
c        Step is between B and C: try it
c        --------------------------------
         if((bx-u)*(u-cx).gt.0.)then
            fu=func(u,iflag)
            if (iflag.lt.0) return
c
c           ---------------------------------------
c           U is a minimum: bracket is B,U,C (exit)
c           ---------------------------------------
            if(fu.lt.fc)then
               ax=bx
               fa=fb
               bx=u
               fb=fu
               go to 1
c
c           -----------------------------------
c           f(U)>f(B): bracket is A,B,U (exit)
c           -----------------------------------
            else if(fu.gt.fb)then
               cx=u
               fc=fu
               go to 1
            endif
c
c           ------------------------------------------------------
c           Have f(A) > f(B) > f(U) > f(C). Try stepping further.
c           ------------------------------------------------------
            u=cx+GOLD*(cx-bx)
            fu=func(u,iflag)
            if (iflag.lt.0) return
c
c        ------------------------------------------------
c        Step is between C and its allowed limit: try it
c        ------------------------------------------------
         else if((cx-u)*(u-ulim).gt.0.)then
            fu=func(u,iflag)
            if (iflag.lt.0) return
c
c           -------------------------------------------------------
c           Have f(A) > f(B) > f(C) > f(U): reset upper bound to B 
c           -------------------------------------------------------
            if(fu.lt.fc)then
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fb=fc
               fc=fu
               fu=func(u,iflag)
               if (iflag.lt.0) return
            endif
c
c        -----------------------------------------------------------
c        Step went beyond the limiting value: try function at limit 
c        -----------------------------------------------------------
         else if((u-ulim)*(ulim-cx).ge.0.)then
            u=ulim
            fu=func(u,iflag)
            if (iflag.lt.0) return
c
c        ------------------------------------------------------------
c        Reject parabolic step and use golden section magnification
c        ------------------------------------------------------------
         else
            u=cx+GOLD*(cx-bx)
            if ((u-ulim)*(ulim-cx).ge.0.) u=ulim
            fu=func(u,iflag)
            if (iflag.lt.0) return
         endif
c
c        ---------------------------------
c        Discard oldest point and continue
c        ---------------------------------
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         go to 1
      endif
c
      iflag=0
      return
      end
