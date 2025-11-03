      subroutine mnbrak(ax,bx,cx,fa,fb,fc,func,iflag,bound)
!
      use errmsg
!
      implicit none
      double precision ax,bx,cx,fa,fb,fc,func,bound
      integer iflag
      external func
!
      double precision dum,fu,plim,r,q,u,ulim
!
      double precision GLIMIT,GOLD,TINY
      parameter (GOLD=1.618034D0,GLIMIT=100.0D0,TINY=1.0D-20)
!
!  -- Evaluate function at points a and b
!
      fa=func(ax,iflag)
      if (iflag.lt.0) return
!
      fb=func(bx,iflag)
      if (iflag.lt.0) return
!
!  -- Rearrange the two points so that a to b is always "downhill"
!     Set parameter search boundary along this direction
!
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
!
!     Choose a third point using golden ratio
!
      cx=bx+GOLD*(bx-ax)
!
      fc=func(cx,iflag)
      if (iflag.lt.0) return
!
!==================================================
!     Main loop: do while f(B) > f(C)
!                (if f(B) < f(C), we're done!)
!==================================================
!
 1    if (fb.ge.fc) then
!
!        ----------------------------------------------------
!        Check that the point with the lowest function value
!        does not go beyond the initial step bound
!        ----------------------------------------------------
         if ((cx-ax)*(plim-cx).le.0.) then
            iflag=NOMIN
            return
         end if
!
!        -----------------------------------------------
!        Take parabolic step based on three known points
!        -----------------------------------------------
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
!
!        ---------------------------------------------------------------
!        Choose a limit that is the smaller of the parabolic step limit
!        and the step bound from the initial point (ax)
!        ----------------------------------------------------------------
         ulim=bx+GLIMIT*(cx-bx)
         if ((ulim-ax)*(plim-ulim).lt.0.) ulim=plim
!
!        --------------------------------
!        Step is between B and C: try it
!        --------------------------------
         if((bx-u)*(u-cx).gt.0.)then
            fu=func(u,iflag)
            if (iflag.lt.0) return
!
!           ---------------------------------------
!           U is a minimum: bracket is B,U,C (exit)
!           ---------------------------------------
            if(fu.lt.fc)then
               ax=bx
               fa=fb
               bx=u
               fb=fu
               go to 1
!
!           -----------------------------------
!           f(U)>f(B): bracket is A,B,U (exit)
!           -----------------------------------
            else if(fu.gt.fb)then
               cx=u
               fc=fu
               go to 1
            endif
!
!           ------------------------------------------------------
!           Have f(A) > f(B) > f(U) > f(C). Try stepping further.
!           ------------------------------------------------------
            u=cx+GOLD*(cx-bx)
            fu=func(u,iflag)
            if (iflag.lt.0) return
!
!        ------------------------------------------------
!        Step is between C and its allowed limit: try it
!        ------------------------------------------------
         else if((cx-u)*(u-ulim).gt.0.)then
            fu=func(u,iflag)
            if (iflag.lt.0) return
!
!           -------------------------------------------------------
!           Have f(A) > f(B) > f(C) > f(U): reset upper bound to B 
!           -------------------------------------------------------
            if(fu.lt.fc)then
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fb=fc
               fc=fu
               fu=func(u,iflag)
               if (iflag.lt.0) return
            endif
!
!        -----------------------------------------------------------
!        Step went beyond the limiting value: try function at limit 
!        -----------------------------------------------------------
         else if((u-ulim)*(ulim-cx).ge.0.)then
            u=ulim
            fu=func(u,iflag)
            if (iflag.lt.0) return
!
!        ------------------------------------------------------------
!        Reject parabolic step and use golden section magnification
!        ------------------------------------------------------------
         else
            u=cx+GOLD*(cx-bx)
            if ((u-ulim)*(ulim-cx).ge.0.) u=ulim
            fu=func(u,iflag)
            if (iflag.lt.0) return
         endif
!
!        ---------------------------------
!        Discard oldest point and continue
!        ---------------------------------
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         go to 1
      endif
!
      iflag=0
      return
      end
