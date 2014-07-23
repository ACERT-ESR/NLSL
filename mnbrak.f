c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                         =======================
c                            subroutine MNBRAK
c                         =======================
c
c   Modified from the routine in NUMERICAL RECIPES to restrict the search
c   region according to the user-specified range.  The error flag (ierr=1)
c   is returned if it does not find the bracket for the minimum.  Note
c   that the region for the parabolic extrapolation to be valid is reduced
c   (glimit=10.0).
c
c   Error flag of negative value means the function "func" failed.
c   
c   Given a function func, and given distinct initial points ax and bx,
c   routine searches in the downhill direction (defined by the function as
c   evaluated at the initial points) and returns new points ax, bx, cx
c   that bracket a minimum of the function.  Also returned are the function
c   values at the three points, fa, fb, and fc.
c----------------------------------------------------------------------
c
      subroutine mnbrak(ax,bx,cx,fa,fb,fc,func,smin,smax,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'parms.inc'
c      include 'iterat.inc'
c
      double precision ax,bx,cx,fa,fb,fc,smin,smax,func
      external func
c
      integer ierr,idir,i
      double precision dum,fu,r,q,u,ulim,bnd
c
      double precision glimit,gold,tiny
      parameter (gold=1.618034d0,glimit=10.,tiny=1.d-20)
c
c  -- Rearrange the two points so that a to b is always "downhill"
c

      bnd=smax
      idir=1
c
      fa=func(ax,ierr)
      if (ierr.lt.0 .or. ihltcmd.ne.0) return
c loop on sites
c      do 101 i=1,ncomps
c        seteval(i)=.true.	! got eigenvalues
c 101    continue
      fb=func(bx,ierr)
      if (ierr.lt.0 .or. ihltcmd.ne.0) return
      if(fb.gt.fa)then
         bnd=smin
         idir=-1
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
      endif
c
      cx=bx+gold*(bx-ax)
c
      if ( idir*(cx-bnd).gt.0.)            cx=bnd
c
      fc=func(cx,ierr)
      if (ierr.lt.0 .or. ihltcmd.ne.0) return
c
1     if(fb.ge.fc)then
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
         if ( idir*(u-bnd).gt.0.)          u=bnd
         ulim=bx+glimit*(cx-bx)
         if ( idir*(ulim-bnd).gt.0.)       ulim=bnd
         if((bx-u)*(u-cx).gt.0.)then
c                                    ** parabolic u between b and c
            fu=func(u,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) return
            if(fu.lt.fc)then
c                                      * minimum between b and c
               ax=bx
               fa=fb
               bx=u
               fb=fu
               ierr=0
               return
            else if(fu.gt.fb)then
c                                      * minimum between a and u
               cx=u
               fc=fu
               ierr=0
               return
            endif
c                                      * default magnification
            if (cx.eq.bnd) go to 99
            u=cx+gold*(cx-bx)
            if ( idir*(u-bnd).gt.0.)       u=bnd
            fu=func(u,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) return
         else if((cx-u)*(u-ulim).gt.0.)then
c                                    ** parabolic u between c and ulim
            if (cx.eq.bnd) go to 99
            fu=func(u,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) return
            if(fu.lt.fc)then
               bx=cx
               cx=u
               u=cx+gold*(cx-bx)
               if ( idir*(u-bnd).gt.0.)    u=bnd
               fb=fc
               fc=fu
               fu=func(u,ierr)
               if (ierr.lt.0 .or. ihltcmd.ne.0) return
            endif
         else if((u-ulim)*(ulim-cx).ge.0.)then
            if (cx.eq.bnd) go to 99
            u=ulim
            fu=func(u,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) return
         else
            if (cx.eq.bnd) go to 99
            u=cx+gold*(cx-bx)
            if ( idir*(u-bnd).gt.0.)       u=bnd
            fu=func(u,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) return
         endif
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         go to 1
      endif
      ierr=0
      return
 99   ierr=1
      return
      end
