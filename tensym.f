c NLSPMC VERSION (VERSION 1.0) 2/5/99
c----------------------------------------------------------------------
c            =====================================================
c             subroutines TOCART, TOSPHR, and TOAXIL, and TCHECK
c            =====================================================
c
c These subroutines interconvert among the three possible representations
c of the various tensors (in 3-space, including magnetic and diffusion
c tensors) used in EPR calculations, namely:
c   (1) Cartesian representation
c   (2) Spherical representation
c   (3) Axial (actually also in Cartesian coordinates)
c
c The Cartesian representation is the familiar set of (x,y,z) components.
c If the Cartesian components of a tensor M are Mx, My, and Mz, then
c the "spherical" components M1, M2, and M3 are given by
c
c   M1=(Mx+My+Mz)/3   (isotropic component)
c   M2=Mz-(Mx+My)/2   (axial component)
c   M3=Mx-My          (rhombic component)
c
c Note that these components differ from the "true" second-order
c spherical tensor components, M(0,0), M(2,0), and M(2,2), respectively,
c by constant factors. They are expressed in this form so that they may
c be more directly correlated with experimental ESR spectra. These 
c transformations are carried out by the TOSPHR routine.
c
c The third representation, "axial" is really in Cartesian coordinates
c as well, but affords the user the option of maintaining axial tensor
c symmetry while varying the Cartesian components.
c
c The transformation from spherical to cartesian components, the inverse
c of that given above, is carried out by TOCART and given by
c
c   Mx=M1 - M2/3 + M3/2
c   My=M1 - M2/3 - M3/2
c   Mz=M1 + 2*M2/3
c----------------------------------------------------------------------
      subroutine tocart( t, iflg )
      implicit none
      double precision two,three,t(3),t1,t2,t3
      integer x,y,z,iflg,CARTESIAN,SPHERICAL,AXIAL
c
      parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
      parameter (two=2.0D0,three=3.0D0,x=1,y=2,z=3)
c
      if (iflg.eq.CARTESIAN) return
      if (iflg.eq.AXIAL) then
         t1=t(x)
         t2=t(y)
         t3=t(z)
         t(x)=t1
         t(y)=t1
         t(z)=t3
      else if (iflg.eq.SPHERICAL) then
         t1=t(x)
         t2=t(y)
         t3=t(z)
         t(x)=t1 - t2/three + t3/two
         t(y)=t1 - t2/three - t3/two
         t(z)=t1 + two*t2/three
      end if
      return
      end


      subroutine tosphr( t, iflg )
      implicit none
      double precision two,three,t(3),tx,ty,tz,t1,t2,t3
      integer iflg,CARTESIAN,SPHERICAL,AXIAL
c
      parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
      parameter (two=2.0D0,three=3.0D0)
c
      if (iflg.eq.SPHERICAL) return
      if (iflg.eq.AXIAL) then
         tx=t(1)
         ty=t(1)
         tz=t(3)
      else if (iflg.eq.CARTESIAN) then
         tx=t(1)
         ty=t(2)
         tz=t(3)
      end if
      t(1)=(tx+ty+tz)/three
      t(2)=tz-(tx+ty)/two
      t(3)=(tx-ty)
      return
      end


      subroutine toaxil( t, iflg )
      implicit none
      include 'stdio.inc'
      include 'rndoff.inc'
      double precision two,t(3),tmp
      integer x,y,z,iflg,CARTESIAN,SPHERICAL,AXIAL
c
      parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
      parameter (two=2.0D0)
c
      if (iflg.eq.AXIAL) return
      if (iflg.eq.SPHERICAL) call tocart( t,iflg )
      if (abs(t(1)-t(2)) .gt. rndoff) write(luout,*)'caution ',
     #	'information lost, tensor was not axial before conversion'
      tmp=t(2)
      t(1)=(t(1)+t(2))/two
      t(2)=0.0D0
      return
      end
c
c----------------------------------------------------------------------
c                     =========================
c                        function TCHECK
c                     =========================
c
c  Checks whether the component of the g, A, or R (T) tensors specified
c  by the ix index is consistent with the previous setting of the tensor
c  mode flags. (igflg, iaflg, and irflg) TCHECK returns .true. if
c  it is consistent (or if ix does not specify one of these tensors)
c  and .false. if an inconsistency is detected.  We require all sites 
c  and spectra to have the same mode.  
c  If a nonzero logical unit number is specified, any error messages will
c  be directed to that unit.
c
c         ix specifies tensor symmetry in the following way:
c             ix < -100: axial mode
c      -100 < ix < 0   : spherical mode
c             ix > 0   : cartesian mode     
c 
c         Mode flag =  0 indicates that mode has not yet been set
c                      1 indicates that Cartesian mode has been set
c                      2 indicates that spherical mode has been set
c                      3 indicates that axial mode has been set
c
c----------------------------------------------------------------------
c
      function tcheck(ix,token,lumsg)
      implicit none
      logical tcheck
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'parms.inc'
      include 'lpnam.inc'
      include 'miscel.inc'
      include 'stdio.inc'
c
      integer ispect,ix,ixa,ixf,lumsg,mflag,CARTESIAN,SPHERICAL,AXIAL
      integer isite,ispec
      parameter(CARTESIAN=1,SPHERICAL=2,AXIAL=3)
c
      character token*(*)
c
c......................................................................
c
      if (ix.gt.0) ispect=CARTESIAN
      if (ix.lt.-100) ispect=AXIAL
      if (ix.gt.-100 .and. ix.lt.0) ispect=SPHERICAL
c
      ixa=abs(mod(ix,100))
      ixf=(ixa-IGXX)/3		! 0,1 or 2
c
c     --- Return .true. for all parameters that have no aliases
c
      if (ixa.lt.IGXX .or. ixf.gt.mxsph-1) then
        tcheck=.true.
        return
      end if
c
      mflag=iparm(IIGFLG+ixf,1,1)    ! get the flag for this parameter
c
c Check if mode is not yet set
c
      if (mflag.eq.0) then
         if (lumsg.ne.0) write (lumsg,1003) symstr(ispect),token(1:1)
c
c Check if tensors are specified as Cartesian when another mode is set
c
      else if (ispect.eq.CARTESIAN .and.
     #        (mflag.eq.SPHERICAL .or. mflag.eq.AXIAL)) then
         if (lumsg.ne.0) write(lumsg,1004) token,symstr(mflag)
         tcheck=.false.
         return
c
c Check if tensors are specified as spherical when another mode is set
c
      else if (ispect.eq.SPHERICAL .and. 
     #        (mflag.eq.CARTESIAN .or. mflag.eq.AXIAL)) then
          if (lumsg.ne.0) write(lumsg,1004) token,symstr(mflag)
          tcheck=.false.
          return
c
c Check if tensors are specified as axial when another mode is set
c
      else if (ispect.eq.AXIAL .and. 
     #         (mflag.eq.CARTESIAN .or. mflag.eq.SPHERICAL)) then
         if (lumsg.ne.0) write(lumsg,1004) token,symstr(mflag)
         tcheck=.false.
         return
      end if
c
c  Set tensor mode flags according to type of tensor that has been 
c  specified in all sites and spectra.
c
      do 30 isite=1,ncomps
        do 30 ispec=1,nspectra
          iparm(IIGFLG+ixf,ispec,isite)=ispect
 30   continue
      tcheck=.true.
      return
c
c ### format statements ########################################
c 
 1003 format(' *** ',a,' form assumed for ',a,' tensor ***')
 1004 format(' *** Cannot modify ''',a,''': ',a,' form',
     #       ' has been specified ***')
      end
