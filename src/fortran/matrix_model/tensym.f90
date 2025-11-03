! Version 1.3.2  5/15/94
!----------------------------------------------------------------------
!        =====================================================
!          subroutines TOCART, TOSPHR, and TOAXIL, and TCHECK
!         =====================================================
!
! These subroutines interconvert among the three possible representations
! of the various tensors (in 3-space, including magnetic and diffusion
! tensors) used in EPR calculations, namely:
!   (1) Cartesian representation
!   (2) Spherical representation
!   (3) Axial (actually also in Cartesian coordinates)
!
! The Cartesian representation is the familiar set of (x,y,z) components.
! If the Cartesian components of a tensor M are Mx, My, and Mz, then
! the "spherical" components M1, M2, and M3 are given by
!
!   M1=(Mx+My+Mz)/3   (isotropic component)
!   M2=Mz-(Mx+My)/2   (axial component)
!   M3=Mx-My          (rhombic component)
!
! Note that these components differ from the "true" second-order
! spherical tensor components, M(0,0), M(2,0), and M(2,2), respectively,
! by constant factors. They are expressed in this form so that they may
! be more directly correlated with experimental ESR spectra. These 
! transformations are carried out by the TOSPHR routine.
!
! The third representation, "axial" is really in Cartesian coordinates
! as well, but affords the user the option of maintaining axial tensor
! symmetry while varying the Cartesian components.
!
! The transformation from spherical to cartesian components, the inverse
! of that given above, is carried out by TOCART and given by
!
!   Mx=M1 - M2/3 + M3/2
!   My=M1 - M2/3 - M3/2
!   Mz=M1 + 2*M2/3
!----------------------------------------------------------------------
      subroutine tocart( t, iflg )
!
      use symdef
!
      implicit none
      double precision TWO,THREE,t(3),t1,t2,t3
      integer iflg
      integer, parameter :: x=1,y=2,z=3
!
      parameter (TWO=2.0D0,THREE=3.0D0)
!
      if (iflg.eq.CARTESIAN) return
      if (iflg.eq.AXIAL) then
         t(y)=t(x)
      else if (iflg.eq.SPHERICAL) then
         t1=t(x)
         t2=t(y)
         t3=t(z)
         t(x)=t1 - t2/THREE + t3/TWO
         t(y)=t1 - t2/THREE - t3/TWO
         t(z)=t1 + TWO*t2/THREE
      end if
      return
      end
!
!- tosphr -----------------------------------------------------------------
!
      subroutine tosphr( t, iflg )
!
      use symdef
!
      implicit none
      double precision TWO,THREE,t(3),tx,ty,tz
      integer iflg
!
      parameter (TWO=2.0D0,THREE=3.0D0)
!
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
      t(1)=(tx+ty+tz)/THREE
      t(2)=tz-(tx+ty)/TWO
      t(3)=(tx-ty)
      return
      end
!
!- toaxil -----------------------------------------------------------------
!
      subroutine toaxil( t, iflg )
!
      use symdef
!
      implicit none
      double precision t(3)
      integer x,y,z,iflg
!
      double precision TWO,ZERO
      parameter (TWO=2.0D0,ZERO=0.0D0)
!
      if (iflg.eq.AXIAL) return
      if (iflg.eq.SPHERICAL) call tocart( t,iflg )
      t(1)=(t(1)+t(2))/TWO
      t(2)=ZERO
      return
      end
!
!----------------------------------------------------------------------
!                     =========================
!                        function TCHECK
!                     =========================
!
!  Checks whether the component of the g, A, or R (T) tensors specified
!  by the ix index is consistent with the previous setting of the tensor
!  mode flags. (igflg, iaflg, and irflg) TCHECK returns .true. if
!  it is consistent (or if ix does not specify one of these tensors)
!  and .false. if an inconsistency is detected.
!
!  If a nonzero logical unit number is specified, any error messages will
!  be directed to that unit.
!
!         ix specifies tensor symmetry in the following way:
!             ix < -100: axial mode
!      -100 < ix < 0   : spherical mode
!             ix > 0   : cartesian mode     
! 
!         Mode flag =  0 indicates that mode has not yet been set
!                      1 indicates that Cartesian mode has been set
!                      2 indicates that spherical mode has been set
!                      3 indicates that axial mode has been set
!
!----------------------------------------------------------------------
!
      function tcheck(ix,ix2,ident,lumsg)
!
      use nlsdim
      use eprprm
      use parcom
      use lpnam
      use stdio
      use symdef
!
      implicit none
      logical tcheck
!
      integer ispec,ix,ix2,ixa,ixf,jx,jx1,jx2,lu1,lu2,lumsg,mflag
      character ident*9
!
      logical symmsg
      external symmsg
!......................................................................
!
!    lu1 and lu2 are used to avoid redundant error messages when
!    ix2 specifies a range of sites 
!
      lu1=lumsg
      lu2=lumsg
!
      if (ix .gt.0) ispec=CARTESIAN
      if (ix.gt.-100 .and. ix.lt.0) ispec=SPHERICAL
      if (ix.lt.-100) ispec=AXIAL
!
      ixa=abs(mod(ix,100))
      ixf=(ixa-IWXX)/3
!
!     --- Return .true. for all parameters that have no aliases
!         and for index out-of-bounds
!
      tcheck = .true.
      if (ixa.lt.IWXX .or. ixf.gt.MXSPH-1 .or. ix2.gt.MXSITE) return
!
      if (ix2.le.0) then
         jx1=1
         jx2=MXSITE
      else
         jx1=ix2
         jx2=ix2
      end if
!
      do jx=jx1,jx2
         mflag=iparm(IIWFLG+ixf,jx)
!
!        ------------------------------
!        Check if mode is not yet set
!        ------------------------------
!
         if (mflag.eq.0) then
            if (lu1.ne.0) then
               write (lu1,1003) symstr(ispec),ident(1:1)
               lu1=0
            end if
!
!     --------------------------------------------------------------------
!     Check if tensors are specified as Cartesian when another mode is set
!     --------------------------------------------------------------------
!
         else if (ispec.eq.CARTESIAN .and.
     #           (mflag.eq.SPHERICAL .or. mflag.eq.AXIAL)) then
            tcheck=symmsg(lu2,ident,jx,symstr(mflag))
!
!     --------------------------------------------------------------------
!     Check if tensors are specified as spherical when another mode is set
!     --------------------------------------------------------------------
!
         else if (ispec.eq.SPHERICAL .and. 
     #           (mflag.eq.CARTESIAN .or. mflag.eq.AXIAL)) then
            tcheck=symmsg(lu2,ident,jx,symstr(mflag))

!     --------------------------------------------------------------------
!     Check if tensors are specified as axial when another mode is set
!     --------------------------------------------------------------------
!
         else if (ispec.eq.AXIAL .and. 
     #           (mflag.eq.CARTESIAN .or. mflag.eq.SPHERICAL)) then
            tcheck=symmsg(lu2,ident,jx,symstr(mflag))
         end if
!
!     --------------------------------------------------------------------
!     Set tensor mode flags according to type of tensor that has been 
!     specified
!     --------------------------------------------------------------------
!
         iparm(IIWFLG+ixf,jx)=ispec
!
      end do
      return
!
! ### format statements ########################################
! 
 1003 format(' *** ',a,'form assumed for ',a,' tensor ***')
      end



      function symmsg( lu, ident, jx, form )
!
      use stdio
!
      implicit none
      integer lu,jx
      character ident*9,form*10
      logical symmsg
!
      integer itrim
      external itrim
!
      if (lu.ne.0) then
         write (lu,1004) ident(:itrim(ident)),jx,form
         if (lu.ne.luttyo) write (luttyo,1004) ident(:itrim(ident)),
     #                                         jx,form
      end if
      lu=0
      symmsg=.false.
      return
 1004 format(' *** Cannot modify ''',a,'(',i1,')'': ',a,'form',
     #       ' has been specified ***')
      end



      function isaxial( t, iflg )
!
      use rnddbl
      use symdef
!
      implicit none
      double precision t(3)
      integer iflg
      logical isaxial
!
      isaxial = (iflg.eq.AXIAL) .or.
     #          (iflg.eq.SPHERICAL .and. abs(t(3)).lt.RNDOFF) .or.
     #          (iflg.eq.CARTESIAN .and. abs(t(1)-t(2)).lt.RNDOFF)   
      return
      end


      function isrhomb( t, iflg )
!
      use rnddbl
      use symdef
!
      implicit none
      double precision t(3)
      integer iflg
      logical isrhomb
!
      isrhomb = (iflg.eq.CARTESIAN .and. 
     #             abs(t(1)-t(2)).lt.RNDOFF .and.
     #             abs(t(2)-t(3)).lt.RNDOFF .and. 
     #             abs(t(1)-t(3)).lt.RNDOFF       ) .or.
     #          (iflg.eq.SPHERICAL .and. 
     #             abs(t(2)).gt.RNDOFF .and.
     #             abs(t(3)).gt.RNDOFF      )
      return
      end

      function isisotr( t, iflg )
!
      use rnddbl
      use symdef
!
      implicit none
      double precision t(3)
      integer iflg
      logical isisotr
!
      isisotr = (iflg.eq.CARTESIAN .and. 
     #              abs(t(1)-t(2)).lt.RNDOFF .and.
     #              abs(t(2)-t(3)).lt.RNDOFF      ) .or.
     #          (iflg.eq.SPHERICAL .and.
     #              abs(t(2)).lt.RNDOFF .and.
     #              abs(t(3)).lt.RNDOFF           ) .or.
     #          (iflg.eq.AXIAL .and. 
     #              abs(t(1)-t(3)).lt.RNDOFF)
      return
      end

