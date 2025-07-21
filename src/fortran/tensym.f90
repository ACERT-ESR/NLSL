c Version 1.3.2  5/15/94
c----------------------------------------------------------------------
c        =====================================================
c          subroutines TOCART, TOSPHR, and TOAXIL, and TCHECK
c         =====================================================
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
c
      use symdef
c
      implicit none
      double precision TWO,THREE,t(3),t1,t2,t3
      integer x,y,z,iflg
c
      parameter (TWO=2.0D0,THREE=3.0D0,x=1,y=2,z=3)
c
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
c
c- tosphr -----------------------------------------------------------------
c
      subroutine tosphr( t, iflg )
c
      use symdef
c
      implicit none
      double precision TWO,THREE,t(3),tx,ty,tz
      integer iflg
c
      parameter (TWO=2.0D0,THREE=3.0D0)
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
      t(1)=(tx+ty+tz)/THREE
      t(2)=tz-(tx+ty)/TWO
      t(3)=(tx-ty)
      return
      end
c
c- toaxil -----------------------------------------------------------------
c
      subroutine toaxil( t, iflg )
c
      use symdef
c
      implicit none
      double precision t(3)
      integer x,y,z,iflg
c
      double precision TWO,ZERO
      parameter (TWO=2.0D0,ZERO=0.0D0)
c
      if (iflg.eq.AXIAL) return
      if (iflg.eq.SPHERICAL) call tocart( t,iflg )
      t(1)=(t(1)+t(2))/TWO
      t(2)=ZERO
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
c  and .false. if an inconsistency is detected.
c
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
      function tcheck(ix,ix2,ident,lumsg)
c
      use nlsdim
      use eprprm
      use parcom
      use lpnam
      use stdio
      use symdef
c
      implicit none
      logical tcheck
c
      integer ispec,ix,ix2,ixa,ixf,jx,jx1,jx2,lu1,lu2,lumsg,mflag
      character ident*9
c
      logical symmsg
      external symmsg
c......................................................................
c
c    lu1 and lu2 are used to avoid redundant error messages when
c    ix2 specifies a range of sites 
c
      lu1=lumsg
      lu2=lumsg
c
      if (ix .gt.0) ispec=CARTESIAN
      if (ix.gt.-100 .and. ix.lt.0) ispec=SPHERICAL
      if (ix.lt.-100) ispec=AXIAL
c
      ixa=abs(mod(ix,100))
      ixf=(ixa-IWXX)/3
c
c     --- Return .true. for all parameters that have no aliases
c         and for index out-of-bounds
c
      tcheck = .true.
      if (ixa.lt.IWXX .or. ixf.gt.MXSPH-1 .or. ix2.gt.MXSITE) return
c
      if (ix2.le.0) then
         jx1=1
         jx2=MXSITE
      else
         jx1=ix2
         jx2=ix2
      end if
c
      do jx=jx1,jx2
         mflag=iparm(IIWFLG+ixf,jx)
c
c        ------------------------------
c        Check if mode is not yet set
c        ------------------------------
c
         if (mflag.eq.0) then
            if (lu1.ne.0) then
               write (lu1,1003) symstr(ispec),ident(1:1)
               lu1=0
            end if
c
c     --------------------------------------------------------------------
c     Check if tensors are specified as Cartesian when another mode is set
c     --------------------------------------------------------------------
c
         else if (ispec.eq.CARTESIAN .and.
     #           (mflag.eq.SPHERICAL .or. mflag.eq.AXIAL)) then
            tcheck=symmsg(lu2,ident,jx,symstr(mflag))
c
c     --------------------------------------------------------------------
c     Check if tensors are specified as spherical when another mode is set
c     --------------------------------------------------------------------
c
         else if (ispec.eq.SPHERICAL .and. 
     #           (mflag.eq.CARTESIAN .or. mflag.eq.AXIAL)) then
            tcheck=symmsg(lu2,ident,jx,symstr(mflag))

c     --------------------------------------------------------------------
c     Check if tensors are specified as axial when another mode is set
c     --------------------------------------------------------------------
c
         else if (ispec.eq.AXIAL .and. 
     #           (mflag.eq.CARTESIAN .or. mflag.eq.SPHERICAL)) then
            tcheck=symmsg(lu2,ident,jx,symstr(mflag))
         end if
c
c     --------------------------------------------------------------------
c     Set tensor mode flags according to type of tensor that has been 
c     specified
c     --------------------------------------------------------------------
c
         iparm(IIWFLG+ixf,jx)=ispec
c
      end do
      return
c
c ### format statements ########################################
c 
 1003 format(' *** ',a,'form assumed for ',a,' tensor ***')
      end



      function symmsg( lu, ident, jx, form )
c
      use stdio
c
      implicit none
      integer lu,jx
      character ident*9,form*10
      logical symmsg
c
      integer itrim
      external itrim
c
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
c
      use rnddbl
      use symdef
c
      implicit none
      double precision t(3)
      integer iflg
      logical isaxial
c
      isaxial = (iflg.eq.AXIAL) .or.
     #          (iflg.eq.SPHERICAL .and. abs(t(3)).lt.RNDOFF) .or.
     #          (iflg.eq.CARTESIAN .and. abs(t(1)-t(2)).lt.RNDOFF)   
      return
      end


      function isrhomb( t, iflg )
c
      use rnddbl
      use symdef
c
      implicit none
      double precision t(3)
      integer iflg
      logical isrhomb
c
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
c
      use rnddbl
      use symdef
c
      implicit none
      double precision t(3)
      integer iflg
      logical isisotr
c
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

