! NLSL Version 1.5.1 beta 11/25/95
!*********************************************************************
!                       ==================
!                        subroutine:EPRLS
!                       ==================
!
!>@brief      Subroutine version of EPRLL family of programs by D. Schneider.
!> This routine is intended for use with nonlinear least-squares
!> applications. The routine calculates the tridiagonal matrix
!> for the parameters passed in the fparm and iparm arrays using
!> the Conjugate Gradients version of the Lanczos algorithm.
!>
!> Calculation parameters are input in the arrays fparm and iparm.
!> They are expected in the same order as that found in common 
!> /eprprm/. The input flag icalc is specified as a three-digit
!> number ijk, where
!>
!>                 j.ne.0  => Perform matrix and CG calculations (matrll/cscg) 
!>                 k.ne.0  => Perform starting vector and CG calculations (stvect/cscg)
!>
!> The continued-fractions calculation is always performed.
!
!         Subroutine arguments are output as follows:
!
!>          @param  al(ndone)    Diagonal of tridiagonal matrix for spectrum
!>          @param  be(ndone)    Off-diagonal of tridiagonal matrix for spectrum
!>          @param  ndone        Number of CG steps taken. Zero for Lanczos
!>                             calculations, and set negative if CG did
!>                             not converge
!>          @param  ierr         Error flag. Meaning of values assigned to ierr
!>                         are defined in 'errmsg.inc'.
!>
!>      Written by DEB, 26-Sep-91 based on programs from DJS 
!
!       Uses :
!               pmatrl.f     Liouville matrix calculation (pruning version)
!               pstvec.f     Starting vector calculation (pruning version)
!               cscg.f       Conjugate gradients tridiagonalization
!               cfs.f        Continued-fraction spectral calculation
!
!*********************************************************************
!
      subroutine eprls( icalc,al,be,bss,iprune,spectr,nft,ndone,ierr)
!
      use nlsdim
      use eprprm
      use eprmat
      use tridag
      use ftwork
!      use prmeqv
      use errmsg
      use pidef
      use stdio
!
!     The parameter TOLEXP is used to determine the smallest number for
!     which the Gaussian convolution function needs to be calculated,
!     exp(-2*TOLEXP).  
!
      implicit none
      double precision EIGHT,ONE,THIRD,TWO,ZERO,EPS
      double complex CZERO
      parameter (ZERO=0.0D0,ONE=1.0D0,THIRD=0.333333333333333D0,
     #           TWO=2.0D0,EIGHT=8.0D0,EPS=1.0D-3,CZERO=(0.0D0,0.0D0))
!
      integer bss(5,MXDIM),ndone,icalc,ierr,iprune,nft
!     double precision spectr(iepr(INFLD))
      double precision spectr(*)
      double complex al(MXDIM),be(MXDIM)
      logical matrix,start
!
      double precision cgerr,wline,df,f,g,gib,gnorm
      integer i,m,n,no2,nstp
!
      double precision dblint
      external dblint
!
!#####################################################################
!
      matrix = mod(icalc,100)/10 .ne. 0
      start  = mod(icalc,10) .ne. 0
!
      if (matrix) then
         if (iprune.ne.0) then
            call pmatrl(bss,ierr)
         else
            call matrll(ierr)
         end if
!
         call wpoll
         if (hltcmd.ne.0 .or. hltfit.ne.0 .or. ierr.gt.FATAL) return
      end if
!
!----------------------------------------------------------------------
!        Generate starting vector
!----------------------------------------------------------------------
      if (start) then
         if (iprune.ne.0) then
            call pstvec(bss,stv,ierr)
         else
            call stvect(stv,ierr)
         endif
      end if
!
      if (start .or. matrix) then
!
!---------------------------------------------------------------------
!       Call conjugate gradients tridiagonalization routine
!---------------------------------------------------------------------
!
         do i=1,ndim
            y(i)=CZERO
         end do
!
!       -----------------------------------------------------
!        Set the number of steps for CG tridiagonalization
!        If nstep is unspecified, use the matrix dimension or
!        the number of steps previously taken, whichever is
!        smaller
!       -----------------------------------------------------
         nstp=nstep
         if (nstep.le.0) nstp=min(ndim,ndone)
!
         call cscg( stv,ndim,nstp,cgtol,dcmplx(shiftr,shifti),
     #              y,al,be,ndone,cgerr )
!
!         Check for user input or halt before proceeding
!
         call wpoll
         if (hltcmd.ne.0 .or. hltfit.ne.0) then
            ierr=CGHLT
            return
         end if
!
!                                        *** CG did not converge
         if (ndone.lt.0) then
            ndone=-ndone
            ierr=NOCONVRG
            write(eprerr(ierr)(27:31),'(i5)') ndone
         end if
!
      end if
!
!     ---------------------------------------------
!     Continued-fraction calculation of spectrum
!     ---------------------------------------------
      call cfs( al,be,ndone,b0-fldi,-dfld,nfld,w0+lb,ideriv,phase,
     #          spectr )
!
!     ------------------------------------
!     Convolution with Gaussian lineshape
!     ------------------------------------
      gib=gib0+gib2*dsin( RADIAN*psi )**2
      wline=dsqrt(gib*gib+lb*lb)
      call gconvl(spectr,wline,dfld,nfld,nft)
!
!     ------------------
!     Normalize spectrum (PI comes from Lorentzian shape)
!     ------------------
      do i=1,nfld
         spectr(i)=spectr(i)/PI
      end do
!
      return
!
      end
