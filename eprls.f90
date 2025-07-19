c NLSL Version 1.5.1 beta 11/25/95
c*********************************************************************
c                       ==================
c                        subroutine:EPRLS
c                       ==================
c
c>@brief      Subroutine version of EPRLL family of programs by D. Schneider.
c> This routine is intended for use with nonlinear least-squares
c> applications. The routine calculates the tridiagonal matrix
c> for the parameters passed in the fparm and iparm arrays using
c> the Conjugate Gradients version of the Lanczos algorithm.
c>
c> Calculation parameters are input in the arrays fparm and iparm.
c> They are expected in the same order as that found in common 
c> /eprprm/. The input flag icalc is specified as a three-digit
c> number ijk, where
c>
c>                 j.ne.0  => Perform matrix and CG calculations (matrll/cscg) 
c>                 k.ne.0  => Perform starting vector and CG calculations (stvect/cscg)
c>
c> The continued-fractions calculation is always performed.
c
c         Subroutine arguments are output as follows:
c
c>          @param  al(ndone)    Diagonal of tridiagonal matrix for spectrum
c>          @param  be(ndone)    Off-diagonal of tridiagonal matrix for spectrum
c>          @param  ndone        Number of CG steps taken. Zero for Lanczos
c>                             calculations, and set negative if CG did
c>                             not converge
c>          @param  ierr         Error flag. Meaning of values assigned to ierr
c>                         are defined in 'errmsg.inc'.
c>
c>      Written by DEB, 26-Sep-91 based on programs from DJS 
c
c       Uses :
c               pmatrl.f     Liouville matrix calculation (pruning version)
c               pstvec.f     Starting vector calculation (pruning version)
c               cscg.f       Conjugate gradients tridiagonalization
c               cfs.f        Continued-fraction spectral calculation
c
c*********************************************************************
c
      subroutine eprls( icalc,al,be,bss,iprune,spectr,nft,ndone,ierr)
c
      use nlsdim
      use eprprm
      use eprmat
      use tridag
      use ftwork
c      use prmeqv
      use errmsg
      use pidef
      use stdio
c
c     The parameter TOLEXP is used to determine the smallest number for
c     which the Gaussian convolution function needs to be calculated,
c     exp(-2*TOLEXP).  
c
      implicit none
      double complex CZERO
      parameter (CZERO=(0.0D0,0.0D0))
c
      integer bss(5,MXDIM),ndone,icalc,ierr,iprune,nft
c     double precision spectr(iepr(INFLD))
      double precision spectr(*)
      double complex al(MXDIM),be(MXDIM)
      logical matrix,start
c
      double precision cgerr,wline,gib
      integer i,m,n,no2,nstp
c
      double precision dblint
      external dblint
c
c#####################################################################
c
      matrix = mod(icalc,100)/10 .ne. 0
      start  = mod(icalc,10) .ne. 0
c
      if (matrix) then
         if (iprune.ne.0) then
            call pmatrl(bss,ierr)
         else
            call matrll(ierr)
         end if
c
         call wpoll
         if (hltcmd.ne.0 .or. hltfit.ne.0 .or. ierr.gt.FATAL) return
      end if
c
c----------------------------------------------------------------------
c        Generate starting vector
c----------------------------------------------------------------------
      if (start) then
         if (iprune.ne.0) then
            call pstvec(bss,stv,ierr)
         else
            call stvect(stv,ierr)
         endif
      end if
c
      if (start .or. matrix) then
c
c---------------------------------------------------------------------
c       Call conjugate gradients tridiagonalization routine
c---------------------------------------------------------------------
c
         do i=1,ndim
            y(i)=CZERO
         end do
c
c       -----------------------------------------------------
c        Set the number of steps for CG tridiagonalization
c        If nstep is unspecified, use the matrix dimension or
c        the number of steps previously taken, whichever is
c        smaller
c       -----------------------------------------------------
         nstp=nstep
         if (nstep.le.0) nstp=min(ndim,ndone)
c
         call cscg( stv,ndim,nstp,cgtol,dcmplx(shiftr,shifti),
     #              y,al,be,ndone,cgerr )
c
c         Check for user input or halt before proceeding
c
         call wpoll
         if (hltcmd.ne.0 .or. hltfit.ne.0) then
            ierr=CGHLT
            return
         end if
c
c                                        *** CG did not converge
         if (ndone.lt.0) then
            ndone=-ndone
            ierr=NOCONVRG
            write(eprerr(ierr)(27:31),'(i5)') ndone
         end if
c
      end if
c
c     ---------------------------------------------
c     Continued-fraction calculation of spectrum
c     ---------------------------------------------
      call cfs( al,be,ndone,b0-fldi,-dfld,nfld,w0+lb,ideriv,phase,
     #          spectr )
c
c     ------------------------------------
c     Convolution with Gaussian lineshape
c     ------------------------------------
      gib=gib0+gib2*dsin( RADIAN*psi )**2
      wline=dsqrt(gib*gib+lb*lb)
      call gconvl(spectr,wline,dfld,nfld,nft)
c
c     ------------------
c     Normalize spectrum (PI comes from Lorentzian shape)
c     ------------------
      do i=1,nfld
         spectr(i)=spectr(i)/PI
      end do
c
      return
c
      end
