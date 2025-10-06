! NLSL Version 1.5 beta 11/23/95

!*********************************************************************
!
!                       ==================
!                        subroutine:MOMDLS
!                       ==================
!
!>@brief Subroutine version of EPRLL family of programs.
!> This routine is intended for use with nonlinear least-squares
!> applications. The routine calculates a MOMD spectrum by
!> averaging spectra for a specified number of orientation using
!> the parameters passed in the fparm and iparm arrays. Spectra
!> are calculated using the Conjugate Gradients version of the Lanczos
!> algorithm.
!>
!> If the parameter specifying the number of MOMD orientations
!> is not greater than 1, or if there is no potential defined, 
!> a single spectral calculation is carried out.
!
!          Subroutine arguments are exactly output as follows:
!
!>@param al(mxdim)    Diagonal of tridiagonal matrix for spectrum
!>@param be(mxdim)    Off-diagonal of tridiagonal matrix for spectrum
!>
!>@param cgerr          Residual error in CG solution vector
!>
!>@param ntotal         Number of CG steps taken (total number for
!> all orientations in MOMD calculations)
!>
!>@param ierr           Error flag
!                            Return values are defined in 'errmsg.f'
!
!
!>@author Written by DEB, 26-Sep-91 based on programs from DJS 
!
!       Includes :
!               nlsdim.inc
!               eprmat.inc
!               eprprm.inc
!               pidef.inc
!               stdio.inc
!               rndoff.inc
!
!       Uses :
!               eprls.f
!
!*********************************************************************
!
      subroutine momdls( fparmi,iparmi,icalc,al,be,bss,iprune,
     #                   spectr,work,nft,ntotal,ierr )
!
      use eprprm
      use nlsdim
      use eprmat
      use tridag
!      use prmeqv
      use errmsg
      use pidef
      use stdio
      use rnddbl
      use dfunc
!
      implicit none
      integer bss(5,MXDIM),iparmi(NIPRM),ntotal,icalc,ierr,iprune,
     #        nft
!     double precision fparmi(NFPRM),spectr(iparmi(INFLD)),
!    #                 work(iparmi(INFLD)),cgerr
      double precision fparmi(NFPRM),spectr(*),
     #                 work(*)
      double complex al(ntotal),be(ntotal)
!
      integer i,id,ixt,j,nptr,nstp,nevl
      double precision acc,cospsi,cosxi,dcos,dnorm,dwt,onorm,sinxi
      logical init
!
      double precision ZERO,ONE,SMALL,D180
      parameter(ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-16,D180=180.0D0)
!
      double precision fu20,fu20phi
      external fu20,fu20phi

      cosxi=ZERO
      sinxi=ZERO
!
!#####################################################################
!
      init = mod(icalc,1000)/100 .ne. 0
!
!     ----------------------------------------------
!     Load values from parameter array into /eprprm/
!     ----------------------------------------------
!
!     *** In the future this could be done by select_site(isi) ***
!      call select_site(isi)
!
      do i=1,NFPRM
         fepr(i)=fparmi(i)
      end do

      do i=1,NIPRM
         iepr(i)=iparmi(i)
      end do
!
      if (init) then
!                                     *** Fatal error in parameters
         call lcheck(ierr)
         if (ierr.ge.FATAL) return
!
!     -------------------------------------------------------
!       If lcheck was not called, we still need to:
!         (1) convert W tensor to Cartesian form and find w0
!         (2) count potential coefficients
!         (3) set director tilt flag
!     -------------------------------------------------------
      else
         call tocart( fepr(IWXX), iwflg )
         w0=(wxx+wyy+wzz)/3.0d0
         ipt=0
         if (dabs(c20).gt.RNDOFF) ipt=ipt+1
         if (dabs(c22).gt.RNDOFF) ipt=ipt+1
         if (dabs(c40).gt.RNDOFF) ipt=ipt+1
         if (dabs(c42).gt.RNDOFF) ipt=ipt+1
         if (dabs(c44).gt.RNDOFF) ipt=ipt+1
         ipsi0=0
         if( ((abs(psi).gt.RNDOFF).and.(abs(psi)-D180.gt.RNDOFF))
     #   .or. nort.gt.1) ipsi0=1
      end if
!
      do j=1,nfld
         spectr(j)=ZERO
      end do
!
      if (nort.le.1 .or. ipt.le.0) then
!
!       --------------------------
!       --- Single orientation ---
!       --------------------------
!
         call eprls( icalc,al,be,bss,iprune,spectr,nft,ntotal,ierr)
!
!      --------------------------
!      ---- MOMD calculation ----
!      --------------------------
!
      else
         ixt=1
         nptr=0
         onorm=ONE/dfloat(nort-1)
!
!        ------------------------------------------------------
!        Partial MOMD
!        If there is ordering of the director, calculate 
!        normalization factor for the pseudopotential function
!        In the case of partial director ordering, the parameter
!        psi actually represents the angle between the director 
!        ordering axis and the field, whereas the calculation
!        subroutine always interprets psi as the director tilt.
!        Save the given psi value as the angle xi. 
!        ------------------------------------------------------
!
         if (abs(dc20).gt.RNDOFF) then
            xi=abs(mod(psi,D180))
            cosxi=cos(RADIAN*xi)
            sinxi=sin(RADIAN*xi)
!
            acc=1.0d-4
            dlam=dc20
            call ccrint(ZERO,ONE,acc,SMALL,dnorm,nevl,fu20,id)
         else
            dlam=ZERO
            dnorm=ONE
         end if
!
!        -----------------------
!        Loop over orientations
!        -----------------------
!
         cospsi=ZERO
         dcos=ONE/(nort-1)
         
         do i=1,nort
            cospsi=min(cospsi,ONE)
            fepr(IPSI)=acos(cospsi)/RADIAN
!
            if (icalc.ne.0) then
               nstp=nstep
               if (nstep.le.0) nstp=min0(ndim,ntotal/nort)
!
!           --------------------------------------------------
!           Repeat continued-fractions for this orientation
!           with the tridiagonal matrix: find number of steps
!           by locating next zero in beta array
!           --------------------------------------------------
!
            else
               nstp=0
 3             nstp=nstp+1
               nptr=nptr+1
               if (nptr.gt.ntotal) then 
                  ierr=TDGERR
                  return
               end if
               if ( abs(be(nptr)).gt.RNDOFF ) go to 3
            end if
!
            call eprls( icalc,al(ixt),be(ixt),bss,iprune,work,
     #                  nft,nstp,ierr )
            if (ierr.ge.FATAL) return
            ixt=ixt+nstp
!
!           --------------------------------------------------
!           If there is ordering of the director, calculate the
!           weighting coefficient for this tilt angle and
!           weight the spectrum
!           --------------------------------------------------
!
            if (abs(dc20).gt.RNDOFF) then
               if (xi.gt.RNDOFF) then
                  ss=sinxi*sqrt(ONE-cospsi*cospsi)
                  cc=cosxi*cospsi
                  call ccrint(ZERO,PI,acc,SMALL,dwt,nevl,fu20phi,id)
                  dwt=dwt/(dnorm*PI)
               else
                  dwt=fu20(cospsi)/dnorm
               end if
!                  
               do j=1,nfld
                  work(j)=work(j)*dwt
               end do
            end if
!
!          --------------------------------------------------------
!          Trapezoidal rule: first and last pts. get factor of 1/2
!          --------------------------------------------------------
!
            if (i.eq.1 .or. i.eq.nort) then
               do j=1,nfld
                  work(j)=work(j)*0.5d0
               end do
            end if
!
!        ------------------------------------------------
!        Add spectrum for current orientation into total 
!        ------------------------------------------------
!
            do j=1,nfld
               spectr(j)=spectr(j)+work(j)
            end do
!
            cospsi=cospsi+dcos
         end do
!
!        -------------------
!        Normalize spectrum
!        -------------------
!
         do j=1,nfld
            spectr(j)=spectr(j)*onorm
         end do
!
!        -------------------------------------------
!        Return total number of steps actually taken
!        -------------------------------------------
! 
         if (icalc.ne.0) ntotal=ixt-1
!
      end if
!
      return
      end

