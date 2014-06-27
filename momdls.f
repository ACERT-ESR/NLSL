c NLSL Version 1.5 beta 11/23/95

c*********************************************************************
c
c                       ==================
c                        subroutine:MOMDLS
c                       ==================
c
c>@brief Subroutine version of EPRLL family of programs.
c> This routine is intended for use with nonlinear least-squares
c> applications. The routine calculates a MOMD spectrum by
c> averaging spectra for a specified number of orientation using
c> the parameters passed in the fparm and iparm arrays. Spectra
c> are calculated using the Conjugate Gradients version of the Lanczos
c> algorithm.
c>
c> If the parameter specifying the number of MOMD orientations
c> is not greater than 1, or if there is no potential defined, 
c> a single spectral calculation is carried out.
c
c          Subroutine arguments are exactly output as follows:
c
c>@param al(mxdim)    Diagonal of tridiagonal matrix for spectrum
c>@param be(mxdim)    Off-diagonal of tridiagonal matrix for spectrum
c>
c>@param cgerr          Residual error in CG solution vector
c>
c>@param ntotal         Number of CG steps taken (total number for
c> all orientations in MOMD calculations)
c>
c>@param ierr           Error flag
c                            Return values are defined in 'errmsg.f'
c
c
c>@author Written by DEB, 26-Sep-91 based on programs from DJS 
c
c       Includes :
c               nlsdim.inc
c               eprmat.inc
c               eprprm.inc
c               pidef.inc
c               stdio.inc
c               rndoff.inc
c
c       Uses :
c               eprls.f
c
c*********************************************************************
c
      subroutine momdls( fparm,iparm,icalc,al,be,bss,iprune,spectr,
     #                   work,nft,ntotal,ierr )
      implicit none
c
      include 'eprprm.inc'
      include 'nlsdim.inc'
c
      integer bss(5,MXDIM),iparm(NIPRM),ntotal,icalc,ierr,iprune,
     #        nft
c     double precision fparm(NFPRM),spectr(iparm(INFLD)),
c    #                 work(iparm(INFLD)),cgerr
      double precision fparm(NFPRM),spectr(*),
     #                 work(*),cgerr
      double complex al(ntotal),be(ntotal)
c
      include 'eprmat.inc'
      include 'tridag.inc'
      include 'prmeqv.inc'
      include 'errmsg.inc'
      include 'pidef.inc'
      include 'stdio.inc'
      include 'rndoff.inc'
      include 'dfunc.inc'
c
      integer i,id,ixt,j,nptr,nstp,nevl
      double precision acc,cospsi,cosxi,dcos,dnorm,dwt,onorm,sinxi
      logical init
c
      double precision ZERO,ONE,SMALL,D180
      parameter(ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-16,D180=180.0D0)
c
      double precision fu20,fu20phi
      external fu20,fu20phi
c
c#####################################################################
c
      init = mod(icalc,1000)/100 .ne. 0
c
c     ----------------------------------------------
c     Load values from parameter array into /eprprm/
c     ----------------------------------------------
c
      do i=1,NFPRM
         fepr(i)=fparm(i)
      end do
c
      do i=1,NIPRM
         iepr(i)=iparm(i)
      end do
c
      if (init) then
c                                     *** Fatal error in parameters
         call lcheck(ierr)
         if (ierr.ge.FATAL) return
c
c     -------------------------------------------------------
c       If lcheck was not called, we still need to:
c         (1) convert W tensor to Cartesian form and find w0
c         (2) count potential coefficients
c         (3) set director tilt flag
c     -------------------------------------------------------
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
c
      do j=1,nfld
         spectr(j)=ZERO
      end do
c
      if (nort.le.1 .or. ipt.le.0) then
c
c       --------------------------
c       --- Single orientation ---
c       --------------------------
c
         call eprls( icalc,al,be,bss,iprune,spectr,nft,ntotal,ierr)
c
c      --------------------------
c      ---- MOMD calculation ----
c      --------------------------
c
      else
         ixt=1
         nptr=0
         onorm=ONE/dfloat(nort-1)
c
c        ------------------------------------------------------
c        Partial MOMD
c        If there is ordering of the director, calculate 
c        normalization factor for the pseudopotential function
c        In the case of partial director ordering, the parameter
c        psi actually represents the angle between the director 
c        ordering axis and the field, whereas the calculation
c        subroutine always interprets psi as the director tilt.
c        Save the given psi value as the angle xi. 
c        ------------------------------------------------------
c
         if (abs(dc20).gt.RNDOFF) then
            xi=abs(mod(psi,D180))
            cosxi=cos(RADIAN*xi)
            sinxi=sin(RADIAN*xi)
c
            acc=1.0d-4
            dlam=dc20
            call ccrint(ZERO,ONE,acc,SMALL,dnorm,nevl,fu20,id)
         else
            dlam=ZERO
            dnorm=ONE
         end if
c
c        -----------------------
c        Loop over orientations
c        -----------------------
c
         cospsi=ZERO
         dcos=ONE/(nort-1)
         
         do i=1,nort
            cospsi=min(cospsi,ONE)
            fepr(IPSI)=acos(cospsi)/RADIAN
c
            if (icalc.ne.0) then
               nstp=nstep
               if (nstep.le.0) nstp=min0(ndim,ntotal/nort)
c
c           --------------------------------------------------
c           Repeat continued-fractions for this orientation
c           with the tridiagonal matrix: find number of steps
c           by locating next zero in beta array
c           --------------------------------------------------
c
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
c
            call eprls( icalc,al(ixt),be(ixt),bss,iprune,work,
     #                  nft,nstp,ierr )
            if (ierr.ge.FATAL) return
            ixt=ixt+nstp
c
c           --------------------------------------------------
c           If there is ordering of the director, calculate the
c           weighting coefficient for this tilt angle and
c           weight the spectrum
c           --------------------------------------------------
c
            if (abs(dc20).gt.RNDOFF) then
               if (xi.gt.RNDOFF) then
                  ss=sinxi*sqrt(ONE-cospsi*cospsi)
                  cc=cosxi*cospsi
                  call ccrint(ZERO,PI,acc,SMALL,dwt,nevl,fu20phi,id)
                  dwt=dwt/(dnorm*PI)
               else
                  dwt=fu20(cospsi)/dnorm
               end if
c                  
               do j=1,nfld
                  work(j)=work(j)*dwt
               end do
            end if
c
c          --------------------------------------------------------
c          Trapezoidal rule: first and last pts. get factor of 1/2
c          --------------------------------------------------------
c
            if (i.eq.1 .or. i.eq.nort) then
               do j=1,nfld
                  work(j)=work(j)*0.5d0
               end do
            end if
c
c        ------------------------------------------------
c        Add spectrum for current orientation into total 
c        ------------------------------------------------
c
            do j=1,nfld
               spectr(j)=spectr(j)+work(j)
            end do
c
            cospsi=cospsi+dcos
         end do
c
c        -------------------
c        Normalize spectrum
c        -------------------
c
         do j=1,nfld
            spectr(j)=spectr(j)*onorm
         end do
c
c        -------------------------------------------
c        Return total number of steps actually taken
c        -------------------------------------------
c 
         if (icalc.ne.0) ntotal=ixt-1
c
      end if
c
      return
      end

