c Version 1.6 8/12/94
c*********************************************************************
c
c                       ================
c                       PROGRAM : MOMDL
c                       ================
c
c       This program calculates the MOMD (Microscopic Order Macroscopic
c       Disorder) spectrum for a randomly oriented dispersion of domains
c       of a locally ordered fluid. The calculation utilizes the EPRCGL
c       program to calculate a spectrum for each effective director tilt
c       angle and averages the spectra according to an isotropic 
c       distribution of tilt angles. 
c
c       written by DJS 24-OCT-87
c       modified by DEB 13-MAR-92 to allow longer filenames
c                       22-OCT-92 Handles case of version incompatibility
c
c       Uses :
c               getids.f
c               setnam.f
c               rddat.f
c               wrdat.f
c               matrll.f
c               stvect.f
c               cscg.f
c
c*********************************************************************
c
      program momdl
c
      implicit none
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'eprmat.inc'
      include 'tridag.inc'
      include 'vectrs.inc'
      include 'spectr.inc'
      include 'fnames.inc'
      include 'pidef.inc'
      include 'timer.inc'
c
      double precision ZERO
      parameter (ZERO=0.0D0)
c
      logical fexist
c
      double precision bsweep,db,gib
      double precision temp(2,MXPT),ascale,cospsi
c
      integer i,iang,j,nfiles,ncalc,ierr,ndone
      double precision error,cputot,cpumat,cpuvec,cputri
      character*30 sname
c
      include 'version.inc'
c
c#####################################################################
c
      write (luttyo,1000)
      write (luttyo,1010) version,vdate
c
      call getids(nfiles,mxcalc)
      if (nfiles.eq.0) goto 9999
c
c======================================================================
c     loop over calculations
c======================================================================
c
      cputot=dtime( tarray )
      do 110 ncalc=1,nfiles
c
        fileid = files(ncalc)
c
c---------------------------------------------------------------------
c     construct file names
c---------------------------------------------------------------------
c
        call setnam
        call setflg( flags(ncalc) )
c     
c----------------------------------------------------------------------
c     open log file and write out header  
c----------------------------------------------------------------------
c
	open (unit=lulog,file=rlname(:namlth),access='sequential',
     #        status='unknown',form='formatted')
	rewind (unit=lulog)
c
        write (luttyo,1000)
	write (lulog,1000)
        write (luttyo,5000) ncalc,prname(:namlth)
    	write (lulog,5000) ncalc,prname(:namlth)
c
c----------------------------------------------------------------------
c     read in parameters from input parameter file
c----------------------------------------------------------------------
c
	inquire(file=prname(:namlth),exist=fexist)
c
	if (fexist) then 
           call rddat(prname,namlth,ierr)
           if (ierr.eq.-100) then
              write (luttyo,1201) prname(:namlth),version
              write (lulog,1201) prname(:namlth),version
              go to 110
           else if (ierr.ne.0) then
              write (luttyo,1200) prname(:namlth)
              write (lulog,1200) prname(:namlth)
              go to 110
           end if
        else
	   write (luttyo,1100) prname(:namlth)
           write (lulog,1100) prname(:namlth)
           go to 110
        endif
c
c----------------------------------------------------------------------
c     check to make sure that MOMD calculation desired
c----------------------------------------------------------------------
c
        if (itype.ne.1) then
           write (luttyo,5100) prname(:namlth)
           write (lulog,5100)  prname(:namlth)
           go to 110
        end if
c
c----------------------------------------------------------------------
c Initialize field array and spectrum
c     Note that we are calculating an energy spectrum but
c     want to represent a field-swept spectrum. To calculate the
c     energy spectrum, the energy sweep variable must be swept in the 
c     opposite sense to that specified for the field sweep (i.e.
c     higher transition energy=lower resonance field for a given
c     spectrometer frequency). Thus, instead of sweeping up from 
c     fieldi to fieldf, we sweep down from -fieldi to -fieldf. 
c     To save space, the energy sweep values are kept in the "field" 
c     array 
c----------------------------------------------------------------------
c
        db=(fieldf-fieldi)/dble(nfield-1)
        do i=1,nfield
           spect(1,i)=ZERO
           spect(2,i)=ZERO
        end do
c
      write (luttyo,2000)
      write (lulog,2000)
c
c----------------------------------------------------------------------
c    Loop over angles
c----------------------------------------------------------------------
c
      call lbasix( ixname,namlth,0 )
c
      do 100 iang=1,nort
         if (nort.gt.1) then
            if (angflg) then
               psi=(iang-1)*90.0d0/(nort-1)
               cospsi=cos( PI*psi/180.0d0 )
               ascale=dsqrt( 1.0d0 - cospsi*cospsi )
            else
               cospsi=dble(iang-1)/dble(nort-1)
               psi=180.0d0*acos( cospsi )/pi
               ascale=1.0d0
               if (iang.eq.1 .or. iang.eq.nort) ascale=0.5d0
            end if
            ascale=ascale/(PI*(nort-1))
         end if
c
c---------------------------------------------------------------------
c     create matrix
c---------------------------------------------------------------------
c
           call matrll(ierr)
           cpumat=dtime( tarray )
c
        if (ierr.ne.0) then
           write (luttyo,1070) iang,cpumat
           write (lulog,1070) iang,cpumat
           if (ierr.eq.1) then
              write (luttyo,1080) mxel
              write (lulog,1080) mxel
           else if (ierr.eq.2) then
              write (luttyo,1090) mxdim
              write (lulog,1090) mxdim
           end if
           go to 110
        end if
c
c---------------------------------------------------------------------
c     starting vector calculation
c---------------------------------------------------------------------
c
        do i=1,MXDIM
           y(1,i)=ZERO
           y(2,i)=ZERO
        end do
c
        call stvect( x )
        cpuvec=dtime( tarray )
c
c---------------------------------------------------------------------
c    call conjugate gradients routine
c---------------------------------------------------------------------
c
        call cscg(x,ndim,nstep,cgtol,dcmplx(shiftr,shifti),
     #       y,alpha,beta,ndone,error)
        cputri=dtime( tarray )
c
        if (ndone.lt.0) then
           write(luttyo,3020) nstep,-ndone,cgtol,error,cputri
           write(lulog,3020) nstep,-ndone,cgtol,error,cputri
           ndone=-ndone
        end if
c
        do j=1,nfield
           temp(1,j)=ZERO
           temp(2,j)=ZERO
        end do
c
        gib=gib0
        if (ipt.ne.0) gib=gib0+gib2*(1.0D0-cospsi*cospsi)
c
c----------------------------------------------------------------------
c       Call continued-fractions routine to calculate spectrum
c----------------------------------------------------------------------
c
        call cfvd(alpha,beta,ndone,-fieldi,-db,nfield,w0,
     #            temp,ideriv)
c
c----------------------------------------------------------------------
c      Convolute spectrum with Gaussian lineshape if specified
c----------------------------------------------------------------------
c
        call gconvl(temp,gib,db,nfield)
c
c----------------------------------------------------------------------
c      Write individual spectrum to disk if needed
c      On output, add B0 offset to field and transpose spectum around B0 
c      (by transposing calculated field value relative to B0) to 
c      represent a field-swept spectrum
c----------------------------------------------------------------------
c 
        if (savflg.and.nort.gt.1) then
           sname=spname
           write(sname(namlth-2:),'(i3.3)') iang
           open ( unit=ludisk, file=sname(:namlth), status='unknown',
     #          form='formatted' )
           write (ludisk,1120)(b0+fieldi+db*(i-1),temp(1,i),temp(2,i),
     #          i=1,nfield)
           close (unit=ludisk)
        end if
c 
c----------------------------------------------------------------------
c  Integrate over orientation using appropriate weighting
c  (Or normalize a single-orientation spectrum)
c----------------------------------------------------------------------
c
        if (nort.gt.1) then
           do j=1,nfield
              spect(1,j)=spect(1,j)+ascale*temp(1,j)
              spect(2,j)=spect(2,j)+ascale*temp(2,j)
           end do
        else
           do j=1,nfield
              spect(1,j)=spect(1,j)*ascale
              spect(2,j)=spect(2,j)*ascale
           end do
        end if
c
        write (lulog,2010) iang,psi,neltot,nelv,ndim,ndone,
     #       cpumat+cpuvec+cputri
        write (luttyo,2010) iang,psi,neltot,nelv,ndim,ndone,
     #       cpumat+cpuvec+cputri
c
        cputot=cputot+cpumat+cpuvec+cputri
c
 100  continue
c
c
c----------------------------------------------------------------------
c     write accumulated MOMD spectrum to disk (calculate field as
c     noted above for field-swept spectrum)
c----------------------------------------------------------------------
c
          open (unit=ludisk,file=spname(:namlth),status='unknown',
     #         form='formatted')
          write (ludisk,1120)(b0+fieldi+db*(i-1),spect(1,i),spect(2,i),
     #         i=1,nfield)
          close (unit=ludisk)
c
c
c----------------------------------------------------------------------
c       End of loop over input files
c----------------------------------------------------------------------
          write (luttyo,6000) cputot
          cputot=ZERO
 110	continue
c
c---------------------------------------------------------------------
c     print out trailer and close log file
c---------------------------------------------------------------------
c
	write (luttyo,1000)
	write (lulog,1000)
c
	close (unit=lulog)
c
 9999   stop 
c
c=====================================================================
c     format statements
c=====================================================================
c
 1000 format(2x,70('#'))
 1010 format(25x,'program MOMDL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
c
 1070 format(10x,
     #     '*** Serious error encountered in matrix generation ***'
     #     /20x,'Step ',i3,'; CPU time (seconds) : ',f10.2)
 1080 format(10x,
     #  '*** Number of off-diagonal matrix elements exceeds ',i8,' ***')
 1090 format(20x,'*** Matrix dimension exceeds ',i6,' ***')
 1100 format(20x,'*** File ''',a,''' not found ***')
 1120 format(3(g14.7,','))
 1200 format(10x,'*** Error reading file ''',a,''' ***')
 1201 format(10x,'*** ''',a,''' is not an EPRLL Version ',a,' file',
     #' ***')
 2000 format(2x,70('-'),/,4x,'i',3x,'angle',4x,'neltot',3x,
     #      'nelv',3x,'ndim',3x,'nstep',2x,'cpu',/,2x,70('-'))  
 2010 format(2x,i3,3x,f5.2,3x,i7,3x,i3,4x,i4,3x,i4,3x,f6.2)
 3020 format(10x,'*** CG calculation ',i3,' did not converge ***',
     #     /,10x,'*** Maximum number of CG steps  : ',i3,' ***',
     #     /,10x,'*** Number of CG steps executed : ',i3,' ***',
     #     /,10x,'*** Maximum CG error allowed    : ',g14.7,' ***',
     #     /,10x,'*** Final CG error              : ',g14.7,' ***',
     #     /,10x,'*** True error                  : ',g14.7,' ***',
     #     /,10x,'*** CPU time (seconds)          : ',g14.2,' ***')
c
 4000 format(20x,'*** Writing ',a,' ***')
 5000 format(15x,38('='),/,
     #     15x,'Calculation no. ',I2,' using file ',a,/,
     #     15x,38('='))
 5100 format(10x,
     #     '*** ''',a,''' is not set to run MOMDL (itype=1) ***')
 6000 format(2x,'Total CPU time (seconds) : ',f10.2)
c
c=====================================================================
c
      end
