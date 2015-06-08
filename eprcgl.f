c Version 1.6  8/12/94
c*********************************************************************
c
c                          ================
c                          PROGRAM : EPRCGL
c                          ================
c                    *** Pruned basis version ***
c
c       This program does the EPR spectral calculation using the
c       complex symmetric conjugate gradients algorithm and the 
c       reconstruction of the "Lanczos" tridiagonal matrix.
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
c               dtime.c
c
c*********************************************************************
c
      program eprcgl
c
      implicit none
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'eprmat.inc'
      include 'tridag.inc'
      include 'vectrs.inc'
      include 'fnames.inc'
      include 'timer.inc'
c
      double precision ZERO
      parameter (ZERO=0.0D0)
c
      logical fexist
c
      integer i,j,nfiles,ncalc,ierr,ndone
      double precision error
c
      real cputot,cpumat,cpuvec,cpucg
c
      include 'version.inc'
c
c#####################################################################
c
      write (luttyo,1000)
      write (luttyo,1010) version,vdate
c
      call getids(nfiles,MXCALC)
      if (nfiles.eq.0) goto 9999
c
c======================================================================
c     loop over calculations
c======================================================================
c
      cputot=dtime( tarray )
      do 100 ncalc=1,nfiles
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
     #         status='unknown',form='formatted')
         rewind (unit=lulog)
c
         write (luttyo,1000)
	 write (lulog,1000)
         write (luttyo,5000) ncalc,prname(:namlth)
         write (lulog,5000) ncalc,prname(:namlth)
c
c----------------------------------------------------------------------
c     read in parameters from matrix
c----------------------------------------------------------------------
c
        inquire( file=prname(:namlth),exist=fexist )
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
c     check to make sure that CG calculation desired
c----------------------------------------------------------------------
c
        if (itype.ne.1) then
           write (luttyo,5100) prname(:namlth)
           write (lulog,5100)  prname(:namlth)
           go to 110
        end if
c
c---------------------------------------------------------------------
c    print out input parameter set
c---------------------------------------------------------------------     
c		
	call wrparm( prname,namlth,lulog )
c
c----------------------------------------------------------------------
c      Determine basis set indices (0: read .ind file if present)
c----------------------------------------------------------------------
c
        call lbasix( ixname,namlth,0 )
c
c---------------------------------------------------------------------
c     create matrix
c---------------------------------------------------------------------
c
        write(luttyo,1040)
     	write(lulog,1040) 
c
        call matrll(ierr)
        cpumat=dtime( tarray )
c
        if (ierr.eq.0) then
           write (luttyo,1050) ndim,nelre,nelim,neltot,cpumat
	   write (lulog,1050) ndim,nelre,nelim,neltot,cpumat
        else 
           write (luttyo,1070) cpumat
           write (lulog,1070) cpumat
           if (ierr.eq.1) then
              write (luttyo,1080) MXEL
              write (lulog,1080) MXEL
           else if (ierr.eq.2) then
              write (luttyo,1090) MXDIM
	      write (lulog,1090) MXDIM
           end if
           go to 110
        end if
c
c----------------------------------------------------------------------
c     Output matrix if desired
c----------------------------------------------------------------------
 
        if (stmflg) then
           open (unit=ludisk,file=mtname(:namlth),status='unknown',
     #          access='sequential',form='unformatted')
           rewind (unit=ludisk)
           write (luttyo,4000) mtname(:namlth)
           write (ludisk) (zdiag(1,i),zdiag(2,i),i=1,ndim)
           write (ludisk) (jzmat(i),i=1,ndim+1)
           write (ludisk) (izmat(i),zmat(i),i=1,nelim)
           write (ludisk) (kzmat(i),i=1,ndim+1)
           write (ludisk) (izmat(mxel-i+1),zmat(mxel-i+1),i=1,nelre)
           close (unit=ludisk)
        end if
c
c---------------------------------------------------------------------
c     starting vector calculation
c---------------------------------------------------------------------
c
       write (luttyo,2000) 
       write (lulog,2000)
c
       do i=1,mxdim
          y(1,i)=ZERO
          y(2,i)=ZERO
       end do
c
       call stvect( x )
       cpuvec=dtime( tarray )
c
       write (luttyo,2010) nelv,cpuvec
       write (lulog,2010) nelv,cpuvec
c
c----------------------------------------------------------------------
c     remove comments to store starting vector
c----------------------------------------------------------------------
c
       if (stvflg) then
          open (unit=ludisk,file=vtname,status='unknown',
     #         access='sequential',form='unformatted')
          rewind (unit=ludisk)
          write (ludisk) (x(1,i),x(2,i),i=1,ndim)
          close (unit=ludisk)
       end if
c
c---------------------------------------------------------------------
c     call conjugate gradients routine
c---------------------------------------------------------------------
c
       write (luttyo,3000) 
       write (lulog,3000) 
c
       call cscg(x,ndim,nstep,cgtol,dcmplx(shiftr,shifti),
     #           y,alpha,beta,ndone,error)
       cpucg = dtime( tarray )
c
       if (ndone.gt.0) then
          write(luttyo,3010) ndone,error,cpucg
          write(lulog,3010) ndone,error,cpucg
       else
          write(luttyo,3020) nstep,-ndone,cgtol,error,cpucg
          write(lulog,3020) nstep,-ndone,cgtol,error,cpucg
       end if
c
c----------------------------------------------------------------------
c     set nstep equal to the number of CG steps done
c----------------------------------------------------------------------
c
       nstep=abs(ndone)
c
c---------------------------------------------------------------------
c     update parameter file
c---------------------------------------------------------------------
c
       write (luttyo,4000) prname(:namlth)
       write (lulog,4000) prname(:namlth)
c
       call wrdat(prname,namlth)
c
c---------------------------------------------------------------------
c     store tridiagonal matrix
c---------------------------------------------------------------------
c
       write (luttyo,4000) tdname(:namlth)
       write (lulog,4000) tdname(:namlth)
c
       open (unit=ludisk,file=tdname,status='unknown',
     #       access='sequential',form='unformatted')
       rewind (unit=ludisk)
       write(ludisk) (alpha(1,j),alpha(2,j),beta(1,j),beta(2,j),
     #                j=1,nstep)
       close (unit=ludisk)
c
c----------------------------------------------------------------------
c       End of loop over input files
c----------------------------------------------------------------------
 110   continue
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
c======================================================================
c     end of loop over calculations
c======================================================================
c
      cputot=cputot+cpumat+cpuvec+cpucg
 100  continue
c
c---------------------------------------------------------------------
c     exit program
c---------------------------------------------------------------------
c
 9999 write (luttyo,6000) cputot
      stop
c
c=====================================================================
c     format statements
c=====================================================================
c
 1000 format(2x,70('#'))
 1010 format(25x,'program EPRCGL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
c
 1040 format(/15x,'*** Calculating matrix elements ***')
 1050 format(2x,'Dimension of the matrix : ',i6,
     #     /,2x,'Number of off-diagonal real elements: ',i8,
     #     /,2x,'Number of off-diagonal imaginary elements: ',i8,
     #     /,2x,'Total number of off-diagonal elements: ',i8,
     #     /,2x,'CPU time (seconds) : ',f10.2)
 1070 format(10x,
     #     '*** Serious error encountered in matrix generation ***'
     #     /20x,'CPU time (seconds) : ',f10.2)
 1080 format(10x,
     #  '*** Number of off-diagonal matrix elements exceeds ',i8,' ***')
 1090 format(20x,'*** Matrix dimension exceeds ',i6,' ***')
 1100 format(20x,'*** File ''',a,''' not found ***')
 1200 format(10x,'*** Error reading file ''',a,''' ***')
 1201 format(10x,'*** ''',a,''' is not an EPRLL Version ',a,' file',
     #' ***')
 2000 format(/15x,'*** Calculating starting vector ***')
 2010 format(2x,'Number of non-zero elements in ',
     #            'starting vector : ',i5,
     #       /,2x,'CPU time (seconds) : ',f10.2)
c
 3000 format(/15x,'*** Calculating tridiagonal matrix ***')
 3010 format(2x,'Number of CG steps executed : ',i5,
     #       2x,'Error in solution vector : ',g14.7,
     #       2x,'CPU time (seconds) : ',f10.2)
 3020 format(10x,'*** CG calculation did not converge ***'
     #     /,10x,'*** Maximum number of CG steps : ',i3,' ***',
     #     /,10x,'*** Number of CG steps executed : ',i3,' ***',
     #     /,10x,'*** Maximum error allowed : ',g14.7,' ***',
     #     /,10x,'*** Final error  : ',g14.7,' ***', 
     #     /,10x,'*** CPU time (seconds) : ',f10.2,' ***')
 4000 format(/20x,'*** Writing ',a,' ***')
 5000 format(15x,38('='),/,
     #     15x,'Calculation no. ',I2,' using file ',a,/,
     #     15x,38('='))
 5100 format(10x,
     #     '*** ''',a,''' is not set to run EPRCGL (itype=1) ***')
 6000 format(2x,'Total CPU time (seconds) : ',f10.2)
c
c=====================================================================
c
      end
