c Version 1.6  5/2/94
c**********************************************************************
c
c                       ==============
c                       PROGRAM : TDLL
c                       ==============
c
c       This program processes the tridiagonal matrix generated 
c       by the Lanczos or conjugate gradients algorithms.
c
c       written by DJS 29-AUG-87
c       Updates 1.1 through 1.5 by David E. Budil
c
c
c       Includes:
c               stddim.inc
c               stdio.inc
c               eprdat.inc
c               fnames.inc
c               pidef.inc
c               rndoff.inc
c
c       Uses:
c               rddat.f
c               getids.f
c               setnam.f
c               prmpt.f
c               cmtqli.f
c               cfvd.f
c
c**********************************************************************
c
      program tdll
c
c      implicit none
c
      include 'stddim.inc'
      include 'stdio.inc'
      include 'eprdat.inc'
      include 'fnames.inc'
      include 'tridag.inc'
      include 'spectr.inc'
      include 'pidef.inc'
      include 'rndoff.inc'
c
      double precision ZERO,UNITY,TWO
      parameter (ZERO=0.0D0,UNITY=1.0D0,TWO=2.0D0)
c
      logical fprmpt
c
      integer i,ifile,nstepl,ierr,isort,nfiles,ndiag,npts
c
      double precision bi,bf,db,tr,ti,gib
      double precision temp1(2,MXSTEP),temp2(2,MXSTEP),w(2,MXSTEP)
c
      logical fexist
c
      character option*1,line*80
c
      integer itoken
      logical prmpt
      double precision ftoken
      external prmpt,itoken,ftoken
c
      include 'version.inc'
c
c######################################################################
c
      write (luttyo,1000)
      write (luttyo,1010) version,vdate
c
      ifile=0
      call getids(nfiles,mxcalc)
      if (nfiles.eq.0) goto 9999
c
 10   ifile=ifile+1
      if (ifile.gt.nfiles) go to 9999
      fileid=files(ifile)
      call setnam
c
c----------------------------------------------------------------------
c     read in data from parameter file
c----------------------------------------------------------------------
c
      inquire(file=prname(:namlth),exist=fexist)
c
      if (fexist) then
         call rddat(prname,namlth,ierr)
         if (ierr.ne.0) then
            write (luttyo,5011) prname(:namlth)
            go to 10
         end if
      else
         write (luttyo,5010) prname(:namlth)
         go to 10
      endif
c
c----------------------------------------------------------------------
c     read in tridiagonal matrix
c----------------------------------------------------------------------
c
      open (unit=ludisk,file=tdname(:namlth),status='old',
     #      access='sequential',form='unformatted',iostat=ierr)
      if (ierr.eq.-100) then
         write (luttyo,5012) prname(:namlth)
         go to 10
      else if (ierr.ne.0) then
         write(luttyo,5010) tdname(:namlth)
         go to 10
      end if
      read (ludisk,iostat=ierr) (alpha(1,i),alpha(2,i),
     #                          beta(1,i),beta(2,i),i=1,nstep)
      if (ierr.ne.0) then
         write(luttyo,5011) tdname(:namlth)
         goto 10
      end if
c
      close (unit=ludisk)
c
c======================================================================
c======================================================================
c     prompt user for processing option
c======================================================================
c======================================================================
c
 25   write (luttyo,1040) tdname(:namlth)
      read (luttyi,1050,err=25) option
      if (option.eq.' ') goto 10
      if (option.eq.'5') goto 9999
c
c----------------------------------------------------------------------
c     list matrix elements on screen
c----------------------------------------------------------------------
c
      if (option.eq.'1') then
         write (luttyo,1060) tdname(:namlth)
         do i=1,nstep
            write (luttyo,1070) i,alpha(1,i),alpha(2,i),
     #                          beta(1,i),beta(2,i)
         end do
c
c----------------------------------------------------------------------
c     write formatted tridiagonal matrix file
c----------------------------------------------------------------------
c
      else if (option.eq.'2') then
         open (unit=ludisk,file=tfname(:namlth),status='unknown',
     #         access='sequential',form='formatted')
         write (ludisk,1080) (dble(i),alpha(1,i),alpha(2,i), 
     #                        beta(1,i),beta(2,i),i=1,nstep)
         close (unit=ludisk)
c
c----------------------------------------------------------------------
c     diagonalize tridiagonal matrix 
c----------------------------------------------------------------------
c
      else if (option.eq.'3') then
         ndiag=nstep
         write(luttyo,1150) ndiag
         fprmpt=prmpt(line)
         if (fprmpt) ndiag=itoken(line)
c
         write(luttyo,1151)
         read(luttyi,1152) line
         isort=1
         if(line(1:1).eq.'w' .or. line(1:1).eq.'W') isort=2
c
         if (ndiag.gt.nstep) ndiag=nstep
         do i=1,ndiag
            temp1(1,i)=alpha(1,i)
            temp1(2,i)=alpha(2,i)
         end do
c       
         do i=1,ndiag
            temp2(1,i)=beta(1,i)
            temp2(2,i)=beta(2,i)
         end do
c       
         temp2(1,ndiag)=ZERO
         temp2(2,ndiag)=ZERO
c
         call cmtqli(ndiag,temp1,temp2,w,ierr,isort)
c
         if (ierr.ne.0) write(luttyo,1140) ierr
c
         do i=1,ndiag
            tr=w(1,i)
            ti=w(2,i)
            w(1,i)=tr*tr-ti*ti
            if (abs(w(1,i)).lt.RNDOFF) w(1,i)=ZERO
            w(2,i)=two*tr*ti
            if (abs(w(2,i)).lt.RNDOFF) w(2,i)=ZERO
         end do
c       
c     *** Output the eigenvalues
c 
         open (unit=ludisk,file=egname(:namlth),status='unknown',
     #        access='sequential',form='formatted')
         write (ludisk,1080) (dble(i),temp1(1,i),temp1(2,i),
     #                       w(1,i),w(2,i),i=1,ndiag)
         close (unit=ludisk)
c
c----------------------------------------------------------------------
c     calculate spectrum from continued fraction
c----------------------------------------------------------------------
c
      else if (option.eq.'4') then
         bi=-50.0D0
         bf=50.0D0
         write(luttyo,1090) b0,bi,bf
         fprmpt=prmpt(line)
         if (fprmpt) then
            bi=ftoken(line)
            bf=ftoken(line)
         end if
c
         gib0=ZERO
         gib2=ZERO
         if (ipsi0.ne.0) then
            write(luttyo,1100) gib0,gib2
            fprmpt=prmpt(line)
            if (fprmpt) then
               gib0=ftoken(line)
               gib2=ftoken(line)
            end if
         else
            write (luttyo,1101) gib0
            fprmpt=prmpt(line)
            if (fprmpt) gib0=ftoken(line)
         end if
c
         gib0=dabs(gib0)
         if (dabs(gib2).gt.gib0) gib2=dsign(gib0,gib2)
c
         gib=gib0
         if (ipt.ne.0) gib=gib0+gib2*(dsin(psi*PI/180.0D0))**2
c
         nstepl=nstep
         write(luttyo,1110) nstep
         fprmpt=prmpt(line)
         if (fprmpt) nstepl=itoken(line)
c
         ideriv = 1
         write(luttyo,1115) ideriv
         fprmpt=prmpt(line)
         if (fprmpt) ideriv=itoken(line)
c
         npts=512
         write(luttyo,1116) MXPT,npts
         fprmpt=prmpt(line)
         if (fprmpt) npts=itoken(line)
         if(npts.lt.1) npts=1
         if(npts.gt.mxpt) npts=MXPT
c
c----------------------------------------------------------------------
c Initialize field array and spectrum
c     Note that we are calculating an energy spectrum but
c     want to represent a field-swept spectrum. To calculate the
c     energy spectrum, the energy sweep variable must be swept in 
c     the opposite sense to that specified for the field sweep (i.e.
c     higher transition energy=lower resonance field for a given
c     spectrometer frequency). Thus, instead of sweeping up from 
c     bi to bf, we sweep down from -bi to -bf.
c----------------------------------------------------------------------
         db=(bf-bi)/dble(npts-1)
         do i=1,npts
            spect(1,i)=ZERO
            spect(2,i)=ZERO
         end do
c
         call cfvd(alpha,beta,nstep,-bi,-db,npts,w0,spect,ideriv)
c
         do i=1,npts
            spect(1,i)=spect(1,i)/PI
            spect(2,i)=spect(2,i)/PI
         end do
c----------------------------------------------------------------------
c      Convolute spectrum with Gaussian lineshape if specified
c----------------------------------------------------------------------
c
         call gconvl(spect,gib,db,npts)
c
c----------------------------------------------------------------------
c  Add B0 offset to field and transpose spectrum around B0 in order
c  to represent a field-swept spectrum
c     write spectrum to disk
c----------------------------------------------------------------------
         open (unit=ludisk,file=spname(:namlth),status='unknown',
     #         form='formatted')
         write (ludisk,1120) (b0+bi+(i-1)*db,spect(1,i),
     #                        spect(2,i),i=1,npts)
         close (unit=ludisk)
c
      end if
c
      go to 25
c
c----------------------------------------------------------------------
c     exit program
c----------------------------------------------------------------------
c
 9999 write(luttyo,1000)
      stop
c
c======================================================================
c     format statements
c======================================================================
c
 1000 format(/70('#')/)
 1010 format(25x,'program TDLL'/22x,'Version ',a,1x,a/
     #       22x,'----------------')
 1040 format(//,25x,'file : ',a,/,
     #       /,5x,'1 - list elements of tridiagonal matrix on screen',
     #       /,5x,'2 - write tridiagonal matrix into formatted file',
     #       /,5x,'3 - diagonalize tridiagonal matrix',
     #       /,5x,'4 - generate cw spectrum from continued fraction',
     #       /,5x,'5 - exit program',
     #       /,5x,'<ENTER> - next file',/
     #       /,10x,'please select option : ')
 1050 format(a)
 1060 format(//,25x,'file : ',a,//,
     #       ' no.',16x,'alpha',32x,'beta',/,70('-'))
 1070 format(i3,2(4x,2(2x,g14.7)))
 1080 format(5(g14.7,1x))
 1090 format(//,2x,'Please enter initial and final fields',
     #       ' relative to B0 = ',f8.2,/,2x,
     #       '[ bi = ',f7.2,', bf = ',f7.2,' ]')
 1100 format(/,2x,'Please enter inhomogeneous linewidths',/,2x,
     #       '[ gib0 = ',g10.5,', gib2 = ',g10.5,' ]')
 1101 format(/,2x,'Please enter inhomogeneous linewidth',/,2x,
     #       '[ gib0 = ',g10.5,' ]')
 1110 format(/,2x,'Please enter number of Lanczos steps',/2x,
     #       ' for calculation',/,2x,'[ nstep =',i4,' ]')
 1115 format(/,2x,'Please specify 0th or 1st derivative spectrum',
     #       /2x,'[ ideriv =',i2, ' ]' )
 1116 format(/,2x,'Please specify number of points in spectrum ',
     #       '[max ',i5,']'/,2x,'[ npts =',i5,' ]')
 1120 format(3(g14.7,1x))
 1140 format(/,10x,'*** Error in diagonalization: ierr=',i3,' ***')
 1150 format(/,2x,'Dimension of matrix to be diagonalized ',
     #       '[max=',i3,'] : ')
 1151 format('Sort eigenvalues by (F)ield or (W)eight? ',
     # ' [F] : ')
 1152 format(a)
 5010 format(/,10x,'*** File ',a,' not found ***',/)
 5011 format(/,10x,'*** Error reading file ',a,' ***',/)
 5012 format(/,10x,'*** ''',a,''' is not an EPRLL Version 1.5 file',
     #' ***')
c
c======================================================================
c
        end
