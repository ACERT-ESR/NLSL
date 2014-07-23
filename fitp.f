c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                         SUBROUTINE FITP
c                    =========================
c     This routine carries out nonlinear least-squares analysis of
c  slow-motional 2D spectra with EPRESA calculation.   It uses a 
c  modification of the Levenberg-Marquardt trust-region algorithm
c  available in MINPACK. 
c
c  Includes:
c     nlsdim.inc
c     stdio.inc
c     nlsnam.inc
c     eprprm.inc 
c     expdat.inc 
c     parcom.inc
c     tdspec.inc
c     lmcom.inc
c     lmtxt.inc
c     iterat.inc
c
c  Uses
c     setdat
c----------------------------------------------------------------------
c
      subroutine fitp
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'names.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'tdspec.inc'
      include 'lmcomm.inc'
      include 'miscel.inc'
      include 'lmtxt.inc'
c      include 'iterat.inc'
c
      logical fexist
      integer i,ifile,iflag,iret,ix,j,k,l,mrgn,nfev,njev,
     #        npt1,npt2,nstot,iparmptr,isite,ispec
      double precision field,rmsdv,sigma,denom,inif1,res1,
     #                 inif2,res2,amax,amin,scfac
      character line*132,chr*2
c
      double precision enorm,ten
      external pfun,enorm
      parameter(ten=10.0D0)
c
c######################################################################
c
c set setbas from basinfo
c
      setbas=.true.
      do 1 isite=1,ncomps
c        seteval(isite)=.false.	! make no assumption about matrix
        do 1 ispec=1,nspectra
	  setbas=setbas .and. (basinfo(1,ispec,isite) .ne. 0)
 1    continue
      if (.not.setbas) then
         write (luout,1001)
         if (luout.ne.luttyo) write (luttyo,1001)
 1001    format(/,'** Basis Index set is not properly set **')
         return
      end if
c
c----------------------------------------------------------------------
c     Open trace file (if desired)
c----------------------------------------------------------------------
      if (itrace.ne.0) then
         inquire(file=trname(:lthfnm),exist=fexist)
         open(lutrc,file=trname(:lthfnm),status='unknown',
     #        access='sequential',form='formatted')
         if (fexist) then
c                              * append the trace file
 5          read (lutrc,'(a2)',end=7) chr
            go to 5
         end if
 7       itrace=lutrc
      end if
c
      if (nspc.ne.nser) write(luttyo,1040)
c
c  Put correct values for the unknowns into the xxx vector.
c
      do 10 j=1,nprm	! for each parameter to be varied
        iparmptr=ixpr(j)	! identify variable
c for all sites, spectra, find one which is being varied:
        do 20 isite=1,ncomps
          do 20 ispec=1,nspectra
c check matching sites with requested site
            if (ixx(iparmptr,ispec,isite).eq.j) then
              xxx(j)=fparm(iparmptr,ispec,isite)
              go to 21
            end if
 20     continue
        write(*,*)'error in fitp, never should get here',j
	stop
 21     continue
 10   continue
c
c----------------------------------------------------------------------
c     Call pfun once if only a single calculation is specified
c     or if there is no data file, implying a single calculation
c----------------------------------------------------------------------
      if (maxitr.le.1 .or. maxev.le.1 .or. nprm.le.0 .or.
     #       .not.dataOK) then
c         seteval=.false.	! make no assumption about matrix
         if (.not.dataOK) then
            nstot=nptot
            nptot=0
            do 11 j=1,nser
 11            nptot=nptot+ndata(j)
         end if
         iflag=1
         call pfun(nptot,nprm,xxx,fvec,fjac,mxpt,iflag)
         if (iflag.lt.0 .or. ihltcmd.ne.0) return
c
         if (dataOK) then
            rmsdv=enorm(nptot,fvec)/dsqrt(dfloat(nptot))
            write(luout,1046) rmsdv
            if (luout.ne.luttyo) write(luttyo,1046) rmsdv
         else
            nptot=nstot
         end if
c
      else
c
c======================================================================
c     Call Marquardt-Levenberg nonlinear least squares routine
c======================================================================
c----------------------------------------------------------------------
c     Search parameter list and set flags for different types of parameters
c     : g-tensor, a-tensor, motional parameter, orienting potential
c----------------------------------------------------------------------
         gflag=.false.
         aflag=.false.
         rotflg=.false.
         potflg=.false.
c         seteval=.false.	! make no assumption about matrix
         do 12 i=1,nprm
            ix=ixpr(i)
            gflag=gflag .or. (ix.ge.IGXX .and. ix.le.IGZZ)
            aflag=aflag .or. (ix.ge.IAXX .and. ix.le.IAZZ)
            rotflg=rotflg .or. (ix.ge.IDX .and. ix.le.IPMZZ)
     #           .or. (ix.eq.IDJF .or. ix.eq.IDJFPRP)
     #           .or. (ix.eq.IOSS)
            potflg=potflg .or. (ix.ge.IC20 .and. ix.le.IC44)
 12      continue
c
         nprint=1
c     
         call lmnls( pfun,nptot,nprm,xxx,fvec,fjac,mxpt,ftol,
     1          xtol,gtol,maxev,maxitr,diag,prscl,factor,nprint,
     2          itrace,jacobi,info,nfev,njev,ipvt,qtf,gnvec,gradf,
     3          work1,work2,work3,work4 )
c
c----------------------------------------------------------------------
c     calculate covariance matrix and estimate errors assuming
c     normal distribution
c----------------------------------------------------------------------
c
	 if (ihltcmd.ne.0) return
c                                         *** MINPACK error return
         if (info.eq.0 .or. info.gt.3) then
            write (luout,2034) minerr(info)
            if (luout.ne.luttyo) write (luttyo,2034) minerr(info)
            if (itrace.ne.0)  write (lutrc,2034) minerr(info)
            return
c                                         *** Normal return
	 else
            fnorm=enorm(nptot,fvec)
            rmsdv=fnorm/dsqrt( dfloat(nptot) )
            sigma=fnorm/dsqrt( dfloat(nptot-nprm) )
            write(luout,1000)
	    write(luout,2035) minerr(info),nfev,nprm,nptot
            write(luout,2036) rmsdv,fnorm*fnorm
            write(luout,1000)
c
	    if (luout.ne.luttyo) then
               write(luttyo,2035) minerr(info),nfev,nprm,nptot
               write(luttyo,2036) rmsdv,fnorm*fnorm
            end if   
            if (itrace.ne.0) then
               write(lutrc,2035) minerr(info),nfev,nprm,nptot
               write(lutrc,2036) rmsdv,fnorm*fnorm
            end if
         end if
c
c-----------------------------------------------------------------------
c     Calculate residuals and covariance matrix
c     COVAR calculates covariance matrix from the R matrix of QR 
c     factorization, which lmder stored in upper triangular form in fjac 
c     ipvt is permutation array returned by lmder
c------------------------------------------------------------------------
c
         call covar( nprm,fjac,mxpt,ipvt,xtol*xtol,work1 )
c
c----------------------------------------------------------------------
c     Calculate correlation matrix from covariance and output it
c----------------------------------------------------------------------
c
         line=' '
         mrgn=max(1,36-4*nprm)
         if (mrgn.lt.1) mrgn=1
	 write (luout,2000) 
         write(line(mrgn:),2002) (tag(i),i=1,nprm)
         write(luout,2001) line(:8*(nprm+1)+mrgn)
c
	 do 23 i=1,nprm
            line=' '
            do 15 j=i,nprm
	       if (i.eq.j) xerr(i)=sigma*dsqrt( fjac(i,i) )
               denom = dsqrt( fjac(i,i) * fjac(j,j) )
               if (denom .ne. 0.0D0) corr(i,j) = fjac(i,j)/denom
               k=8*(j-1)+mrgn
               write(line(k:),2003) corr(i,j)
 15	    continue
            write (luout,2001) line(:8*(nprm+1)+mrgn)
 23	 continue
c
c----------------------------------------------------------------------
c     output final fit of parameters
c----------------------------------------------------------------------
         if (.not.(info.eq.0.or.info.gt.3)) then
            write (luout,1000)
            write (luout,2004) 
            write (luout,2007) (tag(i),xxx(i),xerr(i),i=1,nprm)
            write (luout,*)
         end if
c
c----------------------------------------------------------------------
c     For fit of spherical tensor components, report Cartesian tensor
c     and propagate uncertainty estimates
c----------------------------------------------------------------------
c
c       g-tensor
c
c        if (gflag.and.igflg.ne.1) then
c           write(luout,1060) gxx,gyy,gzz
c        end if
c
c       A-tensor
c
c        if (aflag.and.iaflg.ne.1) then
c           write(luout,1061) axx,ayy,azz
c        end if
c
c       R-tensor
c
c       Calculate uncertainty in each rate from uncertainty in log(rate)
c
c        if (rotflg) then
c           
c           if (irflg.ne.1) write(luout,1062) dx,dy,dz
c
c      end if
c
c         if (potflg)  call ordrpr(c20,ordpar)
c
      end if
c
c----------------------------------------------------------------------
c     Output fit (or calculated spectrum) for each datafile
c----------------------------------------------------------------------
      do 40 i=1,nser
c     
c *** 2D-output format **
c
c---------------------------------------------------------------------
c     Preparation of graph (Frequencies in MHz)
c---------------------------------------------------------------------
c
         amin=0.0d0
         amax=0.0d0
c                                   * unscale the spectrum
         scfac=1.0d0
         if (dataOK) then
	   scfac=1.0d0/sratio(i)
           do 45 j=1,ndata(i)
             k=ixsp(i)+j-1
             totspec(k)=(fvec(k)/datwgt(i)+data(k))*scfac
             if ( totspec(k).gt.amax ) amax=totspec(k)
 45        continue
	 else
c totspec already calculated, use it
	   do 46 j=1,ndata(i)
             k=ixsp(i)+j-1
             if ( totspec(k).gt.amax ) amax=totspec(k)
 46        continue
         end if
c
c         npt1=snpt1(i)
c         npt2=snpt2(i)
c         inif2=-5.0d2*dble(snpt2(i)-1)/sstept2(i)/dble(snpt2(i))
c         res2=-2.0d0*inif2/(snpt2(i)-1)
c	 inif2=-5.0d2/sstept2(i)
c         inif1=-5.0d2*dble(snpt1(i)-1)/sstept1(i)/dble(snpt1(i))
c         res1=-2.0d0*inif1/(snpt1(i)-1)
c	 inif1=-5.0d2/sstept1(i)
c                              * obtain output file name
         call setdat (dataid(i))
c
c         call wrfit( totspec(ixsp(i)),npt1,npt2,inif2,res2,
c     #               inif1,res1,amin,amax,ftname,lthdnm )
         call wrfit( totspec(ixsp(i)),i,amin,amax,ftname )

c
 40   continue
c
c----------------------------------------------------------------------
c     Close files and return to calling program
c----------------------------------------------------------------------
      if (itrace.ne.0) close (lutrc)
      return
c   
c# format statements ##################################################
c
 1000 format(/,2x,70('-'),/)
 1040 format(/10x,'*** Some datafiles are not read yet. ***',
     #       /10x,'*** Fit is independent of the data. ***')
 1046 format(10x,'Single calculation : RMS deviation =',g14.7)
c 1060 format('Cartesian g-tensor = (',2(f10.7,','),f10.7,')')
c 1061 format('Cartesian A-tensor = (',2(f9.3,','),f9.3,')')
c 1062 format(10x,'Rx =',e13.6,'   Ry =',e13.6,'   Rz =',e13.6)
 2000 format(/24x,'*** Correlation Matrix ***'/)
 2001 format(a)
 2002 format(8(1x,a6,1x))
 2003 format(f8.4)
 2004 format(/9x,'*** Final Parameters ***',/,
     1       7x,' Parameter      Value'/7x,26('-')/)
 2007 format(10x,a6,' = ',g13.7,' +/- ',g13.7)
 2034 format(/2x,'MINPACK could not find a solution: ', a )
 2035 format(/9x,'MINPACK completed: ', a/
     #       12x,'Function evaluations: ',i4/
     #       12x,'There were ',i2,' parameters and ',i6,' data points')
 2036 format(12x,'Rms Deviation = ',g13.6/12x,'chi-squared = ',g13.6/)
c
      end
