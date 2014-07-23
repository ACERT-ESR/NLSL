c  VERSION 1.0  (NLSPMC version)   2/5/99
c*********************************************************************
c 
c                       ==================
c                       SUBROUTINE : SIM2D
c                       ==================
c
c     Subroutine version of MOMD program
c     
c     This routine is intended for use with nonlinear least-squares
c     applications.  The routine calculates 2D-spectrum for the 
c     parameters given in /eprprm/ using the eigenvalues and eigenvectors
c     already obtained by evcgf routine.
c
c     On Entry :
c
c        icalc : = 0  recalculate the spectrum using previous xspec.
c                     (valid when varying gib or lib parameters)
c                = 1  recalculate the spectrum using the previous
c                        eigenvalues and eigenvectors.
c                     (valid for series of 2D spectra with different
c                        experimental type, combination, mixing time
c                        and so on
c                      also valid when varying hwid parameter)
c                = 2  recalculate the eigenvalues and the spectrum
c        xspec : complex time-domain spectrum (used for icalc=0 case)
c
c     On Exit  :
c
c        xspec : complex time-domain spectrum (to be used with icalc=0)
c        cspec  : complex frequency spectrum
c        ierr  : error flag
c                      = 0   normal return
c                      < 0   critical error in one of
c                            pcheck ; wrong parameter specification
c                            matrx? ; wrong matrix generation
c                            cmtqli ; failure in QL-decomposition
c     Includes :
c               nlsdim.inc
c               stdio.inc
c               eprprm.inc
c               prmeqv.inc
c               parcom.inc
c               stvcom.inc
c               egvcom.inc
c
c     Uses :
c               evcgf.f
c               spcalc.f
c               convft.f
c
c*********************************************************************
c
      subroutine sim2d(icalc,xspec,cspec,npt1,npt2,nvar,
     #                 iflag,ixi,delx,ispec,isite,ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'stdio.inc'
      include 'simparm.inc'
      include 'parmequ.inc'
      include 'parms.inc'
      include 'stvcom.inc'
      include 'basis.inc'
      include 'egvcom.inc'
      include 'miscel.inc'
c
      double precision zero,scal,pi
      complex*16 czero,ci
      parameter (zero=0.0D0,czero=(0.0D0,0.0D0),ci=(0.0d0,1.0d0))
c
      integer npt1,npt2,icalc,i,j,k,l,iort,bdw,
     #        ntmp,nvar,iflag,ixi,ierr,ispec,isite
      integer tt1,tt2,tt3,nbindx
      double precision t1,t2,cspsi,amin,amax,delx
c
      complex*16 xspec(npt1,npt2),cspec(npt1,npt2)
c      real*8 spec(npt1,npt2)
c
      logical evalOK
      external evalOK
c
c#####################################################################
c

      ierr=0
      pi=4.0d0*datan(1.0d0)
c
c----------------------------------------------------------------------
c     Load values from parameter array into /eprprm/
c----------------------------------------------------------------------
c
10    continue
c
c Move the current ispec, isite parameters into the fepr and iepr 
c arrays.
c
      do 2 i=1,nfprm
         fepr(i)=fparm(i,ispec,isite)
c         write(*,*)'i,fepr ',i,fepr(i)
 2    continue
c
c only copy part of eprprm, NIPRM=23 ???  Rest must not matter?
c same above too...
c
      do 4 i=1,niprm
         iepr(i)=iparm(i,ispec,isite)
c         write(*,*)'i,iepr ',i,iepr(i)
 4    continue
c record site and specturm being calculated:
c      specnum=ispec
c      sitenum=isite
c
c Move the basis set into position for this simulation if necessary:
c
      if (basinfo(1,ispec,isite).ne.basisptr) then
        nbindx=basinfo(1,ispec,isite)
        basisptr=nbindx
        ndimo=ndimoss(nbindx)
        ndimd=ndimdss(nbindx)
        j=0
        tt1=pidptr(nbindx)
        tt2=pidptr(nbindx)+ndimoss(nbindx)-1
         do 41 i=tt1,tt2
c        do 41 i=pidptr(nbindx),pidptr(nbindx)+ndimoss(nbindx)-1
          j=j+1
          jqe1(j)=mjqe1(i)
          pi1(j)=mpi1(i)
          qi1(j)=mqi1(i)
          l1(j)=ml1(i)
          jk1(j)=mjk1(i)
          k1(j)=mk1(i)
          jm1(j)=mjm1(i)
          m1(j)=mm1(i)
 41     continue
        j=0
        tt1=dpidptr(nbindx)
        tt2=dpidptr(nbindx)+ndimdss(nbindx)-1
        do 42 i=tt1,tt2
c        do 42 i=1+dpidptr(nbindx),dpidptr(nbindx)+ndimdss(nbindx)-1
          j=j+1
          djqe1(j)=mdjqe1(i)
          dpi1(j)=mdpi1(i)
          dqi1(j)=mdqi1(i)
          dl1(j)=mdl1(i)
          djk1(j)=mdjk1(i)
          dk1(j)=mdk1(i)
          djm1(j)=mdjm1(i)
          dm1(j)=mdm1(i)
          pp(j)=mpp(i)
 42     continue
        do 420 i=1,ndimoss(nbindx)
          pid(i)=mpid(i+pidptr(nbindx)-1)
 420    continue
      end if
      if (icalc.ge.1) then
c                                    *** test validity of parameters
	if (basinfo(1,ispec,isite) .eq. 0) then
          write(luout,1063)
          if (luout.ne.luttyo) write(luttyo,1063)
 1063	  format(/5x,'** ERROR, basis not set in SIM2D **')
	  return
        else
	  setbas=.true.
	end if
        call pcheck(luttyo,ierr)
        if (ierr.lt.0) then
          write(luout,1000)
          if (luout.ne.luttyo) write(luttyo,1000)
          write(luout,1100)
          if (luout.ne.luttyo) write(luttyo,1100)
          return
        end if
c                                    *** off-diagonal starting vector
c get stvo
c
        call stveco
c
      end if
c
c---------------------------------------------------------------------
c     calculate eigenvalues & eigenvectors, egval and egvec
c     Only do this once for each site.
c---------------------------------------------------------------------
c
      if (icalc.ge.2) then
c                        ===============================
c                     === loop over orientations (MOMD) ===
c                        ===============================
         do 100 iort=1,nort
c
            if (nort.gt.1) then
               cspsi=dfloat(iort-1)/dfloat(nort-1)
               psi=dacos(cspsi)*1.8D2/pi
            end if
c
            if (idebug.ne.0) then
               write (ludeb,2000) iort,psi
               write (ludeb,2001)
            end if
c                     *** solve for eigenvalues & vectors for pS=1
            call evcgf(1,stvo,nevo,egvalx(1,iort,isite),
     #        egvecx(1,1,iort,isite),ierr)	! nevo returned is # of imp't evals.
c
            if (ierr.lt.0) then
               write(luout,1100)
               if (luout.ne.luttyo) write(luttyo,1100)
               return
            end if
c 
	    nev1(iort,isite)=nevo
            do 105 i=1,nevo
 105           egvalx(i,iort,isite)=egvalx(i,iort,isite)+ci*b0
c     #		egvalx(2,1,isite),egvecx(1,1,1,isite)
c	stop
c
            if ( .not.evalOK(1,egvalx(1,iort,isite),nvar,luout) ) then
c
c   ###                                                            ###
c   ### This is a temporary fix for negative eigenvalues occuring  ###
c   ### in function evaluation.  If this problem happens in OFF-   ###
c   ### diagonal space, it USUALLY means the basis set being used  ###
c   ### is not properly pruned yet and the program issues an error ###
c   ### message and exit.  If the eigenvalues of the diagonal sub- ###
c   ### space have negative real part, the first component of the  ###
c   ### R tensor is increased by 1.0d-4 and the spectra are re-    ###
c   ### calculated.  This seems to be related with the pruning     ###
c   ### scheme for the diagonal basis and further development for  ###
c   ### a new scheme is necessary to fix this problem.             ###
c   ###                     July 19, 1993       Sanghyuk Lee       ###
c   ###                                                            ###
c
               write (luout,1012)
               if (luout.ne.luttyo) write (luttyo,1012)
               ierr=-1
               return
c
            end if
c                     *** solve for eigenvalues & vectors for pS=0
            if (diaflg.or.(iexp.gt.2)) then
c
               if (idebug.ne.0) write (ludeb,2011)
c
               call stvecd(egvalx(1,iort,isite),egvecx(1,1,iort,isite))
               call evcgf(0,stvd,nevd,egvalz(1,iort,isite),
     #                    egvecz(1,1,iort,isite),ierr)
               if (ierr.lt.0) then
                  write(luout,1100)
                  if (luout.ne.luttyo) write(luttyo,1100)
                  return
               end if
c
               nev0(iort,isite)=nevd
c these conditions also eliminated by forcing negative eigenvalues
c to zero:
               if(.not.evalOK(0,egvalz(1,iort,isite),nvar,luout))then
                  write(luout,*)'evalOK error in egvalz, stopping'
                  if (.true.) stop
                  if (iflag.eq.1) then
                     fparm(IDX,ispec,isite)=fparm(IDX,ispec,isite)
     #				+1.0d-4
                     write (luout,1014)
                     if (luout.ne.luttyo) write (luttyo,1014)
                  else
                     write (luout,1016)
                     fparm(ixi,ispec,isite)=fparm(ixi,ispec,isite)
     #			+9.0d0*delx
                     delx=delx*10.0d0
                     if (luout.ne.luttyo) write (luttyo,1016)
                  end if
                  go to 10
               end if
c
            end if
c
 100     continue
c
      end if
c
c---------------------------------------------------------------------
c     calculate complex 2D spectrum in time domain
c---------------------------------------------------------------------
c
      if (icalc.ge.1) then
c
         do 200 j=1,npt2
         do 200 i=1,npt1
 200        xspec(i,j)=czero
c                        ===============================
c                     === loop over orientations (MOMD) ===
c                        ===============================
         do 210 iort=1,nort
            scal=1.0d0
            if ((nort.gt.1).and.((iort.eq.1).or.(iort.eq.nort)))
     #               scal=0.5d0
            nevo=nev1(iort,isite)
            nevd=nev0(iort,isite)
            call spcalc(scal,egvalx(1,iort,isite),egvalz(1,iort,isite),
     #          egvecx(1,1,iort,isite),egvecz(1,1,iort,isite),
     #		xspec,npt1,npt2)
 210     continue
c	stop

c
      end if
c
c---------------------------------------------------------------------
c     Calculate magnitude 2D spectrum with Gaussian and/or Lorentzian
c     inhomogeneous Broadening.  Don't modify xspec.
c---------------------------------------------------------------------
c 
       call convft(xspec,cspec,npt1,npt2)
      return
c
c=====================================================================
c     format statements
c=====================================================================
c
 1000 format(/5x,'** Calculation attempted with illegal ',
     #       'parameters **')
 1012 format(/5x,'** Basis set with looser pruning criterion is ',
     #       'recommended **')
 1014 format(/8x,'First diffusion component is increased by 1e-4',/,
     # 'This should not happen')
 1016 format(/8x,'Forward step size for current Jacobian evaluation ',
     #       'is multiplied by 10'/)
 1100 format(/5x,'** Critical ERROR in SIM2D routine **')
c
 2000 format(/70('#'),/5x,'Orientation ',i2,'  : psi = ',f5.2,
     #       /70('#'))
 2001 format(/'### OFF-DIAGONAL ###')
 2011 format(/'### DIAGONAL ###')
c
c=====================================================================
c
      end
