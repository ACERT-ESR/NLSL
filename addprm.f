c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                       subroutine ADDPRM
c                    =========================
c
c Add a parameter to list of parameters being varied for nonlinear
c least-squares. Also maintain the lists of 
c    (1) ibnd    : boundary flag for variable parameter
c    (2) prmin   : minimum for variable parameter
c    (3) prmax   : maximum for variable parameter
c    (4) prscl   : desired accuracy for given parameter
c    (5) xfdstp  : Size of forward-differences step
c    (6) tag     : Name of the parameter (character*6)
c    (7) ixx     : index of each element in the fparm array into the x vector
c
c The values of these parameters for the given variable parameter
c are passed to ADDPRM.
c
c ix - returned by ipfindc
c ibd - ibnd   : Flag for boundaries imposed on each parameter:
c          0=none, 1=minimum, 2=maximum, 3=both
c prmn,prmx,prsc,step - min, max, accuracy, fd-step
c prmID - parameter prefix
c ixrspec,ixrsite - specify spectra and sites to vary via indtkn.
c
c----------------------------------------------------------------------
      subroutine addprm(ixr,ibd,prmn,prmx,prsc,step,prmID,
     #		ixrspec,ixrsite)
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'lmcomm.inc'
      include 'stdio.inc'
      include 'rndoff.inc'
      include 'miscel.inc'
c
      integer ixr,ibd,ispec,isite,specptr,siteptr,iparmptr,
     #	ixrspec,ixrsite
      double precision prmn,prmx,prsc,step
      character*6 prmID
      integer i,ixabs,j,k,k1,nmov
c
      logical tcheck
      external tcheck
c
c######################################################################
c
c----------------------------------------------------------------------
c See if there is enough room in the parameter list
c----------------------------------------------------------------------
      if (nprm+1.gt.MXVAR) then
         write(luout,1001) MXVAR
         if (luout.ne.luttyo) write(luttyo,1001) MXVAR
         return
      end if
c Check axis specification for consistency
      if (.not.tcheck(ixr,prmID,luttyo)) then
         write(*,*)'axis information inconsistent in addprm'
         return
      end if
c
c----------------------------------------------------------------------
c Search parameter list to see if parameter is already there
c Must check specific site/spectrum information here.  
c----------------------------------------------------------------------
c 
c      ixabs=abs(mod(ixr,100))
      specptr=1                 ! get the old information
      siteptr=1
      if (ixrspec.ne.0) specptr=ixrspec         ! get first spec
      if (ixrsite.ne.0) siteptr=ixrsite         ! get first site
      iparmptr=iabs(mod(ixr,100))
c consider all sites, spec:
      do 30 isite=1,ncomps
        do 30 ispec=1,nspectra
c check matching sites with requested site
          if ((isite.eq.ixrsite .or. ixrsite.eq.0) .and. 
     #       (ispec.eq.ixrspec .or. ixrspec.eq.0)) then
c if vary previously selected
            if (ixx(iparmptr,ispec,isite) .ne. 0) then
               write(luttyo,1000) prmID
 1000          format('*** Parameter ''',a,''' already varied ***')
	       return
            end if
          end if
 30   continue
c so add parameter to list
c We use ixpr to keep track of the parameters in the list of 
c variable paremeters. 
c
c      do 10 i=1,nprm
c                                      *** parameter already in list
c         if (ixpr(i).eq.ixabs) then
c            write(luout,1000) prmID
c            if (luout.ne.luttyo) write(luttyo,1000) prmID
c 1000       format('*** Parameter ''',a,''' already being varied ***')
c            return
c         else if (ixpr(i).gt.ixabs) then
c            go to 14
c         end if
c 10   continue
c
c----------------------------------------------------------------------
c Reached end of list: parameter may be added to the end
c----------------------------------------------------------------------
      nprm=nprm+1
c
c----------------------------------------------------------------------
c Found proper location for new parameter: move top of list up
c----------------------------------------------------------------------
c 14   nmov=nprm-i+1
c      do 16 j=1,nmov
c         k=nprm-j+1
c         k1=k+1
c         x(k1)=x(k)
c         ixpr(k1)=ixpr(k)
c         prmin(k1)=prmin(k)
c         prmax(k1)=prmax(k)
c         prscl(k1)=prscl(k)
c         xfdstp(k1)=xfdstp(k)
c         ibnd(k1)=ibnd(k)
c         tag(k1)=tag(k)
c         ixx(ixpr(k1))=ixx(ixpr(k1))+1
c 16   continue
c      nprm=nprm+1
c
c----------------------------------------------------------------------
c  Add new parameter(s) to proper location list
c----------------------------------------------------------------------
 18   xxx(nprm)=fparm(iparmptr,specptr,siteptr)
      ixpr(nprm)=iparmptr
      do 20 isite=1,ncomps
        do 20 ispec=1,nspectra
c check matching sites with requested site
          if ((isite.eq.ixrsite .or. ixrsite.eq.0) .and. 
     #       (ispec.eq.ixrspec .or. ixrspec.eq.0)) then
            ixx(iparmptr,ispec,isite)=nprm
c check that all fparm values are equal.
            if ((abs(fparm(iparmptr,specptr,siteptr)-fparm(iparmptr,
     #		ispec,isite))).gt.0.00001*abs(fparm(
     #		iparmptr,specptr,siteptr)))then
            write(luttyo,1110) prmID
	 write(*,*)'values ',fparm(iparmptr,specptr,siteptr),
     #		fparm(iparmptr,ispec,isite)
 1110       format('*** Parameter ''',a,''' not consistent ***')
            end if
          end if
 20   continue
c
      prmin(nprm)=prmn
      prmax(nprm)=prmx
      prscl(nprm)=prsc
      xfdstp(nprm)=step*xxx(nprm)
      if(dabs(xfdstp(nprm)).lt.1000.0d0*rndoff) xfdstp(nprm)=step
      ibnd(nprm)=ibd
      tag(nprm)=prmID
c
      write(luout,1003) nprm
      if (luout.ne.luttyo) write(luttyo,1003) nprm
      return
c
c #### format statements ###########################################
c
 1001 format('*** Limit of',i2,' variable parameters exceeded ***')
 1003 format(i3,' parameters now being varied')
      end
