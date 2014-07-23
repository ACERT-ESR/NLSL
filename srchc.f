c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                       subroutine SRCHC
c                    =========================
c
c   search  <parameter> {spectrum,site} { tol <tol> step <step>
c                         min <smn> max <smx> start <val> }
c
c    spectrum,site : specify which spectra or sites to apply parameter to
c    tol     : Convergence tolerance for fitting parameter
c    step    : Size of first step to take away from initial parameter value  
c    smn,smx : minimum and maximum that defines the search range
c    val     : starting value of the Brent's search procedure
c              (Note that the bracketing procedure is not performed
c               if this option is used.  Therefore the used have to
c               make sure that there exists a minimum chi-squared within
c               the range.)
c
c----------------------------------------------------------------------
      subroutine srchc( line )
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'parms.inc'
      include 'lpnam.inc'
      include 'stdio.inc'
      include 'miscel.inc'
c      include 'iterat.inc'
c
      logical vchange,gotit
      external vchange
      character line*(LINELG)
      integer ixp1p,i,itmp,ierr,lth,ibrckt,ixxsv(NFPRM,MXSPEC,MXSITE)
      double precision ax,bx,cx,fa,fb,fc,fval,xmin,dummy
      character token*(WORDLG),prmID*6
c
      integer nkeywd,ipprm
      parameter(nkeywd=5)
c
      character*8 keywrd(nkeywd)
      data keywrd / 'TOL','STEP','MIN','MAX','START'/
c
      logical ftoken,prmsOK,tcheck
      integer isite,ispec,specinfo,siteinfo
      integer nprmsv,varspec
      double precision brent,p1pfun
      external ftoken,prmsOK,tcheck,brent,p1pfun
c
c######################################################################
c
c test all basis sets here.
      setbas=.true.
      gotit=.true.
      do 1 isite=1,ncomps
c        seteval(isite)=.false.	 ! make no assumption about having matrix
        do 1 ispec=1,nspectra
          setbas=setbas .and. (basinfo(1,ispec,isite) .ne. 0)
          if (basinfo(1,1,isite).ne.basinfo(1,ispec,isite)) then
            write(luttyo,*)'Error, multiple basis sets not allowed'
            write(luttyo,*)'except for different sites.'
            return
          end if
 1    continue
      if (.not.setbas) then
         write (luout,1200)
         if (luout.ne.luttyo) write (luttyo,1200)
 1200    format(/,'** Basis sets not properly set **')
         return
      end if
c
c  -- Check whether any starting parameters have been input before 
c     proceeding with the FIT command
c
      if (.not.prmsOK(luttyo) ) then	! this checks all site/spect
        write(luttyo,1040)
        return
      end if
c
      ierr=0
      ibrckt=0
c
c----------------------------------------------------------------------
c Get the name of the parameter
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
      call touppr(token,lth)
c
c      ixp1p=ipfind(token,lth)
      call ipfind(ixp1p,token,lth,varspec)
c
c----------------------------------------------------------------------
c     Check whether parameter may be varied 
c----------------------------------------------------------------------
      if (ixp1p.eq.0 .or. ixp1p.gt.NFPRM) then
         write(luttyo,1002) token(:lth)
         return
      end if
c
      if (ixp1p.lt.-100) then 
        prmID=alias2( -99-(IGXX+ixp1p) )
      else if (ixp1p.lt.0) then
        prmID=alias1( 1-(IGXX+ixp1p) )
      else
        prmID=parnam(ixp1p)
      end if
c get spec/site information:
      call indtkn(line,specinfo,siteinfo)
      isite=siteinfo
      ispec=specinfo
      if (specinfo.eq.0) ispec=1
      if (siteinfo.eq.0) isite=1
      if (specinfo.ne.0 .and. nspectra.gt.1) then
        write(luttyo,*)'Error, search must apply to all spectra.'
        write(luttyo,*)'Search may apply to one site.'
        return
      end if

c
c  --- Check whether proper symmetry has been specified for
c      tensor components - all spec,sites have same symmetry.
c
      if (.not.tcheck(ixp1p,prmID,luttyo)) return
      ixp1p=iabs( mod(ixp1p,100) )
c
c Set default search range
c
      smn=fparm(ixp1p,ispec,isite)-1.0d0	! ispec must = 1
      smx=fparm(ixp1p,ispec,isite)+1.0d0
      stval=fparm(ixp1p,ispec,isite)
c					! stol and sstep set in nlsinit.
      if (ixp1p.eq.IC20) smn=1.0d-2
c
c----------------------------------------------------------------------
c  Look for a keyword
c----------------------------------------------------------------------
c
   14 call gettkn(line,token,lth)
      lth=min(lth,8)
c
c------------------------------------------------
c************************************************
c No more keywords: call the line search routine
c************************************************
c------------------------------------------------
      if (lth.eq.0) then
c
c         seteval=.false.
c
c Set two initial parameter values for the parameter to be searched
c
         ax=fparm(ixp1p,ispec,isite)
         bx=fparm(ixp1p,ispec,isite)+sstep
c
         if (smn.gt.ax) smn=ax-1.0d0
         if (smx.lt.ax) smx=ax+1.0d0
c
c  Swap indices for the search parameter into the first elements 
c  of the parameter for the NLS x-vector
c  (This is so x can be used by the pfun routine as in a NLS procedure)
c
c proceed as follows: copy ixx, the index identifying which parameters
c are in xxx, into temporary array (don't modify xxx), set all ixx 
c entries = 0, input search variable id into ixx(1) at all applicable 
c sites and spectra.  Check all sites in fvec are equal, then continue.  
c At end of search undo these steps.  Be sure to reset this on return.
c
         nprmsv=nprm
         itmp=ixpr(1)
         ixpr(1)=ixp1p
         nprm=1
         do 10 ispec=1,nspectra
           do 10 isite=1,ncomps
c loop over all parms
             do 10 ipprm=1,NFPRM
	       ixxsv(ipprm,ispec,isite)=ixx(ipprm,ispec,isite)
	       ixx(ipprm,ispec,isite)=0
c if this parmeter, site and spectra are to be searched, set it as xxx(1):
               if (ixp1p.eq.ipprm .and.
     #		 (isite.eq.siteinfo .or. siteinfo.eq.0) .and.
     #		 (ispec.eq.specinfo .or. specinfo.eq.0)) then
                 ixx(ipprm,ispec,isite) = 1	! set this variable
               end if
10       continue
c
         if (ibrckt.eq.0) then
            write (luout,1100) smn,ax,smx
            if (luout.ne.luttyo) write (luttyo,1100) smn,ax,smx
            call mnbrak(ax,bx,cx,fa,fb,fc,p1pfun,smn,smx,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) go to 202	! bad pfun return
            write (luout,1105) cx,bx,ax,fc,fb,fa
            if (luout.ne.luttyo) write (luttyo,1105) cx,bx,ax,fc,fb,fa
            if(ierr.ne.0) then
              write(luout,*)'error returned from mnbrak ',ierr
            if (luout.ne.luttyo) write (luttyo,*)
     #		'error returned from mnbrak ',ierr
              go to 202
            end if
         else
            ax=smn
            bx=stval
            cx=smx
            write (luout,1100) smn,stval,smx
            if (luout.ne.luttyo) write (luttyo,1100) smn,stval,smx
         end if
c
         if (ierr.eq.0) then
            dummy=brent(ax,bx,cx,p1pfun,stol,xmin,prmID,ierr)
            if (ierr.lt.0 .or. ihltcmd.ne.0) go to 202
         else
            gotit=.false.
            write (luout,1110)
            if (luout.ne.luttyo) write (luttyo,1110)
            if (ax.gt.bx) then
               write (luout,1120) cx,ax,fc,fa
               if (luout.ne.luttyo) write (luttyo,1120) cx,ax,fc,fa
            else
               write (luout,1120) ax,cx,fa,fc
               if (luout.ne.luttyo) write (luttyo,1120) ax,cx,fa,fc
            end if
         end if
c put the parameters back, put the new found value into fparm
         nprm=nprmsv
	 ixpr(1)=itmp
         do 20 ispec=1,nspectra
           do 20 isite=1,ncomps
             do 20 ipprm=1,NFPRM
	       ixx(ipprm,ispec,isite)=ixxsv(ipprm,ispec,isite)
               if (ixp1p.eq.ipprm .and.
     #		 (isite.eq.siteinfo .or. siteinfo.eq.0) .and.
     #		 (ispec.eq.specinfo .or. specinfo.eq.0)) then
c put back new found minimum, No need to change basinfo, was set when
c spectrum was calculated
c                   if (basinfo(2,ispec,isite) .ge. varspec .and.
c     #               (abs(fparm(iparmptr,ispec,isite)-fval)) .gt.
c     #               1.d-8*abs(fval)) then
c                     basinfo(2,ispec,isite)= varspec-1
c                   end if
c
c if there is a parameter change, set matrixok for new calculation
c
                 if (gotit.and.vchange(xmin,
     #			fparm(ixp1p,ispec,isite))) then
     	  	   matrxok(isite)=min(matrxok(isite),specid(ixp1p))
                   fparm(ixp1p,ispec,isite)=xmin
                 end if
               end if
20       continue
         return
c error return, reset ixx first:
 202     nprm=nprmsv
	 ixpr(1)=itmp
         do 201 ispec=1,nspectra
           do 201 isite=1,ncomps
             do 201 ipprm=1,NFPRM
	       ixx(ipprm,ispec,isite)=ixxsv(ipprm,ispec,isite)
 201	 continue
         return
      end if
c
c------------------------------
c Check keyword list (lth .ne. 0)
c------------------------------
      call touppr(token,lth)
      do 15 i=1,nkeywd
         if (token(:lth).eq.keywrd(i)(:lth)) goto 16
 15   continue
c                                        *** Unrecognized keyword
      write (luttyo,1000) token(:lth)
      return
c
c----------------------------------------------------------------------
c  Keyword found: for keywords requiring an argument, convert 
c  next token and assign appropriate value
c----------------------------------------------------------------------
 16   call gettkn(line,token,lth)
c                                        *** No value given
      if (lth.eq.0) then
         write(luttyo,1003) keywrd(i)
         return
      end if
c
      if (ftoken(token,lth,fval)) then
c                                          *** TOL keyword
         if (i.eq.1) then
            stol=fval
c                                          *** STEP keyword
         else if (i.eq.2) then 
            sstep=fval
c                                          *** MIN keyword
         else if (i.eq.3) then
            smn=fval
c                                          *** MIN keyword
         else if (i.eq.4) then
            smx=fval
c                                          *** START keyword
         else if (i.eq.5) then
            ibrckt=1
            stval=fval
c                                      *** Illegal numeric value
         end if
      else
         write(luttyo,1001) token(:lth)
      end if
      go to 14
c
c ##### Format statements ########################################
c
 1000 format('*** Unrecognized SEARCH keyword: ''',a,''' ***')
 1001 format('*** Numeric value expected: ''',a,''' ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
 1003 format('*** No value given for ''',a,''' ***')
 1040 format('*** Initial parameter values are required before',
     #      ' a fit ***')
 1100 format(/4x,'Initial Range  : ',g13.6,2x,g13.6,2x,g13.6)
 1105 format(4x,'Search Bracket : ',g13.6,2x,g13.6,2x,g13.6,
     #      /4x,'Rms Deviation  : ',g13.6,2x,g13.6,2x,g13.6)
 1110 format('*** No minimum found within the valid range ***')
 1120 format(4x,'Last Bracket   : ',g13.6,2x,g13.6,
     #      /4x,'Rms Deviation  : ',g13.6,2x,g13.6)
      end
