c  VERSION 1.0  (NLSPMC version)   2/5/99
c**********************************************************************
c                    =========================
c                          file CMDS
c                    =========================
c
c  Contains the following routines:
c    varyc(line)    Interprets given line as an NLS "vary" command
c    fixc(line)     Interprets given line as an NLS "fix" command
c    basis(line)    Interprets given line as an NLS "basis" command
c    series(line)   Interprets given line as an NLS "series" command
c    fitc(line)     Interprets given line as an NLS "fit" command
c    parc(line,lu)  Interprets line as NLS "parms" command and prints
c                    parameter values to the logical unit specified
c    statc(line,lu) Interprets line as NLS "status" command
c
c    Includes
c      limits.inc
c      simparm.inc
c      datas.inc
c      parcom.inc
c      lpnam.inc
c      lmcomm.inc
c      stdio.inc
c
c    Uses
c      ipfind
c      gettkn
c      ftoken
c      itoken
c
c----------------------------------------------------------------------
c                    =========================
c                      subroutine VARYC
c                    =========================
c
c vary <name> {  {sprectrum,site} minimum <minval> maximum <maxval> 
c                scale <scl> fdstep <step> }
c
c      name    : name of the parameter to be varied
c      sprectrum,site : parameters or * or absent, identifying same. 
c      minval  : Minimum permissible value for parameter
c      maxval  : Maximum permissible value for paramete
c      scale   : Factor by which to scale search vector for this parameter
c      step    : Relative step size for forward-differences approximation
c
c----------------------------------------------------------------------
      subroutine varyc(line)
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'parmequ.inc'
      include 'parms.inc'
      include 'lpnam.inc'
      include 'stdio.inc'
      include 'miscel.inc'

c
      character*(LINELG) line
      integer i,ibd,ix,lth,specinfo,siteinfo
      double precision fval,prmn,prmx,prsc,step
      character token*(WORDLG),prmID*6
c
      integer nkeywd
      parameter(nkeywd=4)
      character*8 keywrd(nkeywd)
      data keywrd /'MINIMUM','MAXIMUM','SCALE','FDSTEP'/
c
      double precision one,zero
      parameter (one=1.0d0,zero=0.0d0)
c
      integer varspec
      logical ftoken
      external ftoken
c
c----------------------------------------------------------------------
c Get the name of the parameter
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
      call touppr(token,lth)
c
      call ipfind(ix,token,lth,varspec)
c
c----------------------------------------------------------------------
c     Check whether parameter may be varied 
c----------------------------------------------------------------------
      if (ix.eq.0 .or. ix.gt.nvprm) then
         write(luttyo,1002) token(:lth)
         return
      end if
c
      if (ix.lt.-100) then 
         prmID=alias2( -99-(IGXX+ix) )
      else if (ix.lt.0) then
         prmID=alias1( 1-(IGXX+ix) )
      else
         prmID=parnam(ix)
      end if
c
c get parenthisized information
c
      call  indtkn(line,specinfo,siteinfo)
c
c set defaults
c
      prmn=zero
      prmx=zero
      prsc=one
      step=1.0D-5
      ibd=0
c
c----------------------------------------------------------------------
c  Look for a keyword
c----------------------------------------------------------------------
   14 call gettkn(line,token,lth)
      lth=min(lth,8)
c
c------------------------------------------------
c No more keywords: add the parameter and exit
c------------------------------------------------
      if (lth.eq.0) then
c         call addprm(ix,ibd,prmn,prmx,prsc,step,prmID)
      if (specinfo.ne.0 .and. nspectra.gt.1) then
        write(luttyo,*)'Vary parameter must apply to all spectra.'
        write(luttyo,*)'but may apply to only one site if you wish.'
        return
      end if
c
        call addprm(ix,ibd,prmn,prmx,prsc,step,prmID,specinfo,siteinfo)
        return
      end if
c
c------------------------------
c Check keyword list
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
c  Keyword found: convert next token and assign appropriate value
c----------------------------------------------------------------------
 16   call gettkn(line,token,lth)
c                                        *** No value given
      if (lth.eq.0) then
         write(luttyo,1003) keywrd(i)
         return
      end if
c
      if (ftoken(token,lth,fval)) then
c                                        *** MINIMUM keyword
         if (i.eq.1) then
            prmn=fval
            if (mod(ibd,2).eq.0) ibd=ibd+1
c                                        *** MAXIMUM keyword
         else if (i.eq.2) then
            prmx=fval
            ibd=ibd+2
c                                        *** SCALE keyword
         else if (i.eq.3) then
            prsc=fval
c                                        *** FDSTEP keyword
         else if (i.eq.4) then
            step=fval
         end if
c                               *** Illegal real value
      else
         write(luttyo,1001) token(:lth)
      end if
c
      go to 14
c
c ###### format statements ########################################
c
 1000 format('*** Unrecognized VARY keyword: ''',a,''' ***')
 1001 format('*** Real value expected: ''',a,''' ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
 1003 format('*** No value given for ''',a,''' ***')

      end


c----------------------------------------------------------------------
c                    =========================
c                      subroutine FIXC
c                    =========================
c
c  Removes the given parameter name from the list of variable parameters.
c
c fix <name>
c
c      name    : name of the parameter to be removed
c----------------------------------------------------------------------
      subroutine fixc(line)
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'parmequ.inc'
      include 'parms.inc'
c      include 'lpnam.inc'
      include 'stdio.inc'
      include 'miscel.inc'
c
      character*(LINELG) line
      integer i,j,lth,specinfo,siteinfo,isite,ispec,specptr,siteptr
      character token*(LINELG),prmID*6
c
      integer ixr,ixabs,iprm,varspec,ident
      logical ftoken,okparm
      external ftoken
c
c----------------------------------------------------------------------
c Get the name of the parameter
c----------------------------------------------------------------------
      call gettkn(line,token,lth)
      lth=min(lth,6)
      if (lth.le.0) return
c
      call touppr(token,lth)
c
      if (token(:lth).eq.'ALL') then
c                                           ** Fix all parameters
        nprm=0
	do 1 iprm=1,nvprm
         do 1 isite=1,ncomps
          do 1 ispec=1,nspectra
c zero parameter being fixed, decrement location for moved parms
            ixx(iprm,ispec,isite)=0
 1      continue
        return
      end if
c
      call ipfind(ixr,token,lth,varspec)
      ixabs=abs(mod(ixr,100))	! token id without axis info
c
c----------------------------------------------------------------------
c     Check whether parameter may be varied 
c----------------------------------------------------------------------
      if (ixabs.eq.0 .or. ixabs.gt.nvprm) then
         write(luttyo,1002) token(:lth)
         return
      end if
c
c get parenthisized information
c
      call  indtkn(line,specinfo,siteinfo)
c
c fix all relevent sites
c
      okparm=.false.
      specptr=1                 ! get the old information
      siteptr=1
      if (specinfo.ne.0) specptr=specinfo         ! get first spec
      if (siteinfo.ne.0) siteptr=siteinfo         ! get first site
c consider all sites, spec:
      do 30 isite=1,ncomps
        do 30 ispec=1,nspectra
c check matching sites with requested site
          if ((isite.eq.siteinfo .or. siteinfo.eq.0) .and.
     #       (ispec.eq.specinfo .or. specinfo.eq.0)) then
c if vary previously selected
               if (ixx(ixabs,ispec,isite) .ne. 0) then
                 ident=ixx(ixabs,ispec,isite)
                 okparm=.true.
                 call rmvprm(ident)
c all sites/spectra should be fixed with the first call to rmvprm.
               end if
          end if
 30   continue
      if (.not.okparm) then
        if (specinfo.eq.0 .and. siteinfo.eq.0) then
          write(luttyo,1003) token(:lth)
        else if (specinfo.eq.0) then
          write(luttyo,1004) token(:lth),siteinfo
        else if (siteinfo.eq.0) then
          write(luttyo,1005) token(:lth),specinfo
        else 
          write(luttyo,1006) token(:lth),specinfo,siteinfo
        end if
      end if
      return
c
c ###### format statements ########################################
c
 1002 format('*** ''',a,''' is not a variable parameter ***')
 1003 format('*** ''',a,'(*,*)',''' is not a variable parameter ***')
 1004 format('*** ''',a,' *,',i1,''' is not a variable parameter ***')
 1005 format('*** ''',a,i2,',*',''' is not a variable parameter ***')
 1006 format('*** ''',a,'(',i2,',',i1,')',
     #	''' is not a variable parm ***')
c
      end


c----------------------------------------------------------------------
c                    =========================
c                      subroutine BASIS
c                    =========================
c
c Generates a basis set to be used in simulation.
c The program just reads the index file if it exists.
c If the index file does not exist or the filename is not specified,
c the full basis vectors within the maximum allowed indices (lemx,
c lomx,kmx,mmx,ipnmx) are generated.  The eigenvector calculation
c with the full basis vectors is very costly especially with the
c Rutishauser algorithm, and the pruned basis set should be used
c as much as possible.
c  
c basis {<name> {site,spectrum} print }
c
c      name    : name (in full) of the basis index file
c      print   : prints out the indices on the screen
c
c----------------------------------------------------------------------
      subroutine basis(line,lu)
      implicit none
c
      include 'limits.inc'
      include 'names.inc'
      include 'simparm.inc'
      include 'parms.inc'
      include 'parmequ.inc'
      include 'basis.inc'
      include 'miscel.inc'
      include 'stdio.inc'
      character*(LINELG) line
c
      integer i,n,lth,lu,ierr
      character token*(WORDLG)
c
      integer ipar,iprt,specinfo,siteinfo,isite,ispec,specptr,siteptr
      external ipar
c
      iprt=0
c
c get parenthisized information, if any:
c
      call  indtkn(line,specinfo,siteinfo)
c
c store basis id # in basinfo
c
      if (nbasis.lt.nspectra*ncomps) then
        nbasis=nbasis+1
      else		! Too many basis sets
        if (lu.ne.0) write (lu,1000)
        if (lu.ne.luttyo) write (luttyo,1000)
	return
      end if
c
      specptr=1                 ! get the old information
      siteptr=1
      if (specinfo.ne.0) specptr=specinfo         ! get first spec
      if (siteinfo.ne.0) siteptr=siteinfo         ! get first site
c
c Set up lemx etc for calculation from this site:
c (Presumes this has been entered before the basis command)
c
      do 28 i=1,NFPRM
        fepr(i)=fparm(i,specptr,siteptr)
 28   continue
      do 4 i=1,NIPRM
        iepr(i)=iparm(i,specptr,siteptr)
 4    continue
c
c----------------------------------------------------------------------
c Get the name of the basis set if any:
c----------------------------------------------------------------------
c
 10   call gettkn(line,ixname,lth)
c
      if (lth.ne.0) then
        if ((ixname(:lth).eq.'print').or.(ixname(:lth).eq.'PRINT')) then
           iprt=1
           go to 10
        end if
      end if
c set up mjqe1 etc:
        if (lth.ne.0) then
          call fbasis(ixname,0,ierr)	! read file
        else
          call fbasis(ixname,1,ierr)	! full basis
          if (lu.ne.0) write (lu,1001)
          if (lu.ne.luttyo) write (luttyo,1001)
        end if
c
      if (ierr.ne.0) then
         if (ierr.eq.1) then
            if (lu.ne.0) then
               write (lu,1002) ixname(:lth)	! file not exist
               write (lu,1001)
            end if
            if (lu.ne.luttyo) then
               write (luttyo,1002) ixname(:lth)
               write (luttyo,1001)
            end if
         else if (ierr.eq.2) then 	! dimension exceeded
            if (lu.ne.0) then
               write (lu,1003) mxdim
               write (lu,1100) ndimoss(nbasis),ndimdss(nbasis)
            end if
            if (lu.ne.luttyo) then
               write (luttyo,1003) mxdim
               write (luttyo,1100) ndimoss(nbasis),ndimdss(nbasis)
            end if
            return
         else if (ierr.eq.3) then	! mag parms not set
            if (lu.ne.0) write (lu,1004)
            if (lu.ne.luttyo) write (luttyo,1004)
            return
         end if
      end if
c
c----------------------------------------------------------------------
c     Obtain MTS index for this basis
c----------------------------------------------------------------------
c
      lemx=-1
      lomx=-1
      kmx=-1
      mmx=-1
      ipnmx=-1
c loop over this basis set:
      do 210 i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1
         n=ml1(i)
         if (ipar(n).eq.1) then
            if (n.gt.lemx) lemx=n
         else
            if (n.gt.lomx) lomx=n
         end if
         n=mk1(i)
         if (n.gt.kmx) kmx=n
         n=abs(mm1(i))
         if (n.gt.mmx) mmx=n
         n=abs(mpi1(i))
         if (n.gt.ipnmx) ipnmx=n
 210  continue
c
      if (lomx.eq.-1) lomx=0
      if (kmx.eq.-1) kmx=0
      if (mmx.eq.-1) mmx=0
      if (ipnmx.eq.-1) ipnmx=0
c
      write (lu,1100) ndimoss(nbasis),ndimdss(nbasis)
      write (lu,1110) lemx,lomx,kmx,mmx,ipnmx
      if (lu.ne.luttyo) then
         write (luttyo,1100) ndimoss(nbasis),ndimdss(nbasis)
         write (luttyo,1110) lemx,lomx,kmx,mmx,ipnmx
      end if
c
c store the lemx etc in relevent sites/spectra
c
      do 40 isite=1,ncomps
        do 40 ispec=1,nspectra
c Fill matching sites with correct contents:
          if ((isite.eq.siteinfo .or. siteinfo.eq.0) .and.
     #       (ispec.eq.specinfo .or. specinfo.eq.0)) then
            iparm(ILEMX,ispec,isite)=lemx
            iparm(ILEMX+1,ispec,isite)=lomx
            iparm(ILEMX+2,ispec,isite)=kmx
            iparm(ILEMX+3,ispec,isite)=mmx
            iparm(ILEMX+4,ispec,isite)=ipnmx
            matrxok(isite)=0	! need new matrix!
c store basis info for this site/spec
	    basinfo(1,ispec,isite)=nbasis
c store basis size information
            iparm(INDIMO,ispec,isite)=ndimoss(nbasis)
            iparm(INDIMD,ispec,isite)=ndimdss(nbasis)
c
          end if
 40   continue
c
      call gettkn(line,token,lth)
      if ((iprt.ne.1).and.(lth.le.0)) return
c
      write(lu,1200)
      do 100 i=pidptr(nbasis),pidptr(nbasis)+ndimoss(nbasis)-1
         write(lu,1210) i-pidptr(nbasis)+1,ml1(i),mjk1(i),mk1(i),
     #                 mjm1(i),mm1(i),mpi1(i),mqi1(i)
 100  continue
c
c      write(lu,1220)
c      do 110 i=dpidptr(nbasis),dpidptr(nbasis)+ndimdss(nbasis)-1
c         write(lu,1210) i-dpidptr(nbasis)+1,mdl1(i),mdjk1(i),mdk1(i),
c     #                  mdjqe1(i),mdm1(i),mdpi1(i),mdqi1(i)
c 110  continue
c
      return
c
c ###### format statements ############################################
c
 1000 format('** Too many basis sets requested **')
 1001 format('** Full basis set is generated **')
 1002 format(' Index file : ',a,' does not exist')
 1003 format(' Dimension of matrix too large : maximum dimension = ',
     #     i4)
 1004 format(' magnetic parameters are not properly set')
 1100 format(' Dimension of matrix: ndimoss = ',i4,' ndimdss = ',i4)
 1110 format(' MTS lemx,lomx,kmx,mmx,ipnmx = ',i3,4(',',i2))
 1200 format(/,15x,'BASIS SET : off-diagonal space',//,2x,
     #     'element',7x,'L  jK   K  jqM  M  pI  qI',/,2x,42('-'))
 1210 format(4x,i5,4x,'|',6(i3,','),i3,'  >')
 1220 format(/,15x,'BASIS SET : diagonal space',//,2x,'element',7x,
     #     'L  jK   K  jqM  M  pI  qI',/,2x,42('-'))
c
      end


c----------------------------------------------------------------------
c                    =========================
c                       subroutine SERIES
c                    =========================
c
c Interpret a "series" command from the given line.
c The command is of the form:
c
c    series <name> {=} <value>, <value>, ...
c 
c For now allowed series parameters apply to all sites.
c
c         name  : Name of a variable parameter
c         value : List of values assumed by the named parameter
c                 for each spectrum in the series
c    Note that nspectra (number of series calculations)
c    values are expected, if fewer, the last is repeated. 
c
c----------------------------------------------------------------------
      subroutine series( line )
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'miscel.inc'
      include 'stdio.inc'
c
      character line*(LINELG)
      logical serprm
      integer itmp,lth,lprm,ix,i
      double precision fval,serval(MXSPEC)
      character token*(WORDLG),prmID*(6),defnm*(6)
c
      logical ftoken
      external ftoken
c
c----------------------------------------------------------------------
c     Get name of series parameter
c----------------------------------------------------------------------
      call gettkn(line,prmID,lth)
      lth=min(lth,6)
c                               *** No series variable specified
      if (lth.eq.0) then
         write(luttyo,1002)
         return
      end if
c
c-----------------------------------------------------
c Check if the specified parameter is series variable.
c-----------------------------------------------------
      call touppr(prmID,lth)
c
      serprm = (prmID(:lth).eq.'IEXP').or.(prmID(:lth).eq.'ICOMB')
     #    .or.(prmID(:lth).eq.'NPT1').or.(prmID(:lth).eq.'NPT2')
     #    .or.(prmID(:lth).eq.'INIT1').or.(prmID(:lth).eq.'INIT2')
     #    .or.(prmID(:lth).eq.'STEPT1').or.(prmID(:lth).eq.'STEPT2')
     #    .or.(prmID(:lth).eq.'TFIX').or.(prmID(:lth).eq.'WEIGHT')
c
      lprm=lth
      if (.not.serprm) then
         write (luttyo,1003) prmID(:lth)
         return
      end if
c
      itmp=0
c
c----------------------------------------------------------------------
c         Get a list of values for the series parameter
c----------------------------------------------------------------------
   10 call gettkn(line,token,lth)
      if (lth.ne.0) then		! have a value or =
c
        if (token(:lth).eq.'=') goto 10
c 
        if (ftoken(token,lth,fval)) then	! if is a fp #
          if (itmp.ge.nspectra) then
            write (luttyo,1005) nspectra
 1005       format('*** SERIES may not have more than',i2,
     #       ' values: did you set nspectra? ***')

            return
          end if
          itmp=itmp+1
          serval(itmp)=fval
c                                *** Illegal real number
        else
          write(luttyo,1006) token(:lth)
 1006     format('*** ',a,' : Not a real value ***')
        end if
c
        go to 10   
      else
c
c------------------------------
c No more values--map them.
c------------------------------
c                                   ** Set nser here **
c if this is the first call to series, set nser
        if (nsparm.eq.0) nser=nspectra
c
c if later calls, fill in if list is short
        if (itmp.lt.nser) then
            write(luttyo,1004)
            do 11 ix=itmp+1,nspectra
               serval(ix)=serval(itmp)
 11         continue
        end if
c
c--------------------------------------------------
c Store current series variables into proper array
c--------------------------------------------------
c
         if (prmID(:lprm).eq.'IEXP') then
            do 20 ix=1,nser
c diaflg indicates ANY of the experiments are type greater than 
c or equal to 3.  
               if (serval(ix).ge.3) diaflg=.true.
 20            siexp(ix)=nint(serval(ix))
         else if (prmID(:lprm).eq.'ICOMB') then
            do 21 ix=1,nser
 21            sicomb(ix)=nint(serval(ix))
         else if (prmID(:lprm).eq.'NPT1') then
            do 22 ix=1,nser
 22            snpt1(ix)=nint(serval(ix))
         else if (prmID(:lprm).eq.'NPT2') then
            do 23 ix=1,nser
 23            snpt2(ix)=nint(serval(ix))
         else if (prmID(:lprm).eq.'INIT1') then
            do 24 ix=1,nser
 24            sinit1(ix)=serval(ix)
         else if (prmID(:lprm).eq.'INIT2') then
            do 25 ix=1,nser
 25            sinit2(ix)=serval(ix)
         else if (prmID(:lprm).eq.'STEPT1') then
            do 26 ix=1,nser
 26            sstept1(ix)=serval(ix)
         else if (prmID(:lprm).eq.'STEPT2') then
            do 27 ix=1,nser
 27            sstept2(ix)=serval(ix)
         else if (prmID(:lprm).eq.'TFIX') then
            do 28 ix=1,nser
 28            stfix(ix)=serval(ix)
         else if (prmID(:lprm).eq.'WEIGHT') then
            do 29 ix=1,nser
 29            datwgt(ix)=serval(ix)
         end if
c                * update # of series parameters set so far *
         nsparm=nsparm+1
c                * initialize some arrays, after have snpt1, snpt2
c This assumes we will be reading experimental data sets.  This is not 
c necessary in the case of no data, but we don't know that at this time.
c 
         if (nsparm.ge.9) then		! got all series parms
            ix=1
            defnm='simu'
            do 30 i=1,nspectra
               ixsp(i)=ix
               ndata(i)=snpt1(i)*snpt2(i)
               ix=ix+ndata(i)
               if (i.lt.10)write(defnm(5:5),'(i1)') i
               if (i.gt.10)write(defnm(5:6),'(i2)') i
               dataid(i)=defnm
 30         continue
         end if
c
         return
c
      end if
c
c
c ###### format statements ######################################
c
 1001 format('*** nser NOT set: nser=1 assumed ***')
 1002 format('*** No SERIES variable specified ***')
 1003 format('*** ',a,' : Not a SERIES variable ***')
 1004 format('*** Not enough SERIES values : last repeated ***')
c
      end


c----------------------------------------------------------------------
c                    =========================
c                       subroutine FITC
c                    =========================
c
c   fit      { trace xtol <xtol> ftol <ftol> gtol <ftol> 
c              maxfun <mxf> maxitr <mxi> maxfit <mxfit> bound <factor> }
c
c          trace:  Specifies that a ".trc" file should be produced for
c                  the fit
c          xtol:   Convergence tolerance for scaled fitting parameters
c          ftol:   Convergence tolerance for chi-squared
c          gtol:   Convergence tolerance for gradient of chi-squared with
c                  respect to the fitting parameters
c          mxf:    Maximum number of function calls allowed
c          mxi:    Maximum number of iterations allowed
c          mxfit:  Maximum number of fit commands allowed until convergence
c          factor: Factor defining initial step bound used in parameter search
c
c----------------------------------------------------------------------
      subroutine fitc( line )
      implicit none
c
      include 'limits.inc'
      include 'lmcomm.inc'
      include 'parms.inc'
c      include 'iterat.inc'
      include 'miscel.inc'
      include 'stdio.inc'
c
      character line*(LINELG)
c
      logical ftoken
      external ftoken
      integer i,lth,ifit,ispec,isite
      double precision fval
      character*(WORDLG) token
c
      integer nkeywd
      parameter(nkeywd=9)
c 
      character*8 keywrd(nkeywd)
      data keywrd / 'FTOL', 'GTOL', 'XTOL', 'BOUND', 'MAXFUN',
     #              'MAXITR', 'MAXFIT', 'TRACE', 'JACOBI' /
c
      logical prmsOK
      external prmsOK
c
c######################################################################
c
c  -- Check whether any starting parameters have been input before 
c     proceeding with the FIT command
c
      if (.not.prmsOK(luttyo) ) then
        write(luttyo,1040)
 1040   format('*** Initial parameter values are required before',
     #      ' a fit ***')
        return
      end if
c  check basis set settings:
      do 1 isite=1,ncomps
        do 1 ispec=1,nspectra
          if (basinfo(1,1,isite).ne.basinfo(1,ispec,isite)) then
            write(luttyo,*)'Error, multiple basis sets not allowed'
            write(luttyo,*)'except for different sites.'
            return
          end if
 1    continue
c
      ifit=0
c----------------------------------------------------------------------
c  Look for a keyword
c----------------------------------------------------------------------
c
 14   call gettkn(line,token,lth)
      lth=min(lth,8)
c
c   ***** THIS IS IT!!! ************
c------------------------------------------------
c No more keywords: call the NLS fitting routine
c------------------------------------------------
      if (lth.eq.0) then
c
c        call catchc( ihltcmd )
c
 15     call fitp
c
        if (ihltcmd.ne.0) then
c          call uncatchc( ihltcmd )
          return
        else        
          ifit=ifit+1
          if ((ifit.lt.maxfit).and.((info.eq.0).or.(info.gt.3)))
     #                   go to 15
        end if
c
        return
      end if
c
c------------------------------
c Check keyword list
c------------------------------
      call touppr(token,lth)
      do 16 i=1,nkeywd
        if (token(:lth).eq.keywrd(i)(:lth)) goto 17
 16   continue
c                                        *** Unrecognized keyword
      write (luttyo,1000) token(:lth)
      return
c
c----------------------------------------------------------------------
c  Keyword found: for keywords requiring an argument, convert 
c  next token and assign appropriate value
c----------------------------------------------------------------------
 17   if (i.ge.1 .and.i.le.7) then
        call gettkn(line,token,lth)
c                                        *** No value given
        if (lth.eq.0) then
          write(luttyo,1003) keywrd(i)
          return
        end if
c
        if (ftoken(token,lth,fval)) then
c                                          *** FTOL keyword
          if (i.eq.1) then
            ftol=fval
c                                          *** GTOL keyword
          else if (i.eq.2) then
            gtol=fval
c                                          *** XTOL keyword
          else if (i.eq.3) then
            xtol=fval
c                                          *** BOUND keyword
          else if (i.eq.4) then
            factor=fval
c                                          *** MAXFUN keyword
          else if (i.eq.5) then
            maxev=int(fval)
c                                          *** MAXITR keyword
          else if (i.eq.6) then
            maxitr=int(fval)
c                                          *** MAXFIT keyword
          else if (i.eq.7) then
            maxfit=int(fval)
          end if
c                                          *** Illegal numeric value
        else
          write(luttyo,1001) token(:lth)
        end if
c                                          *** TRACE keyword
      else if (i.eq.8) then
        if (luout.eq.luttyo) then
          write (luttyo,1050)
          itrace=0
        else
          itrace=1
        end if
c                                          *** JACOBI keyword
      else if (i.eq.9) then
        jacobi=1
      end if
c
      go to 14
c
 1000 format('*** Unrecognized FIT keyword: ''',a,''' ***')
 1001 format('*** Numeric value expected: ''',a,''' ***')
 1003 format('*** No value given for ''',a,''' ***')
 1050 format('*** A log file must be opened before using TRACE ***')
      end


c----------------------------------------------------------------------
c                    =========================
c                       subroutine PARC 
c                    =========================
c
c  Prints out a list of parameter values on the given logical unit
c  number.
c  Add lulog to do parms conversion and output parms into log file. 
c----------------------------------------------------------------------
      subroutine parc( line )
      implicit none
c
      include 'limits.inc'
      include 'names.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'stdio.inc'
      include 'miscel.inc'
      include 'rndoff.inc'
c
      character line*(LINELG),hdr*7,shdr*4,fileID*(WORDLG),
     #	flname*(WORDLG)
      integer i,ioerr,jx,lth,lu,npotn
      double precision totwgt
c
      character*2 nonbr(3)
      data nonbr/ 'l', 'xy', 'zz'/
c
c parms outputs varaiables if they are spec/site = 1/1 or if they vary
c from those valuse for other spectra and sites.  This test may still
c output copies that are equal but it should limit the output.
c
      integer idxx,isp,istx,ispec,isite
      logical test
      test(idxx,isp,istx)=((isp.eq.1 .and. istx.eq.1) .or. ((isp.eq.1)
     #  .and. (abs(fparm(idxx,isp,istx)-fparm(idxx,1,1)) .gt. 
     #  0.00001*fparm(idxx,1,1))) .or. ((isp.ne.1) .and. 
     #  (abs(fparm(idxx,isp,istx)-fparm(idxx,1,istx)) .gt.
     #  0.00001*fparm(idxx,1,istx))))
c
c......................................................................
c
      call gettkn(line,fileID,lth)
c
c----------------------------------------------------------------------
c     No name specified: output to terminal 
c----------------------------------------------------------------------
      if (lth.eq.0) then
         lu=luttyo
         hdr=' '
         shdr=' '
      else
c
c----------------------------------------------------------------------
c     Set for output of "let" commands to specified file
c----------------------------------------------------------------------
         flname=fileID(:lth)//'.run'
         open(ludisk,file=flname,status='unknown',
     #        access='sequential',form='formatted',iostat=ioerr)
         if (ioerr.ne.0) then
            write (luout,3000) flname
            return
         end if
         lu=ludisk
c the idea of writing a file that can be read does not work
c  for multiple spectra/sites.
         hdr='let '
         shdr='c '
         write (lu,1022) fileID(:lth)
         write (lu,1024) ixname
      end if
c
c Loop over spectra and sites, outputting only when information 
c differes from spectra 1
c
      write(lulog,4007)
      write(lulog,4008)
      do 10 ispec=1,nspectra
	do 10 isite=1,ncomps
          write(lu,101)ispec,isite
 101      format(/,' For Spectrum ',i3,' Site ',i3,/)
          write(lulog,101)ispec,isite
	  if(test(IB0,ispec,isite)) 
     #       write (lu,1020) hdr,fparm(IB0,ispec,isite)
c----------------------------------------------------------------------
c     g-tensor
c----------------------------------------------------------------------
          if(test(IGXX,ispec,isite)) then
            if (iparm(IIGFLG,ispec,isite).eq.3) then
              write (lu,2001) hdr,(fparm(IGXX+i,ispec,isite),i=0,2)
            else if (iparm(IIGFLG,ispec,isite).eq.2) then
              write (lu,1001) hdr,(fparm(IGXX+i,ispec,isite),i=0,2)
            else
              write (lu,1000) hdr,(fparm(IGXX+i,ispec,isite),i=0,2)
            end if
          end if
c
c----------------------------------------------------------------------
c     Nuclear spin and hyperfine tensor
c----------------------------------------------------------------------
          if(test(IAXX,ispec,isite)) then
            write (lu,1002) hdr,iparm(IIN2,ispec,isite)
          end if
c
          if (iparm(IIN2,ispec,isite).ne.0) then
            if(test(IAXX,ispec,isite)) then
              if (iparm(IIGFLG+1,ispec,isite).eq.3) then
                write(lu,2004) hdr,(fparm(IAXX+i,ispec,isite),i=0,2)
              else if (iparm(IIGFLG+1,ispec,isite).eq.2) then
                write(lu,1004) hdr,(fparm(IAXX+i,ispec,isite),i=0,2)
              else
                write(lu,1003) hdr,(fparm(IAXX+i,ispec,isite),i=0,2)
              end if
            end if
          end if
c
c----------------------------------------------------------------------
c     Diffusion tensor
c----------------------------------------------------------------------
          if(test(IDX,ispec,isite)) then 
            if (iparm(IIGFLG+2,ispec,isite).eq.3) then
              write(lu,2006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
              write(lulog,2006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
              write(lulog,4001) hdr,(10.d0**fparm(IDX+i,ispec,isite)
     #        ,i=0,2)
            else if (iparm(IIGFLG+2,ispec,isite).eq.2) then
              write(lu,1006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
              write(lulog,1006) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
              write(lulog,4003) hdr,(10.d0**fparm(IDX+i,ispec,isite)
     #        ,i=0,2)
            else
              write(lu,1005) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
              write(lulog,1005) hdr,(fparm(IDX+i,ispec,isite),i=0,2)
              write(lulog,4002) hdr,(10.d0**fparm(IDX+i,ispec,isite)
     #        ,i=0,2)
            end if
          end if
c
c----------------------------------------------------------------------
c     Non-Brownian rotational model parameters
c----------------------------------------------------------------------
          if(test(IDJF,ispec,isite) .or. test(IPML,ispec,isite)) then
            write (lu,1007) hdr,iparm(IIPDF,ispec,isite)
            if (iparm(IIPDF,ispec,isite).eq.1) then
              do 2 i=0,2
                write(lu,1008) hdr,nonbr(i+1),nonbr(i+1),
     #		  iparm(IML+i,ispec,isite),fparm(IPML+i,ispec,isite)
 2            continue
c
c   -- Anisotropic viscosity
c
            else if (iparm(IIPDF,ispec,isite).eq.2) then
              write (lu,1009) hdr,fparm(IDJF,ispec,isite),
     #		fparm(IDJFPRP,ispec,isite)
            end if
c
c   --- Discrete jump motion
            if (iparm(IIPDF,ispec,isite).ne.2 .and. 
     #	      iparm(IIST,ispec,isite).ne.0 .and. 
     #        npotn.gt.0) write (lu,1010) hdr,iparm(IIST,ispec,isite),
     #		fparm(IDJF,ispec,isite)
          end if
c
c----------------------------------------------------------------------
c   Orienting potential
c----------------------------------------------------------------------
          if(test(IC20,ispec,isite) .or. test(IC20+1,ispec,isite)
     #       .or. test(IPSI,ispec,isite)) then
            npotn=0
            do 3 i=0,4
              if (dabs(fparm(IC20+i,ispec,isite)).gt.rndoff) npotn=i+1
 3          continue
            if (npotn.gt.0) then
              write (lu,1011) hdr,iparm(INORT,ispec,isite),
     #		fparm(IPSI,ispec,isite)
              write (lu,1012) hdr,(fparm(IC20+i,ispec,isite),i=0,4)
              write (lulog,1011) hdr,iparm(INORT,ispec,isite),
     #          fparm(IPSI,ispec,isite)
              write (lulog,1012) hdr,(fparm(IC20+i,ispec,isite),i=0,4)
            end if
          end if
c
c   --- Heisenberg exchange
c
          if(test(IOSS,ispec,isite)) then
            if (dabs(fparm(IOSS,ispec,isite)).gt.rndoff) then 
              write (lu,1013) hdr,fparm(IOSS,ispec,isite)
              write (lulog,1013) hdr,fparm(IOSS,ispec,isite)
              write (lulog,4004) hdr,10.d0**fparm(IOSS,ispec,isite)
            end if
          end if
c
c----------------------------------------------------------------------
c     Diffusion tilt & molecular tilt
c----------------------------------------------------------------------
          if(test(IALD,ispec,isite).or.test(IALM,ispec,isite)) then
            if ((dabs(fparm(IALD,ispec,isite)).gt.rndoff).or.
     #        (dabs(fparm(IALD+1,ispec,isite)).gt.rndoff).or.
     #        (dabs(fparm(IALD+2,ispec,isite)).gt.rndoff))
     #        write (lu,1014) hdr,(fparm(IALD+i,ispec,isite),i=0,2)
C
            if ((dabs(fparm(IALM,ispec,isite)).gt.rndoff).or.
     #        (dabs(fparm(IALM+1,ispec,isite)).gt.rndoff).or.
     #        (dabs(fparm(IALM+2,ispec,isite)).gt.rndoff))
     #        write (lu,1015) hdr,(fparm(IALM+i,ispec,isite),i=0,2)
          end if
c
c----------------------------------------------------------------------
c  Relaxation constants
c----------------------------------------------------------------------
          if(test(IT2EDI,ispec,isite).or.test(IT2EDI+1,ispec,isite).or.
     #        test(IT2EDI+2,ispec,isite).or.test(IT2EDI+3,ispec,isite)
     #        .or. test(IT2EDI+4,ispec,isite)) then
            write (lu,1016) hdr,(fparm(IT2EDI+i,ispec,isite),i=0,4)
            write (lulog,1016) hdr,(fparm(IT2EDI+i,ispec,isite),i=0,4)
            write (lulog,4005) hdr,(10.d0**fparm(IT2EDI+i,ispec,
     #      isite),i=0,3,3)
            write (lulog,4006) hdr,(1.d0/(10.d0**fparm(IT2EDI+i,
     #      ispec,isite))*1.d9,i=0,3,3)
          end if
c
c----------------------------------------------------------------------
c  Siteweight
c----------------------------------------------------------------------
          if(test(ISWGT,ispec,isite))then
            totwgt=0.0d0
            do 31 i=1,ncomps
             totwgt=totwgt+fparm(ISWGT,ispec,i)
 31         continue
cc
            write(lu,1116) hdr,fparm(ISWGT,ispec,isite)/totwgt
            write(lulog,1116) hdr,fparm(ISWGT,ispec,isite)/totwgt
          end if
c
c----------------------------------------------------------------------
c  Broadening width
c----------------------------------------------------------------------
          if(test(IMWID,ispec,isite).or.test(IGIB,ispec,isite) .or.
     #        test(ILIB,ispec,isite).or.test(IHWID,ispec,isite)) then
            write (lu,1017) hdr,fparm(IMWID,ispec,isite),
     #		fparm(IGIB,ispec,isite),fparm(ILIB,ispec,isite),
     #                fparm(IHWID,ispec,isite)
            write (lulog,1017) hdr,fparm(IMWID,ispec,isite),
     #          fparm(IGIB,ispec,isite),fparm(ILIB,ispec,isite),
     #                fparm(IHWID,ispec,isite)
          end if
c
c----------------------------------------------------------------------
c  Basis set and pruning tolerance
c----------------------------------------------------------------------
c
          if(test(ICGTOL,ispec,isite).and.nbasis.ne.0) then
            write (lu,1018) hdr,(iparm(ILEMX+i,ispec,isite),i=0,4),
     #		hdr,ndimoss(basinfo(1,ispec,isite)),
     #		ndimdss(basinfo(1,ispec,isite))
            write (lu,1019) hdr,fparm(ICGTOL,ispec,isite),
     #		iparm(INSTEP,ispec,isite)
          end if
c
c----------------------------------------------------------------------
c  Series parameters
c----------------------------------------------------------------------
c
          if(ispec.eq.1.and.isite.eq.1) then  ! parms apply to all sites
            if (lu.eq.ludisk) hdr='series '
            write (lu,1030) shdr,nser
            write (lu,1031) hdr,(siexp(i),i=1,nser)
            write (lu,1032) hdr,(sicomb(i),i=1,nser)
            write (lu,1033) hdr,(snpt1(i),i=1,nser)
            write (lu,1034) hdr,(snpt2(i),i=1,nser)
            write (lu,1035) hdr,(sinit1(i),i=1,nser)
            write (lu,1036) hdr,(sstept1(i),i=1,nser)
            write (lu,1037) hdr,(sinit2(i),i=1,nser)
            write (lu,1038) hdr,(sstept2(i),i=1,nser)
            write (lu,1039) hdr,(stfix(i),i=1,nser)
            write (lu,1040) hdr,(sfac(i),i=1,nser)
            write (lu,1140) hdr,(sratio(i),i=1,nser)
            if (lu.eq.luttyo) write (lu,1041) hdr,(idepnd(i),i=1,nser)
          end if
10    continue
      write(lulog,4009)
      write(lulog,4007)
c
      return
c
c######################################################################
c   
 1000 format(a,'gxx,gyy,gzz = ',2(f9.6,','),f9.6)
 1001 format(a,'g1,g2,g3 = ',2(f9.6,','),f9.6)
 2001 format(a,'gprp,grhm,gpll = ',2(f9.6,','),f9.6)
 1002 format(a,'in2 =',i2)
 1003 format(a,'Axx,Ayy,Azz = ',2(f9.4,','),f9.4)
 1004 format(a,'A1,A2,A3 = ',2(f9.4,','),f9.4)
 2004 format(a,'Aprp,Arhm,Apll = ',2(f9.4,','),f9.4)
 1005 format(a,'log(Rx,Ry,Rz) = ',2(f9.4,','),f9.4)
 1006 format(a,'log(Rbar,N,Nxy) = ',2(f9.4,','),f9.4)
 2006 format(a,'log(Rprp,Rrhm,Rpll) = ',2(f9.4,','),f9.4)
 1007 format(a,'ipdf = ',i1)
 1008 format(a,'m',a2,', pm',a2,' = ',i2,',',g10.3)
 1009 format(a,'djf,djfprp = ',g10.3,',',g10.3)
 1010 format(a,'ist,djf =',i3,',',g10.3)
 1011 format(a,'nort,psi = ',i3,',',f8.3)
 1012 format(a,'c20,c22,c40,c42,c44 = ',4(f7.4,','),f7.4)
 1013 format(a,'oss = ',g11.4)
 1014 format(a,'ald,bed,gad =',f7.2,',',f7.2,',',f7.2)
 1015 format(a,'alm,bem,gam =',f7.2,',',f7.2,',',f7.2)
 1016 format(a,'t2edi,t2ndi,t2efi,t1edi,t1ndi = ',4(f6.3,','),f6.3)
 1116 format(a,'Siteweight = ',f10.6)
 1017 format(a,'mwid,gib,lib,hwid = ',3(f8.4,','),f8.4)
 1018 format(a,'lemx,lomx,kmx,mmx,ipnmx = ',4(i3,','),i3,
     #       /a,'ndimoss,ndimdss = ',i4,i5)
 1019 format(a,'cgtol,nstep = ',g10.3,i5)
 1020 format(a,'B0 =',f10.3)
 1022 format('log    ',a)
 1024 format('basis  ',a)
 1030 format(a,'*** Parameters for ',i2,' spectra in series ***')
 1031 format(a,'iexp = ',9(i3,','),i3)
 1032 format(a,'icomb = ',9(i3,','),i3)
 1033 format(a,'npt1 = ',9(i4,','),i4)
 1034 format(a,'npt2 = ',9(i4,','),i4)
 1035 format(a,'init1  = ',9(f6.1,','),f6.1)
 1036 format(a,'stept1 = ',9(f6.1,','),f6.1)
 1037 format(a,'init2  = ',9(f6.1,','),f6.1)
 1038 format(a,'stept2 = ',9(f6.1,','),f6.1)
 1039 format(a,'tfix = ',9(f7.1,','),f7.1)
 1040 format(a,'spect weight = ',9(f10.6,','),f10.6)
 1140 format(a,'data scale factor = ',9(f10.6,','),f10.6)
 1041 format(a,'idepnd = ',9(i3,','),i3)
 3000 format('*** Unable to open file ',a,' for parameter output ***')
 4001 format(a,'Rprp,Rrhm,Rpll [1/sec]= ',2(e9.4,','),e9.4)
 4002 format(a,'Rx,Ry,Rz [1/sec] = ',2(e9.4,','),e9.4)
 4003 format(a,'Rbar[1/sec],N,Nxy = ',2(e9.4,','),e9.4)
 4004 format(a,'oss [1/sec]= ',e9.4)
 4005 format(a,'t2edi,t1edi [1/sec]= ',1(e9.4,','),e9.4)
 4006 format(a,'T2,T1 [nsec]= ',1(f9.2,','),f9.2)
 4007 format('=======================================================')
 4008 format('************** PARAMETER conversion *************')
 4009 format('************** PARAMETER conversion END *************')
c
      end


c----------------------------------------------------------------------
c                    =========================
c                       subroutine STATC
c                    =========================
c
c  Prints out the information on variables on the given logical unit
c  number.
c
c----------------------------------------------------------------------
      subroutine statc( line,lu )
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'lmcomm.inc'
      include 'miscel.inc'
c      include 'iterat.inc'
c
      character line*(LINELG)
      integer lu,i,isite,ispec
c
c----------------------------------------------------------------------
c     parameters for fit command
c----------------------------------------------------------------------
c
c only output spec 1 site 1:
      ispec=1
      isite=1
      write(lu,1000)
      write(lu,1001) xtol,ftol,gtol
      write(lu,1002) maxitr
      write(lu,1003) maxev
      write(lu,1004) factor
c
c----------------------------------------------------------------------
c     parameters for search command
c----------------------------------------------------------------------
      write(lu,1006)
      write(lu,1007) stol,sstep,smn,smx
c
c----------------------------------------------------------------------
c     variable informations
c----------------------------------------------------------------------
      if (nprm.ge.1) then
         write(lu,1010) nprm
         write(lu,1011)
         do 10 i=1,nprm
            write(lu,1012) tag(i),fparm(ixpr(i),ispec,isite),
     #		prscl(i),xfdstp(i)
 10      continue
      else
         write(lu,1020)
      end if
c
      write(lu,1030) lemx,lomx,kmx,mmx,ipnmx
      write(lu,1032) ndimoss(basinfo(1,ispec,isite)),
     #	ndimdss(basinfo(1,ispec,isite))
      write(lu,1040) fnorm/dsqrt(dfloat(nptot))
c
      return
c
c######################################################################
c   
 1000 format(/,2x,'FIT parameters : ')
 1001 format(5x,'xtol = ',1p,e8.1,3x,'ftol = ',1p,e8.1,3x,'gtol = ',
     #       1p,e8.1)
 1002 format(5x,'maximum iterations = ',i2)
 1003 format(5x,'maximum function evaluations = ',i3)
 1004 format(5x,'initial step bound in parameter search = ',f8.1)
c
 1006 format(/,2x,'SEARCH parameters : ')
 1007 format(5x,'tol = ',1p,e8.1,3x,'step = ',1p,e8.1,3x,'min = ',
     #       1p,e8.1,3x,'max = ',1p,e8.1)
c
 1010 format(/,2x,'VARY parameters : ',i1,' variables')
 1011 format(2x,56('-'),/,6x,'VARIABLE',8x,'VALUE',10x,'SCALE',5x,
     #       'JC STEP',/,2x,56('-'))
 1012 format(7x,a6,5x,g14.7,3x,f6.1,5x,1p,e8.1)
c
 1020 format(/,2x,'No variables are being varied.')
 1030 format(/,2x,'MTS lemx,lomx,kmx,mmx,ipnmx = ',i3,4(',',i2))
 1032 format(2x,'Dimension of matrix: ndimoss = ',i4,' ndimdss = ',i4)
 1040 format(/,2x,'Rms Deviation = ',g14.7,/)
c
      end
