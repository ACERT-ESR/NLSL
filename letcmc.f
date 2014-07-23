c  VERSION 1.0  NLSPMC   2/5/99
c----------------------------------------------------------------------
c                    =========================
c                      subroutine LETCMC
c                    =========================
c Interprete the command:
c let <name>{(index)} {, <name>,...} = <value> {, <value> ... }
c
c let only allows setting of variables that apply to a specific site
c and all spectra.  To input variables applying to one or more spectra,
c in a multiple series, use the series command.  In the case of only 
c one of a series, the let may be used or the series command.
c
c      name    : name of the parameter to be assigned a value
c      index   : site or spectral index of parameter
c      value   : value to be assigned to the given parameter
c      
c   Up to 10 variables may be assigned values in a single let command
c
c Parse the input line into arrays (leftid, leftspect and leftsite) 
c containing the name and site information in the form: 0=all 
c sites/spectra, else=specific site/spectrum.
c
c Then parse into the final arrays fparm(variable,spectrum,site) and
c iparm(variable,spectrum,site).  Based on the CW version with mods.
c
c Later for each spectrum,site, we copy the current values into xxx 
c for fitting, or into fepr and iepr in sim2d for simulation.
c
c Errors in various checks for consistency result in return without
c further processing of line.
c
c----------------------------------------------------------------------
      subroutine letcmc(line)
      implicit none
c
c      include 'limits.inc'
c	integer mxdim
c	parameter mxdim=1200

      include 'limits.inc'
      include 'simparm.inc'
      include 'datas.inc'
      include 'parms.inc'
      include 'lpnam.inc'
      include 'miscel.inc'
      include 'stdio.inc'
c
c  maximum number of variables in command line:
      character*(LINELG) line
      integer mxlft,specptr,siteptr,varyinfo,iparmptr,isite,ispec
      parameter (mxlft=10)
      logical vchange
      external vchange
c
      integer ixs,ival,instp,vspc
c      integer ileft,irt,ival,ixf,ixr,ixs,ixri,ixt,lth,left(mxlft)
c  ileft,irt = counters for parameter names and values
c  leftid,leftspec,leftsite = hold identification # for parameter,
c   spectrum and site info (0=all spectra or all sites).  
c
      integer leftid(mxlft),leftspec(mxlft),leftsite(mxlft)
c      integer ileft,irt,lth,left(mxlft),leftx(mxlft),llth
      integer ileft,irt,lth,llth,ixr,ixrspec,ixrsite
      double precision fval
c
      character token*(WORDLG),prmID*6
c      logical rhs,logflg,serprm
      logical rhs,logflg
c
      integer isfind,varspec(mxlft)
      logical ftoken,itoken,tcheck
      external ftoken,itoken,tcheck,isfind
c
c######################################################################
c
      rhs=.false.
      logflg=.false.
c change starting point for irt and ileft counters.
      ileft=0
      irt=0
c
c----------------------------------------------------------------------
c Get the next parameter name (or value for right-hand-side of statement)
c line = line of text to be parsed, token = returned token, lth = length
c----------------------------------------------------------------------
 1    call gettkn(line,token,lth)
      lth=min(lth,6)
      call touppr(token,lth)
c
c----------------------------------------------------------------------
c  Check for error conditions
c----------------------------------------------------------------------
      if (lth.eq.0) then		! expect a value or parameter
	write (luttyo,1000)
 1000   format('*** not enough values ***')
	return
      end if
c
      if (token.eq.')' .or. token.eq.'(') then
c
c If (, not expected at this point.
c

        write(luttyo,1014)
 1014   format('*** Unexpected ( or ) ***')
	return
      end if
c
      if (token.eq.'LOG') then 
	 if(rhs) then
           logflg=.true.
           go to 1
         else
           write(luttyo,1015)
 1015      format('*** Log before =? ***')
         end if
         return
      end if
c
c----------------------------------------------------------------------
c  Check for '=' sign in assignment
c----------------------------------------------------------------------
      if (token.eq.'=') then
         rhs=.true.
         go to 1
      end if
c
c----------------------------------------------------------------------
c  Left-hand side: assign parameter index and site information
c----------------------------------------------------------------------
c
      if ( .not. rhs) then	! if lefthand side:
c
c get parameter id # and save in leftid.
c
	 if(ileft.gt.mxlft-1) then
	   write (luttyo,1016) token(:lth)
 1016      format('*** Too many parmeters in let''',a, ''' ***')
	   return
	 end if
c get the parameter id # including axis information:
c         leftid(ileft+1)=ipfind(token,lth)	! # coded in eprprm.inc
         call ipfind(leftid(ileft+1),token,lth,varspec(ileft+1))
         if (leftid(ileft+1).eq.0) then
	   write (luttyo,1003) token(:lth)	! name not found
 1003      format('*** ''',a,''' is not a parameter ***')
	   go to 1
   	 end if
c
	 call indtkn(line,leftspec(ileft+1),leftsite(ileft+1))
	 if(leftspec(ileft+1).eq.-1.or.leftsite(ileft+1).eq.-1) then
	   write(luttyo,*)'error in indtkn '
	   return
	 end if
c
         ileft=ileft+1

         go to 1			! go get another
      end if

c
c----------------------------------------------------------------------
c  Right-hand side: assign value to parameter specified in leftid array
c----------------------------------------------------------------------
c
      if (rhs) then
        irt=irt+1
	if(irt.gt.ileft) Then
	  write(luttyo,1017)
 1017	  format('*** too many right hand side values ***')
	  return
	end if
        ixr=leftid(irt)		! pointer to parameter name
	ixrspec=leftspec(irt)		! spectra info, 0=all spectra
	ixrsite=leftsite(irt)		! site info, 0=all site
c
c   --- EPR parameter assignment.  See ipfind.f for meaning of values.
c  Note, axial and spherical tensor elements are pointed to by IGXX,
c  while fp and integer are pointed to by ixr directly.
c
        if (ixr.lt.-100 .and. ixr.gt.-200) then		! axial tensor
           prmID=alias2(-99-(IGXX+ixr))
        else if (ixr.lt.0) then				! spherical
           prmID=alias1(1-(IGXX+ixr))
        else if (ixr.gt.100) then			! integer
           prmID=parnam(ixr-100)
        else if (ixr.gt.0 .and. ixr.lt.100) then	! floating point
           prmID=parnam(ixr)
	else 
	   write(luttyo,1018)ixr
 1018	   format('*** Error, parm unknown, ixr ***''',i8,''' ***')
           return
        end if
c
        if (ixr.lt.100) then		! floating point.
c get fval if fp:
          if (.not. ftoken(token,lth,fval)) then	! assigns fval
            write(luttyo,1004) token(:lth)
 1004       format('*** Real value expected: ''',a,''' ***')
	    return
          end if
c got fval, store it:
          if (logflg) then		! take log if necessary
            if (fval.gt.0.0d0) then
              fval=dlog10(fval)
            else
              write(luttyo,1011) fval
 1011         format(' *** Illegal log argument:',g12.5,' ***')
	      return
            end if
            logflg=.false.
          end if
c           -------------------------------
c           Assign the floating point value
c           -------------------------------
c  Check for consistent mode of tensors that depend of axis system.
c  Set mode if previously undefined.  Compare ixr with mode flag.
c  
          if (tcheck(ixr,prmID,luout)) then	! tensym.f, if valid axis
c
c                 ---------------------------------------------
c                 Issue a warning if:
c                    (1) the parameter is being set for a specific
c                        site/spectrum when it is being varied for
c                        *all* sites/spectra
c                    (2) the parameter is being set for all sites/spectra
c                        when its current value is different for the
c                        currently defined sites/spectra
c                 ---------------------------------------------
c  for now, we only allow variation of parameters that apply to all
c  spectra.  Just check the vary information, don't change it.
c
              specptr=1			! get the old information
	      siteptr=1
	      if (ixrspec.ne.0) specptr=ixrspec		! get first spec
	      if (ixrsite.ne.0) siteptr=ixrsite		! get first site
              iparmptr=iabs(mod(ixr,100))
	      varyinfo=ixx(iparmptr,specptr,siteptr)  ! vary this?
c	      mattrisp=varspec(iparmptr)
c consider all sites, spec:
	      do 30 isite=1,ncomps
		do 30 ispec=1,nspectra
		  if ((isite.eq.ixrsite .or. ixrsite.eq.0) .and. 
     #			(ispec.eq.ixrspec .or. ixrspec.eq.0)) then
		    if (ixx(iparmptr,specptr,siteptr) .ne. 
     #			varyinfo) then
                      write(luttyo,1021) prmID
 1021                 format(' *** Parm variation inconsistent:',
     #				a,' ***')
		    else
c if parameter is changed, update "have matrix" information: 
c Drop basinfo(2 information, replace with matrxok.
c                      if (basinfo(2,ispec,isite) .ge. mattrisp .and. 
c     #			(abs(fparm(iparmptr,ispec,isite)-fval)) .gt. 
c     #			1.d-8*abs(fval)) then
c                        basinfo(2,ispec,isite)= mattrisp-1
c                      end if
                      if (vchange(fval,fparm(iparmptr,ispec,isite)))
     #                  matrxok(isite)=
     #                    min(matrxok(isite),specid(iparmptr))
		      fparm(iparmptr,ispec,isite)=fval
                    end if
                  end if
 30           continue
c
	  else
	    write(luout,*)'Tcheck failure in letcmc'
	    return
	  end if		! end if tcheck
c
c  --- EPR parameter integer assignment 
c
        else if (ixr.gt.100) then		! integer parameter
c do we have a symbol?
          ixs=isfind(token,lth)		! zero if not symbol
          if (ixs.gt.0) then
            ival=symval(ixs)
          else if (.not.itoken(token,lth,ival)) then  ! integer->ival
c error result (itoken=.false.):
            write(luttyo,1005) token(:lth)
 1005       format('*** Integer value expected: ''',a,''' ***')
            return
          end if      
c
c have ival, do some checks here:
c 
c
c check # of CG steps against maximum
c
          call ipfind(instp,'NSTEP',5,vspc)  ! get pointer to 'NSTEP'
          if (ixr .eq. instp) then	! if the keyword is NSTEP:
            if (ival .gt. MXSTEP-2) then
              write(luttyo,*)'nstep too big, set to ',MXSTEP-2
              ival=MXSTEP-2
            end if
          end if
c
c
c
c
c set for specified sites, spectra:
c
          do 40 isite=1,ncomps
            do 40 ispec=1,nspectra
c if this spectra/site is to be changed:
              if ((isite.eq.ixrsite .or. ixrsite.eq.0) .and. 
     #	          (ispec.eq.ixrspec .or. ixrspec.eq.0)) then
                  iparm(ixr-100,ispec,isite)=ival
              end if
 40       continue
        end if
      end if	! end of if(rhs)
c
c    --- Return if all assignments have been made
c
      if(irt.eq.ileft) then
        return
      end if
      go to 1
c ###### format statements ########################################
c
 1001 format('*** variable name expected ***')
 1002 format('*** no ''='' specified ***')
 1012 format(' *** Series parameters in let command : series ',
     #       'assumed ***')
      end
