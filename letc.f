c  Version 1.5.1 beta 2/3/96
c----------------------------------------------------------------------
c                    =========================
c                      subroutine LETC
c                    =========================
c
c let <name>{(index)} {, <name>,...} = <value> {, <value> ... }
c
c      name    : name of the parameter to be assigned a value
c      index   : site or spectral index of parameter
c      value   : value to be assigned to the given parameter
c      
c   Up to 10 variables may be assigned values in a single let command
c
c----------------------------------------------------------------------
      subroutine letc( line )
      implicit none
      character*80 line
c
      include 'nlsdim.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'prmeqv.inc'
      include 'parcom.inc'
      include 'lpnam.inc'
      include 'stdio.inc'
c
      integer MXLFT
      parameter (MXLFT=10)
c
      integer i,ileft,irt,ival,ixf,ixr,ixr2,ixs,ixt,ixxa,
     #        lth,llth,left(MXLFT),leftx(MXLFT)
      double precision fval
      character token*30,prmID*9
      logical rhs,logflg
c
      integer ipfind,isfind,indtkn,itrim,ixlim
      double precision getprm
      logical ftoken,itoken,tcheck,spcpar
      character*7 ptype
      external ftoken,itoken,tcheck,ipfind,isfind,itrim,indtkn,
     #         spcpar,ptype,ixlim,getprm
c
c######################################################################
c
      rhs=.false.
      logflg=.false.
      ileft=1
c
c----------------------------------------------------------------------
c Get the next parameter name (or value for right-hand-side of statement)
c----------------------------------------------------------------------
 1    call gettkn(line,token,lth)
      lth=min(lth,6)
      call touppr(token,lth)
c
c----------------------------------------------------------------------
c  No more tokens: check for error conditions
c----------------------------------------------------------------------
      if (lth.eq.0) then
         if (rhs) then
c                                *** not enough values
            write (luttyo,1000)
         else
            if (ileft.eq.1) then
c                                *** variable name expected
               write (luttyo,1001)
            else
c                                *** no '=' given
               write(luttyo,1002)
            end if
         end if
         return
      end if
c
      if (token.eq.')' .or. (logflg.and.token.eq.'(')) goto 1
c
      if (token.eq.'LOG'.and.rhs) then
         logflg=.true.
         go to 1
      end if
c
c----------------------------------------------------------------------
c  Check for '=' sign in assignment
c----------------------------------------------------------------------
      if (token.eq.'=') then
         rhs=.true.
         irt=1
         go to 1
      end if
c
c----------------------------------------------------------------------
c  Right-hand side: assign value to parameter specified in left array
c----------------------------------------------------------------------
      if (rhs) then
         ixr=left(irt)
         ixr2=leftx(irt)
         irt=irt+1
c
c        ---------------------------------------------------
c        Make assignment for legal parameter names/indices
c        ---------------------------------------------------
         if (ixr.ne.0) then
            if (ixr2.lt.1) ixr2=0
c
c           ---------------------------------------------
c           Set the name of the parameter being assigned
c           ---------------------------------------------
            if (ixr.lt.-100 .and. ixr.gt.-200) then
               prmID=alias2(-99-(IWXX+ixr))
            else if (ixr.lt.0) then
               prmID=alias1(1-(IWXX+ixr))
             else if (ixr.gt.100) then
                prmID=parnam(ixr-100)
             else if (ixr.gt.0 .and. ixr.lt.100) then
               prmID=parnam(ixr)
            end if
c
c           -----------------------------
c           Get the value being assigned
c           -----------------------------
            if (ixr.lt.100) then
c
c              -------------------------------
c              Look for a floating-point token
c              -------------------------------
               if (.not.ftoken(token,lth,fval)) then
                  write(luttyo,1004) token(:lth)
               else
                  if(logflg) then
                     if (fval.gt.0.0d0) then
                        fval=dlog10(fval)
                     else
c                                              *** illegal log arg
                        write(luttyo,1011) fval
                     end if
                     logflg=.false.
                  end if
               end if
c
            else
c
c              --------------------------------------------------------
c              Look for a symbolic value: pre-defined symbol or ID name
c              --------------------------------------------------------
               ixs=isfind(token,lth)
               if (ixs.ne.0) then
                  if (ixs.gt.0) ival=symval(ixs)
                  if (ixs.lt.0) ival=abs(ixs)
c
c              -------------------------
c              Look for an integer token
c              -------------------------
               else if (.not.itoken(token,lth,ival)) then
                  write(luttyo,1005) token(:lth)
               end if
c
c                      if (ixr.lt.100)...else...
            end if
c
c           -------------------------------
c           Assign the floating point value
c           -------------------------------   
            if (ixr.lt.100) then
               if (tcheck(ixr,ixr2,prmID,luout)) then
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
c
c
                  llth=itrim(prmID)
                  ixxa=ixx(iabs(mod(ixr,100)),ixr2+1)
                  if (ixxa.gt.0 .and. (ixr2.ge.1 .and. ixst(ixxa).lt.1)) 
     #               write (luttyo,1012) parnam(ixr)(:llth)
c
                  if (ixr2.lt.1) then
                     do i=2,ixlim(ixr)
                        if (getprm(ixr,i).ne.getprm(ixr,1)) then
                           write(luttyo,1013) parnam(ixr)(:llth),
     #                                        ptype(ixr)
                           go to 2
                        end if
                     end do
                  end if
                        
c
 2                call setprm(ixr,ixr2,fval)
c
               end if
c
            else
c              ------------------------
c              Assign the integer value
c              ------------------------
               call setipr(ixr,ixr2,ival)
c
            end if
c
c                   if (ixr.ne.0) ....
         end if
c
c        --------------------------------------------------
c        Return if all assignments have been made
c        --------------------------------------------------
         if (irt.ge.ileft) return
c
c----------------------------------------------------------------------
c  Left-hand side: build a list of indices of parameters to be assigned
c----------------------------------------------------------------------
c
      else if (ileft.le.mxlft) then
         left(ileft)=ipfind(token,lth)
         if (left(ileft).eq.0) write (luttyo,1003) token(:lth)
         leftx(ileft)=indtkn(line)
         ileft=ileft+1
      end if
      go to 1
c
c ###### format statements ########################################
c
 1000 format('*** Not enough values specified ***')
 1001 format('*** Variable name expected ***')
 1002 format('*** No ''='' specified ***')
 1003 format('*** ''',a,''' is not a parameter ***')
 1004 format('*** Real value expected: ''',a,''' ***')
 1005 format('*** Integer value expected: ''',a,''' ***')
 1011 format('*** Illegal log argument:',g12.5,' ***')
 1012 format('*** Warning: ',a,' is being varied globally ***')
 1013 format('*** Warning: ',a,' is now the same for all ',a,' ***')
      end
