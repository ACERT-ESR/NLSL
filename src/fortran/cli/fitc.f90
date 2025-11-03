! Version 1.3 7/3/93
!----------------------------------------------------------------------
!                    =========================
!                       subroutine FITC
!                    =========================
!
!   fit      { trace xtol <xtol> ftol <ftol> gtol <ftol> 
!              maxfun <mxf> maxitr <mxi> bound <factor> }
!
!          trace:  Specifies that a ".trc" file should be produced for
!                  the fit
!          xtol:   Convergence tolerance for scaled fitting parameters
!          ftol:   Convergence tolerance for chi-squared
!          gtol:   Convergence tolerance for gradient of chi-squared with
!                  respect to the fitting parameters
!          mxf:    Maximum number of function calls allowed
!          mxi:    Maximum number of iterations allowed
!          factor: Factor defining initial step bound used in parameter search
!
!----------------------------------------------------------------------
      subroutine fitc( line )
!
      use nlsdim
      use lmcom
      use parcom
      use iterat
      use stdio
!
      implicit none
      character line*80
!
      logical ftoken
      external ftoken
!
      integer i,lth
      double precision fval
      character*30 token
!
      integer NKEYWD
      parameter(NKEYWD=19)
!
      integer itrim
      external itrim
!
      character*8 keywrd(NKEYWD)
      data keywrd / 'FTOL',   'GTOL', 'XTOL',   'BOUND',  'MAXFUN',
     #              'MAXITR', 'SHIFT','SRANGE', 'TRACE',  'JACOBI', 
     #              'NOSHIFT','NEG',  'NONEG',  'TRIDIAG','ITERATES',
     #              'WRITE', 'WEIGHTED', 'UNWEIGHT', 'CTOL' /
!
!######################################################################
!
!  -- Reset "non-sticky" flags for fitting procedure
!
      output=0
!
!----------------------------------------------------------------------
!  Look for a keyword
!----------------------------------------------------------------------
!
   14 call gettkn(line,token,lth)
      lth=min(lth,8)
!
!---------------------------------------------------------------
!                    ****************************************
! No more keywords:  **** call the NLS fitting routine ******
!                    ****************************************
!---------------------------------------------------------------
      if (lth.eq.0) then
        call fitl
        return
      end if
!
!------------------------------
! Check keyword list
!------------------------------
      call touppr(token,lth)
      do 15 i=1,NKEYWD
        if (token(:lth).eq.keywrd(i)(:lth)) goto 16
   15   continue
!                                        *** Unrecognized keyword
      write (luttyo,1000) token(:lth)
      return
!
!----------------------------------------------------------------------
!  Keyword found: for keywords requiring an argument, convert 
!  next token and assign appropriate value
!----------------------------------------------------------------------
   16   if ((i.ge.1 .and.i.le.8) .or. i.eq.19) then
          call gettkn(line,token,lth)
!                                        *** No value given
          if (lth.eq.0) then
!                                        *** default 0 for SHIFT keyword
             if (i.eq.7) then
                nshift=0
!                                        *** otherwise, error 
             else
                write(luttyo,1003) keywrd(i)(:itrim(keywrd(i)))
                return
             end if
          end if
!
          if (ftoken(token,lth,fval)) then
!                                          *** FTOL keyword
             if (i.eq.1) then
                ftol=fval
!                                          *** GTOL keyword
             else if (i.eq.2) then
                gtol=fval
!                                          *** XTOL keyword
             else if (i.eq.3) then
                xtol=fval
!                                          *** BOUND keyword
             else if (i.eq.4) then
                factor=fval
!                                          *** MAXFUN keyword
             else if (i.eq.5) then
                maxev=int(fval)
!                                          *** MAXITR keyword
             else if (i.eq.6) then
                maxitr=int(fval)
!                                          *** SHIFT keyword
             else if (i.eq.7) then
                nshift=int(fval)
!                                          *** SRANGE keyword
             else if (i.eq.8) then
                srange=fval/1.0D2
                if (srange.gt.1.0d0) srange=1.0d0
                if (srange.lt.0.0d0) srange=0.0d0
             end if
!                                      *** Illegal numeric value
          else
             if (i.eq.7) then
                nshift=0
                call ungett(token,lth,line)
             else
                write(luttyo,1001) token(:lth)
             end if
          end if
!                                          *** TRACE keyword
       else if (i.eq.9) then
          if (luout.eq.luttyo) then
             write (luttyo,1050)
             itrace=0
          else
             itrace=1
          end if
!                                          *** JACOBI keyword
       else if (i.eq.10) then
          jacobi=1
!                                          *** NOSHIFT keyword
       else if (i.eq.11) then
          nshift=-1
!                                          *** NEG keyword
       else if (i.eq.12) then
          noneg=0
!                                          *** NONEG keyword
       else if (i.eq.13) then
          noneg=1
!                                          *** TRIDIAG keyword
       else if (i.eq.14) then
          itridg=1
!                                          *** ITERATES keyword
       else if (i.eq.15) then
          iitrfl=1
!                                          *** WRITE keyword
       else if (i.eq.16) then
          output=1
!                                          *** WEIGHTED keyword
       else if (i.eq.17) then
          iwflag=1
!                                          *** UNWEIGHT keyword
       else if (i.eq.18) then
          iwflag=0
!                                          *** CTOL keyword
       else if (i.eq.19) then
          ctol=fval
       end if
!
       go to 14
!
!######################################################################
!
 1000 format('*** Unrecognized FIT keyword: ''',a,''' ***')
 1001 format('*** Numeric value expected: ''',a,''' ***')
 1003 format('*** No value given for ''',a,''' ***')
 1050 format('*** A log file must be opened before using TRACE ***')
      end
