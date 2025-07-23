c Version 1.3 7/3/93
c----------------------------------------------------------------------
c                    =========================
c                       subroutine FITC
c                    =========================
c
c   fit      { trace xtol <xtol> ftol <ftol> gtol <ftol>
c              maxfun <mxf> maxitr <mxi> bound <bound> }
c
c          trace:  Specifies that a ".trc" file should be produced for
c                  the fit
c          xtol:   Convergence tolerance for scaled fitting parameters
c          ftol:   Convergence tolerance for chi-squared
c          gtol:   Convergence tolerance for gradient of chi-squared with
c                  respect to the fitting parameters
c          mxf:    Maximum number of function calls allowed
c          mxi:    Maximum number of iterations allowed
c          bound:  Factor defining initial step bound used in parameter search
c
c----------------------------------------------------------------------
      subroutine fitc( line )
c
      use nlsdim
      use lmcom
      use parcom
      use iterat
      use stdio
c
      implicit none
      character line*80
c
      logical ftoken
      external ftoken
c
      integer i,lth
      double precision fval
      character*30 token
c
      integer NKEYWD
      parameter(NKEYWD=19)
c
      integer itrim
      external itrim
c
      character*8 keywrd(NKEYWD)
      data keywrd / 'FTOL',   'GTOL', 'XTOL',   'BOUND',  'MAXFUN',
     #              'MAXITR', 'SHIFT','SRANGE', 'TRACE',  'JACOBI', 
     #              'NOSHIFT','NEG',  'NONEG',  'TRIDIAG','ITERATES',
     #              'WRITE', 'WEIGHTED', 'UNWEIGHT', 'CTOL' /
c
c######################################################################
c
c  -- Reset "non-sticky" flags for fitting procedure
c
      write_flag=0
c
c----------------------------------------------------------------------
c  Look for a keyword
c----------------------------------------------------------------------
c
   14 call gettkn(line,token,lth)
      lth=min(lth,8)
c
c---------------------------------------------------------------
c                    ****************************************
c No more keywords:  **** call the NLS fitting routine ******
c                    ****************************************
c---------------------------------------------------------------
      if (lth.eq.0) then
        call fitl
        return
      end if
c
c------------------------------
c Check keyword list
c------------------------------
      call touppr(token,lth)
      do 15 i=1,NKEYWD
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
   16   if ((i.ge.1 .and.i.le.8) .or. i.eq.19) then
          call gettkn(line,token,lth)
c                                        *** No value given
          if (lth.eq.0) then
c                                        *** default 0 for SHIFT keyword
             if (i.eq.7) then
                shift_flag=0
c                                        *** otherwise, error 
             else
                write(luttyo,1003) keywrd(i)(:itrim(keywrd(i)))
                return
             end if
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
                bound=fval
c                                          *** MAXFUN keyword
             else if (i.eq.5) then
                maxfun=int(fval)
c                                          *** MAXITR keyword
             else if (i.eq.6) then
                maxitr=int(fval)
c                                          *** SHIFT keyword
             else if (i.eq.7) then
                shift_flag=int(fval)
c                                          *** SRANGE keyword
             else if (i.eq.8) then
                srange=fval/1.0D2
                if (srange.gt.1.0d0) srange=1.0d0
                if (srange.lt.0.0d0) srange=0.0d0
             end if
c                                      *** Illegal numeric value
          else
             if (i.eq.7) then
                shift_flag=0
                call ungett(token,lth,line)
             else
                write(luttyo,1001) token(:lth)
             end if
          end if
c                                          *** TRACE keyword
       else if (i.eq.9) then
          if (luout.eq.luttyo) then
             write (luttyo,1050)
             trace=0
          else
             trace=1
          end if
c                                          *** JACOBI keyword
       else if (i.eq.10) then
          jacobi=1
c                                          *** NOSHIFT keyword
       else if (i.eq.11) then
          shift_flag=-1
c                                          *** NEG keyword
       else if (i.eq.12) then
          neg_flag=0
c                                          *** NONEG keyword
       else if (i.eq.13) then
          neg_flag=1
c                                          *** TRIDIAG keyword
       else if (i.eq.14) then
          tridiag_flag=1
c                                          *** ITERATES keyword
       else if (i.eq.15) then
          iterates_flag=1
c                                          *** WRITE keyword
       else if (i.eq.16) then
          write_flag=1
c                                          *** WEIGHTED keyword
       else if (i.eq.17) then
          weighted_flag=1
c                                          *** UNWEIGHT keyword
       else if (i.eq.18) then
          weighted_flag=0
c                                          *** CTOL keyword
       else if (i.eq.19) then
          ctol=fval
       end if
c
       go to 14
c
c######################################################################
c
 1000 format('*** Unrecognized FIT keyword: ''',a,''' ***')
 1001 format('*** Numeric value expected: ''',a,''' ***')
 1003 format('*** No value given for ''',a,''' ***')
 1050 format('*** A log file must be opened before using TRACE ***')
      end
