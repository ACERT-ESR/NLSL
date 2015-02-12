c NLSL Version 1.9.0 beta 2/9/15
c----------------------------------------------------------------------
c> @file strutl1.f90
c>             I/O UTILITY ROUTINES FOR NLS COMMAND INTERPRETER
c>
c>  Contains four functions and two subroutines:
c>
c>       getlin  -- issue a prompt, retrieve a line from command stream *
c>
c>       gettkn  -- extract a token from the given line (subroutine)
c>
c>       ungett  -- replace a token at the line's beginning (subroutine)
c>
c>       ftoken  -- returns a real number
c>
c>       itoken  -- returns an integer number
c>
c>       indtkn  -- returns an index, i.e., a parenthesized number *
c>
c>   * Note: two of these routines use module stdio. Subroutine touppr
c>     and function itrim were moved into strutl2.f90 because they have
c>     more general utility and a different set of dependencies.
c>----------------------------------------------------------------------


c----------------------------------------------------------------------
c                    =========================
c                       function GETLIN
c                    =========================
c----------------------------------------------------------------------
      function getlin( line )
c> @brief issue a prompt and retreive a line from command stream
c
      use stdio
      implicit none
c
      logical getlin
c
      integer ioerr
      character*80 line
c
      call wpoll
c
      if(lucmd.eq.luttyi) call lprmpt
      read (lucmd,1001,iostat=ioerr) line
c
c   -- Echo line if required
c
      if (ioerr.eq.0 .and. lucmd.ne.luttyi .and. luecho.ne.0)
     #     write(luecho,1001) line
      if (luecho.ne.0.and.luout.eq.lulog) write(lulog,1001) line
      getlin=ioerr.eq.0
      return
c
 1001 format(a)
      end function getlin

c------------------------------------------------------------------------
c                         ===================
c                          subroutine GETTKN
c                         ===================
c 
c> @brief
c>  Written for free-form input of parameters for slow-motional 
c>  calculations.  Returns a token consisting of nonseparator
c>  characters (separators are space, tab, and ',') with all control
c> @details
c>  characters filtered out. 
c>  Special cases:  '(', ')', '=', '*', and end-of-line, which are
c>  returned as single-character tokens.
c
c ------------------------------------------------------------------------
c
      subroutine gettkn(line,token,lth)
      implicit none
      integer lth
c
      integer LINLTH
      parameter (LINLTH=80)
c
      character line*80,token*30,chr*1
c
      integer i,j,ichr,ichar
c
c *** Function definitions
c
       logical issepr,is1tok,isctrl,istab
       isctrl(chr) = ichar(chr).lt.32
       istab(chr) = ichar(chr).eq.9
       issepr(chr) = chr.eq.' '.or. ichar(chr).eq.9 .or. chr.eq.','
     #               .or.chr.eq.';'
       is1tok(chr) = chr.eq.'('.or.chr.eq.')'.or.chr.eq.'*'
     #               .or.chr.eq.'='
c ***
c
c------------------------------------------
c  Find the next non-whitespace character 
c------------------------------------------
      i=0
 2    i=i+1
 3    chr=line(i:i)
c
c     -------------------------
c      skip control characters
c     -------------------------
      if (isctrl(chr).and. .not.istab(chr)) then
         line(i:)=line(i+1:)
         go to 3
      end if
c
      if (issepr(chr).and. i.lt.LINLTH) goto 2
c
      
      if (i.ge.LINLTH) then
         lth=0
         token=' '
         return
      end if
c
c     -----------------------------------
c     Check for single-character tokens 
c     -----------------------------------
      if (is1tok(chr)) then
         token=chr
         lth=1
         line=line(i+1:)
         return
      end if
c
c----------------------------------------------------------------
c Place the next continuous string of characters in the token
c (stop at whitespace, punctuation, and single-character tokens)
c Shorten the rest of the line accordingly
c----------------------------------------------------------------
      j=i
 4    j=j+1
 5    chr=line(j:j)
c
c     -----------------------
c     Skip control characters
c     -----------------------
      if (isctrl(chr).and. .not.istab(chr)) then
         line(j:)=line(j+1:)
         go to 5
      end if
c
      if ( issepr(chr) .or. is1tok(chr) ) then
         token=line(i:j-1)
         lth=j-i
         line=line(j:)
         return
      else
         go to 4
      end if
      end subroutine gettkn

c----------------------------------------------------------------------
c                     =========================
c                         subroutine UNGETT
c                     =========================
c  Replaces given token at the beginning of the given line
c  (Oh for the string functions of C..)
c----------------------------------------------------------------------
      subroutine ungett(token,lth,line)
      implicit none
      character line*80,tmplin*80,token*30
      integer lth
      if (lth.gt.0.and.lth.lt.80) then
        tmplin=line
        line=token(:lth) // tmplin
      end if
      return
      end subroutine ungett

c----------------------------------------------------------------------
c                    =========================
c                         function FTOKEN
c                    =========================
c  Decodes a token into a floating point number
c----------------------------------------------------------------------
      function ftoken( token,lth,val )
      implicit none
      character token*30,tkn*30,tmptkn*30,chr*1
      integer i,lth,idot,ibrk
      double precision val
      logical ftoken
c     
c *** Function definitions --- these don't work with every FORTRAN 
c     implementation.
c
      logical isdot,isexp,isdig
      isdot(chr)=chr.eq.'.'
      isexp(chr)=(chr .eq. 'd' .or. chr .eq. 'D') .or.
     1           (chr .eq. 'e' .or. chr .eq. 'E')
      isdig(chr)=chr.ge.'0' .and. chr.le.'9'
c ***
c----------------------------------------------------------------------
c
       tkn=token
       idot=0
       ibrk=0
c
c----------------------------------------------------------------------
c     Find where a '.' is needed in the string
c     (this is to overcome the implied decimal used by some compilers)
c
c     Also, check for illegal characters (this is for FORTRAN compilers
c     that don't return an error when non-numeric characters 
c     are encountered in the read.)
c----------------------------------------------------------------------
       do 10 i=1,lth
          chr=tkn(i:i)
          if (isdot(chr))  then
             idot=i
          else if (isexp(chr)) then
             ibrk=i
          else if (.not. isdig(chr).and.chr.ne.'-') then
             go to 13
          end if
 10    continue
c
       if (idot.eq.0) then
          if (ibrk .eq. 0) then
             tkn=tkn(:lth)//'.'
          else
             tmptkn=tkn(ibrk:)
             tkn=tkn(:ibrk-1) // '.' // tmptkn
          end if
          lth=lth+1
       end if
 
       read(tkn,1000,err=13) val
       ftoken=.true.
       return
c     
 13    ftoken=.false.
       return
c
 1000  format(bn,f20.10)
       end function ftoken

c----------------------------------------------------------------------
c                    =========================
c                       function ITOKEN
c                    =========================
c----------------------------------------------------------------------
      function itoken( token,lth,ival )
      implicit none
      character*30 token
      integer lth,ival
      logical itoken
c
c.................................................. 
c      read (token,1000,err=13) ival
c 12   itoken=.true.
c      return
c
c 13   itoken=.false.
c      return
c
c 1000 format(bn,i20)
c..................................................
c
      double precision fval
      logical ftoken
      external ftoken
      itoken=ftoken(token,lth,fval)
      ival=fval
      return
      end function itoken

c----------------------------------------------------------------------
c                    =========================
c                      function INDTKN
c                    =========================
c
c     Looks for a secondary index specified by the series of tokens
c     '(' <n> { ')' } or '(' '*' { ')' }. Returns n if found, -1 if 
c     '*' was specified, and 0 if no index was specified
c----------------------------------------------------------------------
      function indtkn(line)
c
      use stdio
      implicit none
c
      integer indtkn
      character line*80,token*30
c
      integer ival,lth
      logical wldcrd
c
      logical itoken
      external itoken
c
c######################################################################
c
      call gettkn(line,token,lth)
c
c----------------------------------------------------------------------
c     Look for parenthesis indicating a second index will be specified
c----------------------------------------------------------------------
      if (token.eq.'(') then
         call gettkn(line,token,lth)
c
c----------------------------------------------------------------------
c     Check for a valid index: '*' or an integer in range
c----------------------------------------------------------------------
         wldcrd=token.eq.'*'
c
c                                   *** Empty '()': index is 0
         if (token.eq.')') then
            indtkn=0
         else
c                                   *** Wildcard: return -1
            if (wldcrd) then
               indtkn=-1
            else
c                                   *** Check for legal number
c 
               if (itoken(token,lth,ival)) then
                  indtkn=ival
               else
c                                            *** Illegal index
c
                  write(luttyo,1000) token(:lth)
                  indtkn=0
               end if
            end if
            call gettkn(line,token,lth)
         end if
c
c----------------------------------------------------------------------
c      Check for a closing parenthesis (actually, this is optional)
c----------------------------------------------------------------------
         if (token.eq.')') call gettkn(line,token,lth)
c
c----------------------------------------------------------------------
c     No '(' found: index is 0
c----------------------------------------------------------------------
      else
         indtkn=0
      end if
c
c----------------------------------------------------------------------
c    Restore last token taken from the line
c----------------------------------------------------------------------
      call ungett(token,lth,line)
      return
c
 1000 format('*** Illegal index: ''',a,''' ***')
      end function indtkn

