! NLSL Version 1.9.0 beta 2/9/15
!----------------------------------------------------------------------
!> @file strutl1.f90
!>             I/O UTILITY ROUTINES FOR NLS COMMAND INTERPRETER
!>
!>  Contains four functions and two subroutines:
!>
!>       getlin  -- issue a prompt, retrieve a line from command stream *
!>
!>       gettkn  -- extract a token from the given line (subroutine)
!>
!>       ungett  -- replace a token at the line's beginning (subroutine)
!>
!>       ftoken  -- returns a real number
!>
!>       itoken  -- returns an integer number
!>
!>       indtkn  -- returns an index, i.e., a parenthesized number *
!>
!>   * Note: two of these routines use module stdio. Subroutine touppr
!>     and function itrim were moved into strutl2.f90 because they have
!>     more general utility and a different set of dependencies.
!>----------------------------------------------------------------------


!----------------------------------------------------------------------
!                    =========================
!                       function GETLIN
!                    =========================
!----------------------------------------------------------------------
      function getlin( line )
!> @brief issue a prompt and retreive a line from command stream
!
      use stdio
      implicit none
!
      logical getlin
!
      integer ioerr
      character*80 line
!
      call wpoll
!
      if(lucmd.eq.luttyi) call lprmpt
      read (lucmd,1001,iostat=ioerr) line
!
!   -- Echo line if required
!
      if (ioerr.eq.0 .and. lucmd.ne.luttyi .and. luecho.ne.0)
     #     write(luecho,1001) line
      if (luecho.ne.0.and.luout.eq.lulog) write(lulog,1001) line
      getlin=ioerr.eq.0
      return
!
 1001 format(a)
      end function getlin

!------------------------------------------------------------------------
!                         ===================
!                          subroutine GETTKN
!                         ===================
! 
!> @brief
!>  Written for free-form input of parameters for slow-motional 
!>  calculations.  Returns a token consisting of nonseparator
!>  characters (separators are space, tab, and ',') with all control
!> @details
!>  characters filtered out. 
!>  Special cases:  '(', ')', '=', '*', and end-of-line, which are
!>  returned as single-character tokens.
!
! ------------------------------------------------------------------------
!
      subroutine gettkn(line,token,lth)
      implicit none
      integer lth
!
      integer LINLTH
      parameter (LINLTH=80)
!
      character line*80,token*30,chr*1
!
      integer i,j,ichr,ichar
!
! *** Function definitions
!
       logical issepr,is1tok,isctrl,istab
       isctrl(chr) = ichar(chr).lt.32
       istab(chr) = ichar(chr).eq.9
       issepr(chr) = chr.eq.' '.or. ichar(chr).eq.9 .or. chr.eq.','
     #               .or.chr.eq.';'
       is1tok(chr) = chr.eq.'('.or.chr.eq.')'.or.chr.eq.'*'
     #               .or.chr.eq.'='
! ***
!
!------------------------------------------
!  Find the next non-whitespace character 
!------------------------------------------
      i=0
 2    i=i+1
 3    chr=line(i:i)
!
!     -------------------------
!      skip control characters
!     -------------------------
      if (isctrl(chr).and. .not.istab(chr)) then
         line(i:)=line(i+1:)
         go to 3
      end if
!
      if (issepr(chr).and. i.lt.LINLTH) goto 2
!
      
      if (i.ge.LINLTH) then
         lth=0
         token=' '
         return
      end if
!
!     -----------------------------------
!     Check for single-character tokens 
!     -----------------------------------
      if (is1tok(chr)) then
         token=chr
         lth=1
         line=line(i+1:)
         return
      end if
!
!----------------------------------------------------------------
! Place the next continuous string of characters in the token
! (stop at whitespace, punctuation, and single-character tokens)
! Shorten the rest of the line accordingly
!----------------------------------------------------------------
      j=i
 4    j=j+1
 5    chr=line(j:j)
!
!     -----------------------
!     Skip control characters
!     -----------------------
      if (isctrl(chr).and. .not.istab(chr)) then
         line(j:)=line(j+1:)
         go to 5
      end if
!
      if ( issepr(chr) .or. is1tok(chr) ) then
         token=line(i:j-1)
         lth=j-i
         line=line(j:)
         return
      else
         go to 4
      end if
      end subroutine gettkn

!----------------------------------------------------------------------
!                     =========================
!                         subroutine UNGETT
!                     =========================
!  Replaces given token at the beginning of the given line
!  (Oh for the string functions of C..)
!----------------------------------------------------------------------
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

!----------------------------------------------------------------------
!                    =========================
!                         function FTOKEN
!                    =========================
!  Decodes a token into a floating point number
!----------------------------------------------------------------------
      function ftoken( token,lth,val )
      implicit none
      character token*30,tkn*30,tmptkn*30,chr*1
      integer i,lth,idot,ibrk
      double precision val
      logical ftoken
!     
! *** Function definitions --- these don't work with every FORTRAN 
!     implementation.
!
      logical isdot,isexp,isdig
      isdot(chr)=chr.eq.'.'
      isexp(chr)=(chr .eq. 'd' .or. chr .eq. 'D') .or.
     1           (chr .eq. 'e' .or. chr .eq. 'E')
      isdig(chr)=chr.ge.'0' .and. chr.le.'9'
! ***
!----------------------------------------------------------------------
!
       tkn=token
       idot=0
       ibrk=0
!
!----------------------------------------------------------------------
!     Find where a '.' is needed in the string
!     (this is to overcome the implied decimal used by some compilers)
!
!     Also, check for illegal characters (this is for FORTRAN compilers
!     that don't return an error when non-numeric characters 
!     are encountered in the read.)
!----------------------------------------------------------------------
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
!
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
!     
 13    ftoken=.false.
       return
!
 1000  format(bn,f20.10)
       end function ftoken

!----------------------------------------------------------------------
!                    =========================
!                       function ITOKEN
!                    =========================
!----------------------------------------------------------------------
      function itoken( token,lth,ival )
      implicit none
      character*30 token
      integer lth,ival
      logical itoken
!
!.................................................. 
!      read (token,1000,err=13) ival
! 12   itoken=.true.
!      return
!
! 13   itoken=.false.
!      return
!
! 1000 format(bn,i20)
!..................................................
!
      double precision fval
      logical ftoken
      external ftoken
      itoken=ftoken(token,lth,fval)
      ival=fval
      return
      end function itoken

!----------------------------------------------------------------------
!                    =========================
!                      function INDTKN
!                    =========================
!
!     Looks for a secondary index specified by the series of tokens
!     '(' <n> { ')' } or '(' '*' { ')' }. Returns n if found, -1 if 
!     '*' was specified, and 0 if no index was specified
!----------------------------------------------------------------------
      function indtkn(line)
!
      use stdio
      implicit none
!
      integer indtkn
      character line*80,token*30
!
      integer ival,lth
      logical wldcrd
!
      logical itoken
      external itoken
!
!######################################################################
!
      call gettkn(line,token,lth)
!
!----------------------------------------------------------------------
!     Look for parenthesis indicating a second index will be specified
!----------------------------------------------------------------------
      if (token.eq.'(') then
         call gettkn(line,token,lth)
!
!----------------------------------------------------------------------
!     Check for a valid index: '*' or an integer in range
!----------------------------------------------------------------------
         wldcrd=token.eq.'*'
!
!                                   *** Empty '()': index is 0
         if (token.eq.')') then
            indtkn=0
         else
!                                   *** Wildcard: return -1
            if (wldcrd) then
               indtkn=-1
            else
!                                   *** Check for legal number
! 
               if (itoken(token,lth,ival)) then
                  indtkn=ival
               else
!                                            *** Illegal index
!
                  write(luttyo,1000) token(:lth)
                  indtkn=0
               end if
            end if
            call gettkn(line,token,lth)
         end if
!
!----------------------------------------------------------------------
!      Check for a closing parenthesis (actually, this is optional)
!----------------------------------------------------------------------
         if (token.eq.')') call gettkn(line,token,lth)
!
!----------------------------------------------------------------------
!     No '(' found: index is 0
!----------------------------------------------------------------------
      else
         indtkn=0
      end if
!
!----------------------------------------------------------------------
!    Restore last token taken from the line
!----------------------------------------------------------------------
      call ungett(token,lth,line)
      return
!
 1000 format('*** Illegal index: ''',a,''' ***')
      end function indtkn

