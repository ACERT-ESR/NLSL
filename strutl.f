c NLSPMC Version 1.0 2/5/99
c----------------------------------------------------------------------
c             I/O UTILITY ROUTINES FOR NLS COMMAND INTERPRETER
c
c  Contains the subroutines:
c       getlin  -- issue a prompt and retreive a line from command stream
c       gettkn  -- return a token from the given line
c       touppr  -- converts a string to uppercase
c       ftoken  -- returns a real number
c       itoken  -- returns an integer number
c       indtkn  -- returns an index (parenthesized number)
c----------------------------------------------------------------------


c----------------------------------------------------------------------
c                    =========================
c                       function GETLIN
c                    =========================
c----------------------------------------------------------------------
      function getlin( line )
      implicit none
      logical getlin
c
      include 'stdio.inc'
      include 'limits.inc'
c
      integer ioerr
      character*(LINELG) line
c
c      call wpoll
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
      end

c------------------------------------------------------------------------
c                         ===================
c                          subroutine GETTKN
c                         ===================
c 
c  Written for free-form input of parameters for slow-motional 
c  calculations.  Returns a token consisting of nonseparator
c  characters (separators are space, tab, and ',') with all control
c  characters filtered out. 
c  Special cases:  '(', ')', '=', '*', and end-of-line, which are
c  returned as single-character tokens.
c
c ------------------------------------------------------------------------
c
      subroutine gettkn(line,token,lth)
      implicit none
      include 'limits.inc'
      integer lth
c
c
      character line*(LINELG),token*(WORDLG),chr*1
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
      if (issepr(chr).and. i.lt.LINELG) goto 2
c
      
      if (i.ge.LINELG) then
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
      end

c----------------------------------------------------------------------
c                     =========================
c                       subroutine UNGETT
c                     =========================
c  Replaces given token at the beginning of the given line
c  (Oh for the string functions of C..)
c----------------------------------------------------------------------
      subroutine ungett(token,lth,line)
      implicit none
      include 'limits.inc'
      character line*(LINELG),tmplin*(LINELG),token*(WORDLG)
      integer lth
      if (lth.gt.0.and.lth.lt.LINELG) then
        tmplin=line
        line=token(:lth) // tmplin
      end if
      return
      end
      
c----------------------------------------------------------------------
c                    =========================
c                       subroutine TOUPPR
c                    =========================
c Converts string argument to uppercase
c----------------------------------------------------------------------

      subroutine touppr(string,lth)
      implicit none
      include 'limits.inc'
      character string*(WORDLG),chr
      integer i,ich,ichar,lth
c
c *** Function definition
c      logical islc
c      islc(chr) = ichar(chr).ge.97 .and. ichar(chr).le.122
c ***
c----------------------------------------------------------------------
c
      do i=1,lth
         chr=string(i:i)
         ich=ichar(chr)
         if (ich.ge.97 .and. ich.le.122) string(i:i)=char( ich-32 )
      end do
      return
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine FTOKEN
c                    =========================
c  Decodes a token into a floating point number
c----------------------------------------------------------------------
      function ftoken( token,lth,val )
      implicit none
      include 'limits.inc'
      character token*(WORDLG),tkn*(WORDLG),tmptkn*(WORDLG),chr*1
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
       end


c----------------------------------------------------------------------
c                    =========================
c                       function ITOKEN
c                    =========================
c----------------------------------------------------------------------
      function itoken( token,lth,ival )
      implicit none
      include 'limits.inc'
      character*(WORDLG) token
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
      end

c----------------------------------------------------------------------
c                       =========================
c                          function ITRIM
c                       =========================
c
c   Returns position of first blank in a string
c----------------------------------------------------------------------

      function itrim( string )
      implicit none
      integer itrim,j,lth
      character*(*) string
      lth=len(string)
      do 1 j=1,lth
         if (string(j:j).eq.' ') go to 2
 1    continue
      itrim=lth
      return
c
 2    itrim=j-1
      return
      end

c----------------------------------------------------------------------
c                    =========================
c                      subroutine INDTKN
c                    =========================
c
c     Looks for a secondary index specified by the series of tokens
c     '(' <n,m> { ')' } or '(' '*,*' { ')' }.  The information is 
c     returned in the variables spectid and compid.  In the case of 
c     no parameters or *, zeros are returned meaning all sites or 
c     spectra.  Return -1 on error.  
c     
c----------------------------------------------------------------------
      subroutine indtkn(line,spectid,compid)
      implicit none
      include 'limits.inc'
      include 'stdio.inc'
      include 'miscel.inc'
c
      character line*(LINELG),token*(WORDLG)
      integer ival,lth,spectid,compid
      logical wldcrd
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
      if ( token .ne. '(' ) then		! got ()?
        spectid=0	! no
        compid=0
        call ungett(token,lth,line)	! restore last token to line
        return
      else	! yes, parse the info
         call gettkn(line,token,lth)
c
c----------------------------------------------------------------------
c     Check for a valid index: '*' or an integer in range
c----------------------------------------------------------------------
         wldcrd=token.eq.'*'
c
c                                   *** Empty '()': indices are 0
         if (token.eq.')') then
            spectid=0
	    compid=0
            return
         else
c                                   *** Wildcard: return 0
            if (wldcrd) then
               spectid=0
            else
c                                   *** Check for legal number
c 
               if (itoken(token,lth,ival)) then
                  spectid=ival
               else
c                                            *** Illegal index
c
                  write(luttyo,1000) token(:lth)
                  spectid=-1
                  compid=-1
                  return
               end if
            end if
c now go for second index:
            call gettkn(line,token,lth)
            wldcrd=token.eq.'*'
            if (wldcrd) then
              compid=0		! second entry is *
	      call gettkn(line,token,lth)
	      if (token.eq.')') return
	      call ungett(token,lth,line)
	      return	! if not ), put it back and return
	    else
              if (itoken(token,lth,ival)) then
                compid=ival
                call gettkn(line,token,lth)
	        if (token.eq.')') return
	        call ungett(token,lth,line)
                return  ! if not ), put it back and return
              else
                if (token.eq.')') then
                  compid=0
	          return
	        else
	          write(luttyo,1000) token(:lth)
                  spectid=-1
                  compid=-1
                  return
                end if
              end if
            end if
         end if
      end if
c
 1000 format('*** Illegal index: ''',a,''' ***')
      end

