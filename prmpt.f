c Version 1.5 5/2/95
c----------------------------------------------------------------------
c                    =========================
c                       function PRMPT
c                    =========================
c
c  Issues a prompt for new values and reads an input line from the
c  user
c----------------------------------------------------------------------
      function prmpt( line )
      implicit none
c
      logical prmpt
      character*80 line
c
      include 'stdio.inc'
c
      write(luttyo,1000)
      read (luttyi,1001) line
      prmpt= .not. (line .eq. ' ')
      return
 1000 format(' New value(s) > ',$)
 1001 format(a)
      end

c----------------------------------------------------------------------
c                    =========================
c                        function FTOKEN
c                    =========================
c
c Interprets a string token as a floating point number
c
c----------------------------------------------------------------------
      function ftoken(line)
      implicit none
      double precision ftoken,ftemp
      character line*80,fstr*20,chr*1
      integer i1,i2,ixexp,ixdot
c
      integer MXLINE
      parameter(MXLINE=80)
c
      logical isnum,isdot,isexp
      isnum(chr)=((chr .ge. '0' .and. chr .le. '9') .or.
     1            (chr .eq. '+' .or. chr .eq. '-') )
      isdot(chr)=chr.eq.'.'
      isexp(chr)=(chr .eq. 'd' .or. chr .eq. 'D') .or.
     1           (chr .eq. 'e' .or. chr .eq. 'E')
c
c------------------------------------------------------------------------
c Skip all non-numeric characters before the beginning of next token
c------------------------------------------------------------------------
c
        i1=1
    1   chr = line(i1:i1)
        if ( i1.le.MXLINE .and. .not.
     !  (isnum(chr).or.isdot(chr).or.isexp(chr)) ) then
          i1=i1+1
          goto 1
        endif
c
        if (i1.gt.MXLINE) then
          ftoken=0.0d0
          return
        endif
c
c----------------------------------------------------------------------
c Include all numeric characters in the token
c----------------------------------------------------------------------
        ixdot=0
        ixexp=0
        if (isdot(chr)) ixdot=i1
        i2=i1+1
2       chr=line(i2:i2)
        if (isdot(chr)) ixdot=i2
        if (isexp(chr)) ixexp=i2
        if (i2.le.MXLINE .and.
     #     (isnum(chr).or.isdot(chr).or.isexp(chr)) ) then
          i2=i2+1
          goto 2
        endif
c
       i2=i2-1
c
c----------------------------------------------------------------------
c If user did not specify a decimal point, insert it into the string
c where appropriate. This complete kludge is necessary to make sure 
c internal FORTRAN 77 floating-point reads give the correct number. 
c Specifically, if neither '.' nor an exponent is given, append the '.'
c at the end of the string. If exponent but no '.' is given, insert
c the '.' before the exponent. Where is C when you need it?
c----------------------------------------------------------------------
       if (ixdot.eq.0) then
          if (ixexp.eq.0) then
            fstr = line(i1:i2)//'.'
          else
            fstr = line(i1:ixexp-1)//'.'//line(ixexp:i2)
          end if
       else
          fstr = line(i1:i2)
       end if
c
       line=line(i2+1:)
c
       read(fstr,1000,err=3) ftemp
       ftoken=ftemp
       return
c
 3    ftoken=0.0d0    
      return
 1000 format(bn,f20.10)
      end


c----------------------------------------------------------------------
c                    =========================
c                       function ITOKEN
c                    =========================
c
c  Interprets a string token as an integer number
c
c----------------------------------------------------------------------
      function itoken(line)
      implicit none
c
      character line*80,istr*20,chr*1
      logical isnum,isdot,dotfnd
      integer i1,i2,itoken,itemp
c
      integer MXLINE
      parameter(MXLINE=80)
c
      isnum(chr) = ( (chr .ge. '0' .and. chr .le. '9') .or.
     1       (chr .eq. '+' .or. chr .eq. '-') )
c
c------------------------------------------------------------------------
c Skip all non-numeric characters before the beginning of next token
c------------------------------------------------------------------------
c
        i1=1
    1   chr=line(i1:i1)
        if ( i1.le.MXLINE .and. .not. isnum(chr) ) then
          i1=i1+1
          goto 1
        endif
c
        if (i1.gt.MXLINE) then
          itoken=0
          return
        endif
c
c----------------------------------------------------------------------
c Include all numeric characters in the token
c----------------------------------------------------------------------
        i2=i1+1
2       chr=line(i2:i2)
        if (i2.le.MXLINE .and. isnum(chr) ) then
          i2=i2+1
          goto 2
        endif
c
       i2=i2-1
       istr=line(i1:i2)
       line=line(i2+1:)
       read(istr,1000,err=3) itemp
       itoken=itemp
       return
c
 3    itoken=0
      return
c
 1000 format(bn,i10)
      end
