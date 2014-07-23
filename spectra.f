c NLSPMC Version 1.0 2/5/99
c----------------------------------------------------------------------
c                    =========================
c                      subroutine SPECTRA
c                    =========================
c
c spec <nspectra>
c
c      nspectra  : number of spectra
c
c  NOTE: defaults to 1
c
c      
c----------------------------------------------------------------------
      subroutine spectra ( line )
      implicit none
      character*80 line
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'parmequ.inc'
      include 'parms.inc'
      include 'lpnam.inc'
      include 'stdio.inc'
      include 'miscel.inc'
c
      integer ival,lth
      character token*30
      logical itoken
c
      call gettkn(line,token,lth)
c                                        *** No value given
      if (lth.eq.0) then
         write(luttyo,*) 'Number of spectra expected'
         return
      end if
c
      if (itoken(token,lth,ival)) then
        if (ival.gt.MXSPEC) then
	  write(luttyo,1002) MXSPEC
 1002     format('*** Maximum of ''',i2,''' spectra allowed ***')
	  nspectra = MXSPEC
	return
        else
	  nspectra = ival
	return
        end if
      else
        write(luttyo,1001) token(:lth)
	return
      end if
 1001 format('*** number of spectra expected: ''',a,''' ***')
      end
