c NLSPMC Version 1.0 2/5/99
c----------------------------------------------------------------------
c                    =========================
c                      subroutine COMPS
c                    =========================
c
c comp <ncomps>
c
c      ncomps  : number of components for all spectra
c
c  NOTE: must apply to all spectra, ncomps defaults to 1
c
c      
c Special rules:
c
c----------------------------------------------------------------------
      subroutine comps ( line )
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
         write(luttyo,*) 'Number of components expected'
         return
      end if

      if (itoken(token,lth,ival)) then
* check max ival
	ncomps=ival
	return
      else
        write(luttyo,1001) token(:lth)
	return
      end if
 1001 format('*** number of components expected: ''',a,''' ***')
      end
