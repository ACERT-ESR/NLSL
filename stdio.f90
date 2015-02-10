c NLSL Version 1.9.0 beta 2/10/15
c*********************************************************************
c
c       STDIO: declarations of standard I/O parameters 
c
c       Notes:
c               1) The parameters define here are used as
c
c                       luttyi : logical unit for terminal input
c                       luttyo : logical unit for terminal output
c                       lulog  : logical unit for log file
c                       lutrc  : logical unit for trace file
c                       ludisk : logical unit for other disk files
c                       
c                       NOTE: The units that will be used for additional
c                       disk input depend on the mxfile parameter defined
c                       in "nlsdim.inc". If mxfile is greater than 1, 
c                       units <ludisk+1>, <ludisk+2>, ... <ludisk+mxfile>
c                       are assumed to be available.
c
c               2) The appropriate values for these parameters are
c                  highly operating system and compiler dependent.
c
c               3) Fortran 90/95 does not have intrinsic functions for
c                  converting a string to uppercase, or for finding the
c                  index of the first blank character in a string.
c                  Therefore, subroutine touppr and function itrim are
c                  being supplied in this I/O-related module. They were
c                  previously included in the strutl.f collection.
c
c*********************************************************************
c
      module stdio
      implicit none
c
      integer, parameter :: lulog=4, luttyi=5, luttyo=6,
     #                      lutrc=7, ludisk=8
c
      integer, save :: lucmd, luecho, luout, hltcomd, hltfit, itrace
      logical :: warn
c
      contains
c
c----------------------------------------------------------------------
c                    =========================
c                       subroutine TOUPPR
c                    =========================
c
c Converts string argument to uppercase
c----------------------------------------------------------------------
c
      subroutine touppr(string,lth)
      implicit none
      character(len=*) :: string
      character(len=1) :: chr
      integer :: i, ich, ichar, lth
c
      do i=1,lth
         chr=string(i:i)
         ich=iachar(chr)
         if (ich.ge.97 .and. ich.le.122) string(i:i)=achar(ich-32)
      end do
      return
      end subroutine touppr
c
c----------------------------------------------------------------------
c                    =========================
c                         function ITRIM
c                    =========================
c
c   Returns position of first blank in a string
c----------------------------------------------------------------------
c
      function itrim( string )
      implicit none
      integer :: itrim, j, lth
      character(len=*) :: string
      lth=len(string)
      do j=1,lth
         if (string(j:j).eq.' ') then
            itrim=j-1
            return
         end if
      end do
c
c     If here: no blanks in string, return full length
c     This differs from F77 index() intrinsic, which would return 0 
      itrim=lth
      return
      end function itrim
c
      end module stdio
