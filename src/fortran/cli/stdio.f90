! NLSL Version 1.9.0 beta 2/11/15
!*********************************************************************
!
!       STDIO: declarations of standard I/O parameters 
!
!       Notes:
!               1) The parameters define here are used as
!
!                       luttyi : logical unit for terminal input
!                       luttyo : logical unit for terminal output
!                       lulog  : logical unit for log file
!                       lutrc  : logical unit for trace file
!                       ludisk : logical unit for other disk files
!                       exitprog : flag to exit the program
!                       
!                       NOTE: The units that will be used for additional
!                       disk input depend on the mxfile parameter defined
!                       in "nlsdim.inc". If mxfile is greater than 1, 
!                       units <ludisk+1>, <ludisk+2>, ... <ludisk+mxfile>
!                       are assumed to be available.
!
!               2) The appropriate values for these parameters are
!                  highly operating system and compiler dependent.
!
!*********************************************************************
!
      module stdio
      implicit none
!
      integer, parameter :: lulog=4, luttyi=5, luttyo=6,
     #                      lutrc=7, ludisk=8
!
      integer, save :: lucmd, luecho, luout, hltcmd, hltfit, itrace
      logical, save :: warn, exitprog
!
      end module stdio
