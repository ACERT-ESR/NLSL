c NLSL Version 1.9.0 beta 2/11/15
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
c*********************************************************************
c
      module stdio
      implicit none
c
      integer, parameter :: lulog=4, luttyi=5, luttyo=6,
     #                      lutrc=7, ludisk=8
c
      integer, save :: lucmd, luecho, luout, hltcmd, hltfit, itrace
      logical, save :: warn
c
      end module stdio
