! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                   =========================
!                         module NLSNAM
!                   =========================
!
! Module for holding all the filenames associated with a given
! fitting parameter file ID, or datafile ID.
!
!      prname : <fileid>.PAR   Parameter file (formatted)
!      lgname : <fileid>.LOG   Log file for fitting results
!      trname : <fileid>.TRC   Trace of NLS steps
!
!      dtname : <dataid>.DAT   Datafile
!      spname : <dataid>.SPC   Splined experimental spectrum + fit
!
! NOTE: this module uses the nlsdim and stdio modules
! 
!----------------------------------------------------------------------
!
      module nlsnam
      use nlsdim
      use stdio
      implicit none
!
!     These were previously in common block nlsnam
      integer, save :: lthfnm, lthdnm
      character*30, save :: prname, lgname, trname, dtname, spname

!     This character*30 variable was unused - inlist
!
!     These were previously in common block filcom
      integer, save :: nfiles
      character*30, save :: files(MXFILE)
!
      end module nlsnam
