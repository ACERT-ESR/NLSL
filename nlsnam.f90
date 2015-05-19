c NLSL Version 1.9.0 beta 2/12/15
c----------------------------------------------------------------------
c                   =========================
c                         module NLSNAM
c                   =========================
c
c Module for holding all the filenames associated with a given
c fitting parameter file ID, or datafile ID.
c
c      prname : <fileid>.PAR   Parameter file (formatted)
c      lgname : <fileid>.LOG   Log file for fitting results
c      trname : <fileid>.TRC   Trace of NLS steps
c
c      dtname : <dataid>.DAT   Datafile
c      spname : <dataid>.SPC   Splined experimental spectrum + fit
c
c NOTE: this module uses the nlsdim and stdio modules
c 
c----------------------------------------------------------------------
c
      module nlsnam
      use nlsdim
      use stdio
      implicit none
c
c     These were previously in common block nlsnam
      integer, save :: lthfnm, lthdnm
      character*30, save :: prname, lgname, trname, dtname, spname

c     This character*30 variable was unused - inlist
c
c     These were previously in common block filcom
      integer, save :: nfiles
      character*30, save :: files(MXFILE)
c
      end module nlsnam
