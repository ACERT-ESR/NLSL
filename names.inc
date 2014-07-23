c Version 1.4
c----------------------------------------------------------------------
c                   =========================
c                      common block NLSNAM
c                   =========================
c
c Common block for all the filenames associated with a given
c fitting parameter file ID, or datafile ID.
c
c      prname : <fileid>.PAR   Parameter file (formatted)
c      lgname : <fileid>.LOG   Log file for fitting results
c      trname : <fileid>.TRC   Trace of NLS steps
c
c      dtname : <dataid>.DAT   Datafile
c      spname : <dataid>.SPC   Splined experimental spectrum + fit
c
c NOTE: Files nlsdim.inc AND stdio.inc must be included BEFORE this file 
c in the source code
c----------------------------------------------------------------------
c
      character*30 prname,lgname,spname,dtname,trname,inlist
      integer lthfnm,lthdnm
c
      common /nlsnam/ lthfnm,lthdnm,prname,lgname,trname,dtname,
     #                spname
c
      integer nfiles
      character*30 files
      common /filcom/ nfiles,files(MXFILE)
