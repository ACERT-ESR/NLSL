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
c      dbname : <fileid>.DEB   Debugging aids
c      dtname : <dataid>.SPC   Input Data file
c      nsname : <dataid>.NSP   2D splined Data file (used in NLSP fit)
c      ftname : <dataid>.FIT   NLSP FIT
c      ixname : <index fileid> Specified in basid command
c
c NOTE: Files nlsdim.inc must be included BEFORE this file 
c in the source code
c----------------------------------------------------------------------
c
      character*(NAMELG) prname,lgname,dtname,ftname,nsname,trname,
     #             dbname,ixname
      integer lthfnm,lthdnm
c
      common /nlsnam/ lthfnm,lthdnm,prname,lgname,trname,dtname,
     #                ftname,nsname,ixname,dbname
c
      integer nfiles
      character*(NAMELG) files
      common /filcom/ nfiles,files(MXSPEC)
c
