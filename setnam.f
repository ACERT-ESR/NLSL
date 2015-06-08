c     Version 1.5    5/2/94
c----------------------------------------------------------------------
c                    =====================
c                     subroutine SETNAM
c                    =====================
c
c Given a generic name in <fileid>, this subroutine removes any
c file extension from <fileid> and determines the names for the
c various possible output files associated with EPRLL/EPRLF
c slow-motional calculations by appending a file extension
c in the form '.xxx'. The extensions are
c
c      prname : <fileid>.PAR   Parameter file (unformatted)
c      fmname : <fileid>.FMT   Parameter listing (formatted)
c      bsname : <fileid>.BSS   Basis set indices
c      spname : <fileid>.SPC   Calculated spectrum
c      rlname : <fileid>.LOG   Log file from last program run
c      tdname : <fileid>.TDL   Tridiagonal SLE matrix (unformatted)
c      tfname : <fileid>.TDF   Tridiagonal SLE matrix (formatted)
c      tsname : <fileid>.MTS   MTS database file
c      mtname : <fileid>.MTX   SLE matrix (unformatted)
c      mfname : <fileid>.MTF   SLE matrix (formatted)
c      egname : <fileid>.EGV   SLE eigenvalues/weights
c      vtname : <fileid>.STV   Starting vector
c      ixname : <fileid>.IND   Basis set index list
c
c  Includes
c     stddim.inc
c     fnames.inc   Definition of names in common /fnames/
c
c  NB: Although the code in namdef.inc and this routine are standard
c      FORTRAN 77, they may not be portable to installations that
c      do not fully implement the standard.
c
c----------------------------------------------------------------------
      subroutine setnam
c
      include 'stddim.inc'
      include 'fnames.inc'
c
      integer i,iend
      character*1 chr
c
      i=0
    1 i=i+1
      chr=fileid(i:i)
      if(chr.ne.'.' .and. chr.ne.' '.and.i.lt.mxlth) goto 1
c
      iend=i-1
c
      if (iend .lt. 1) then
         fileid = 'noname'
         iend = 6
      endif
c
      prname=fileid(:iend)//'.par'
      fmname=fileid(:iend)//'.fmt'
      bsname=fileid(:iend)//'.bss'
      spname=fileid(:iend)//'.spc'
      rlname=fileid(:iend)//'.log'
      tdname=fileid(:iend)//'.tdl'
      tfname=fileid(:iend)//'.tdf'
      tsname=fileid(:iend)//'.mts'
      mtname=fileid(:iend)//'.mtx'
      mfname=fileid(:iend)//'.mtf'
      egname=fileid(:iend)//'.egv'
      vtname=fileid(:iend)//'.stv'
      ixname=fileid(:iend)//'.ind'
      nameOK=.true.
      namlth=iend+4
      return
      end
