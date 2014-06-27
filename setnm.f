c----------------------------------------------------------------------
c                    =====================
c                     subroutine SETFIL
c                    =====================
c
c Given a generic name in <fileid>, this subroutine removes any
c file extension from <fileid> and determines the names for the
c various possible output files associated with the NLSL slow-
c motional fitting calculations by appending a file extension
c in the form '.xxx'. The extensions are
c
c      prname : <fileid>.PAR   Parameter file (formatted)
c      lgname : <fileid>.LOG   Log file for fitting results
c      trname : <fileid>.TRC   Trace of NLS steps
c
c  Includes
c     nlsdim.inc
c     nlsnam.inc   Definition of names in common /nlsnam/
c
c----------------------------------------------------------------------
      subroutine setfil( fileid )
      implicit none
c
      include 'nlsdim.inc'
      include 'nlsnam.inc'
c
      character*30 fileid
c
      integer iend,iroot
      character*1 chr
      external iroot
c
c######################################################################
c
      iend = iroot( fileid )
      if (iend .lt. 1) then
         fileid = 'noname'
         iend = 6
      endif
c
      prname=fileid(:iend)//'.par'
      lgname=fileid(:iend)//'.log'
      trname=fileid(:iend)//'.trc'
      lthfnm=iend+4
      return
      end
c----------------------------------------------------------------------
c                    =====================
c                     subroutine SETDAT
c                    =====================
c
c Given a generic name in <dataid>, this subroutine removes any
c file extension from <dataid> and determines the names for the
c various input and output files associated with the NLSL slow-
c motional fitting calculations by appending a file extension
c in the form '.xxx'. The extensions are
c
c      dtname : <dataid>.DAT   Datafile
c      spname : <dataid>.SPC   Calculated spectrum and fit
c
c  Includes
c     nlsdim.inc
c     nlsnam.inc   Definition of names in common /nlsnam/
c
c----------------------------------------------------------------------
      subroutine setdat( dataid )
      implicit none
c
      include 'nlsdim.inc'
      include 'nlsnam.inc'
c
      character*30 dataid
c
      integer iend,iroot
      character*1 chr
      external iroot
c
c######################################################################
c
      iend = iroot( dataid )
c
      if (iend .lt. 1) then
         dataid = 'noname'
         iend = 6
      endif
c
      dtname=dataid(:iend)//'.dat'
      spname=dataid(:iend)//'.spc'
      lthdnm=iend+4
      return
      end


c----------------------------------------------------------------------
c                    =========================
c                       function IROOT
c                    =========================
c  Strips the dot and trailing file extension from a file name
c  and returns the length of the root name
c----------------------------------------------------------------------
      function iroot( fileid )
      implicit none
c
      include 'nlsdim.inc'
      integer i,iroot
      character fileid*30,chr*1
c
      i=0
    1 i=i+1
      chr=fileid(i:i)
      if(chr.ne.'.' .and. chr.ne.' '.and.i.lt.30) goto 1
c
      fileid=fileid(:i-1)
      iroot=i-1
      return
      end
