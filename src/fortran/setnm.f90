!----------------------------------------------------------------------
!                    =====================
!                     subroutine SETFIL
!                    =====================
!
! Given a generic name in <fileid>, this subroutine removes any
! file extension from <fileid> and determines the names for the
! various possible output files associated with the NLSL slow-
! motional fitting calculations by appending a file extension
! in the form '.xxx'. The extensions are
!
!      prname : <fileid>.PAR   Parameter file (formatted)
!      lgname : <fileid>.LOG   Log file for fitting results
!      trname : <fileid>.TRC   Trace of NLS steps
!
!  Includes
!     nlsdim.inc
!     nlsnam.inc   Definition of names in common /nlsnam/
!
!----------------------------------------------------------------------
      subroutine setfil( fileid )
!
      use nlsdim
      use nlsnam
!
      implicit none
      character*30 fileid
!
      integer iend,iroot
      character*1 chr
      external iroot
!
!######################################################################
!
      iend = iroot( fileid )
      if (iend .lt. 1) then
         fileid = 'noname'
         iend = 6
      endif
!
      prname=fileid(:iend)//'.par'
      lgname=fileid(:iend)//'.log'
      trname=fileid(:iend)//'.trc'
      lthfnm=iend+4
      return
      end
!----------------------------------------------------------------------
!                    =====================
!                     subroutine SETDAT
!                    =====================
!
! Given a generic name in <dataid>, this subroutine removes any
! file extension from <dataid> and determines the names for the
! various input and output files associated with the NLSL slow-
! motional fitting calculations by appending a file extension
! in the form '.xxx'. The extensions are
!
!      dtname : <dataid>.DAT   Datafile
!      spname : <dataid>.SPC   Calculated spectrum and fit
!
!  Includes
!     nlsdim.inc
!     nlsnam.inc   Definition of names in common /nlsnam/
!
!----------------------------------------------------------------------
      subroutine setdat( dataid )
!
      use nlsdim
      use nlsnam
!
      implicit none
      character*30 dataid
!
      integer iend,iroot
      character*1 chr
      external iroot
!
!######################################################################
!
      iend = iroot( dataid )
!
      if (iend .lt. 1) then
         dataid = 'noname'
         iend = 6
      endif
!
      dtname=dataid(:iend)//'.dat'
      spname=dataid(:iend)//'.spc'
      lthdnm=iend+4
      return
      end


!----------------------------------------------------------------------
!                    =========================
!                       function IROOT
!                    =========================
!  Strips the dot and trailing file extension from a file name
!  and returns the length of the root name
!----------------------------------------------------------------------
      function iroot( fileid )
!
      use nlsdim
!
      implicit none
      integer i,iroot
      character fileid*30,chr*1
!
      i=0
    1 i=i+1
      chr=fileid(i:i)
      if(chr.ne.'.' .and. chr.ne.' '.and.i.lt.30) goto 1
!
      fileid=fileid(:i-1)
      iroot=i-1
      return
      end
