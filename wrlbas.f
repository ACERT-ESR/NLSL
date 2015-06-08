c    Version 1.5  5/2/94
c**********************************************************************
c
c                       ===================
c                       SUBROUTINE : WRLBAS
c                       ===================
c
c       This routine lists the basis set indices contained in common
c       /indexl/ to the specified logical output unit. It is assumed
c       that subroutine lbasix has been called prior to calling this
c       routine, which is primarily designed to be called by lbll.f.
c
c       Modified by DEB 13-MAR-93 from program EPRBAS by DJS
c
c       Includes:
c               stddim.inc
c               eprdat.inc
c               indexl.inc
c
c       Uses:
c               ipar.f
c
c**********************************************************************
c
      subroutine wrlbas(prname,lth,lu)
c
      include 'stddim.inc'
      include 'eprdat.inc'
      include 'indexl.inc'
c
      character*30 prname
      integer lu,lth,i
c
c######################################################################
c
      write(lu,1050) prname(:lth)
      write(lu,1060)
      do i = 1,ndim
         write(lu,1070) i,l1(i),k1(i),m1(i),pi1(i),qi1(i)
      end do
c
      write(lu,1080) ndim
      write(lu,1000)
c
      return
c
c======================================================================
c     format statements
c======================================================================
c
 1000 format(//,2x,70('#'),//)
 1050 format(//,2x,70('#'),///,28x,'file : ',a,/,
     #     28x,15('-'),//)
 1060 format(/,15x,'BASIS SET',//,2x,'element',7x,'L',3x,
     #     'K',3x,'M',2x,'pI',2x,'qI',/,2x,32('-'))
 1070 format(4x,i5,4x,'|',4(i3,','),i3,'  >')
 1080 format(/,2x,'The dimension of the matrix is : ',i5)
c
c======================================================================
c
      end
