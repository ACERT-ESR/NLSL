! Version 1.5.1 beta 2/3/96
!----------------------------------------------------------------------
!                     ====================
!                      subroutine WRITEC
!                     ====================
!
!  Process "write" command. Forces writing of calculated spectra
!  having no corresponding data
!
!----------------------------------------------------------------------
      subroutine writec( line )
!
      use nlsdim
      use expdat
      use parcom
      use lmcom
      use mspctr
      use stdio
!
      implicit none
      character*80 line
!
      integer i,j,k,l,lth
      double precision field
      character*40 oname
!
      if (written.eq.0) call wrspc()
      if (nspc.lt.nser) then
         do i = nspc+1,nser
            call gettkn( line, oname, lth )
            if (lth.eq.0) then
               write(luout,1000) i
               return
            end if
!
           open(ludisk,file=oname(:lth),status='unknown',
     #        access='sequential',form='formatted',err=10)
           field=sbi(i)
           do j=1,npts(i)
              k=ixsp(i)+j-1
!                                *** Multiple sites
              if (nsite.gt.1) then
                 write(ludisk,1045,err=10) field,-fvec(k),
     #              (sfac(l,i)*spectr(k,l),l=1,nsite)
!
!                                *** Single site
              else
                 write(ludisk,1045,err=10) field,-fvec(k)
              end if
              field=field+sdb(i)
           end do
           close (ludisk)
        end do
      end if
      return
!
 10   write (luout,1001) oname(:lth)
      if (luout.ne.luttyo) write (luttyo,1001) oname(:lth)
      close(ludisk)
      return
!
!
 1000 format('*** Output filename required for spectrum',i2,' ***')
 1001 format('*** Error opening or writing file ''',a,''' ***')
 1045 format(f10.3,6(' ',g14.7))
      end
            
