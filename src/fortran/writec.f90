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
      use pylog_mod, only: log_enabled, log_buffer, ensure_log_buffer, flush_log_buffer
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
               if (log_enabled) then
                  call ensure_log_buffer(log_buffer)
                  write(log_buffer,1000) i
                  call flush_log_buffer()
               end if
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
 10   if (log_enabled) then
         call ensure_log_buffer(log_buffer)
         write(log_buffer,1001) oname(:lth)
         call flush_log_buffer()
      end if
      if (luout.ne.luttyo) write (luttyo,1001) oname(:lth)
      close(ludisk)
      return
!
!
 1000 format('*** Output filename required for spectrum',i2,' ***')
 1001 format('*** Error opening or writing file ''',a,''' ***')
 1045 format(f10.3,6(' ',g14.7))
      end
            
