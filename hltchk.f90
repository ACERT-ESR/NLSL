c Version 1.5.2 b 7/15/96 
c----------------------------------------------------------------------
c                    =========================
c                      subroutine HLTCHK
c                    =========================
c
c     Check whether a user halt (control-C) or other error 
c     has occurred during the spectral calculation and set the
c     iflag error flag for LMNLS accordingly. Returns .true.
c     for fatal errors or user halts, and .false. otherwise.
c
c----------------------------------------------------------------------
      function hltchk(ierr,isite,ispec,iflag)
      implicit none
      integer ierr,iflag,isite,ispec
      logical hltchk
c
      use stdio
      use errmsg
c
      hltchk=.false.
      if (ierr.eq.0) return
c
c      ------------------------------------
c      Non-fatal errors: warn if requested    
c      ------------------------------------
      if (ierr.lt.FATAL) then
         if (warn) then
            write (luout,1000) 'Warning',isite,ispec,eprerr(ierr)
            if (luout.ne.luttyo) write (luttyo,1000) 'Warning',
     #                                   isite,ispec,eprerr(ierr)
            warn=.false.
         end if
c
c        ------------------------------
c        Fatal errors (non-user related)
c        ------------------------------
      else if (ierr.ge.FATAL .and. ierr.ne.MTXHLT
     #        .and. ierr.ne.CGHLT) then 
               
         write (luout,1000) 'Fatal err',isite,ispec,eprerr(ierr)
         if (luout.ne.luttyo) then
            write (luttyo,1000) 'Fatal err',isite,ispec,eprerr(ierr)
         end if
         hltchk=.true.
         iflag=-1
c
c        ------------------------------
c        Other (user halt)
c        ------------------------------
      else
         write (luout,1001) eprerr(ierr)
         if (luout.ne.luttyo) write (luttyo,1001) eprerr(ierr)
         hltchk=.true.
         iflag=-1
      end if
c
      return
c
c----------------------------------------------------------------------
c
 1000 format(/2x,'*** ',a,' site',i2,' spctrm',i2,': ',a)
 1001 format(/2x,'*** ',a)
      end
