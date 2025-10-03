! NLSL Version 1.5 beta 11/25/95
!----------------------------------------------------------------------
!                    =========================
!                      subroutine FIXC
!                    =========================
!
!  Removes the named parameter(s) from the list of variable parameters.
!
! fix all | <name> { ( <index> | * ) } {, <name>...}
!
!      name    : name of the parameter to be removed
!      index   : optional site/series index
!      *       : specifies all indices 
!----------------------------------------------------------------------
      subroutine fixc(line)
!
      use nlsdim
      use eprprm
!      use prmeqv
      use parcom
      use lpnam
      use stdio
      use pylog_mod, only: log_enabled, log_buffer, ensure_log_buffer, flush_log_buffer
!
      implicit none
      character*80 line
!
      integer i,ix,ix2,j,lth
      character token*30,prmID*9
!
      integer ipfind,indtkn
      logical ftoken
      external ftoken,ipfind,indtkn
!
!----------------------------------------------------------------------
! Get the name of the parameter
!----------------------------------------------------------------------
 1    call gettkn(line,token,lth)
      lth=min(lth,6)
      if (lth.le.0) return
!
      call touppr(token,lth)
      if (token(:lth).eq.'ALL') then
         if (nprm.gt.0) then
            if (log_enabled) then
               call ensure_log_buffer(log_buffer)
               write(log_buffer,1001)
               call flush_log_buffer()
            end if
            if (luout.ne.luttyo) write(luttyo,1001)
         end if
!
         do i=1,NFPRM
            do j=1,MXSITE
               ixx(i,j)=0
            end do
         end do
         nprm=0
         return
      end if
!
      ix=ipfind(token,lth)
!
!----------------------------------------------------------------------
!     Check whether parameter may be varied 
!----------------------------------------------------------------------
      if (ix.eq.0 .or. ix.gt.nvprm) then
         write(luttyo,1002) token(:lth)
         return
      end if
!
      if (ix.lt.-100) then 
        prmID=alias2( -99-(IWXX+ix) )
      else if (ix.lt.0) then 
        prmID=alias1( 1-(IWXX+ix) )
      else
        prmID=parnam(ix)
      end if
!
!     --- Get secondary index
!
      ix2=indtkn( line )
!
      call rmvprm(ix,ix2,prmID)
      go to 1
!
! ###### format statements ########################################
!
 1001 format('*** All variable parameters have been fixed ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
      end
