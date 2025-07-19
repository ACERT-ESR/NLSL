c NLSL Version 1.5 beta 11/25/95
c----------------------------------------------------------------------
c                    =========================
c                      subroutine FIXC
c                    =========================
c
c  Removes the named parameter(s) from the list of variable parameters.
c
c fix all | <name> { ( <index> | * ) } {, <name>...}
c
c      name    : name of the parameter to be removed
c      index   : optional site/series index
c      *       : specifies all indices 
c----------------------------------------------------------------------
      subroutine fixc(line)
c
      use nlsdim
      use eprprm
c      use prmeqv
      use parcom
      use lpnam
      use stdio
c
      implicit none
      character*80 line
c
      integer i,ibd,ix,ix2,j,lth
      double precision prmn,prmx,prsc,step
      character token*30,prmID*9
c
      integer ipfind,indtkn
      logical ftoken
      external ftoken,ipfind,indtkn
c
c----------------------------------------------------------------------
c Get the name of the parameter
c----------------------------------------------------------------------
 1    call gettkn(line,token,lth)
      lth=min(lth,6)
      if (lth.le.0) return
c
      call touppr(token,lth)
      if (token(:lth).eq.'ALL') then
         if (nprm.gt.0) then
            write(luout,1001)
            if (luout.ne.luttyo) write(luttyo,1001)
         end if
c
         do i=1,NFPRM
            do j=1,MXSITE
               ixx(i,j)=0
            end do
         end do
         nprm=0
         return
      end if
c
      ix=ipfind(token,lth)
c
c----------------------------------------------------------------------
c     Check whether parameter may be varied 
c----------------------------------------------------------------------
      if (ix.eq.0 .or. ix.gt.nvprm) then
         write(luttyo,1002) token(:lth)
         return
      end if
c
      if (ix.lt.-100) then 
        prmID=alias2( -99-(IWXX+ix) )
      else if (ix.lt.0) then 
        prmID=alias1( 1-(IWXX+ix) )
      else
        prmID=parnam(ix)
      end if
c
c     --- Get secondary index
c
      ix2=indtkn( line )
c
      call rmvprm(ix,ix2,prmID)
      go to 1
c
c ###### format statements ########################################
c
 1001 format('*** All variable parameters have been fixed ***')
 1002 format('*** ''',a,''' is not a variable parameter ***')
      end
