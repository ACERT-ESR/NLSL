! NLSL Version 1.9.0 beta 2/7/15
!----------------------------------------------------------------------
!                    =========================
!                         function IPFIND
!                    =========================
!
!  Search for the given token in one of the lists of the names of EPRLL
!  spectral calculation parameters, which are defined in module lpnam. 
!  Function return values are as follows:
!
! -200 < ipfind < -100   Specified name was the (axial tensor) alias 
!                        of a floating point parameter. 
!                        Returned value is the negative of the index
!                        of the parameter in the fepr array, minus 100. 
!
! -100 < ipfind < 0      Specified name was the (spherical tensor) alias 
!                        of a floating point parameter. 
!                        Returned value is the negative of the index of
!                        the parameter in the fepr array. 
!
!    0 < ipfind < 100    Specified name was a floating point parameter.
!                        Returned value is the index of the parameter
!                        in the fepr array
!
!        ipfind > 100    Specified name was an integer parameter.
!                        Returned value minus 100 is the index of the
!                        parameter in the iepr array
!                       
!
!  NOTE: The order of names in this list MUST correspond to the order
!   of parameters in modules eprprm and lpnam for this routine to work 
!   properly.
! 
!----------------------------------------------------------------------
      function ipfind(token,lth)

      use nlsdim
      use parcom
      use eprprm
      use lpnam

      implicit none
      integer :: ipfind,lth
      character :: token*30
!
      integer :: i
!
!----------------------------------------------------------------------
!     Search the list of floating point parameter names
!----------------------------------------------------------------------
      do i=1,nfprm
        if (token(:lth).eq.parnam(i)(:lth)) then
          ipfind=i
          return
        end if
      end do
!
!----------------------------------------------------------------------
!     Search both lists of floating point parameter aliases 
!     (these are all names of spherical tensor components)
!     Return negative index if found in alias1
!     Return negative index-100 if found in alias2
!----------------------------------------------------------------------
      do i=1,nalias
        if (token(:lth).eq.alias1(i)(:lth)) then
           ipfind=1-(IWXX+i)
           return
        else if (token(:lth).eq.alias2(i)(:lth)) then
           ipfind=-(99+IWXX+i)
           return
        end if 
      end do
!
!----------------------------------------------------------------------
!     Search the list of integer parameter names
!----------------------------------------------------------------------
      do i=1,niprm
        if (token(:lth).eq.iprnam(i)(:lth)) then
          ipfind=100+i
          return
        end if
      end do
!
!----------------------------------------------------------------------
!     Token was not found
!----------------------------------------------------------------------
      ipfind=0
      return
      end function ipfind

!----------------------------------------------------------------------
!                    =========================
!                        function ISFIND
!                    =========================
!
!  Search for the given token in one of several lists of strings.
!  Potential matches are stored in the arrays that hold the datafile,
!  basisid, and symbolic names. Datafile and basisid names are defined
!  via I/O within the datac and basisc routines. Symbolic names are
!  defined in the initialization segment of module lpnam--formerly
!  these were encoded in the block data section of nlstxt.f.
!
!  If a match is found, ISFIND returns an index into the array which 
!  contains the match. The index is positive if the match is in the
!  symbol array. The index is multipled by -1 if the match is in the
!  datafile or basisid array. Otherwise, 0 is returned.
!----------------------------------------------------------------------
      function isfind(token,lth)
!
      use nlsdim
      use expdat
      use basis
      use lpnam
      use stdio
      implicit none
!
      integer :: isfind,lth
      character :: token*30,tmpstr*30
!
      integer :: i,nlth
      integer, external :: itrim
!
!----------------------------------------------------------------------
!     Search the list of datafile names
!----------------------------------------------------------------------
      do i=1,nspc
         tmpstr=dataid(i)
         nlth=itrim(tmpstr)
         if (nlth.ge.lth) then
            if (token(:lth).eq.tmpstr(:lth)) then
               isfind=-i
               return
            end if
         end if
      end do
!
!----------------------------------------------------------------------
!     Search the list of basis set names
!----------------------------------------------------------------------
      do i=1,nbas
         tmpstr=basisid(i)
         nlth=itrim(tmpstr)
         if (nlth.ge.lth) then
            if (token(:lth).eq.tmpstr(:lth)) then
               isfind=-i
               return
            end if
         end if
      end do
!
!
!----------------------------------------------------------------------
!     Search the list of symbolic names
!----------------------------------------------------------------------
      tmpstr=token
      call touppr(tmpstr,lth)
      do i=1,nsymbl
         if (tmpstr(:lth).eq.symbol(i)(:lth)) then
            isfind=i
            return
         end if
      end do
!
      isfind=0
      return
      end
