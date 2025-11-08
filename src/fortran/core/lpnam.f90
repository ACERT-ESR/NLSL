! NLSL Version 1.9.0 beta 2/3/15
!----------------------------------------------------------------------
!                    =========================
!                           module LPNAM
!                    =========================
!
!   Declaration and initialization of character arrays containing 
!   names of parameters, aliases, and various symbols.
!
!   This module replaces lpnam.inc and the first half of nlstxt.f.
!
!   NOTE: The order of names in this list MUST correspond to the
!   order of parameters in eprprm.f90 as well as the function ipfind
!   for the combined code to work properly.
!
!----------------------------------------------------------------------

      module lpnam
      use nlsdim
      implicit none

! ########## Definition of names of EPRLL parameters ###############
!
!  Note that blanks are specified for parameters that are not 
!  intended to be user-accessible, but must still be passed to EPRLS
!
      integer, save :: parnam_strlen
      data parnam_strlen /6/
      character*6 parnam(NFPRM)
      save parnam
      data parnam /'PHASE ', 'GIB0  ', 'GIB2  ', 'WXX   ', 'WYY   ',
     #            'WZZ   ', 'GXX   ', 'GYY   ', 'GZZ   ', 'AXX   ',
     #            'AYY   ', 'AZZ   ', 'RX    ', 'RY    ', 'RZ    ', 
     #            'PML   ', 'PMXY  ', 'PMZZ  ', 'DJF   ', 'DJFPRP',
     #            'OSS   ', 'PSI   ', 'ALPHAD', 'BETAD ', 'GAMMAD',  
     #            'ALPHAM', 'BETAM ', 'GAMMAM', 'C20   ', 'C22   ',
     #            'C40   ', 'C42   ', 'C44   ', 'LB    ', 'DC20  ',
     #            'B0    ', 'GAMMAN', 'CGTOL ', 'SHIFTR', 'SHIFTI',
     #            'RANGE ', '      ', '      ' /
!
      integer, save :: iprnam_strlen
      data iprnam_strlen /6/
      character*6 iprnam(NIPRM)
      save iprnam
      data iprnam /'IN2   ', 'IPDF  ', 'IST   ', 'ML    ', 'MXY   ',
     #            'MZZ   ', 'LEMX  ', 'LOMX  ', 'KMN   ', 'KMX   ',   
     #            'MMN   ', 'MMX   ', 'IPNMX ', 'NORT  ', 'NSTEP ',
     #            'NFIELD', 'IDERIV', '      ', '      ', '      ',
     #            '      ', '      ', '      ', '      ' /
!
      integer, save :: alias_strlen
      data alias_strlen /6/
      character*6 alias1(NALIAS)
      save alias1
      data alias1 /  'W1    ', 'W2    ', 'W3    ', 
     #              'G1    ', 'G2    ', 'G3    ',  
     #              'A1    ', 'A2    ', 'A3    ', 
     #              'RBAR  ', 'N     ', 'NXY   ' /
!
      character*6 alias2(NALIAS)
      save alias2
      data alias2 /  'WPRP  ', '      ', 'WPLL  ',  
     #              'GPRP  ', '      ', 'GPLL  ',
     #              'APRP  ', '      ', 'APLL  ',
     #              'RPRP  ', '      ', 'RPLL  ' /
!
      integer, save :: symstr_strlen
      data symstr_strlen /10/
      character*10, dimension(NSYMTR), save, target ::
     # symstr = (/'CARTESIAN ', 'SPHERICAL ', 'AXIAL     '/)
!
      integer, save :: symbol_strlen
      data symbol_strlen /10/
      character*10 symbol(NSYMBL)
      save symbol
      data symbol /'BROWNIAN  ', 'NONBROWNIA', 'ANISOVISCO',
     #            'FREE      ', 'JUMP      '/
!
      integer, dimension(NSYMBL), save, target ::
     # symval = (/ 0,            1,            2,
     #             1,            0          /)
      end module lpnam
