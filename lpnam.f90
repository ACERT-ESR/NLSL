c NLSL Version 1.9.0 2/3/15
c----------------------------------------------------------------------
c                    =========================
c                           module LPNAM
c                    =========================
c
c   Declaration and initialization of character arrays containing 
c   names of parameters, aliases, and various symbols.
c
c   This module replaces lpnam.inc and the first half of nlstxt.f.
c
c   NOTE: The order of names in this list MUST correspond to the
c   order of parameters in eprprm.f90 as well as the function ipfind
c   for the combined code to work properly.
c
c----------------------------------------------------------------------

      module lpnam
      use nlsdim
      implicit none

c ########## Definition of names of EPRLL parameters ###############
c
c  Note that blanks are specified for parameters that are not 
c  intended to be user-accessible, but must still be passed to EPRLS
c
      character*6, dimension(NFPRM), save ::
     # parnam = (/'PHASE ', 'GIB0  ', 'GIB2  ', 'WXX   ', 'WYY   ',
     #            'WZZ   ', 'GXX   ', 'GYY   ', 'GZZ   ', 'AXX   ',
     #            'AYY   ', 'AZZ   ', 'RX    ', 'RY    ', 'RZ    ', 
     #            'PML   ', 'PMXY  ', 'PMZZ  ', 'DJF   ', 'DJFPRP',
     #            'OSS   ', 'PSI   ', 'ALPHAD', 'BETAD ', 'GAMMAD',  
     #            'ALPHAM', 'BETAM ', 'GAMMAM', 'C20   ', 'C22   ',
     #            'C40   ', 'C42   ', 'C44   ', 'LB    ', 'DC20  ',
     #            'B0    ', 'GAMMAN', 'CGTOL ', 'SHIFTR', 'SHIFTI',
     #            'RANGE ', '      ', '      ' /)
c
      character*6, dimension(NIPRM), save ::
     # iprnam = (/'IN2   ', 'IPDF  ', 'IST   ', 'ML    ', 'MXY   ',
     #            'MZZ   ', 'LEMX  ', 'LOMX  ', 'KMN   ', 'KMX   ',   
     #            'MMN   ', 'MMX   ', 'IPNMX ', 'NORT  ', 'NSTEP ',
     #            'NFIELD', 'IDERIV', '      ', '      ', '      ',
     #            '      ', '      ', '      ', '      ' /)
c
      character*6, dimension(NALIAS), save ::
     # alias1 = (/  'W1    ', 'W2    ', 'W3    ', 
     #              'G1    ', 'G2    ', 'G3    ',  
     #              'A1    ', 'A2    ', 'A3    ', 
     #              'RBAR  ', 'N     ', 'NXY   ' /)
c
      character*6, dimension(NALIAS), save ::
     # alias2 = (/  'WPRP  ', '      ', 'WPLL  ',  
     #              'GPRP  ', '      ', 'GPLL  ',
     #              'APRP  ', '      ', 'APLL  ',
     #              'RPRP  ', '      ', 'RPLL  ' /)
c
      character*10, dimension(NSYMTR), save ::
     # symstr = (/'CARTESIAN ', 'SPHERICAL ', 'AXIAL     '/)
c
      character*10, dimension(NSYMBL), save ::
     # symbol = (/'BROWNIAN  ', 'NONBROWNIA', 'ANISOVISCO',
     #            'FREE      ', 'JUMP      '/)
c
      integer, dimension(NSYMBL), save ::
     # symval = (/ 0,            1,            2,
     #             1,            0          /)

      end module lpnam
