c Version 1.5.1 beta 2/3/96
c----------------------------------------------------------------------
c                    =========================
c                       block data NLSTXT
c                    =========================
c
c     Initialization of the strings used for pre-defined parameter
c     names for the EPRLL family of programs.
c
c     NOTE: The order of names in this list MUST correspond to the
c      order of parameters in EPRPRM.INC for this routine to work 
c      properly.
c----------------------------------------------------------------------
      block data
c
      use nlsdim
      use eprprm
      use lpnam
      use errmsg
c
c----------------------------------------------------------------------
c
c ########## Definition of names of EPRLL parameters ###############
c
c  Note that blanks are specified for parameters that are not 
c  intended to be user-accessible, but must still be passed to EPRLS
c
      data parnam/'PHASE', 'GIB0',  'GIB2',  'WXX',   'WYY',
     #            'WZZ',   'GXX',   'GYY',   'GZZ',   'AXX',
     #            'AYY',   'AZZ',   'RX',    'RY',    'RZ', 
     #            'PML',   'PMXY',  'PMZZ',  'DJF',   'DJFPRP',
     #            'OSS',   'PSI',   'ALPHAD','BETAD', 'GAMMAD',  
     #            'ALPHAM','BETAM', 'GAMMAM','C20',   'C22',
     #            'C40',   'C42',   'C44',   'LB',    'DC20',
     #            'B0',    'GAMMAN','CGTOL', 'SHIFTR','SHIFTI',
     #            'RANGE', ' ',     ' '/
c
      data iprnam/'IN2',   'IPDF',  'IST',   'ML',    'MXY',
     #            'MZZ',   'LEMX',  'LOMX',  'KMN',   'KMX',   
     #            'MMN',   'MMX',   'IPNMX', 'NORT',  'NSTEP',
     #            'NFIELD','IDERIV',' ',     ' ',     ' ',
     #            ' ',     ' ',     ' ',     ' ' /
c
      data alias1/  'W1',   'W2',  'W3', 
     #              'G1',   'G2',  'G3',  
     #              'A1',   'A2',  'A3', 
     #              'RBAR', 'N',   'NXY' /
c
      data alias2/  'WPRP',  ' ',  'WPLL',  
     #              'GPRP',  ' ',  'GPLL',
     #              'APRP',  ' ',  'APLL',
     #              'RPRP',  ' ',  'RPLL' /
c
      data symstr/'CARTESIAN','SPHERICAL','AXIAL'/
c
      data symbol/'BROWNIAN','NONBROWNIA','ANISOVISCO','FREE','JUMP'/
      data symval/ 0,         1,           2,           1,     0     /
c  
      data eprerr /
     1     'Zero nuclear spin assumed',
     2     'Positive B0 assumed',
     3     'High-field approx. may not apply',
     4     'Nonzero potential with jump/free: ipdf set to 0',
     5     'Zero potential with aniso. viscos.: ipdf set to 0',
     6     'Discrete jumps with aniso. viscos.: ipdf set 0',
     7     'lemx must be 48 or less with potential',
     8     'Basis set adjusted to (iii,iii,iii,iii,iii,iii,ii)',
     9     'NSTEP too large: adjusted to iiiii',
     A     'CGTOL too small: adjusted to xxxxxxxxx',
     1     'SHIFTR too small: adjusted to xxxxxxxxx',
     2     'CG did not converge after iiiii steps',
     3     'Angle alpham not needed w/ axial g-matrix',
     4     'Angle alphad not needed w/ axial diffusion tensor',
     5     'Taylor series for Bessel function did not converge',
     6     'Zero B0',
     7     'Zero g-tensor element',
     8     'Error in tridiagonal matrix storage',
     9     'Too many matrix elements',
     A     'Matrix dimension too large',
     1     'User halt during matrix calculation',
     2     'User halt during CG tridiagonalization',
     3     'Too many tridiagonal matrix elements',
     4     'Bracketing did not find minimum within step bound',
     5     'Number of points in spectrum is zero',
     6     'Field range of spectrum is zero' /
c
       data minerr / 'Bad input parameters',
     1               'Chi-squared convergence',
     2               'X vector convergence',
     3               'Chi-squared/X vector convergence',
     4               'Gradient convergence',
     5               'Max function evaluations reached',
     5               'Max iterations reached',
     6               'FTOL too small',
     7               'XTOL too small',
     8               'Zero Jacobian or GTOL too small',
     9               'Terminated internally',
     A               'Single calculation completed'/
      end
