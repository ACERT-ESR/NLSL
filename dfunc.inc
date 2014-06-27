c Version 1.5.1 beta 3/30/96
c----------------------------------------------------------------------
c                    =========================
c                       common block DFUNC
c                    =========================
c
c Contains storage for parameters related to the partial MOMD model.
c   dlam     The C20 coefficient for the director ordering potential
c   xi       Angle (in degrees) between the director ordering axis 
c            and B0
c   cc       The product cos(psi)*cos(xi) where psi is the angle betwwen
c            the local director and B0
c   ss       The product sin(psi)*sin(xi)
c----------------------------------------------------------------------
c
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
