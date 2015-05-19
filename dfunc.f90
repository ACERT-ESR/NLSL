c NLSL Version 1.9.0 beta 2/12/15
c----------------------------------------------------------------------
c                    =========================
c                           module DFUNC
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
      module dfunc
      implicit none
c
      double precision, save :: dlam, xi, cc, ss
c
      end module dfunc
