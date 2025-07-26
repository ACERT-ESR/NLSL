! NLSL Version 1.9.0 beta 2/12/15
!----------------------------------------------------------------------
!                    =========================
!                           module DFUNC
!                    =========================
!
! Contains storage for parameters related to the partial MOMD model.
!   dlam     The C20 coefficient for the director ordering potential
!   xi       Angle (in degrees) between the director ordering axis 
!            and B0
!   cc       The product cos(psi)*cos(xi) where psi is the angle betwwen
!            the local director and B0
!   ss       The product sin(psi)*sin(xi)
!----------------------------------------------------------------------
!
      module dfunc
      implicit none
!
      double precision, save :: dlam, xi, cc, ss
!
      end module dfunc
