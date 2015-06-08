c   Version 1.5   5/2/94
c---------------------------------------------------------------------
c                      =========================
c                            program ORDER
c                      =========================
c This program calculates the molecular order parameters corresponding 
c to a user-specified axial orienting potential. The potential is 
c expressed as
c
c                                L   L
c    U(phi,theta,0) = Sum   Sum c   D   (phi,theta,0)
c                      L     K   K   0K
c 
c where D^L_{0K} is a Wigner rotation matrix element (equivalent to
c the spherical harmonic Y^L_K), c^L_K are user-specified coefficients,
c and the sums are restricted to L=0 to 4 and K=0 to L (for even L,K only)
c
c The coefficients are specified in a list separated by spaces in 
c the order c20, c22, c40, c42, c44. The user may specify only as 
c many coefficients as needed; however, in order to specify a potential
c with a nonzero c40, for example, it is necessary to specify zeroes
c for c20 and c22. 
c
c The program calculates an order parameter for each nonzero coefficient
c specified. The order parameter is the average of <D^L_{0K}(phi,theta)> 
c weighted by the value of exp(-U(phi,theta)) and integrated over
c angles theta and phi.
c
c A blank input line exits the program.
c
c Includes:
c   rndoff.inc
c 
c Uses:
c    ccrints.f  (Integration routines ccrint ccrin1)
c    ftht       (Function returning theta-dependent part of U)
c    fphi       (Function returning phi-dependent part of U)
c    nxttok     (Function returning next numeric string on input line)
c
c    The last three functions are included in this file.
c    Written by David E. Budil 
c    Cornell Department of Chemistry, March 1992
c 
c----------------------------------------------------------------------

        program order
        implicit real*8(a-h,o-z)
        parameter (zero=0.0D0, one=1.0D0, small=1.0D-16)
        parameter (ncoeff=5)
        character*5 ident(ncoeff)
        character*3 coeff(ncoeff)
        character token*20, nxttkn*20, line*80
        external ftht

        common/param/acc,c(ncoeff),kount
        common/func/d20,d40,bd22,bd42,bd44

        data coeff / 'c20', 'c22','c40','c42','c44' /
        data ident / '<D20>', '<D22>', '<D40>', '<D42>', '<D44>'/

1       ipt = 0
        line = ' '
        write (*,1002) (coeff(j),j=1,5)
1002    format(' Enter coefficients (',4(a,','),a,'): ')
        read (*,'(a)') line

        do 10 i = 1, ncoeff
          token = nxttkn(line)

c----------------------------------------------------------------------
c Blank line stops input (stops program if there has been no input)
c----------------------------------------------------------------------
          if (token .eq. ' ') then
            if (i .eq. 1) stop
            goto 15
          endif

c----------------------------------------------------------------------
c Append a dot to all tokens which do not have one
c (this kludge is for the fortran internal read of floating-point 
c numbers)
c----------------------------------------------------------------------
          idot=0
          do 5 ix=1,20          
          if (token(ix:ix).eq.'.') idot=i
    5     continue
c          
          iblk=0
          do 6 ix=1,20
          j=21-ix
          if (token(j:j).eq.' ') iblk=j
    6     continue
c
          if (idot.eq.0.and.iblk.ne.0)  token(iblk:iblk) = '.'
          read (token,'(f20.10)',err=12) c(i)
          if (c(i) .ne. zero) ipt = i
10        continue

        goto 15

12      write(*,'('' Error in input line.'')')
        goto 1

15      acc=0.001

        kount=0
        if (ipt .gt. 0) then
c---------------------------------------------
c Integrate unweighted potential function ftht
c---------------------------------------------

          call ccrint(zero,one,acc,small,pfn,nup,ftht,id)

c-----------------------------------------------
c Integrate potential weighted by D20, D22, etc.
c-----------------------------------------------
          write (*,1015) pfn
 1015     format(' Integral of U =',g12.6/' Order parameters:')
          do 25 kount=1,ipt
            call ccrint(zero,one,acc,small,fz,nup,ftht,id)
            pder=fz/pfn
            write(*,'(3x,a,'' = '',f9.5 )') ident(kount),pder
25          continue
            write(*,'(/)')
        endif

        goto 1
        end



        double precision function ftht(ctht)
        implicit real*8(a-h,o-z)
        common/param/acc,c(5),kount
        common/func/d20,d40,bd22,bd42,bd44
        external fphi

c---------------------------------------------------------
c  definition of some constants
c    a22 = sqrt(3/2), a42 = sqrt(5/8),  a44 = sqrt(35/8)
c---------------------------------------------------------
        parameter (pi=3.14159265358979d0, 
     1             a22=1.22474487139159d0,
     2             a42=0.790569415042095d0,
     3             a44=1.04582503316759d0 )

        parameter(one=1.0D0, zero=0.0D0, small=1.0D-16 )

        ctht2 = ctht*ctht
        stht2 = one-ctht2
        d20 = 1.5d0*ctht2-0.5d0
        d40 = ( (4.375d0*ctht2)-3.75d0)*ctht2+0.375d0
        bd22 = a22*stht2
        bd42 = a42*stht2*(7.0d0*ctht2-one)
        bd44 = a44*stht2*stht2
c
        call ccrin1(zero,pi,acc,small,result,nup,fphi,id)
        ftht = result
        return
        end

                               
        double precision function fphi(phi)
c
        implicit real*8(a-h,o-z)
        common/param/acc,c(5),kount
        common/func/d20,d40,bd22,bd42,bd44
c
        parameter (one=1.0D0,two=2.0D0)
c
        c2phi = dcos(phi+phi)
        c4phi = two*c2phi*c2phi - one
c
        fphi = dexp(  c(1)*d20
     2              + c(2)*bd22*c2phi
     3              + c(3)*d40
     4              + c(4)*bd42*c2phi
     5              + c(5)*bd44*c4phi )
c
        if(kount.eq.0) return
        if(kount.eq.1) fphi=d20*fphi
        if(kount.eq.2) fphi=bd22*fphi*c2phi
        if(kount.eq.3) fphi=d40*fphi
        if(kount.eq.4) fphi=bd42*fphi*c2phi
        if(kount.eq.5) fphi=bd44*fphi*c4phi
        return
        end                                                     


        function nxttkn( line )
        character line*80, nxttkn*20
        character*1 chr
        logical isnum

        isnum(chr) = ( (chr .ge. '0' .and. chr .le. '9') .or.
     1       (chr .eq. '+' .or. chr .eq. '-') .or.
     2       (chr .eq. 'd' .or. chr .eq. 'D') .or.
     3       (chr .eq. 'e' .or. chr .eq. 'E') .or.
     4       (chr .eq. '.') )

c------------------------------------------------------------------------
c Skip all non-numeric characters before the beginning of next token
c------------------------------------------------------------------------
        i1 = 1
1       if ( i1 .le. 80 .and. .not. isnum( line(i1:i1) )  ) then
          i1 = i1+1
          goto 1
        endif

        if (i1 .gt. 80) then
          nxttkn = ' '
          return
        endif
c----------------------------------------------------------------------
c Include all numeric characters in the token
c----------------------------------------------------------------------
        i2 = i1+1
2        if (i2 .le. 80 .and. isnum( line(i2:i2) ) ) then
          i2 = i2+1
          goto 2
        endif

        i2 = i2-1
        nxttkn = line(i1:i2)
        line = line(i2+1:)
       return
       end
