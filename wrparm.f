c Version 1.6  5/2/94
c**********************************************************************
c
c                    =========================
c                     subroutine WRPARM
c                    =========================
c
c  Outputs a formatted summary of EPRL calculation parameters to the
c  specified logical unit
c
c  modified from subroutine of same name by DJS
c  DEB 22-APRIL-1992
c**********************************************************************
c
      subroutine wrparm(prname,lth,lu)
c
      include 'rndoff.inc'
      include 'eprdat.inc'
c
      character*30 prname
      integer lu,lth
c
      integer i,j
c
      double precision ZERO
      parameter (ZERO=0.0D0) 
c
c######################################################################
c
      write (lu,3000) prname(:lth)
c
      write (lu,2000) gxx,gyy,gzz
      write (lu,2010) in2
      if (in2.ne.0) write (lu,2020) axx,ayy,azz
      write (lu,2280) wxx,wyy,wzz
      write (lu,2030) b0
      write (lu,3010) ipdf
      write (lu,2050) dx,dy,dz
      if (ipdf.ne.0) then
        write (lu,3020) tl,tkxy,tkzz
        write (lu,3030) pl,pkxy,pkzz
      end if
      if ((ipdf.eq.2).and.(ist.eq.0).and.(abs(tl).lt.rndoff).and.
     #     (abs(tkxy).lt.rndoff).and.(abs(tkzz).lt.rndoff)) then
        write (lu,2072) djfprp,djf
      else if (ist.ne.0) then
        write (lu,3040) ist,djf
      end if
      write (lu,2100) oss
      write (lu,2110) ipt
      if (ipt.ne.0) then
        write (lu,2120)
        do 910 i=1,ipt
          write (lu,2130) i,(ixp(i,j),j=1,2),cxp(i)
 910    continue
      end if
      write (lu,2140) psi
      write (lu,2150) itd
      if (itd.eq.1) write (lu,2160) ald,bed,gad
      write (lu,2150) itm
      if (itm.eq.1) write (lu,2162) alm,bem,gam
      write (lu,2170) lemx,lomx,kmx,mmx,ipnmx
      write (lu,2180) nstep
      write (lu,2190) itype
      if (itype.ge.1.or.nort.gt.1) then
        write (lu,2200) cgtol
        write (lu,2210) shiftr,shifti
      end if
      if (itype.ge.2.or.nort.gt.1) then
        write (lu,3050) fieldi,fieldf
        write (lu,3060) nfield
        write (lu,3065) btol
      end if
      if (nort.gt.1) then
         write(lu,3070) gib0,gib2
         write(lu,3080) nort
         write(lu,3090) ideriv
      end if
c
      write (lu,1000) 
c
      return 
c
c======================================================================
c               format statements
c======================================================================
c
 1000 format(//,2x,70('#'),//)
c
 2000 format(2x,'g-tensor [gxx,gyy,gzz] : ',3(f10.7,2x))
 2010 format(2x,'twice the nuclear spin [in2] : ',i2)
 2020 format(2x,'A-tensor [axx,ayy,azz] (gauss) : ',3(f9.3,2x))
 2030 format(2x,'static field [B0] (gauss) : ',f11.3)
c
 2050 format(2x,'diffusion tensor [dx,dy,dz] (1/sec) =',3(g11.4,2x))
 2072 format(2x,'diffusion tensor in the director frame ',
     #     '[djfprp,djf] (1/sec)',/,5x,2(g11.4,2x))
c
 2100 format(2x,'Heisenberg spin exchange frequency [oss] = ',g11.4)
c
 2110 format(2x,'number of terms in the potential [ipt] = ',i1)
 2120 format(2x,'coefficients of the potential : ')
 2130 format(10x,'ipt = ',i1,2x,'[l,k,coef.] = ',2(i1,2x),g11.4)
 2140 format(2x,'angle between B0 and local director [psi] ',
     #     '(degrees) : ',f7.2)
c
 2150 format(2x,'diffusion tilt flag [itd] = ',i1)
 2152 format(2x,'magnetic tilt flag [itm] = ',i1)
 2160 format(2x,'diffusion tilt [ald,bed,gad] (degrees) = ',
     # 2(f6.1,','),f6.1)
 2162 format(2x,'magnetic tilt [alm,bem,gam] (degrees) = ',
     # 2(f6.1,','),f6.1)
c
 2170 format(2x,'truncation values  [lemx,lomx,kmx,'
     #     ,'mmx,ipnmx] : ',5(I3,2x))
c
 2180 format(2x,'number of Lanczos/CG steps [nstep] : ',i5)
 2190 format(2x,'calculation type (0=L,1=CG,2=FS) [itype] : ',i1)
 2200 format(2x,'error tolerance for residual [cgtol] : ',g11.4)
 2210 format(2x,'origin shift [shiftr,shifti] :',2(2x,g11.4))
 2280 format(2x,'Inhomog. linewidth tensor [wxx,wyy,wzz] (gauss) : ',
     #     3(f9.3,2x))
 3000 format(//,2x,25('#'),3x,'file : ',a,3x,25('#'),//)
 3010 format(2x,'diffusion parameter [ipdf] = ',i1)
 3020 format(2x,'residence times [tl,tkxy,tkzz] : ',
     #     3(g11.4,2x))
 3030 format(2x,'non-Brownian exponents [pl,pkxy,pkzz] : ',/,
     #     10x,'0.5 --> free',/,10x,'1.0 --> jump',/,
     #     2x,3(f6.2,2x))
 3040 format(2x,'discrete jumps parameters [ist,djf] : ',
     #     i2,2x,g11.4)
 3050 format(2x,'field sweep [fieldi,fieldf] ',
     #     '(gauss) : ',f10.3,2x,f10.3)
 3060 format(2x,'number of field positions [nfield] : ',i3)
 3065 format(2x,'pruning tolerance [btol] : ',f7.4)
 3070 format(2x,'Gaussian line-broadening parameters [gib0,gib2] : ',
     #     2x,f8.3,2x,f8.3)
 3080 format(2x,'number of MOMD orientations [nort] : ',i3)
 3090 format(2x,'derivative mode [ideriv] ',i2)
c
c======================================================================
c
      end
