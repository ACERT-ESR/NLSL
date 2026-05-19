! NLSL Version 1.5.1 beta 11/25/95
!**********************************************************************
!                    =========================
!                       subroutine MATRLL
!                    =========================
!
!       Subroutine for calculating matrix elements for EPRLL calculation.
!
!       Updated to improve matrix storage algorithm. This version differs
!       from the published version (Biol. Magnetic Resonance, Vol. 8)
!       in that it only stores half of the matrix. Changes are noted
!       in the code given below.
!
!       Updated to provide for a nonaxial diffusion tensor R
!       (only in the case of Brownian motion).
!
!       Updated to replace input of residence times tl, tkxy, tkzz
!       with the product of residence time and diffusion rate constants:
!       pml, pmzz, pmxy. (Differs from parameters used in the EPRLL/EPRCGL
!       programs in this respect)
!
!       flags
!       -----
!
!             fld :   True if |lr-lc| <=2 , false otherwise.
!                     Used as flag for bandwidth of Liouville operator
!                     matrix elements. The restriction on the 
!                     L-bandwidth due the inclusion of spin operators 
!                     of ranks 0 and 2.
!
!             fldmd : True if |lr-lc| <=2 or mr=mc, false otherwise.
!                     Used as flag for possible non-zero matrix 
!                     elements. The L-bandwidth restriction is
!                     explained above and the diffusion operator 
!                     matrix elements can be found only when mr=mc.
!
!             flk0 :  True if lr=lc and kr=kc, false otherwise.
!                     Used as a flag for determing the existence of 
!                     non-zero Liouville (A-tensor) and diffusion 
!                     (Heisenberg spin exchange) matrix elements.
!
!       written by DJS 29-AUG-87
!
!       modified by JPG 06-JUN-88 to correct the anisotropic 
!                                 viscosity terms in the diffusion 
!                                 operator to conform to Polnaszek
!                                 and Freed 1975 paper
!       modified by DJS JUN-90 to update matrix storage. Comments by DEB.
!
!       modified by DEB: calculation of isotropic part of diffusion
!                        operator moved outside iqnr loop. Fully anisotropic
!                        rotational diffusion tensor implemented.
!
!       Uses:
!               cd2km.f
!               anxlk.f
!               ipar.f
!               w3j.f
!
!**********************************************************************
!
      subroutine matrll(ierr)
!
      use nlsdim
      use rnddbl
      use eprmat
      use eprprm
      use errmsg
      use physcn
      use stdio
      use maxl
!
      implicit none
      integer i,j,l,kdi,ksi,li,ierr
      integer lr,lrprod,iparlr,kr,krmx,mr,mrmx,
     #     ipnr,ipnrmn,ipnrmx,iqnr,iqnrmn,iqnrmx,
     #     lc,lcmax,iparlc,ld,ldabs,kc,kcmn,kcmx,kd,ks,
     #     kdabs,ksabs,mc,mcmn,mcmx,md,ms,mdabs,msabs,
     #     ipnc,ipncmn,ipncmx,iqnc,iqncmn,iqncmx,
     #     ipnd,ipns,ipndab,ipnsab,iqnd,iqns,iqndab,
     #     kip,nrow,ncol,nelr,neli,nel,ipnrsg,ipncsg,
     #     jmc,jmr,krmn,mrmn,jkc,jkr,jkd,jmd,kcsgn,krsgn,mcsgn,mrsgn,
     #     ipn2r,ipn2mn,ipn2mx,iqn2r,iqn2mn,iqn2mx,
     #     ipn2c,ipn2cmn,ipn2cmx,iqn2c,iqn2cmn,iqn2cmx,
     #     ipn2d,ipn2s,ipn2dab,ipn2sab,iqn2d,iqn2s,iqn2dab,
     #     kip2
!
      double precision dsq2,dsq6,dsq15,dsq23
      double precision ct,rpl,rmi,rz,rj,rp,ossc,fii,ra0,
     #     cpmkr,cnpr,c1,c2,c3,c4,
     #     cdfd,cnl,cplkc,cnk,z,ra1,rg1,rw1,ra2,rg2,rw2,
     #     ra,rg,rw,cliou,cgam,cplmc,wliou1,wliou2,cdff,cdffr,
     #     cdffu,cnpc,cnp,cnorm,d1,d2,sa1,sa2,sw1,sw2,clioi1,fki,
     #     clioi2,clioua,sg1,sg2,clioug,ctemp,cgamw
      double precision cjkr,cjkc,cjmc,zeen
      double precision d2km(2,5,5)
      double precision fii2,ra02,zeen2,cnp2r,cnp2c,cnp2
      double precision sa2b1,sa2b2,fki2,clioi2b1,clioi2b2,ra2b
      double precision ra2b1,ra2b2,d1b,d2b
      logical fdpqi2
!
      logical flk0,fld,fldmd,fdpqi,fdjkm
!
      double precision ZERO,ONE
      parameter (ZERO=0.0D0,ONE=1.0D0)
!
      integer ipar 
      double precision w3j
      external ipar,w3j
!
!######################################################################
!
!.... Debugging purposes only!
!      open (unit=20,file='nlmat.tst',status='unknown',
!     #     access='sequential',form='formatted')
!............................
!
!----------------------------------------------------------------------
!     Scale diffusion constants
!----------------------------------------------------------------------
!
      ct=g0*BETAE/HBAR
      rpl=0.5d0*(dx+dy)/ct
      rmi=0.25*(dx-dy)/ct
      rz=dz/ct
      rj=djf/ct
      rp=djfprp/ct
      ossc=oss/ct
!
!----------------------------------------------------------------------
!     Define constants
!----------------------------------------------------------------------
!
      dsq2=dsqrt(2.0D0)
      dsq6=dsqrt(6.0D0)
      dsq15=dsqrt(1.5D0)
      dsq23=dsqrt(2.0D0/3.0D0)
!
      fii=0.25D0*dble(in2*(in2+2))
!
      if(ipdf.ne.1) then
        pml=ZERO
        pmzz=ZERO
        pmxy=ZERO
      end if
!
      if ((ipdf.eq.2) .and. (ist.eq.0)) then
        rpl=rpl+rp
        rz=rz+rp
      end if
!
      if(in2.eq.0) then
        ra0=ZERO
        zeen=ZERO
      else
        ra0=-dsq15*a0
        zeen=b0*gamman/GAMMAE
      end if
!
      if(in2b.le.0) then
        ra02=ZERO
        zeen2=ZERO
        fii2=ZERO
      else
        fii2=0.25D0*dble(in2b*(in2b+2))
        ra02=-dsq15*a20
        zeen2=b0*gamman2/GAMMAE
      end if
!
!----------------------------------------------------------------------
!                 2
!     Calculate d   (psi) as in Brink and Satchler
!                k,m
!                          2
!     Note : dlm(i+3,j+3)=d   (psi)
!                          i,j
!----------------------------------------------------------------------
!
      call cd2km(d2km,ZERO,psi,ZERO)
!
!----------------------------------------------------------------------
!                               L
!     Calculate the quantities X  used in potential-dependent
!                               K
!     portion of diffusion operator
!----------------------------------------------------------------------
!
      call anxlk(rpl,rmi,rz)
!
!----------------------------------------------------------------------
!     Initialize counters
!----------------------------------------------------------------------
!
      nrow=0
      ncol=0
      nelr=0
      neli=0
      nel=0
!
!----------------------------------------------------------------------
!     *** loop over lr ***
!----------------------------------------------------------------------
!
      do 300 lr=0,lemx,ldelta
!
        lrprod=lr*(lr+1)
        iparlr=ipar(lr)
!
!      --- Skip odd L greater than odd L maximum
!
        if((lr.gt.lomx).and.(iparlr.ne.1)) go to 300
!
!----------------------------------------------------------------------
!      *** loop over kr ***
!----------------------------------------------------------------------
!
!      --- Find limits for K index for row
!
        krmx = min0(lr,kmx)
        krmn = max0(-lr,kmn)
        do 320 krsgn=krmn,krmx,kdelta
!
!         --- Find symmetry index jK from sign of K
!
          kr=abs(krsgn)
          jkr=isign(1,krsgn)
          if (krsgn.eq.0) jkr=iparlr
          if (jkr.lt.jkmn) goto 320
!
!        --------------
!        Find parities
!        --------------
          cjkr=jkr
!
!
!----------------------------------------------------------------------
!     *** loop over mr ***
!----------------------------------------------------------------------
!
!         --- Find limits for M index for row
!
          mrmx=min0(lr,mmx)
          mrmn=max0(-lr,mmn)
          do 340 mrsgn=mrmn,mrmx
!
!         --- Find symmetry index jM from sign of M
!             (if M=0, jM is specified by sign of pI)
!
             mr=abs(mrsgn)
             jmr=isign(1,mrsgn)
!
!           --------------
!           find parities 
!           --------------
            cpmkr=ipar(mr+kr)
!
!----------------------------------------------------------------------
!     Calculate isotropic part of diffusion operator
!     (Code moved to here since it does not depend on nuclear
!     transition indices--DEB)
!----------------------------------------------------------------------
!
                c1=dble(lrprod)
                if(ipdf.ne.0) then
                  c2=ONE+pml*c1
                  if(dabs(expl).gt.RNDOFF) c2=c2**expl
                  c1=c1/c2
                end if
!
                cdfd=rpl*c1
                c1=dble(kr*kr)
                c2=rz
                c3=rpl
!
                if(ipdf.ne.0) then
                  c4=ONE+pmzz*c1
                  if(dabs(expkzz).gt.RNDOFF) c4=c4**expkzz
                  c2=c2/c4
                  c4=ONE+pmxy*c1
                  if(dabs(expkxy).gt.RNDOFF) c4=c4**expkxy
                  c3=c3/c4
                end if
!
                cdfd=cdfd+c1*(c2-c3)
                if(ipdf.eq.2) cdfd=cdfd+dble(mr*mr)*(rj-rp)
!
                if(ist.ne.0) then
                  i=kr/ist
                  if((i*ist).ne.kr) cdfd=cdfd+rj
                end if
!
!            --------------------------------
!            Find limits for pI index for row
!            --------------------------------
            if(ipsi0.eq.0) then
              ipnrmn=mr
              ipnrmx=mr
            else
              ipnrmx = min(in2,ipnmx)
              if(mr.eq.0 .and. jmmn.eq.1) then
                ipnrmn=0
              else
                ipnrmn=-ipnrmx
              end if
            end if
!
!----------------------------------------------------------------------
!      *** loop over ipnr ***
!----------------------------------------------------------------------
!
            do 350 ipnrsg=ipnrmn,ipnrmx
!
!         --- Find symmetry index jM from sign of pI if undetermined
!
               if (mr.eq.0) then
                  jmr=isign(1,ipnrsg)
                  ipnr=abs(ipnrsg)
                  if (ipnr.eq.0) jmr=iparlr
               else
                  ipnr=ipnrsg
               end if
               if(jmr.lt.jmmn) goto 350
!
!----------------------------------------------------------------------
!     Find scaling factor for matrix elements
!----------------------------------------------------------------------
!
              if((ipnr.eq.0).and.(mr.eq.0)) then
                cnpr=dsq2
              else
                cnpr=ONE
              end if
!
!----------------------------------------------------------------------
!     Find limits for qI index for row
!----------------------------------------------------------------------
!
              iqnrmx=in2-abs(ipnr)
              iqnrmn=-iqnrmx
!
!----------------------------------------------------------------------
!      *** loop over iqnr ***
!----------------------------------------------------------------------
!
              do 360 iqnr=iqnrmn,iqnrmx,2
!
!----------------------------------------------------------------------
!     *** loops over ipn2r, iqn2r (nucleus 2) ***
!----------------------------------------------------------------------
              if (in2b.gt.0) then
                ipn2mn=-in2b
                ipn2mx=in2b
              else
                ipn2mn=0
                ipn2mx=0
              end if
              do 362 ipn2r=ipn2mn,ipn2mx
              if (in2b.gt.0 .and. ipn2r.eq.0) then
                cnp2r=dsq2
              else
                cnp2r=ONE
              end if
              if (in2b.gt.0) then
                iqn2mn=-(in2b-iabs(ipn2r))
                iqn2mx=in2b-iabs(ipn2r)
              else
                iqn2mn=0
                iqn2mx=0
              end if
              do 364 iqn2r=iqn2mn,iqn2mx,2
              if (in2b.le.0 .and. iqn2r.ne.0) goto 364
!
!----------------------------------------------------------------------
!     Increment row counter
!----------------------------------------------------------------------
!
                nrow=nrow+1
                if (nrow.gt.MXDIM) then
                  ierr=DIMBIG
                  go to 9999
                end if
!
                jzmat(nrow)=neli+1
                kzmat(nrow)=nelr+1
!
!----------------------------------------------------------------------
!      Calculate limits for L index for columns
!----------------------------------------------------------------------
!
                lcmax = min0( lr+lband, lemx )
!
!----------------------------------------------------------------------
!     Initialize column counter
!----------------------------------------------------------------------
!
                ncol=nrow
!
!               Check for user input or halt
!
                call wpoll
                if (hltcmd.ne.0 .or. hltfit.ne.0) then
                   ierr=MTXHLT
                   return
                end if
!
!----------------------------------------------------------------------
!      *** loop over lc ***
!----------------------------------------------------------------------
!
                do 400 lc=lr,lcmax,ldelta
!
                  iparlc=ipar(lc)
!
!----------------------------------------------------------------------
!      skip if L is odd and greater than odd L maximum
!----------------------------------------------------------------------
!
                  if((lc.gt.lomx).and.(iparlc.ne.1)) go to 400
!
!----------------------------------------------------------------------
!     find limits for non-zero Liouville operator
!     matrix elements
!----------------------------------------------------------------------
!
                  ld=lr-lc
                  ldabs=abs(ld)
                  fld=abs(ld).le.2
!
!                 --- find scaling factor
!
                  cnl=dsqrt((2.D0*lr+1.D0)*(2.D0*lc+1.D0))
!
!                 --- Find limits of K index for columns

                  kcmx=min(lc,kmx)
!
                  if (nrow.eq.ncol) then
                    kcmn=kr
                  else
                    kcmn=max(-lc,kmn)
                  endif
!
!----------------------------------------------------------------------
!      *** loop over kc ***
!----------------------------------------------------------------------
!
                  do 420 kcsgn=kcmn,kcmx,kdelta
!
!                    --- Find symmetry index of jK from sign of K
!
                     kc=abs(kcsgn)
                     jkc=isign(1,kcsgn)
                     if (kcsgn.eq.0) jkc=iparlc
                     if(jkc.lt.jkmn) go to 420
!
!                   --- find sums, differences and parities
!
                     kd=kr-kc
                     ks=kr+kc
                     jkd=jkr-jkc
                     kdabs=abs(kd)
                     ksabs=abs(ks)
                     cplkc=ipar(lc+kc)
                     cjkc=jkc
!
!                  -- Find scaling factor
!
                    if((kr.eq.0).and.(kc.eq.0)) then
                      cnk=0.5D0
                    else if((kr.ne.0).and.(kc.ne.0)) then
                      cnk=ONE
                    else
                      cnk=ONE/dsq2
                    end if
!
!----------------------------------------------------------------------
!              Calculate rhombic part of diffusion tensor
!----------------------------------------------------------------------
!
                    if (kr+2.eq.kc .and. jkd.eq.0) then
                       cdffr=rmi*sqrt(dble((lr+kc-1)*(lr+kc)*
     #                      (lr-kc+1)*(lr-kc+2) ) )/cnk
                    else if (kr-2.eq.kc .and. jkd.eq.0) then
                       cdffr=rmi*sqrt(dble((lr-kc-1)*(lr-kc)*
     #                      (lr+kc+1)*(lr+kc+2) ) )/cnk
                    else
                       cdffr=ZERO
                    end if
!
!---------------------------------------------------------------------
!       Calculate R(mu,l) as in Meirovitch, Igner, Igner, Moro, & Freed
!       J. Phys. Chem. 77 (1982) p. 3935, Eqs. A42 & A44
!       mu = g, A, W
!---------------------------------------------------------------------
!
                    if((kdabs.le.2).and.fld) then
                      z=w3j(lr,2,lc,kr,-kd,-kc)
                      if (jkd.eq.0) then
                         ra1=fad(1,kd+3)*z
                         rg1=fgd(1,kd+3)*z
                         rw1=fwd(1,kd+3)*z
                      else
                         ra1=fad(2,kd+3)*z*cjkr
                         rg1=fgd(2,kd+3)*z*cjkr
                         rw1=fwd(2,kd+3)*z*cjkr
                      end if
                    else
                      ra1=ZERO
                      rg1=ZERO
                      rw1=ZERO
                    end if
!
                    if((ksabs.le.2).and.fld) then
                       if (jkd.eq.0) then
                          z=w3j(lr,2,lc,kr,-ks,kc)
                          ra2=fad(1,ks+3)*z
                          rg2=fgd(1,ks+3)*z
                          rw2=fwd(1,ks+3)*z
                       else
                          ra2=fad(2,ks+3)*z*cjkr
                          rg2=fgd(2,ks+3)*z*cjkr
                          rw2=fwd(2,ks+3)*z*cjkr
                       end if
                    else
                      ra2=ZERO
                      rg2=ZERO
                      rw2=ZERO
                    end if
!
                    if(in2.ne.0) then
                      ra=ra1+cplkc*cjkc*ra2
                    else
                      ra=ZERO
                    end if
!
!                   --- Nucleus 2 hyperfine angular factor
                    if(in2b.gt.0) then
                      if((kdabs.le.2).and.fld) then
                        z=w3j(lr,2,lc,kr,-kd,-kc)
                        if (jkd.eq.0) then
                          ra2b1=fad2(1,kd+3)*z
                        else
                          ra2b1=fad2(2,kd+3)*z*cjkr
                        end if
                      else
                        ra2b1=ZERO
                      end if
                      if((ksabs.le.2).and.fld) then
                        z=w3j(lr,2,lc,kr,-ks,kc)
                        if (jkd.eq.0) then
                          ra2b2=fad2(1,ks+3)*z
                        else
                          ra2b2=fad2(2,ks+3)*z*cjkr
                        end if
                      else
                        ra2b2=ZERO
                      end if
                      ra2b=ra2b1+cplkc*cjkc*ra2b2
                    else
                      ra2b=ZERO
                    end if
!
                    rg=rg1+cplkc*cjkc*rg2
                    rw=rw1+cplkc*cjkc*rw2
!
!----------------------------------------------------------------------
!     find limits for M index for columns
!----------------------------------------------------------------------
!
                    mcmx=min(lc,mmx)
                    if (nrow.eq.ncol) then
                      mcmn=mr
                    else
                      mcmn=max(-lc,mmn)
                    endif
!
!----------------------------------------------------------------------
!     *** loop over mc ***
!----------------------------------------------------------------------
!
                    do 440 mcsgn=mcmn,mcmx
!
!                  --- Find symmetry index jM from sign of M
!                      (or pI if M is 0)
!
                      mc=abs(mcsgn)
                      jmc=isign(1,mcsgn)
!
!                    --- Find M differences, sums, and parities
!
                      md=mr-mc
                      ms=mr+mc
                      mdabs=abs(md)
                      msabs=abs(ms)
                      cplmc=ipar(lc+mc)
!
!                  --------------------------------------
!                  set flags and zero out matrix elements
!                  if no non-zero elements are possible
!                  --------------------------------------
                      flk0 = (ld.eq.0).and.(kd.eq.0)
                      if(fld.or.(md.eq.0)) then
                        fldmd=.true.
                      else
                        cliou=ZERO
                        cgam=ZERO
                        fldmd=.false.
                      end if
!
!                    -------------------------------------------------
!                    get appropriate 3-J symbols if non-zero Liouville
!                    operator matrix elements are possible
!                    -------------------------------------------------
                      if(fld) then
                        wliou1=w3j(lr,2,lc,mr,-md,-mc)
                        wliou2=w3j(lr,2,lc,mr,-ms,mc)
                      else
                        wliou1=ZERO
                        wliou2=ZERO
                      end if
!
!----------------------------------------------------------------------
!     Calculate off-diagonal terms of the diffusion operator, including
!     the potential-dependent portion and the rhombic component of the
!     isotropic diffusion operator.
!
!     See Meirovitch,Igner,Igner,Moro & Freed, J. Chem. Phys. 77 (1982)
!     pp.3933-3934 especially Eqns. A22-24 and A40 for more information
!----------------------------------------------------------------------
!
!                 --- Rhombic part of isotropic diffusion operator
!
                      if (ld.eq.0 .and. md.eq.0) then
                         cdff=cdffr
                      else
                         cdff=ZERO
                      end if
!
!                 --- Potential-dependent part of diffusion operator
! 
                      if((ipt.ne.0).and.(ldabs.le.lband).and.(md.eq.0)
     #                  .and.((kdabs.le.kband).or.(ksabs.le.kband)).and.
     #                  (ipar(ks).eq.1).and.(jkd.eq.0)) then
!
                      kdi=1+kdabs/2
                      ksi=1+ksabs/2
                      c1=cpmkr*cnl*cnk 
!
                      do l=0,lband,2
                        cdffu=ZERO
                        li=1+l/2
                        if (ksi.le.li .and. abs(xlk(li,ksi)).ge.RNDOFF)
     #                   cdffu=cjkc*cplkc*xlk(li,ksi)
     #                   *w3j(lr,l,lc,kr,-ks,kc)
!
                        if (kdi.le.li .and. abs(xlk(li,kdi)).ge.RNDOFF)
     #                   cdffu=cdffu+xlk(li,kdi)*w3j(lr,l,lc,kr,-kd,-kc)
!
                        if (abs(cdffu).gt.RNDOFF) cdff=cdff+
     #                     w3j(lr,l,lc,mr,0,-mr)*c1*cdffu
                       end do
!
                     end if
!
!----------------------------------------------------------------------
!     Find limits for pI index for columns
!----------------------------------------------------------------------
!
                      if(ipsi0.eq.0) then
                        ipncmn=mc
                        ipncmx=mc
                      else
                        ipncmx=min(ipnmx,in2)
!
                        if (mc.eq.0.and.jmmn.eq.1) then
                          ipncmn=0
                        else
                          ipncmn=-ipncmx
                        end if
                      end if
!
                      if (nrow.eq.ncol) ipncmn=ipnr
!
!----------------------------------------------------------------------
!     *** loop over ipnc ***
!----------------------------------------------------------------------
!
                      do 450 ipncsg=ipncmn,ipncmx
!
                        if(mc.eq.0) then
                           ipnc=abs(ipncsg)
                           jmc=isign(1,ipncsg)
                           if (ipnc.eq.0) jmc=iparlc
                        else
                           ipnc=ipncsg
                        end if
                        if (jmc.lt.jmmn) goto 450
!
!                     ---  Find pI sums and differences
!
                        ipnd=ipnr-ipnc
                        ipns=ipnr+ipnc
                        ipndab=abs(ipnd)
                        ipnsab=abs(ipns)
!
!                      ---find scaling factors
!
                        if((ipnc.eq.0).and.(mc.eq.0)) then
                          cnpc=dsq2
                        else
                          cnpc=ONE
                        end if
!
                        cjmc=jmc
                        jmd=jmr-jmc
                        fdjkm=(jkd.eq.0).and.(jmd.eq.0)
                        cnp=ONE/(cnpr*cnpc)
                        cnorm=cnl*cnk*cnp*cpmkr
!
!----------------------------------------------------------------------
!     get director tilt dependence of Liouville
!     operator matrix elements
!----------------------------------------------------------------------
!
                        if((ipndab.le.2).and.(mdabs.le.2)) then
                          d1=d2km(1,ipnd+3,md+3)*wliou1
                        else
                          d1=ZERO
                        end if
!
                        if((ipnsab.le.2).and.(msabs.le.2)) then
                          d2=d2km(1,ipns+3,ms+3)*wliou2
                        else
                          d2=ZERO
                        end if
!
!----------------------------------------------------------------------
!      find limits to qI index for columns
!----------------------------------------------------------------------
!
                        iqncmx=in2-abs(ipnc)
                        iqncmn=-iqncmx
!
!
!----------------------------------------------------------------------
!      *** loop over iqnc ***
!----------------------------------------------------------------------
!
                        if (nrow.eq.ncol) iqncmn=iqnr
                        do 460 iqnc=iqncmn,iqncmx,2
!
!----------------------------------------------------------------------
!     *** loops over ipn2c, iqn2c (nucleus 2 column) ***
!----------------------------------------------------------------------
                        if (in2b.gt.0) then
                          ipn2cmn=-in2b
                          ipn2cmx=in2b
                        else
                          ipn2cmn=0
                          ipn2cmx=0
                        end if
                        if (nrow.eq.ncol) ipn2cmn=ipn2r
                        do 462 ipn2c=ipn2cmn,ipn2cmx
                        if (in2b.gt.0) then
                          cnp2c=ONE
                          if (ipn2c.eq.0) cnp2c=dsq2
                          iqn2cmn=-(in2b-iabs(ipn2c))
                          iqn2cmx=in2b-iabs(ipn2c)
                        else
                          cnp2c=ONE
                          iqn2cmn=0
                          iqn2cmx=0
                        end if
                        if (nrow.eq.ncol .and. ipn2c.eq.ipn2r)
     #                     iqn2cmn=iqn2r
                        do 464 iqn2c=iqn2cmn,iqn2cmx,2
                        if (in2b.le.0.and.iqn2c.ne.0) goto 470
!
!     --- Compute nucleus 2 index differences
                        ipn2d=ipn2r-ipn2c
                        ipn2s=ipn2r+ipn2c
                        ipn2dab=iabs(ipn2d)
                        ipn2sab=iabs(ipn2s)
                        iqn2d=iqn2r-iqn2c
                        iqn2s=iqn2r+iqn2c
                        iqn2dab=iabs(iqn2d)
                        fdpqi2=(ipn2d.eq.0).and.(iqn2d.eq.0)
!
                        cnp2=ONE/(cnp2r*cnp2c)
!
!     --- Director tilt factors for nucleus 2
                        if (in2b.gt.0) then
                          if((ipn2dab.le.2).and.(mdabs.le.2)) then
                            d1b=d2km(1,ipn2d+3,md+3)*wliou1
                          else
                            d1b=ZERO
                          end if
                          if((ipn2sab.le.2).and.(msabs.le.2)) then
                            d2b=d2km(1,ipn2s+3,ms+3)*wliou2
                          else
                            d2b=ZERO
                          end if
                        else
                          d1b=ZERO
                          d2b=ZERO
                        end if
!
!----------------------------------------------------------------------
!     skip calculation if non-zero elements not possible
!----------------------------------------------------------------------
!
                          if(.not.fldmd .and. .not.fdpqi2) go to 470
                          if(.not.fldmd) then
                            if (in2b.gt.0 .and. .not.fdpqi2) then
                              go to 470
                            end if
                          end if
!
!----------------------------------------------------------------------
!      find sums and differences for nucleus 1
!----------------------------------------------------------------------
!
                          iqnd=iqnr-iqnc
                          iqns=iqnr+iqnc
                          iqndab=abs(iqnd)
                          fdpqi=(ipnd.eq.0).and.(iqnd.eq.0)
!
!----------------------------------------------------------------------
!     Liouville operator
!----------------------------------------------------------------------
!
                          cliou=ZERO
                          cgamw=ZERO
!
                          if(fld) then
!
                            if (jmd.eq.0) then
!
!----------------------------------------------------------------------
!     Nucleus 1 hyperfine: requires delta(pI2,qI2)
!----------------------------------------------------------------------
                              if (fdpqi2) then
!
                              if((ipndab.eq.iqndab).and.
     #                           (ipndab.le.1).and.(in2.ne.0)) then
!
                                if(ipnd.eq.0) then
                                  sa1=iqnr/dsq6
                                  if(flk0.and.md.eq.0.and.jkd.eq.0) then
                                    clioi1=-sa1*ra0
                                  else
                                    clioi1=ZERO
                                  end if
!
                                else
                                  kip=iqnr*iqnd+ipnr*ipnd
                                  kip=kip*(kip-2)
                                  fki=dsqrt(fii-0.25D0*kip)
                                  sa1=-ipnd*fki*0.25D0
                                  clioi1=ZERO
                                end if
!
                              else
                                sa1=ZERO
                                clioi1=ZERO
                              end if
!
                              if((ipnsab.eq.iqndab).and.
     #                           (ipnsab.le.1).and.(in2.ne.0)) then
!
                                if(ipns.eq.0) then
                                  sa2=iqnr/dsq6
                                  if(flk0.and.ms.eq.0.and.jkd.eq.0) then
                                    clioi2=-sa2*ra0
                                  else
                                    clioi2=ZERO
                                  end if
                                else
                                  kip=iqnr*iqnd+ipnr*ipns
                                  kip=kip*(kip-2)
                                  fki=dsqrt(fii-0.25D0*kip)
                                  sa2=-ipns*fki*0.25D0
                                  clioi2=ZERO
                                end if
!
                              else
                                sa2=ZERO
                                clioi2=ZERO
                              end if
!
                              clioua=(sa1*d1+cjmc*cplmc*sa2*d2)
                              cliou=cnorm*(clioua*ra)+
     #                               cnp*(clioi1+clioi2)
!
                              end if
!
!----------------------------------------------------------------------
!     Nucleus 2 hyperfine: requires delta(pI1,qI1)
!----------------------------------------------------------------------
                              if (fdpqi .and. in2b.gt.0) then
!
                              if((ipn2dab.eq.iqn2dab).and.
     #                           (ipn2dab.le.1)) then
!
                                if(ipn2d.eq.0) then
                                  sa2b1=iqn2r/dsq6
                                  if(flk0.and.md.eq.0.and.jkd.eq.0)then
                                    clioi2b1=-sa2b1*ra02
                                  else
                                    clioi2b1=ZERO
                                  end if
!
                                else
                                  kip2=iqn2r*iqn2d+ipn2r*ipn2d
                                  kip2=kip2*(kip2-2)
                                  fki2=dsqrt(fii2-0.25D0*kip2)
                                  sa2b1=-ipn2d*fki2*0.25D0
                                  clioi2b1=ZERO
                                end if
!
                              else
                                sa2b1=ZERO
                                clioi2b1=ZERO
                              end if
!
                              if((ipn2sab.eq.iqn2dab).and.
     #                           (ipn2sab.le.1)) then
!
                                if(ipn2s.eq.0) then
                                  sa2b2=iqn2r/dsq6
                                  if(flk0.and.ms.eq.0.and.jkd.eq.0)then
                                    clioi2b2=-sa2b2*ra02
                                  else
                                    clioi2b2=ZERO
                                  end if
                                else
                                  kip2=iqn2r*iqn2d+ipn2r*ipn2s
                                  kip2=kip2*(kip2-2)
                                  fki2=dsqrt(fii2-0.25D0*kip2)
                                  sa2b2=-ipn2s*fki2*0.25D0
                                  clioi2b2=ZERO
                                end if
!
                              else
                                sa2b2=ZERO
                                clioi2b2=ZERO
                              end if
!
                              clioua=(sa2b1*d1b+cjmc*cplmc*sa2b2*d2b)
                              cliou=cliou+cnorm*(clioua*ra2b)+
     #                               cnp2*(clioi2b1+clioi2b2)
!
                              end if
!
!                              -----------------------
!                              Electronic Zeeman term
!                              (diagonal in ALL nuclear indices)
!                              -----------------------
                              if(fdpqi.and.fdpqi2.and.
     #                           (iqnd.eq.0).and.
     #                           (abs(rg).gt.RNDOFF)) then
!
                                if(ipnd.eq.0) then
                                  sg1=dsq23
                                else
                                  sg1=ZERO
                                end if
!
                                if(ipns.eq.0) then
                                  sg2=dsq23
                                else
                                  sg2=ZERO
                                end if
!
                                clioug=(sg1*d1+cjmc*cplmc*sg2*d2)
                                cliou=cliou+cnorm*(clioug*rg)
!
                              end if
!
!                             --------------------------------------
!                             Orientation-dependent linewidth term
!                             (diagonal in all nuclear indices)
!                             --------------------------------------
                              if(fdpqi.and.fdpqi2.and.
     #                           (iqnd.eq.0).and.(abs(rw).gt.RNDOFF))
     #                        then
                                if(ipnd.eq.0) then
                                  sw1=dsq23
                                else
                                  sw1=ZERO
                                end if
!
                                if(ipns.eq.0) then
                                  sw2=dsq23
                                else
                                  sw2=ZERO
                                end if
                                cgamw=cnorm*(sw1*d1+cjmc*cplmc*sw2*d2)
     #                                     *rw
                              end if
!
!                                if (jmd.eq.0) then...
                           else
!
!                          --------------------
!                          Nuclear Zeeman terms
!                          (diagonal in all indices)
!                          --------------------
                             if(flk0.and.(md.eq.0).and.(jkd.eq.0)
     #                         .and.fdpqi.and.fdpqi2) then
                               cliou=ipnr*zeen+ipn2r*zeen2
                             end if
!
                           end if
!
!                               if (fld) then...
                         else
                           cliou=ZERO
                           cgamw=ZERO
                         end if
!
!----------------------------------------------------------------------
!                      Heisenberg exchange operator
!----------------------------------------------------------------------
!
                         if((abs(ossc).gt.RNDOFF).and.(ipnd.eq.0).and.
     #                     (md.eq.0).and.flk0.and.fdjkm.and.
     #                     fdpqi2) then
!
                            ctemp=ZERO
                            if(iqnd.eq.0) ctemp=ONE
                            if ((ipnr.eq.0).and.(ipn2r.eq.0)
     #                          .and.(lr.eq.0))
     #                       ctemp=ctemp-ONE/dble((in2+1)*(max(in2b,1)))
                            cgam=cgamw+ctemp*ossc
!
                         else
                            cgam=cgamw
                         end if
!
!                        ----------------------------------------
!                         Add in diffusion terms on the diagonal
!                        ----------------------------------------
                         if(fdpqi.and.fdpqi2.and.fdjkm) cgam=cgam+cdff
!
!----------------------------------------------------------------------
!  Store matrix elements
!----------------------------------------------------------------------
!
!                         --- Diagonal elements
!
                          if (nrow.eq.ncol) then
                             cgam=cgam+cdfd
                             zdiag(1,nrow)=cgam
                             zdiag(2,nrow)=-cliou
!
!                         --- Off-diagonal elements
!
                          else
                             if (abs(cliou).gt.RNDOFF) then
                                nel=nel+1
                                if (nel.gt.MXEL) then
                                   ierr=MTXBIG
                                   go to 9999
                                else
                                   neli=neli+1
                                   zmat(neli)=-cliou
                                   izmat(neli)=ncol
                                end if
                             end if
!
                             if (abs(cgam).gt.RNDOFF) then
                                nel=nel+1
                                if (nel.gt.mxel) then
                                   ierr=MTXBIG
                                   go to 9999
                                else
                                   nelr=nelr+1
                                   zmat(mxel-nelr+1)=cgam
                                   izmat(mxel-nelr+1)=ncol
                                end if
                             end if
!
                          end if
!
!                        --------------------------
!                         increment column counter
!                        --------------------------
!
 470                      ncol=ncol+1
!
                          if (ncol.gt.MXDIM) then
                            ierr=DIMBIG
                            go to 9999
                          end if
!
!----------------------------------------------------------------------
!       end loop over columns
!----------------------------------------------------------------------
!
 464                    continue
 462                    continue
 460                    continue
 450                  continue
 440                continue
 420              continue
 400            continue
!
!----------------------------------------------------------------------
!     end loop over rows
!----------------------------------------------------------------------
!
 364          continue
 362          continue
 360          continue
 350        continue
 340      continue
 320    continue
 300  continue
!
      ierr=0
      ndim=nrow
      neltot=nel
      nelre=nelr
      nelim=neli
!
      jzmat(ndim+1)=neli+1
      kzmat(ndim+1)=nelr+1
!       
!.......... Debugging purposes only!
!      close (20)
!...................................
 9999 return
      end
