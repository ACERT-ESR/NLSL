c Version 1.6  10/25/94
c**********************************************************************
c                        =========================
c                           subroutine MATRLL
c                        =========================
c                      *****  Pruning version *****
c
c       Subroutine for calculating matrix elements for EPRLL calculation.
c       Uses expressions given in Meirovitch, Igner, Igner, Moro, & Freed
c       J. Phys. Chem. 77 (1982) 
c
c       Updated to improve matrix storage algorithm. This version differs
c       from the published version (Biol. Magnetic Resonance, Vol. 8)
c       in that it only stores half of the matrix. Changes are noted
c       in the code given below.
c
c       Further updated to provide for a nonaxial diffusion tensor R
c       (only in the case of Brownian motion) and for an inhomogeneous
c       linewidth tensor. (DEB)
c
c       Also extended to permit antisymmetric linear cominations of
c       basis functions with a given K or M (original EPRLL basis was
c       constructed using only symmetric K and M combinations). This
c       provides additional flexibility in diffusion and magnetic tilt
c       angles, and allows the inclusion of the nuclear Zeeman interaction.
c       (DEB)
c
c       flags
c       -----
c
c             fld :   True if |lr-lc| <=2 , false otherwise.
c                     Used as flag for bandwidth of Liouville operator
c                     matrix elements. The restriction on the 
c                     L-bandwidth due the inclusion of spin operators 
c                     of ranks 0 and 2.
c
c             fldmd : True if |lr-lc| <=2 or mr=mc, false otherwise.
c                     Used as flag for possible non-zero matrix 
c                     elements. The L-bandwidth restriction is 
c                     explained above and the diffusion operator 
c                     matrix elements can be found only when mr=mc.
c
c             flk0 :  True if lr=lc and kr=kc, false otherwise.
c                     Used as a flag for determing the existence of 
c                     non-zero Liouville (A-tensor) and diffusion 
c                     (Heisenberg spin exchange) matrix elements.
c
c       ----------------------------------------
c       written by DJS 29-AUG-87
c
c       modified by JPG 06-JUN-88 to correct the anisotropic 
c                                 viscosity terms in the diffusion 
c                                 operator according to Polnaszek
c                                 and Freed, J. Chem. Phys. (1975) 
c
c       modified by DJS JUN-90 to update matrix storage. Comments by DEB.
c
c       modified by DEB: calculation of isotropic part of diffusion
c                        operator moved outside iqnr loop. Fully anisotropic
c                        rotational diffusion tensor implemented.
c                        Antisymmetric K- and M- combinations added to
c                        accommodate extra tilt angles and nuclear Zeeman
c
c       Uses:
c               cd2km.f
c               anxlk.f
c               ipar.f
c               w3j.f
c
c**********************************************************************
c
      subroutine matrll(ierr)
c
      implicit none
c
      include 'stddim.inc'
      include 'rndoff.inc'
      include 'eprmat.inc'
      include 'eprdat.inc'
      include 'indexl.inc'
      include 'physcn.inc'
      include 'maxl.inc'
c
      integer i,l,kdi,ksi,li,ierr
      integer lr,lrprod,kr,krmx,mr,mrmx,
     #     ipnr,ipnrmn,ipnrmx,iqnr,iqnrmn,iqnrmx,
     #     lc,lcmax,ld,kc,kcmn,kcmx,kd,ks,
     #     kdabs,ksabs,mc,mcmn,mcmx,md,ms,mdabs,msabs,
     #     ipnc,ipncmn,ipncmx,iqnc,iqncmn,iqncmx,
     #     ipnd,ipns,ipndab,ipnsab,iqnd,iqns,iqndab,
     #     kip,nrow,ncol,nelr,neli,nel
      integer ioldlr,ioldkr,ioldmr,ioldpnr,ioldlc,ioldkc,ioldmc,
     #        ioldpnc,iparlr,iparlc,jkc,jkr,jkd,jmc,jmr,jmd,
     #        krsgn,mrsgn,kcsgn,mcsgn,ipnrsg,ipncsg
c
      double precision dsq2,dsq6,dsq15,dsq23
      double precision ct,rpl,rmi,rz,rj,rp,ossc,fii,
     #     cpmkr,cnpr,c1,c2,c3,c4,pmz,pmxk,pmxl,
     #     cdfd,cnl,cplkc,cnk,z,ra0,ra,rg,rw,ra1,rg1,rw1,
     #     ra2,rg2,rw2,cliou,cgam,cplmc,wliou1,wliou2,cdff,cdffr,
     #     cdffu,cnpc,cnp,cnorm,d1,d2,sa1,clioi1,fki,cgamw,cjmc,cjkr,
     #     cjkc,sa2,clioi2,clioua,sg1,sg2,sw1,sw2,clioug,ctemp,zeen
      double precision d2km(2,5,5)
c
      logical flk0,fld,fldmd,fdpqi,fdjkm
      logical newlr,newkr,newmr,newpnr,newlc,newkc,newmc,newpnc
c
      double precision ZERO,ONE
      parameter (ZERO=0.0D0,ONE=1.0D0)
c
      integer ipar 
      double precision w3j
      external ipar,w3j
c
c######################################################################
c.... Debugging purposes only!
      open (unit=20,file='lpmat.tst',status='unknown',
     #     access='sequential',form='formatted')
c............................
c
      if (ndim.gt.MXDIM) then
         ierr=2
         go to 9999
      end if
c
c----------------------------------------------------------------------
c     Scale diffusion constants
c----------------------------------------------------------------------
c
      ct=g0*BETAE/HBAR
      rpl=0.5d0*(dx+dy)/ct
      rmi=0.25*(dx-dy)/ct
      rz=dz/ct
      rj=djf/ct
      rp=djfprp/ct
      ossc=oss/ct
c
c----------------------------------------------------------------------
c     Define constants
c----------------------------------------------------------------------
c
      dsq2=dsqrt(2.0D0)
      dsq6=dsqrt(6.0D0)
      dsq15=dsqrt(1.5D0)
      dsq23=dsqrt(2.0D0/3.0D0)
c
      fii=0.25D0*dble(in2*(in2+2))
c                                        *** Nuclear Zeeman ***
      zeen=b0*gamman/GAMMAE
c
      if(ipdf.ne.0) then
         pmxl=dx*tl
         pmz=dz*tkzz
         pmxk=dx*tkxy
      end if
c
      if ((ipdf.eq.2) .and. (ist.eq.0)) then
         rpl=rpl+rp
         rz=rz+rp
      end if
c
      if(in2.eq.0) then
         ra0=ZERO
      else
         ra0=-dsq15*a0
      end if
c
c----------------------------------------------------------------------
c                2
c     Calculate d   (psi) 
c                k,m
c                             2
c     Note : d2km(1,i+3,j+3)=d   (psi)
c                             i,j
c
c     Result should be real-valued; i.e. d2km(2,i+3,j+3)=0
c----------------------------------------------------------------------
c
      call cd2km(d2km,ZERO,psi,ZERO)
c
c----------------------------------------------------------------------
c                               L
c     Calculate the quantities X  used in potential-dependent
c                               K
c     portion of diffusion operator
c----------------------------------------------------------------------
c
      call anxlk(rpl,rmi,rz)
c
c----------------------------------------------------------------------
c     Initialize counters
c----------------------------------------------------------------------
c
      nelr=0
      neli=0
      nel=0
      ioldlr=-1
      ioldkr=kmn-1
      ioldmr=mmn-1
      ioldpnr=-in2-1
c
c----------------------------------------------------------------------
c     **** Loop over rows ***
c----------------------------------------------------------------------
c
      do 400 nrow=1,ndim
         lr=l1(nrow)
         krsgn=k1(nrow)
         mrsgn=m1(nrow)
         ipnrsg=pi1(nrow)
         iqnr=qi1(nrow)
c
         newlr=lr.ne.ioldlr
         newkr=newlr.or.(krsgn.ne.ioldkr)
         newmr=newkr.or.(mrsgn.ne.ioldmr)
         newpnr=newmr.or.(ipnrsg.ne.ioldpnr)

         if (newlr) then
            ioldlr=lr
            lrprod=lr*(lr+1)
            iparlr=ipar(lr)
         end if
c
         if (newkr) then
            ioldkr=krsgn
            kr=abs(krsgn)
            jkr=isign(1,krsgn)
            if (krsgn.eq.0) jkr=iparlr
c            cplkr=ipar(lr+kr)
c            cpkr=ipar(kr)
            cjkr=jkr
         end if
c
         if (newmr) then
            ioldmr=mrsgn
	    jmr=isign(1,mrsgn)
            mr=abs(mrsgn)
            cpmkr=ipar(mr+kr)
c            cplmr=ipar(lr+mr)
c            cpmr=ipar(mr)
c
c----------------------------------------------------------------------
c     Calculate isotropic part of diffusion operator
c----------------------------------------------------------------------
c
            c1=dble(lrprod)
            if(ipdf.ne.0) then
               c2=ONE+pmxl*c1
               if(abs(pl).gt.RNDOFF) c2=c2**pl
               c1=c1/c2
            end if
c
            cdfd=rpl*c1
            c1=dble(kr*kr)
	    c2=rz
	    c3=rpl
c
            if(ipdf.ne.0) then
               c4=ONE+pmz*c1
               if(abs(pkzz).gt.RNDOFF) c4=c4**pkzz
               c2=c2/c4
               c4=ONE+pmxk*c1
               if(abs(pkxy).gt.RNDOFF) c4=c4**pkxy
               c3=c3/c4
            end if
c
            cdfd=cdfd+c1*(c2-c3)
            if(ipdf.eq.2) cdfd=cdfd+dble(mr*mr)*(rj-rp)
c
c    Discrete jump motion
c
            if(ist.ne.0) then
               i=kr/ist
               if((i*ist).ne.kr) cdfd=cdfd+rj
            end if
         end if
c
         if (newpnr) then
            ioldpnr=ipnrsg
            if (mr.eq.0) then
               jmr=isign(1,ipnrsg)
               ipnr=abs(ipnrsg)
               if (ipnr.eq.0) jmr=iparlr
            else
               ipnr=ipnrsg
            end if
c              
            if((ipnr.eq.0).and.(mr.eq.0)) then
               cnpr=dsq2
            else
               cnpr=ONE
            end if
         end if
c
         jzmat(nrow)=neli+1
         kzmat(nrow)=nelr+1
c
c----------------------------------------------------------------------
c   **** Loop over columns (note only upper diagonal) ****
c----------------------------------------------------------------------
         ioldlc=-1
         ioldkc=kmn-1
         ioldmc=mmn-1
         ioldpnc=-in2-1
         lcmax = min0( lr+lband, lemx )

         do 300 ncol=nrow,ndim
            lc=l1(ncol)
            kcsgn=k1(ncol)
            mcsgn=m1(ncol)
            ipncsg=pi1(ncol)
            iqnc=qi1(ncol)
c
            newlc=lc.ne.ioldlc
            newkc=newlc.or.(kcsgn.ne.ioldkc)
            newmc=newkc.or.(mcsgn.ne.ioldmc)
            newpnc=newmc.or.(ipncsg.ne.ioldpnc)
c
            if (newlc) then
               ioldlc=lc
               if (lc.gt.lcmax) goto 350
	       iparlc=ipar(lc)
c
c----------------------------------------------------------------------
c     find limits for non-zero Liouville operator
c     matrix elements and scaling factor
c----------------------------------------------------------------------
c
               ld=lr-lc
c               ls=lr+lc
               fld=abs(ld).le.2
               cnl=dsqrt((2.D0*lr+1.D0)*(2.D0*lc+1.D0))
            end if
c
c----------------------------------------------------------------------
c      find sums, differences, parities and scaling factor
c----------------------------------------------------------------------
c
            if (newkc) then
               ioldkc=kcsgn
               kc=abs(kcsgn)
               jkc=isign(1,kcsgn)
               if (kcsgn.eq.0) jkc=iparlc
c
               jkd=jkr-jkc
               kd=kr-kc
               ks=kr+kc
               kdabs=abs(kd)
               ksabs=abs(ks)
               cplkc=ipar(lc+kc)
               cjkc=jkc
c
c                  -- Find scaling factor
c
               if ((kr.eq.0).and.(kc.eq.0)) then
                  cnk=0.5D0
               else if ((kr.ne.0).and.(kc.ne.0)) then
                  cnk=ONE
               else
                  cnk=ONE/dsq2
               end if
c
c----------------------------------------------------------------------
c              Calculate rhombic part of diffusion tensor
c----------------------------------------------------------------------
c
               if (kr+2.eq.kc .and. jkd.eq.0) then
                  cdffr=rmi*sqrt(dble((lr+kc-1)*(lr+kc)*
     #                 (lr-kc+1)*(lr-kc+2) ) )/cnk
               else if (kr-2.eq.kc .and. jkd.eq.0) then
                  cdffr=rmi*sqrt(dble((lr-kc-1)*(lr-kc)*
     #                 (lr+kc+1)*(lr+kc+2) ) )/cnk
               else
                  cdffr=ZERO
               end if

c---------------------------------------------------------------------
c       Calculate R(mu,l) as in Meirovitch, Igner, Igner, Moro, & Freed
c       J. Phys. Chem. 77 (1982) p. 3935, Eqs. A42 & A44
c       mu = g, A, W
c---------------------------------------------------------------------
c
               if ((kdabs.le.2).and.fld) then
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
c
               if ((ksabs.le.2).and.fld) then
                  z=w3j(lr,2,lc,kr,-ks,kc)
                  if (jkd.eq.0) then
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
c
               if (in2.ne.0) then
                  ra=ra1+cplkc*cjkc*ra2
               else
                  ra=ZERO
               end if
c
               rg=rg1+cplkc*cjkc*rg2
               rw=rw1+cplkc*cjkc*rw2
            end if
c
            if (newmc) then
               ioldmc=mcsgn
c
c            --- Find symmetry index jM from sign of M
c                (or pI if M is 0)
c
               mc=abs(mcsgn)
	       jmc=isign(1,mcsgn)
c
c              ----------------------------------------
c              Find M differences, sums, and parities
c              ----------------------------------------
               md=mr-mc
               ms=mr+mc
               mdabs=abs(md)
               msabs=abs(ms)
               cplmc=ipar(lc+mc)
c
c             ----------------------------------------
c             Set flags and zero out matrix elements
c             if no non-zero elements are possible
c             ----------------------------------------
               flk0=(ld.eq.0).and.(kd.eq.0)
               if(fld.or.(md.eq.0)) then
                  fldmd=.true.
               else
                  cliou=ZERO
                  cgam=ZERO
                  fldmd=.false.
               end if
c
c             -------------------------------------------------
c             get appropriate 3-J symbols if non-zero Liouville
c             operator matrix elements are possible
c             -------------------------------------------------
               if (fld) then
                  wliou1=w3j(lr,2,lc,mr,-md,-mc)
                  wliou2=w3j(lr,2,lc,mr,-ms,mc)
               else
                  wliou1=ZERO
                  wliou2=ZERO
               end if
c
c----------------------------------------------------------------------
c     Calculate off-diagonal terms of the diffusion operator, including
c     the potential-dependent portion and the rhombic component of the
c     isotropic diffusion operator.
c
c     See Meirovitch,Igner,Igner,Moro & Freed, J. Chem. Phys. 77 (1982)
c     pp.3933-3934 especially Eqns. A22-24 and A40 for more information
c----------------------------------------------------------------------
c
c              --- Rhombic part of isotropic diffusion operator
c
               if (ld.eq.0 .and. md.eq.0) then
                  cdff=cdffr
               else
                  cdff=ZERO
               end if
c
c              --- Potential-dependent part of diffusion operator
c 
               if((md.eq.0).and.(ipt.ne.0).and.(ipar(ks).eq.1)
     #              .and. ((kdabs.le.kband).or.(ksabs.le.kband)) 
     #              .and. jkd.eq.0) then
c
                  kdi=1+kdabs/2
                  ksi=1+ksabs/2
                  c1=cpmkr*cnl*cnk 
c     
c                 --- Sum over terms in potential
c
                  do l=0,lband,2
                     cdffu=ZERO
                     li=1+l/2
c
c QUESTION: Do the following statements give bogus results if k{s,d}i.gt.li?
c           (can there be a bad reference in xlk?)
c QUESTION: Why does the following term not have a factor of jK in it?
c
                     if (ksi.le.li.and. abs(xlk(li,ksi)).ge.RNDOFF)
     #                   cdffu=cplkc*xlk(li,ksi)*w3j(lr,l,lc,kr,-ks,kc)
c
                     if (kdi.le.li.and. abs(xlk(li,kdi)).ge.RNDOFF)
     #                  cdffu=cdffu+xlk(li,kdi)*w3j(lr,l,lc,kr,-kd,-kc)
c
                     if (abs(cdffu).gt.RNDOFF) then
                        cdff=cdff+w3j(lr,l,lc,mr,0,-mr)*c1*cdffu
                     end if

c
                  end do
c
               end if
c
            end if
c
c          ---------------------------------------------
c          Find sums and differences of pI spin indices
c          ---------------------------------------------
            if (newpnc) then
               ioldpnc=ipncsg
	       if (mc.eq.0) then
                  ipnc=abs(ipncsg)
                  jmc=isign(1,ipncsg)
                  if (ipnc.eq.0) jmc=iparlc
               else
                 ipnc=ipncsg
               end if
c
c             --- Find pI sums and differences
c
               ipnd=ipnr-ipnc
               ipns=ipnr+ipnc
               ipndab=abs(ipnd)
               ipnsab=abs(ipns)
               if((ipnc.eq.0).and.(mc.eq.0)) then
                  cnpc=dsq2
               else
                  cnpc=ONE
               end if
c
               cjmc=jmc
               jmd=jmr-jmc
               fdjkm=(jkd.eq.0).and.(jmd.eq.0)
               cnp=ONE/(cnpr*cnpc)
               cnorm=cnl*cnk*cnp*cpmkr
c
c              --------------------------------------------------
c              Get director tilt dependence of Liouville
c              operator matrix elements
c              --------------------------------------------------
               if((ipndab.le.2).and.(mdabs.le.2)) then
                  d1=d2km(1,ipnd+3,md+3)*wliou1
               else
                  d1=ZERO
               end if
c
               if((ipnsab.le.2).and.(msabs.le.2)) then
                  d2=d2km(1,ipns+3,ms+3)*wliou2
               else
                  d2=ZERO
               end if
            end if
c
c----------------------------------------------------------------------
c     Calculate the matrix element if a non-zero element is possible
c----------------------------------------------------------------------
c
            if(fldmd) then

c              --------------------------------------------
c              find sums and differences of qI spin indices
c              --------------------------------------------
               iqnd=iqnr-iqnc
               iqns=iqnr+iqnc
               iqndab=abs(iqnd)
c               iqnsab=abs(iqns)
	       fdpqi=(ipnd.eq.0).and.(iqnd.eq.0)
c
c----------------------------------------------------------------------
c              Liouville operator
c
c              NOTE: see Appendix A p.3936-3937 in M,I,I,M, & F for
c              details.  Sa and Sg in program also contain the
c              appropriate Clebsch-Gordon coefficient. fki is
c              Ki in text. 
c     
c              The isotropic contributions from the A tensor are
c              not described in the text. They appear only when
c              L1 (lr) equals L2 (lc). The Clebsch-Gordon coefficients
c              for these terms (+/- 1/sqrt(3)) are included in the
c              definition of a0. The Wigner 3J symbols appropriate for
c              the isotropic terms exactly cancel the normalization
c              constant cnorm. Instead of calculating the 3J symbols,
c              the cnorm factor is simply omitted from the isotropic
c              terms below.
c----------------------------------------------------------------------
               if(fld) then
c
                  if (jmd.eq.0) then
c
                     if((ipndab.eq.iqndab).and.
     #                  (ipndab.le.1).and.(in2.ne.0)) then
c
c                    ----------------------------
c                     Hyperfine interaction term
c                    ----------------------------
                        if(ipnd.eq.0) then
                           sa1=iqnr/dsq6
                           if(flk0.and.md.eq.0.and.jkd.eq.0) then
                              clioi1=-sa1*ra0
                           else
                              clioi1=ZERO
                           end if
                        else
                           kip=iqnr*iqnd+ipnr*ipnd
                           kip=kip*(kip-2)
                           fki=dsqrt(fii-0.25D0*kip)
                           sa1=-ipnd*fki*0.25D0
                           clioi1=ZERO
                        end if
c
                     else
                        sa1=ZERO
                        clioi1=ZERO
                     end if
c
                     if((ipnsab.eq.iqndab).and.
     #                  (ipnsab.le.1).and.(in2.ne.0)) then
c
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
c
                     else
                        sa2=ZERO
                        clioi2=ZERO
                     end if
c
                     clioua=(sa1*d1+cjmc*cplmc*sa2*d2)
c
c                    ------------------------
c                     Electronic Zeeman term
c                    ------------------------
                     if((iqnd.eq.0).and.(abs(rg).gt.RNDOFF)) then
                        if(ipnd.eq.0) then
                           sg1=dsq23
                        else
                           sg1=ZERO
                        end if
c
                        if(ipns.eq.0) then
                           sg2=dsq23
                        else
                           sg2=ZERO
                        end if
                        clioug=(sg1*d1+cjmc*cplmc*sg2*d2)
c
                     else
                        clioug=ZERO
                     end if
c
c                    --------------------------------------
c                     Orientation-dependent linewidth term
c                    --------------------------------------
                     if((iqnd.eq.0).and.(abs(rw).gt.RNDOFF)) then
                        if(ipnd.eq.0) then
                           sw1=dsq23
                        else
                           sw1=ZERO
                        end if
c
                        if(ipns.eq.0) then
                           sw2=dsq23
                        else
                           sw2=ZERO
                        end if
                        cgamw=cnorm*(sw1*d1+cjmc*cplmc*sw2*d2)*rw
c
                     else
                        cgamw=ZERO
                     end if
c
                     cliou=cnorm*(clioua*ra+clioug*rg)+
     #                     cnp*(clioi1+clioi2)
c
c                           if (jmd.eq.0) then...
                  else
c                     
c                    --------------------
c                    Nuclear Zeeman term
c                    --------------------
                     if (flk0.and.(md.eq.0).and.(jkd.eq.0)
     #                  .and.fdpqi) then
                        cliou=ipnr*zeen
                     else
                        cliou=ZERO
                     endif
                     cgamw=ZERO
c
                  end if
c
c                           if (fld) then...
               else
                  cliou=ZERO
                  cgamw=ZERO
               end if
c
c----------------------------------------------------------------------
c              Heisenberg exchange operator
c              Now also includes orientational averaging via HE
c              in slow-motional region with the delta(L1,0) term.
c              (See Lee,Budil,and Freed, JCP, 1994, Eq. B9)
c----------------------------------------------------------------------
c
               if((abs(ossc).gt.RNDOFF).and.(ipnd.eq.0).and.
     #            (md.eq.0).and.flk0.and.fdjkm) then
c
                  ctemp=ZERO
                  if(iqnd.eq.0) ctemp=ONE
                  if((ipnr.eq.0).and.(lr.eq.0)) 
     #               ctemp=ctemp-ONE/dble(in2+1)
c
                  cgam=cgamw+ctemp*ossc
               else
                  cgam=cgamw
               end if
c
c              ----------------------------------------
c              Add in diffusion terms on the diagonal  
c              ----------------------------------------
               if(fdpqi.and.fdjkm) cgam=cgam+cdff
c
c----------------------------------------------------------------------
c  Store matrix elements
c    Real off-diagonal matrix elements (cgam nonzero) are stored 
c    sequentially starting from the end of the space allocated for zmat 
c    and working down.
c
c    Imaginary off-diagonal matrix elements (cliou nonzero) are stored
c    sequentially starting at the beginning of the zmat array.
c----------------------------------------------------------------------
c
c              ------------------
c              Diagonal elements
c              ------------------
               if (nrow.eq.ncol) then
                  cgam=cgam+cdfd
                  zdiag(1,nrow)=cgam
                  zdiag(2,nrow)=-cliou
               else
c
c              ----------------------
c              Off-diagonal elements
c              ----------------------
                  if (abs(cliou).gt.RNDOFF) then
                     nel=nel+1
                     if (nel.gt.MXEL) then
                        ierr=1
                        go to 9999
                     else
                        neli=neli+1
                        zmat(neli)=-cliou
                        izmat(neli)=ncol
                     end if
                  end if
c
                  if (abs(cgam).gt.RNDOFF) then
                     nel=nel+1
                     if (nel.gt.MXEL) then
                        ierr=1
                        go to 9999
                     else
                        nelr=nelr+1
                        zmat(MXEL-nelr+1)=cgam
                        izmat(MXEL-nelr+1)=ncol
                     end if
                  end if
c
c                            if (nrow.eq.ncol)...else...
               end if
c
c..........Debugging purposes only!
               if (abs(cgam).gt.RNDOFF.or.abs(cliou).gt.RNDOFF)
     #             write (20,7000) lr,krsgn,mrsgn,ipnr,iqnr,
     #                             lc,kcsgn,mcsgn,ipnc,iqnc,cgam,-cliou
 7000 format('<',4(i2,','),i2,'|L|',4(i2,','),i2,'> =',2g14.7)
c...................................
c
c                            if (fldmd) ....
            end if
c
c
c----------------------------------------------------------------------
c       end loop over columns
c----------------------------------------------------------------------
 300     continue
c
c----------------------------------------------------------------------
c        Dummy label for exiting loop over columns
c----------------------------------------------------------------------
c
 350     continue
c
c----------------------------------------------------------------------
c     end loop over rows
c----------------------------------------------------------------------
c
 400  continue
c
c     ---------------------------------------------------------
c     Store statistics, mark end of row index arrays, and exit
c     ---------------------------------------------------------
      ierr=0
      neltot=nel
      nelre=nelr
      nelim=neli
c
      jzmat(ndim+1)=neli+1
      kzmat(ndim+1)=nelr+1
c    
c.......... Debugging purposes only!
      close (20)
c...................................
 9999 return
      end
