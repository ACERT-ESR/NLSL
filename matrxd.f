c  VERSION 1.0  (NLSPMC version)   2/5/99
c----------------------------------------------------------------------
c
c                         ===================
c                          subroutine MATRXD
c                         ===================
c
c *** MATRIX STORAGE SCHEME FOR CG ALGORITHM ***
c
c    This routine calculates the matrix elements diagonal subspace
c (pS=0, qS=+/- 1).  Nonsecular terms in the Hamiltonian is ignored
c due to the high field approximation.  This is done by not allowing
c the coupling between the spaces with different pS indices.  The only
c source of the coupling between diagonal and off-diagonal subspaces is
c pulse propagator.
c
c    qS-symmetrization, equivalent to M-symmetrization in off-diagonal
c space is performed for the basis set.  [Change the sign of qS also
c when the signs of pI and M change].  And M-symmetry is not used since
c it would not make the matrix block-diagonal.  In the absence of
c nuclear Zeeman, the block with jqS=1 is completely decoupled from the
c block with jqS=-1.  The pulse propagator does not couple different
c jqS and and jM.  In the presence of nuclear Zeeman, the block with
c jqS=-1 is included in the calculation.
c
c    This routine is based upon EPRLL routines ver.1.3., whose features are
c      1. improved matrix storage algorithm (store only half of the matrix)
c      2. support nonaxial diffusion tensor R (in case of Brownian motion)
c      3. uses pruning scheme by eprbl and tnll.  eprbl and tnll generates
c         an <fileid>.ind file which keeps only the basis elements with
c         significant projection on the starting vector.
c
c    The phenomenological relaxation processes are also included as
c    T1e, T2n, T1n.
c
c      T1e  : electronic T1 relaxation for the allowed ESR transitions
c      T2n  : orientationally invariant nuclear T2 relaxation
c             This defines the homogeneous line-width of NMR transition.
c      T1n  : nuclear T1 relaxation for the allowed NMR transitions
c             This parameter is implicitly calculated since the pseudo-
c             secular term is included in the calculation of the evolution
c             of the diagonal subspace (pS=0)
c   *Note that the sign of these relaxation constants are reversed
c    with respect to the time domain formalism of the relaxation prcesses     
c    since we are calculating the frequency spectrum directly.
c
c      Original program by G. Moro with modifications by D. Budil
c      Further modification by Sanghyuk Lee to calculate the slow-motional
c              2D-FT spectra.
c
c       flags
c       -----
c
c             fld   : True if |lr-lc| <=2.
c                     Used as flag for bandwidth of Liouville operator
c                     matrix elements. The limit on the L-bandwidth is due
c                     to the restriction of spin operator tensors to
c                     ranks 0 and 2.
c
c             fldmd : True if |lr-lc| <=2 or mr=mc, false otherwise.
c                     Used as flag for possible non-zero matrix 
c                     elements. The L-bandwidth restriction is 
c                     explained above and the diffusion operator 
c                     matrix elements can be found only when mr=mc.
c
c             flk0  : True if lr=lc and kr=kc.
c                     Used as a flag for determing the existence of
c                     non-zero Liouville (A-tensor) and diffusion
c                     (Heisenberg spin exchange) matrix elements.
c
c	      fkm0  : True if kr=kc and mr=mc
c	              
c       Includes :
c               nlsdim.inc
c               eprprm.inc
c               indexf.inc
c               eprmat.inc
c               maxl.inc
c               rndoff.inc
c               physcn.inc
c
c       Uses:
c               cd2km.f
c               anxlk.f
c               ipar.f
c               w3j.f
c
c----------------------------------------------------------------------
      subroutine matrxd(ierr)
c
      implicit none
c
      include 'limits.inc'
      include 'simparm.inc'
      include 'basis.inc'
      include 'eprmat.inc'
      include 'maxl.inc'
      include 'rndoff.inc'
      include 'physcn.inc'
c
      logical fnz,fld,flk0,diftrm,sumtrm,fldmd
      logical newlr,newjkr,newkr,newjqer,newmr,newpnr,
     #        newlc,newjkc,newkc,newjqec,newmc,newpnc
      integer ioldlr,ioldjkr,ioldkr,ioldjqer,ioldmr,ioldpnr,
     #        ioldlc,ioldjkc,ioldkc,ioldjqec,ioldmc,ioldpnc
c
      double precision c,ctemp,c1,c2,c3,c4,cdfd,cdfd1,
     #       cdff,cdffu,cga0,cga2,cgam,
     #       cliou,cliou0,clioua,clioug,
     #       cnorm,cnk,cnl,cpmkr,cplkc,cplmc,ct,
     #       d1,dsq2,dsq12,dsq18,dsq23,
     #       fii,fki,ga,gg,ossc,
     #       ra,ra1,ra2,rg,rg1,rg2,rpl,rmi,rz,rj,rp,
     #       sa1,sg1,wliou,z,t1e,t2n,t1n
c
      double precision zero,one
      parameter (zero=0.0d0,one=1.0d0)
c
      double precision d2km(2,5,5)
      integer i,j,l,li,ierr,iper,
     #       jqer,jqec,jqed,ipnr,ipnc,ipnd,ipndab,ipns,
     #       iqnr,iqnc,iqnd,iqndab,
     #       lr,lc,ld,ldabs,lrprod,jkr,jkc,jkd,
     #       kr,kc,kd,ks,kdabs,ksabs,kdi,ksi,kip,
     #       mr,mc,md,ms,mdabs,
     #       nrow,ncol,nel,nelr,neli
c
      integer ipar
      double precision w3j
      external ipar,w3j
c
c######################################################################
c
c      if(abs(zeen).gt.rndoff) fnz=.true.
c
c----------------------------------------------------------------------
c     define constants
c----------------------------------------------------------------------
c
      dsq2=dsqrt(2.0D0)
      dsq12=one/dsq2
      dsq18=one/dsqrt(8.0D0)
      dsq23=dsqrt(2.0D0/3.0D0)
c
c----------------------------------------------------------------------
c     Scale diffusion constants
c----------------------------------------------------------------------
c
      ct=g0*betae/hbar
      rpl=0.5D0*(dx+dy)/ct
      rmi=0.25D0*(dx-dy)/ct
      rz=dz/ct
      rj=djf/ct
      rp=djfprp/ct
      ossc=oss/ct
c        
      t2n=t2ndi/ct
      t1e=t1edi/(2.0D0*ct)
      t1n=t1ndi/(2.0D0*ct)
c
      fii=0.25D0*dble(in2*(in2+2))
c
c    Non-Brownian motional parameters
c
      if(ipdf.ne.1) then
        pml=zero
        pmzz=zero
        pmxy=zero
      end if
c
c   Correction term for anisotropic viscosity
c
      if((ipdf.eq.2).and.(ist.eq.0)) then
        rpl=rpl+rp
        rz=rz+rp
      end if
c
c----------------------------------------------------------------------
c                                 2
c     Get Wigner rotation matrix D  (0,psi,0) for director tilt.
c                                 KM
c
c     Note: resultant matrix will be real-valued and d2km(2,i,j)=zero
c
c                              2
c           d2km(1,i+3,j+3) = d  (psi)
c                              ij
c----------------------------------------------------------------------
c
      call cd2km(d2km,zero,psi,zero)
c
c----------------------------------------------------------------------
c                               L
c     Calculate the quantities X  used in the potential-dependent
c                               K
c     portion of the diffusion operator
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
c
      iper=0
      ioldjqer=0
      ioldlr=-1
      ioldjkr=0
      ioldkr=1000
      ioldmr=1000
      ioldpnr=-in2-1
c
c----------------------------------------------------------------------
c     **** Loop over rows ***
c----------------------------------------------------------------------
c
      do 200 nrow=1,ndimd
        jqer=djqe1(nrow)
        lr=dl1(nrow)
        jkr=djk1(nrow)
        kr=dk1(nrow)
        mr=dm1(nrow)
        ipnr=dpi1(nrow)
        iqnr=dqi1(nrow)
c
        newjqer=jqer.ne.ioldjqer
        newlr=newjqer.or.(lr.ne.ioldlr)
        newjkr=newlr.or.(jkr.ne.ioldjkr)
        newkr=newjkr.or.(kr.ne.ioldkr)
        newmr=newkr.or.(mr.ne.ioldmr)
        newpnr=newmr.or.(ipnr.ne.ioldpnr)
c
        if(newjqer) ioldjqer=jqer
c
        if(newlr) then
          ioldlr=lr
          lrprod=lr*(lr+1)
        end if
c
        if(newjkr) ioldjkr=jkr
c
        if(newkr) then
          ioldkr=kr
c
c         ** Calculate isotropic part of diffusion operator
c
          c1=dble(lrprod)
          if(ipdf.ne.0) then
            c2=one+pml*c1
            if(dabs(pl).gt.rndoff) c2=c2**pl
            c1=c1/c2
          end if
          cdfd=rpl*c1
c
          c1=dble(kr*kr)
          c2=rz
          c3=rpl
          if(ipdf.ne.0) then
            c4=one+pmzz*c1
            if(dabs(pkzz).gt.rndoff) c4=c4**pkzz
            c2=c2/c4
            c4=one+pmxy*c1
            if(dabs(pkxy).gt.rndoff) c4=c4**pkxy
            c3=c3/c4
          end if
          cdfd=cdfd+c1*(c2-c3)
c
c         ** Discrete jump motion
c
          if(ist.ne.0) then
            i=kr/ist
            if((i*ist).ne.kr) cdfd=cdfd+rj
          end if
c
          cdfd1=cdfd
c
        end if
c
        if(newmr) then
          ioldmr=mr
          cpmkr=ipar(mr+kr)
c
c   Anisotropic viscosity term in diffusion operator
c
          if(ipdf.eq.2) cdfd=cdfd1+dble(mr*mr)*(rj-rp)
c
        end if
c
        if(newpnr) ioldpnr=ipnr
c
c----------------------------------------------------------------------
c   **** Loop over columns (note only upper diagonal) ****
c----------------------------------------------------------------------
c
        ioldjqec=0
        ioldlc=-1
        ioldjkc=0
        ioldkc=1000
        ioldmc=1000
        ioldpnc=-in2-1
c
        jzmat(nrow)=neli+1
        kzmat(nrow)=nelr+1
c
        do 300 ncol=nrow,ndimd
          jqec=djqe1(ncol)
          lc=dl1(ncol)
          jkc=djk1(ncol)
          kc=dk1(ncol)
          mc=dm1(ncol)
          ipnc=dpi1(ncol)
          iqnc=dqi1(ncol)
c
          newjqec=jqec.ne.ioldjqec
          newlc=newjqec.or.(lc.ne.ioldlc)
          newjkc=newlc.or.(jkc.ne.ioldjkc)
          newkc=newjkc.or.(kc.ne.ioldkc)
          newmc=newkc.or.(mc.ne.ioldmc)
          newpnc=newmc.or.(ipnc.ne.ioldpnc)
c
          if(newjqec) then
            ioldjqec=jqec
            jqed=jqer-jqec
          end if
c
          if(newlc) then
            ioldlc=lc
            ld=lr-lc
            ldabs=iabs(ld)
            fld=ldabs.le.2
            cnl=dsqrt((2.D0*lr+1.D0)*(2.D0*lc+1.D0))
          end if
c
          if(newjkc) then
            ioldjkc=jkc
            jkd=jkr-jkc
          end if
c
          if(newkc) then
            ioldkc=kc
            kd=kr-kc
            ks=kr+kc
            kdabs=iabs(kd)
            ksabs=iabs(ks)
            cplkc=ipar(lc+kc)
c
            flk0=(ld.eq.0).and.(kd.eq.0)
c
            if((kr.eq.0).and.(kc.eq.0)) then
              cnk=0.5D0
            else if((kr.ne.0).and.(kc.ne.0)) then
              cnk=one
            else
              cnk=dsq12
            end if
c
c---------------------------------------------------------------------
c       Calculate R(mu,l) as in Meirovitch, Igner, Igner, Moro, & Freed
c       J. Phys. Chem. 77 (1982) p. 3935, Eqs. A42 & A44
c---------------------------------------------------------------------
c
            ra=zero
            rg=zero
            if(fld) then
              ra1=zero
              rg1=zero
              if(kdabs.le.2) then
                if(jkd.eq.0) then
                  ga=fad(1,kd+3)
                  gg=fgd(1,kd+3)
                else
                  ga=fad(2,kd+3)*jkr
                  gg=fgd(2,kd+3)*jkr
                end if
                z=w3j(lr,2,lc,kr,-kd,-kc)
                ra1=ga*z
                rg1=gg*z
              end if
c
              ra2=zero

	      rg2=zero
              if(ksabs.le.2) then
                if(jkd.eq.0) then
                  ga=fad(1,ks+3)
                  gg=fgd(1,ks+3)
                else
                  ga=fad(2,ks+3)*jkr
                  gg=fgd(2,ks+3)*jkr
                end if
                z=w3j(lr,2,lc,kr,-ks,kc)
                ra2=ga*z
                rg2=gg*z
              end if
c                                      --- for g tensor
              rg=rg1+cplkc*jkc*rg2

c                                      --- for A tensor
              if(in2.ne.0) ra=ra1+cplkc*jkc*ra2
	    end if
c
c----------------------------------------------------------------------
c       Calculate off-diagonal terms of the diffusion operator,
c       including the potential-dependent portion and the rhombic
c       component of the isotropic diffusion operator.
c
c       See Meirovitch, Igner, Igner, Moro,and Freed,
c       J. Chem. Phys. 77 (1982) pp.3933-3934, and especially
c       Eqns. A22-24 and A40 for more information.
c----------------------------------------------------------------------
c
c           -------- Rhombic part of isotropic diffusion operator
c
            if((ld.eq.0).and.(kr+2.eq.kc)) then
              cdff=rmi*dsqrt(dble((lr+kc-1)*(lr+kc)*
     #             (lr-kc+1)*(lr-kc+2)))/cnk
            else
              cdff=zero
            end if
c
c           -------- Potential-dependent part of diffusion operator
c
            if((ipt.ne.0).and.(ldabs.le.lband).and.(ipar(ks).eq.1)
     #        .and.(jkd.eq.0).and.((kdabs.le.kband)
     #	      .or.(ksabs.le.kband))) then
c
              kdi=1+kdabs/2
              ksi=1+ksabs/2
              c1=cpmkr*cnl*cnk 
c
c             Loop over L and K indices in sum over terms in potential
c
              do 444 l=0,lband,2
                cdffu=0.0D0
                li=1+l/2
                if (ksi.le.li.and. abs(xlk(li,ksi)).ge.rndoff)
     #               cdffu=cplkc*jkc*xlk(li,ksi)*w3j(lr,l,lc,kr,-ks,kc)
                if (kdi.le.li.and. abs(xlk(li,kdi)).ge.rndoff)
     #               cdffu=cdffu+xlk(li,kdi)*w3j(lr,l,lc,kr,-kd,-kc)
c
                if (abs(cdffu).gt.rndoff) cdff=cdff+
     #                           w3j(lr,l,lc,mr,0,-mr)*c1*cdffu
 444          continue
c
            end if
          end if
c
          if(newmc) then
            ioldmc=mc
            md=mr-mc
            ms=mr+mc
            mdabs=iabs(md)
            cplmc=ipar(lc+mc)
c
c----------------------------------------------------------------------
c     set flags and zero out matrix elements
c     if no non-zero elements are possible
c----------------------------------------------------------------------
c
            if(fld.or.(md.eq.0)) then
              fldmd=.true.
            else
              cliou=0.0D0
              cgam=0.0D0
              fldmd=.false.
            end if
c
c----------------------------------------------------------------------
c     get appropriate 3-J symbols if non-zero Liouville
c     operator matrix elements are possible
c----------------------------------------------------------------------
c
            wliou=zero
            if(fld) wliou=w3j(lr,2,lc,mr,-md,-mc)
c
          end if
c
          if(newpnc) then
            ioldpnc=ipnc
            ipnd=ipnr-ipnc
            ipndab=iabs(ipnd)
            ipns=ipnr+ipnc
c
            cnorm=cnl*cnk*cpmkr
c
c----------------------------------------------------------------------
c     get director tilt dependence of Liouville
c     operator matrix elements
c----------------------------------------------------------------------
c
            d1=zero
            if((ipndab.le.2).and.(mdabs.le.2))
     #                           d1=d2km(1,ipnd+3,md+3)*wliou
c
          end if
c
          if(fldmd) then
            iqnd=iqnr-iqnc
            iqndab=iabs(iqnd)
c
c----------------------------------------------------------------------
c    Liouville operator
c
c    Note: see Meirovitch, Igner,Igner,Moro, and Freed
c    J. Chem. Phys. 77 (1982) Appendix B p.3936-3937 for details.
c
c          Sa and Sg in this program also contain the
c          appropriate Clebsch-Gordon coefficient, and fki is
c          Ki in text.
c
c          The manner in which the isotropic contributions from the
c          spin Hamiltonian are included in this program is not clear
c          from the text, and  deserves some comments.
c          When Lc=Lr, both the (script) l=0 and l=2 ISTO's must
c          be included, whereas only l=2 terms appear off the diagonal.
c          The Clebsch-Gordan coefficients for these terms are
c          c(1,0,1,0;0,0) and c(1,+/-1,1,-/+1;0,0), or -/+ 1/sqrt(3),
c          respectively. Here, cga0 is set to +/- 1 instead of -/+
c          1/sqrt(3) because the extra factors are already included in g0
c          and a0 (note sign change). The isotropic terms are not
c          normalized by any factor.
c	   The cnk factor is not needed because the l=0 ISTO components
c	   have not been transformed by the K-symmetrization.
c          The factor cnl*cpmkr is canceled by the product of
c          w3j(lr,mr,0,0,lr,-mr) and w3j(lr,kr,0,0,lr,-kr), which
c          therefore do not need to be calculated.
c
c----------------------------------------------------------------------
c
            cliou=zero
c
            if(fld) then
              clioua=zero
              clioug=zero
              cliou0=zero
c
c        **   Hyperfine interaction   **
c
              if((jqed.eq.0).and.(in2.ne.0).and.(mdabs.le.2).and.
     #          (ipndab.eq.iqndab).and.(ipndab.le.1)) then
c
                cga0=zero
                if(ipnd.ne.0) then
c                                         *** pseudosecular term ***
c					      (delta pI .ne. 0)
                  kip=iqnr*iqnd+ipnr*ipnd
                  kip=kip*(kip-2)
                  fki=dsqrt(fii-0.25D0*kip)
                  sa1=-iqnd*fki*dsq18
                  cga2=dsq12
                else
c                                         *** secular term ***
c					      (delta pS,pI .eq. 0)
                  sa1=ipnr*0.5D0
	          cga0=one
	          cga2=dsq23
                end if
                clioua=sa1*d1*ra*cga2
                if(flk0.and.(jkd.eq.0).and.(md.eq.0))
     #                    cliou0=cliou0+sa1*cga0*a0
c
              end if
c
c        **   Nuclear Zeeman interaction (isotropic)  **
c
c              if(fnz.and.(jqed.ne.0).and.flk0.and.(jkd.eq.0).and.
c     #          (md.eq.0).and.(ipnd.eq.0).and.(iqnd.eq.0)) then
c                cliou0=cliou0+ipnr*zeen
c              end if
c
              cliou=cnorm*(clioua+clioug)+cliou0
c
            end if
c
c----------------------------------------------------------------------
c           Diffusion operator
c----------------------------------------------------------------------
c
            cgam=zero
c
c      **   Heisenberg exchange   **
c
            if((dabs(ossc).gt.rndoff).and.flk0.and.(jkd.eq.0).and.
     #         (jqed.eq.0)) then
              c=zero
              ctemp=zero
              if((md.eq.0).and.(ipnd.eq.0)) then
                if(iqnd.eq.0) c=1.0D0
                if((ipnr.eq.0).and.(lr.eq.0)) c=c-one/dble(in2+1)
              end if
              if((ms.eq.0).and.(ipns.eq.0)) then
                if(iqnd.eq.0) ctemp=1.0D0
                if((ipnr.eq.0).and.(lr.eq.0)) 
     #                               ctemp=ctemp-one/dble(in2+1)
                ctemp=ctemp*jqec*cplmc
              end if
              cgam=ossc*(c+ctemp)*0.5d0
            end if
c
c      **   Potential-dependent terms   **
c
            if((ipnd.eq.0).and.(iqnd.eq.0).and.(md.eq.0).and.
     #        (jkd.eq.0).and.(jqed.eq.0))     cgam=cgam+cdff
c
c      **   Isotropic electronic T1 relaxation (We)  **
c
            if((jqed.eq.0).and.(ld.eq.0).and.(jkd.eq.0).and.
     #        (kd.eq.0).and.(iqnd.eq.0)) then
c
              if((mr.eq.mc).and.(ipnr.eq.ipnc)) cgam=cgam+t1e
              if((mr.eq.-mc).and.(ipnr.eq.-ipnc))
     #                                   cgam=cgam+jqec*cplmc*t1e
            end if
c
c      **   Isotropic nuclear T1 relaxation (Wn)  **
c           (specifically for I=1/2 and I=1 case)
c
            if((jqed.eq.0).and.(ipnd.eq.0).and.(ld.eq.0).and.
     #        (jkd.eq.0).and.(kd.eq.0).and.(md.eq.0)) then
c                                   * I=1/2 case
              if(in2.eq.1) then
                if(ipnr.eq.0) then
                  if(iqnd.eq.0) then
                    cgam=cgam+t1n
                  else if(iqndab.eq.2) then
                    cgam=cgam-t1n
                  end if
                else if((ipnr.eq.1).or.(ipnr.eq.-1)) then
                  cgam=cgam+t1n*7.0d0/6.0d0
                end if
c                                   * I=1 case
              else if(in2.eq.2) then
                if(ipnr.eq.0) then
                  if(iqnd.eq.0) then
                    if((iqnr.eq.in2).or.(iqnr.eq.-in2)) then
                      cgam=cgam+t1n
                    else
                      cgam=cgam+t1n*2.0d0
                    end if
                  else if(iqndab.eq.2) then
                    cgam=cgam-t1n
                  end if
                else if((ipnr.eq.1).or.(ipnr.eq.-1)) then
                  if(iqndab.eq.0) cgam=cgam+t1n*13.0d0/6.0d0
                  if(iqndab.eq.2) cgam=cgam-t1n
                else if((ipnr.eq.2).or.(ipnr.eq.-2)) then
                  cgam=cgam+t1n*11.0d0/3.0d0
                end if
              end if
            end if
c
c      **   T2n relaxation HERE !!!   **
c
          end if
c
c
c----------------------------------------------------------------------
c         Store matrix elements for Rutishauser algorithm
c         ***  stores upper triangular matrix only
c----------------------------------------------------------------------
c
          if(nrow.eq.ncol) then
            cgam=cgam+cdfd
            zdiag(1,nrow)=cgam
            zdiag(2,nrow)=-cliou
          else
c
            if(abs(cliou).gt.rndoff) then
              nel=nel+1
              if(nel.gt.mxel) then
                ierr=1
                return
              else
                neli=neli+1
                zmat(neli)=-cliou
                izmat(neli)=ncol
              end if
            end if
c
            if(abs(cgam).gt.rndoff) then
              nel=nel+1
              if(nel.gt.mxel) then
                ierr=1
                return
              else
                nelr=nelr+1
                zmat(mxel-nelr+1)=cgam
                izmat(mxel-nelr+1)=ncol
              end if
            end if
          end if
c
c----------------------------------------------------------------------
c     Increment column counter
c----------------------------------------------------------------------
c
c  ****test*****
c         if (nrow.eq.1) then
c         if (ncol.eq.1)       write(6,332)
c         write(6,333) ncol,iqec,lc,jkc,kc,mc,ipnc,iqnc
c 332     format(2x,' ncol  iqec   lc    jkc   mc   ipnc  iqnc')
c 333     format(2x,8i6)
c         end if
c
          if(ncol.gt.mxdim) then
            ierr=2
            return
          end if
c----------------------------------------------------------------------
c         end loop over columns
c----------------------------------------------------------------------
 300    continue       
c----------------------------------------------------------------------
c       end loop over rows
c----------------------------------------------------------------------
 200  continue
c
      ierr=0
      neltot=nel
      nelre=nelr
      nelim=neli
c
      jzmat(ndimd+1)=neli+1
      kzmat(ndimd+1)=nelr+1
c
      return
      end
