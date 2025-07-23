c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine COVRPT  
c                    =========================
c
c    Produces a report of the covariance matrix and parameter 
c    uncertainties on the specified logical unit.
c----------------------------------------------------------------------
c
      subroutine covrpt( lu )
c
      use nlsdim
      use eprprm
      use expdat
      use parcom
      use lmcom
      use lpnam
      use mspctr
      use iterat
      use rnddbl
      use stdio
c
      implicit none
      integer lu
c
      integer i,isi,isp,j,k,mrgn
      double precision denom,s(5,MXSITE)
      character*132 line
      logical potflg
c
      integer itrim
      double precision chi,ts,fs
      external itrim,chi,ts,fs
c
c     ---------------------------------------------
c      Append labels for scale factors to tag array
c     ---------------------------------------------
      if (nsite.gt.1) then
         do isi=1,nsite
            write(tag(nprm+isi),1010) isi
            x(nprm+isi)=sfac(isi,1)
         end do
      else
         do isp=1,nser
            write(tag(nprm+isp),1011) isp
            x(nprm+isp)=sfac(1,isp)
         end do
      end if
c
c-----------------------------------------------------------------------
c  Calculate residuals and covariance matrix
c  COVAR calculates covariance matrix from the R matrix of QR 
c  factorization, which lmder stored in upper triangular form in fjac 
c  ipvt is permutation array returned by lmder
c------------------------------------------------------------------------
c
      if (.not.covarOK) then
         call covar( njcol,fjac,MXPT,ipvt,xtol*xtol,work1 )
         covarOK=.true.
c                                    --- weighted residual fit
         if (weighted_flag.ne.0) then
            delchi1=chi(confid,1)
            ch2bnd=chi(confid,2)
            chnbnd=chi(confid,njcol)
c                                    --- unweighted residual fit
         else
            tbound=ts( (1.0D0-confid)/2.0D0, ndatot-nprm )
            f2bnd=2.0D0*fs( 1.0D0-confid, 2, ndatot-2 )*fnorm*fnorm/
     #            float(ndatot-2)
            fnbnd=float(njcol)*fs( 1.0D0-confid, njcol, ndatot-njcol )*
     #            fnorm*fnorm/float(ndatot-njcol)
         end if
c
      end if
c
c     ------------------------------------------------------------ 
c     Calculate correlation matrix from covariance and output it
c     ------------------------------------------------------------ 
      line=' '
      mrgn=35-5*njcol
      if (mrgn.lt.1) mrgn=1
      write (lu,2000) 
      write(line(mrgn:),2002) (tag(i),i=1,njcol)
      write(lu,2001) line(:10*(njcol+1)+mrgn)
c
      do i=1,njcol
         line=' '
         do j=i,njcol
c
c----------------------------------------------------------------------
c Estimate of error bounds:
c     For weighted residuals, use chi-squared statistics
c     For unweighted residuals use t-statistics
c----------------------------------------------------------------------
            if (i.eq.j) then
               if (weighted_flag.ne.0) then
                  xerr(i)=sqrt(delchi1*fjac(i,i))
               else
                  xerr(i)=tbound*fnorm*sqrt(fjac(i,i))
     #                    /sqrt(float(ndatot-njcol))
               end if
            end if

            denom = dsqrt( fjac(i,i) * fjac(j,j) )
            if (denom .ne. 0.0D0) corr(i,j) = fjac(i,j)/denom
            k=10*(j-1)+mrgn
            write(line(k:),2003) corr(i,j)
         end do
         write (lu,2001) line(:10*(njcol+1)+mrgn)
      end do
c     
c    --------------------------------
c     output final fit of parameters 
c    --------------------------------
      write (lu,1000)
      if (weighted_flag.ne.0) then
         write (lu,2004) 'Chi'
      else
         write(lu,2004) 'T'
      end if
c
      write (lu,2007) (tag(i)(:itrim(tag(i))),x(i),xerr(i),i=1,njcol)
c
c                              --- weighted residual fit ---
      if (weighted_flag.ne.0) then
         write (lu,2008) confid,delchi1,ch2bnd,njcol,chnbnd
c
c                              --- weighted residual fit ---
      else
         write (lu,2009) confid,tbound,f2bnd,njcol,fnbnd
      end if
c
c     --------------------------------------------------------
c     Calculate and output order parameters for each site if
c     appropriate
c     --------------------------------------------------------
      do isi=1,nsite
         potflg=.false.
         do j=0,4
            potflg=potflg.or.(abs(fparm(IC20+j,isi)).gt.RNDOFF)
         end do
c      
         if (potflg) then
            call ordrpr( fparm(IC20,isi),s(1,isi) )
            write (lu,1060) isi
            do j=0,4
               if (abs(fparm(IC20+j,isi)).gt.RNDOFF) 
     #              write (lu,1070) parnam(IC20+j)(2:3),s(j+1,isi)
            end do

         end if

      end do
c
      return
c
c##### Formats #########################################################
c
 1000 format(/,2x,70('-'),/)
 1010 format('SITE',i1)
 1011 format('SPCTRM',i1)
 2000 format(/24x,'*** Correlation Matrix ***'/)
 2001 format(a)
 2002 format(12(2x,a8))
 2003 format(f8.4)
 2004 format(/9x,'*** Final Parameters ***',/,
     1       7x,' Parameter      Value',13x,'Uncertainty (',a,
     2       ' statistics)'/7x,46('-')/)
 2007 format(10x,a9,' = ',g13.7,' +/- ',g13.7)
 2008 format(/7x,'Confidence= ',f5.3,3x,'Del-Chi-Sqr= ',g11.5/
     #        7x,'Del-Chi-Sqr for confidence regions:'/
     #        7x,'2 parameters: ',g11.5,4x,i2,' parameters: ',g11.5)
 2009 format(/7x,'Confidence = ',f5.3,3x,'T-bound=',g11.5/
     #        7x,'F-bound for confidence region:'/
     #        7x,'2 parameters: ',g11.5,4x,i2,' parameters: ',g11.5)
 1060 format(/10x,'Order parameters for site',i2,':')
 1070 format(15x,'<D',a2,'> =  ',f10.4,' +/- ')
c
      end


