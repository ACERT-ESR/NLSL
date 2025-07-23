c Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                       WRITR subroutines
c                    =========================
c
c  Outputs various information regarding intermediate spectral 
c  calculations for a specific iteration during a nonlinear least-squares
c  fit.
c
c  Contains the routines:
c    wrfun     Write function (calculated spectra) into <nlsfun.nnn> 
c    wrspc     Write spectra (experimental) and current fit into <dataid>.spc
c    wrjac     Write Jacobian
c    wrtrdg    Write tridiagonal matrices 
c
c----------------------------------------------------------------------      
c                    =========================
c                       subroutine WRFUN
c                    =========================
c
c  Outputs the current calculated spectrum matrix (spectr in commmon
c  /mspctr/ into an ASCII file named <nlsfun.nnn> where nnn is a three-
c  digit iteration number (passed as an argument)
c----------------------------------------------------------------------
c
      subroutine wrfun(iter)
c
      use nlsdim
      use expdat
      use parcom
      use lmcom
      use mspctr
      use stdio
c
      integer iter
      integer i,isp,ixs,j,k
      character tmpnam*40
c
      integer itrim
      external itrim
c
c    --------------------------------------------------------
c     Open file for printing iterates (function evaluations)
c    --------------------------------------------------------
      write(tmpnam,2001) iter
      open(unit=ludisk,file=tmpnam,status='unknown',
     #     access='sequential',form='formatted')
c
c     --------------------
c      Loop over spectra
c     --------------------
      do isp=1,nser
         ixs=ixsp(isp)
c
c     ----------------------------------------------------------
c      Output calculation information and spectra for all sites
c     ----------------------------------------------------------
c
         write (ludisk,3005) (tag(i)(:itrim(tag(i))),x(i),i=1,nprm)
         write (ludisk,3006) isp,sbi(isp)-shft(isp)-tmpshft(isp),
     #                       sdb(isp),(sfac(j,isp),j=1,nsite)
         do i=1,npts(isp)
            k=ixs+i-1
            write (ludisk,2003) (spectr(k,j),j=1,nsite)
         end do
c
      end do
c
      close(ludisk)
      return
c
c######################################################################
c
 2001 format('nlsfun.',i3.3)
 2003 format(10(g15.8,1x))
 3005 format('c ',a9,' = ',f10.5)
 3006 format('c Spectrum ',i2,', fieldi=',f9.3,' step=',f7.4/
     #       'c Scales:'/
     #       'c',10(g15.8,1x))
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine WRSPC
c                    =========================
c
c  Outputs the current experimental spectra and calculated function
c  into the files <dataID1>.spc <dataID2>.spc, ... for each member
c  of the series being fit.
c----------------------------------------------------------------------
      subroutine wrspc()
c
      use nlsdim
      use expdat
      use parcom
      use lmcom
      use iterat
      use nlsnam
      use mspctr
      use stdio
c
      implicit none
      integer i,j,k,l
      double precision field
c
      do i=1,nspc
         call setdat( dataid(i) )
c     
         open(ludisk,file=spname(:lthdnm),status='unknown',
     #        access='sequential',form='formatted')
         field=sbi(i)
         do j=1,npts(i)
            k=ixsp(i)+j-1
c                                *** Multiple sites
            if (nsite.gt.1) then
c                                        *** weighted residuals
               if (weighted_flag.ne.0) then
                  write(ludisk,1045) field,data(k),
     #                 data(k)-rmsn(i)*fvec(k),
     #                 (sfac(l,i)*spectr(k,l),l=1,nsite)

c                                        *** unweighted residuals
               else
                  write(ludisk,1045) field,data(k),
     #                 data(k)-fvec(k),
     #                 (sfac(l,i)*spectr(k,l),l=1,nsite)
               end if

c                                *** Single site
            else
c                                        *** weighted residuals
               if (weighted_flag.ne.0) then
                  write(ludisk,1045) field,data(k),
     #                 data(k)-rmsn(i)*fvec(k)

c                                        *** unweighted residuals
               else

                  write(ludisk,1045) field,data(k),
     #                 data(k)-fvec(k)
               end if

            end if
            field=field+sdb(i)
         end do
         close (ludisk)
c
      end do
      written=1
      return
c
 1045 format(f10.3,6(' ',g14.7))
c
      end

c----------------------------------------------------------------------
c                    =========================
c                       subroutine WRJAC
c                    =========================
c
c  Outputs the current Jacobian matrix in array fjac (commmon /lmcom/
c  into an ASCII file named <nlsjac.nnn> where nnn is a three-digit
c  iteration number (passed as an argument)
c----------------------------------------------------------------------
      subroutine wrjac(iter)
c
      use nlsdim
      use expdat
      use lmcom
      use parcom
      use stdio
c
      implicit none
      integer iter,i,j
      character tmpnam*40
c
      write(tmpnam,1000) iter
      open(unit=ludisk,file=tmpnam,status='unknown',
     #     access='sequential',form='formatted')
c
      write (ludisk,1002) (tag(j),j=1,njcol)
      write (ludisk,1003) (xfdstp(j),j=1,njcol)
c
      do i=1,ndatot
         write (ludisk,1001) (fjac(i,j),j=1,njcol)
      end do
c
      close(ludisk)
      return
c
c######################################################################
c
 1000 format('nlsjac.',i3.3)
 1001 format(10(g15.8,1x))
 1002 format('c   ',10(a9,6x))
 1003 format('c  ',10(g15.7))     
      end



c----------------------------------------------------------------------
c                    =========================
c                       subroutine WRTRDG
c                    =========================
c
c  Outputs the current set of tridiagonal matrices into an ASCII
c  file named <nlstri.nnn> where nnn is a three-digit iteration
c  number (passed as an argument)
c----------------------------------------------------------------------
      subroutine wrtrdg(iter)
c
      use nlsdim
      use expdat
      use tridag
      use parcom
      use stdio
c
      implicit none
      integer iter,isp,isi,ixt,j
      character tmpnam*40
c
      write(tmpnam,1000) iter
      open(unit=ludisk,file=tmpnam,status='unknown',
     #      access='sequential',form='formatted')
c
      do isp=1,nser
         do isi=1,nsite
            ixt=ixtd(isi,isp)
            write (ludisk,1001) isp,isi
            write (ludisk,1002) (j,alpha(ixt+j-1),
     #             beta(ixt+j-1),j=1,ltd(isi,isp))
         end do
      end do
c
      close(ludisk)
      return
c
c######################################################################
c
 1000 format('nlstri.',i3.3) 
 1001 format('c Spectrum ',i2,', site ',i2)
 1002 format(i4,' (',g12.5,',',g12.5,') (',g12.5,',',g12.5,')')
      end

