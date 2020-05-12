	module mrcnsl0
	use msylm
	implicit none
	contains
	

      subroutine rcnsl0(tau,a,alat, dl) ! lmax,nlm,alat,rlat,nkr,dlat,nkd,vol,cy,
	    use esvar
C  reduced structure constants on lattice for e=0 and q=0.
C  result is periodic sum of (2*l-1)!!*ylm/r**(l+1). additive
C  constant for l=0 is chosen so that function averages to zero.
      implicit none
C passed parameters
      !integer, intent(in):: lmax,nlm,nkd,nkr
      double precision, intent(in):: a,alat !,vol
      double precision, intent(in):: tau(3) !, dlat(3,nkd), rlat(3,nkr),cy(nlm)
      double precision, intent(out):: dl(nlm)

C local parameters
      integer ilm
	    dl = 0.0d0
!       external rcsl01,rcsl02
	    !write(*,'(10000f10.2)') tau,a,lmxst,nlm,alat, vol

!       call tcn('k-space Ewald')
      call rcsl01(tau,a,lmxst,nlm,alat,glat,nkg,vol,dl)
!       call tcx('k-space Ewald')

	    !write(*,'(a,10000f10.4)')'k: dl = ',dl

!       call tcn('r-space Ewald')
      call rcsl02(tau,a,lmxst,nlm,alat,dlat,nkd,dl)
!       call tcx('r-space Ewald')

	    !write(*,'(a,10000f10.4)')'r: dl = ',dl

	    !write(*,'(a,10000f10.2)')'dl = ',dl(1:20)


	    !dl(1:nlm) = dl(1:nlm)*cy(1:nlm)
	    !write(*,'(a,10000f10.2)')'dl = ',dl(1:20)
      end


      subroutine rcsl01(tau,a,lmax,nlm,alat,rlat,nkr,vol,dl)
C  k-space part of reduced structure constants for e=0 and q=0
      implicit none
C passed parameters
      integer, intent(in) :: lmax,nkr, nlm
      real(8), intent(in) :: a,alat,vol,tau(3),rlat(3,nkr)
      real(8), intent(out) ::  dl(nlm)
C local parameters
      integer, parameter :: lmaxx = 10
      real(8), parameter :: tpi = 2.0_8*3.14159265358979323846_8
      integer :: ilm,ir,l,m !,tpi
      real(8) :: fpibv,gamma,scalp,tpiba,yyy
      real(8) :: eiphi(0:1),xxx,r2,r(3),yl((lmaxx+1)**2)
      !external sylm

      if (lmax > lmaxx) write(*,*)'rcnsl0: increase lmaxx'
!       tpi=8d0*atan(1d0)
      gamma = .25d0/(a*a)
      fpibv = 2d0*tpi/vol
      tpiba = tpi/alat
      !nlm = (lmax+1)**2
      dl(1:nlm) = 0d0
      dl(1) = -fpibv*gamma
      do  ir = 2, nkr
        r(1:3) = tpi*rlat(1:3,ir) !tpiba*rlat(1:3,ir)
!         print '(a,x, 3(x,f20.12))', " rlat:", rlat(1:3,ir)
!         print '(a,x, 3(x,f20.12))', "r:", r
        scalp = sum(r(1:3)*tau(1:3)) !alat*sum(r(1:3)*tau(1:3))
        eiphi = [cos(scalp), sin(scalp)]
!         eiphi(0) = cos(scalp)
!         eiphi(1) = sin(scalp)
!         print '(a,x, 2(x,f20.12))', "  eiphi:", eiphi
        call sylm(r,yl,lmax,r2)
!         print '(a,x, f20.12)', "  r2:", r2
        yyy = fpibv*exp(-gamma*r2)/r2
        ilm = 0
        do  l = 0, lmax
!           print '(a,x, 2(x,f20.12))', "  eiphi:", eiphi
          do  m = 1, 2*l+1
            ilm = ilm+1
            dl(ilm) = dl(ilm) + yl(ilm)*yyy*eiphi(0)
!             print '(a,x, 5(x,f20.12),3(x,i0))', "dl:", r, dl(ilm),yl(ilm),l,m,lmax
          enddo
          eiphi = [eiphi(1),-eiphi(0)]
        enddo
      enddo

      end


      subroutine rcsl02(tau,a,lmax,nlm,alat,dlat,nkd,dl)
C  real space summation
      implicit none
C passed parameters
      integer, intent(in) :: lmax,nkd, nlm
      double precision, intent(in) :: a,alat
      double precision, intent(in):: tau(3),dlat(3,nkd)
      double precision, intent(inout):: dl(nlm)

C local parameters
      integer ilm,ir,ir1,l,m,lmaxx
      parameter(lmaxx=9)
      double precision a2,cc,gl,r1,srpi,ta2,derfc,ddot,dl0,r2
      double precision r(3),yl((lmaxx+1)**2),chi(0:lmaxx)
      real(8), parameter :: twoinvsqpi = 1.12837916709551257390_8 ! 2/sqrt(pi)
      !external sylm

!       srpi=dsqrt(4d0*datan(1d0))
      ir1 = 2
!       if (tau(1)**2+tau(2)**2+tau(3)**2 > 1d-6) ir1=1
      if (sum(tau*tau) > 1d-6) ir1=1

      if (lmax > 0) then
        a2 = a*a
        ta2 = 2d0*a2
!         cc = 4d0*a2*a/srpi
        cc = ta2*a*twoinvsqpi
        do  ir = ir1, nkd
          r(1:3) = alat*(tau(1:3)-dlat(1:3,ir))
          call sylm(r,yl,lmax,r2)
          r1 = dsqrt(r2)
          chi(0) = derfc(a*r1)/r1
          gl = -cc*dexp(-a2*r2)/ta2
          do  l = 1, lmax
            chi(l) = ((2*l-1)*chi(l-1) - gl)/r2
            gl = ta2*gl
          enddo
          ilm = 0
          do  l = 0, lmax
!             do  m = 1, 2*l+1
!               ilm = ilm+1
!               dl(ilm) = dl(ilm) + yl(ilm)*chi(l)
!             enddo
            m = 2*l+1
	      dl(ilm+1:ilm+m) = dl(ilm+1:ilm+m) + yl(ilm+1:ilm+m)*chi(l)
            ilm = ilm + m
          enddo
        enddo
c...In case lmax = 0 do everything explicitly
      else
        dl0 = 0d0
        do  ir = ir1, nkd
          r(1:3) = tau(1:3)-dlat(1:3,ir)
          r1 = alat*dsqrt(sum(r*r))
          dl0 = dl0 + derfc(a*r1)/r1
        enddo
        dl(1) = dl(1) + dl0
      endif

C --- add dl3 for diagonal sructure constants ------
!       if (ir1 == 2) dl(1) = dl(1) - 2d0*a/srpi
      if (ir1 == 2) dl(1) = dl(1) - a*twoinvsqpi
      end
	end 	module mrcnsl0
