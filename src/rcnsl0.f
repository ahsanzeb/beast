	module mrcnsl0
	use msylm
	implicit none
	contains
	

      subroutine rcnsl0(tau,a,dl) ! lmax,nlm,alat,rlat,nkr,dlat,nkd,vol,cy,
	    use esvar
C  reduced structure constants on lattice for e=0 and q=0.
C  result is periodic sum of (2*l-1)!!*ylm/r**(l+1). additive
C  constant for l=0 is chosen so that function averages to zero.
      implicit none
C passed parameters
      !integer, intent(in):: lmax,nlm,nkd,nkr
      double precision, intent(in):: a !,vol
      double precision, intent(in):: tau(3) !, dlat(3,nkd), rlat(3,nkr),cy(nlm)
      double precision, intent(out):: dl(nlm)

C local parameters
      integer ilm
	    dl = 0.0d0
	    !write(*,'(10000f10.2)') tau,a,lmxst,nlm,alat, vol

!     a = s_lat%awald = 0.32*ewalda/r_cut [estatic]

      !call rcsl01(tau,a,lmxst,nlm,glat,nkg,vol,dl)
	    !write(*,'(a,10000f10.4)')'k: dl = ',dl

      call rcsl02(tau,a,lmxst,nlm,dlat,nkd,dl)

	    !write(*,'(a,10000f10.4)')'r: dl = ',dl

			!write(*,*)'rcnsl0: normalisation disabled.... testing...'
	    dl(1:nlm) = dl(1:nlm)*cy(1:nlm)

      end subroutine rcnsl0




! all vectors have their natural dimensions, nothing rescaled to make dimless.
      subroutine rcsl01(tau,a,lmax,nlm,rlat,nkr,vol,dl)
C  k-space part of reduced structure constants for e=0 and q=0
      implicit none
C passed parameters
      integer, intent(in) :: lmax,nkr, nlm
      real(8), intent(in) :: a,vol,tau(3),rlat(3,nkr)
      real(8), intent(out) ::  dl(nlm)
C local parameters
      integer, parameter :: lmaxx = 10
      real(8), parameter :: tpi = 2.0_8*3.14159265358979323846_8
      integer :: ilm,ir,l,m
      real(8) :: fpibv,gamma,scalp,yyy,yyy0
      real(8) :: eiphi(0:1),xxx,r2,r(3),yl((lmaxx+1)**2)
      !external sylm

      if (lmax > lmaxx) write(*,*)'rcnsl0: increase lmaxx'
      gamma = .25d0/(a*a)
      fpibv = 2d0*tpi/vol
      dl(1:nlm) = 0d0
      dl(1) = -fpibv*gamma
      do  ir = 2, nkr
        scalp = sum(rlat(1:3,ir)*tau)
        eiphi = [cos(scalp), sin(scalp)]
        call sylm(rlat(1:3,ir),yl,lmax,r2)
        yyy = fpibv*exp(-gamma*r2)/r2
        ilm = 0
        do  l = 0, lmax
!           print '(a,x, 2(x,f20.12))', "  eiphi:", eiphi
	        yyy0 = yyy*eiphi(0)
          do  m = 1, 2*l+1
            ilm = ilm+1
            dl(ilm) = dl(ilm) + yl(ilm) *yyy0
!             print '(a,x, 5(x,f20.12),3(x,i0))', "dl:", r, dl(ilm),yl(ilm),l,m,lmax
          enddo
          eiphi = [eiphi(1),-eiphi(0)]
        enddo
      enddo

      end subroutine rcsl01
!----------------------------------------------------------------------
! maz, ref: W. Smith, CCP5 Newsletter No. 46, 1998
! 
!----------------------------------------------------------------------
      subroutine rcsl02(tau,a,lmax,nlm,dlat,nkd,dl)
C  real space summation
      implicit none
C passed parameters
      integer, intent(in) :: lmax,nkd, nlm
      double precision, intent(in) :: a
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
      if (sum(tau*tau) > 1d-6) ir1=1

      if (lmax > 0) then
        a2 = a*a
        ta2 = 2d0*a2
        ! cc: bug in tb routines fixed by matching our code 
        !   with Eq 33 in Smiths noes--- notes r bugged, tb code correct!!!!!!
        cc = twoinvsqpi*a; !/(2.0d0*a);
        do  ir = ir1, nkd
          r(1:3) = tau(1:3)-dlat(1:3,ir)
          call sylm(r,yl,lmax,r2)
          r1 = dsqrt(r2)
          chi(0) = derfc(a*r1)/r1 ! Eq 25 in Smith's notes, Chi(0): B_0
          gl = cc*dexp(-a2*r2) ! e^{-alpha^2 * r^2}/(alpha*sqrt[pi])
          do  l = 1, lmax
	          gl = ta2*gl ! (2*alpha)^{l} * e^{-alpha^2 * r^2}/(alpha*sqrt[pi])
            chi(l) = ((2*l-1)*chi(l-1) + gl)/r2 ! Eq. 33 in Smith's notes
          enddo
          ilm = 0
          do  l = 0, lmax
            m = 2*l+1
	      dl(ilm+1:ilm+m) = dl(ilm+1:ilm+m) + yl(ilm+1:ilm+m)*chi(l)
            ilm = ilm + m
          enddo
        enddo
c...In case lmax = 0 do everything explicitly
      else
      !write(*,*)'========> lmax = 0'
        dl0 = 0d0
        do  ir = ir1, nkd
          r(1:3) = tau(1:3)-dlat(1:3,ir)
          r1 = norm2(r)
          !if(r1 < 50.0d0) then ! testing......... maz... remove this condition.
          dl0 = dl0 + derfc(a*r1)/r1
          !endif
        enddo
        dl(1) = dl(1) + dl0
      endif

C --- add dl3 for diagonal sructure constants ------
      if (ir1 == 2) dl(1) = dl(1) - a*twoinvsqpi
      
      end subroutine rcsl02
      
	end 	module mrcnsl0
