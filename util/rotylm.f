
	! Ahsan Zeb, Aug 2023
	! Rotation matrix calculation of real spherical harmonics using the algorithm described in:
	! Joseph Ivanic and Klaus Ruedenberg, J. Phys. Chem. 1996, 100, 6342-6347

	! 1. uses m=-l:l order of bases for calculation of R^l as ivanic's paper require
	! 2. Starts with usual R3 rotation matrix in x,y,z order; shuffle it to feed it to ivanic algorithm.
	! 3. Shuffle the resulting full 0:lmax space matrix Rlmax to match Beast/TB bases order.

	!module rotylm
	program rotylm
	implicit none

	integer, parameter :: lmax = 4
	integer, parameter :: nlm = (1+lmax)**2
	double precision, dimension(1:nlm,1:nlm) :: Rlmax ! rot matrix in full l=0-lmax space

	! l-sector rotation matrices in ylm space
	type :: ylmrotmatrix
			double precision, allocatable, dimension(:,:) :: R ! rot matrix in ylm space
	end type
	type(ylmrotmatrix), dimension(0:lmax) :: Rl
	double precision, dimension(-1:1,-1:1) :: R

	! use ishm to make m=-l:l order
	! use ish to make TB/Beast order
	integer, parameter, dimension(9) :: ish=(/1,4,2,3,5,6,8,9,7/);
	integer, parameter, dimension(9) :: ishm=(/1,3,4,2,5,6,9,7,8/)
	integer, parameter, dimension(-1:1) :: ishm1= (/2,3,1/) !(/3,4,2/) - 1; 
	integer, parameter, dimension(5) :: ish2=(/-2,-1,1,2,0/);! (/5,6,8,9,7/) mapped to -2:2
	
	double precision :: theta, phi
	integer :: l, m, mp
	integer :: ilm, it



	!write(*,*) 'Enter theta:'
	!read(*,*) theta
	!phi = theta

	do l=0,lmax
	 allocate(Rl(l)%R(-l:l,-l:l))
	enddo

	open(1,file='Rlmax.dat',action='write')

	do it=1, 9

	 theta = (it-1)*(-2.d0)! 2 degree increment
	 phi = theta; ! for the time being

	 Rlmax = 0.0d0;
	 call getRlmax(theta, phi)	

	 do ilm=1,nlm
	  write(1,*) Rlmax(ilm,:)
	 end do

	end do ! ot
	close(1)

	do l=0,lmax
	 deallocate(Rl(l)%R)
	enddo

	write(*,'(a)')'rotylm: done!'


	contains

	subroutine getRlmax(theta, phi)
	implicit none
	double precision, intent(in) :: theta, phi
	integer :: l, ilm, i1,i2, im, imp, ii,jj,i,j
	double precision, dimension(3,3) ::	R3

	! initialise
	Rlmax = 0.0d0;
	R3 = 0.0d0;
	
	! calc real space rotation matrix
	call getRotMatComb(theta, phi, R3) ! R3 order: x,y,z

	! change order to: y,z,x; ==> R order
	do i=-1,1
	 do j=-1,1
	  ii=ishm1(i); ! ishm1=(/2,3,1/) : y ,z, x
	  jj=ishm1(j);
	  R(i,j) = R3(ii,jj)
	 end do
	end do

	write(*,'(a)')'......................................'
	write(*,'(a,3f10.5)') 'angle theta, phi = ',theta, phi
	write(*,'(3f10.5)') R3

	! l=0 set manually
	l=0;
	!allocate(Rl(0)%R(0:0,0:0))
	Rl(0)%R(0,0) = 1.0d0;
	Rlmax(1,1) = Rl(0)%R(0,0);
	i1=1;

	! l=1 also set manually, ivanic algorithm does not apply to this special case; they miss Eq. 7.3ab case.
	Rl(1)%R = R;
	i2=i1+3;
	Rlmax(i1+1:i2,i1+1:i2) = R3 !R3 has order, x,y,z. !Rl(1)%R has order: y,z,x
	i1=i2; !

	
	do l=2, lmax
		!allocate(Rl(l)%R(-l:l,-l:l))
		Rl(l)%R = 0.0d0;
		i2 = i1 + (2*l + 1)
		
	 do im=-l,l
	  do imp=-l,l
		  if(l <= 2) then ! shuffle ilm index to get beast/TB basis
		  		! ish(/1,4,2,3,5,6,8,9,7/);! 
				m = ish((l+1)**2-l + im) - ((l+1)**2-l) 
				mp = ish((l+1)**2-l + imp) - ((l+1)**2-l)				
		  else ! same index
				m=im; mp=imp;
		  endif
	    Rl(l)%R(m,mp) = cu(l,m,mp)*U(l,m,mp) + cv(l,m,mp)*V(l,m,mp)
     .                 + cw(l,m,mp)*W(l,m,mp) 	    
	  end do
	 end do
		
	!Rlmax(i1+1:i2,i1+1:i2) = Rl(l)%R
	! shuffle indices to make TB/beast code l<=2 order for Rlmax.
	! l loop l=2,lmax; l=1 not included.
	 if(l == 2) then ! 5x5 matrix
	  do i=1,5
	  do j=1,5
	   ii=ish2(i);
	   jj=ish2(j)
	   Rlmax(i1+i,i1+j) = Rl(l)%R(ii,jj) ! Rl(2)%R(-2:2,-2:2)
	  end do
	  end do
	 else ! m=-l:l order
	  Rlmax(i1+1:i2,i1+1:i2) = Rl(l)%R
	 endif

		i1 = i2
	end do ! l

	return
	end subroutine
!------------------------------------------
	double precision function P(i,l,mu,mp)
	implicit none
	integer, intent(in) :: i,l,mu,mp

	if(iabs(mp) < l) then
		P = R(i,0) * Rl(l-1)%r(mu,mp)
	elseif(mp == l) then
		P = R(i,1) * Rl(l-1)%r(mu,l-1) - R(i,-1) * Rl(l-1)%r(mu,-l+1)
	else !(mp == -l) then
		P = R(i,1) * Rl(l-1)%r(mu,-l+1) + R(i,-1) * Rl(l-1)%r(mu,l-1)
	endif

	return
	end function
!------------------------------------------


!------------------------------------------
	double precision function U(l,m,mp)
	implicit none
	integer, intent(in) :: l,m,mp
			U = P(0,l,m,mp)
	end function
!------------------------------------------
	double precision function V(l,m,mp)
	implicit none
	integer, intent(in) :: l,m,mp
	if(m == 0) then
	 V = P(1,l,1,mp) + P(-1,l,-1,mp)
	elseif(m > 0) then
	 V = P(1,l,m-1,mp) * dsqrt(dble(1+del(m,1))) 
	 V = V - P(-1,l,-m+1,mp) * dsqrt(dble(1-del(m,1))) 
	else !(m < 0)
	 V = P(1,l,m+1,mp) * dsqrt(dble(1-del(m,-1))) 
	 V = V + P(-1,l,-m-1,mp) * dsqrt(dble(1+del(m,-1))) !dsqrt term (1-del) in ivanic1996-correction
	endif
	end function
!------------------------------------------
	double precision function W(l,m,mp)
	implicit none
	integer, intent(in) :: l,m,mp
	if(m == 0) then
	 W = 0.0d0 ! at m=0: w^l_mm' =0; so W can be set to any finite number
	elseif(m > 0) then
	 W = P(1,l,m+1,mp) + P(-1,l,-m-1,mp)	
	else !(m < 0)
	 W = P(1,l,m-1,mp) - P(-1,l,-m+1,mp)	
	endif
	end function
!------------------------------------------




!------------------------------------------
	double precision function cu(l,m,mp)
	implicit none
	integer, intent(in) :: l,m,mp
	if(iabs(mp) < l) then
		cu = (l + m)*(l-m) 
		cu = cu/((l + mp)*(l-mp)) ! cu already double precision
		cu = dsqrt(cu)
	else 
		cu = (l + m)*(l-m) 
		cu = cu/((2*l)*(2*l-1)) 
		cu = dsqrt(cu)
	endif
	end function
!------------------------------------------
	double precision function cv(l,m,mp)
	implicit none
	integer, intent(in) :: l,m,mp
	if(iabs(mp) < l) then
	 cv = (1+del(m,0))*(l+iabs(m)-1)*(l+iabs(m))
	 cv = cv/((l+mp)*(l-mp))
	 cv = 0.5d0*dsqrt(cv)*(1-2*del(m,0))
	else 
	 cv = (1+del(m,0))*(l+iabs(m)-1)*(l+iabs(m))
	 cv = cv/((2*l)*(2*l-1))
	 cv = 0.5d0*dsqrt(cv)*(1-2*del(m,0))	! 
	endif
	end function
!------------------------------------------
	double precision function cw(l,m,mp)
	implicit none
	integer, intent(in) :: l,m,mp
	if(iabs(mp) < l) then
	 cw = (l-iabs(m)-1)*(l-iabs(m))
	 cw = cw/((l+mp)*(l-mp))
	 cw = -0.5d0*dsqrt(cw)*(1-del(m,0))
	else 
	 cw = (l-iabs(m)-1)*(l-iabs(m))
	 cw = cw/((2*l)*(2*l-1))
	 cw = -0.5d0*dsqrt(cw)*(1-del(m,0))
	endif
	end function
!------------------------------------------

	double precision function del(m,mp)
	implicit none
	integer, intent(in) :: m,mp
	 del = 0.0d0
	 if(m == mp) del = 1.0d0
	end function
!------------------------------------------







!..............................................................
! getRotMat() and getRotMatComb() from beast code
!..............................................................
	subroutine getRotMatComb(theta,phi,Rmat) !(il,io)
	implicit none
	double precision, intent(in) :: theta, phi !,il,io
	double precision, dimension(3,3), intent(out) ::  Rmat

	! rotation axis
	double precision, dimension(3) :: u1,u2,u3
	double precision, dimension(3,3) :: R1, R2

	u1 = (/0.0d0,0.0d0,1.0d0/); ! z axis
	call getRotMat(u1,phi,R1);

	u2 = (/1.0d0,1.0d0,0.0d0/)/dsqrt(2.0d0);
	call r3mv(R1,u2,u3) ! u3 = rotated 110
	call getRotMat(u3,theta,R2);

	! combined effect: Rmat = R2.R1
	Rmat = matmul(R2,R1);

	return
	end subroutine getRotMatComb
	!..............................................................

	subroutine getRotMat(u,th,Rmat)
	implicit none
	double precision, dimension(3), intent(in) :: u ! unit vector/axis
	double precision, intent(in) :: th ! theta/angle
	double precision, dimension(3,3), intent(out) ::	Rmat
	integer :: i,j
	! rotation axis
	double precision, dimension(3,3) :: uu, ux,id
	double precision :: tt
	
	! https://en.wikipedia.org/wiki/Rotation_matrix 
	! tensor product of u with u:
	do i=1,3
	 do j=1,3
	  uu(i,j) = u(i)*u(j);
	 end do
	end do
	! cross product matrix of u: 
	ux = 0.0d0;
	ux(1,2) = -u(3);
	ux(1,3) =  u(2);
	ux(2,1) =  u(3);
	ux(2,3) = -u(1);
	ux(3,1) = -u(2);
	ux(3,2) =  u(1);	

	!	identity matrix
	id = 0.0d0;
	id(1,1) = 1.0d0
	id(2,2) = 1.0d0	
	id(3,3) = 1.0d0

	! pi/180 = datan(1.0d0)/45.0d0 = 0.017453292519943295769d0
	tt = th*0.017453292519943295769d0; ! degree to radians
	
	! rotation matrix about the unit vector u for an angle th:
	Rmat = dcos(tt) * id + dsin(tt)*ux + (1.0d0-dcos(tt))*uu;

	return
	end subroutine getRotMat
	!..............................................................
	! copeid from elk-6.2.8
	subroutine r3mv(a,x,y)
	implicit none
	real(8), intent(in) :: a(3,3),x(3)
	real(8), intent(out) :: y(3)
	y(1)=a(1,1)*x(1)+a(1,2)*x(2)+a(1,3)*x(3)
	y(2)=a(2,1)*x(1)+a(2,2)*x(2)+a(2,3)*x(3)
	y(3)=a(3,1)*x(1)+a(3,2)*x(2)+a(3,3)*x(3)
	return
	end subroutine
	!..............................................

	end program !module rotylm
