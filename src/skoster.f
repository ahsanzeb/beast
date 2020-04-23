

	module skoster
	use modmain
	implicit none

	contains
!------------------------------------------------------------------------
	subroutine realHij()
	implicit none
	
	integer :: il,io,is, js, ii, k
	double precision, dimension(3,norbtm) :: h	
	double precision, dimension(3):: dr

	!...................................................................
	! TM-O  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	  is = tm(il,io)%is; ! species index of TM atom
	  do k = 1,6 ! 1st nns O
	   allocate(tm(il,io)%nn1(k)%h(norbtm,norbo))
	   ! slatkospd computed (p,d): to get (d,p), r ===> -r and aux h as out
		 dr = tm(il,io)%r - tm(il,io)%nn1(k)%r;
		 tm(il,io)%nn1(k)%dr = dr
	   call slatkospd(-1.d0*dr, skbo(is,:), h)
	   tm(il,io)%nn1(k)%h = transpose(h);
	  end do ! k
	 end do ! io
	end do ! il

	write(*,*)'-------------------- a: tmnn2=',tmnn2

	!...................................................................
	! 	TM-TM (2nd neighbours) 
	!...................................................................
	if(tmnn2) then
	do il=1,nlayers
	 do io=1, noctl
	  	is = 1; !tm(il,io)%is;
	  !	write(*,*)'is,norbtm = ',is, norbtm
	  do k = 1,6 ! 2nd nns TM
	   allocate(tm(il,io)%nn2(k)%h(norbtm,norbtm))
		 js = 1; !tm(il,io)%nn2(k)%is;
		 dr = tm(il,io)%r - tm(il,io)%nn2(k)%r;
		 tm(il,io)%nn2(k)%dr = dr
	   call slatkosdd(dr, skbb(is,js,:),tm(il,io)%nn2(k)%h)

	   !write(*,*)'tm(il,io)%nn2(k)%h = ',tm(il,io)%nn2(k)%h
	  end do ! k
	 end do ! io
	end do ! il
	endif
	!...................................................................

	write(*,*)'-------------------- b'


	!...................................................................
	! O-TM  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  !ia = ox(il,io,ii)%ia;
	  do k = 1,2 ! 1st nns TM
	   allocate(ox(il,io,ii)%nn1(k)%h(norbo,norbtm))
		 js = 1; !ox(il,io,ii)%nn1(k)%is
		 dr = ox(il,io,ii)%r - ox(il,io,ii)%nn1(k)%r
		 ox(il,io,ii)%nn1(k)%dr = dr
	   call slatkospd(dr, skbo(js,:),ox(il,io,ii)%nn1(k)%h)
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	!...................................................................
	write(*,*)'-------------------- c'

	!...................................................................
	! O-O  (2nd neighbours)
	!...................................................................
	if(oxnn2) then
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  !ia = ox(il,io,ii)%ia;
	  do k = 1,8 ! 2nd nns O
	   allocate(ox(il,io,ii)%nn2(k)%h(norbo,norbo))
		 dr = ox(il,io,ii)%r - ox(il,io,ii)%nn2(k)%r
		 ox(il,io,ii)%nn2(k)%dr = dr
	   call slatkospp(dr, skoo, ox(il,io,ii)%nn2(k)%h)
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	endif
	!...................................................................

	write(*,*)'-------------------- d'


	return
	end 	subroutine realHij
!------------------------------------------------------------------------


	! direction cosines, l,m,n for SK method
	function getlmn(x) result (lmn)
	implicit none
	double precision, dimension(3), intent(in) :: x
	double precision, dimension(3) :: lmn
	! local
	double precision :: r
	integer :: i

	r = norm2(x);
	do i=1,3
	 lmn(i) = x(i)/r
	end do
	
	end function getlmn
	!.....................................................
	! Hamiltonian matrix elements between p-p orbitals using SK method
	subroutine slatkospp(r, sk, h)
	implicit none
	double precision, dimension(3), intent(in) :: r
	double precision, dimension(2), intent(in) :: sk
	double precision, dimension(3,3), intent(out) :: h
	! local
	double precision :: rr
	double precision :: l,m,n ! direction cosines; real not integers.
	double precision :: sp, s, p
	double precision :: l2,m2,n2
	
	s = sk(1); p = sk(2); ! s=sigma_pp, p=pi_pp
	sp = s-p;

	rr = norm2(r);
	! direction cosines; r is from first to the second atom
	l = r(1)/rr;
	m = r(2)/rr;
	n = r(3)/rr;
	
	l2=l**2; m2=m**2; n2=n**2;
	
	! (x,x), (x,y), (x,z)
	h(1,1) = l2 *s + (1-l2) *p
	h(1,2) = l*m*sp
	h(1,3) = l*n*sp
	! (y,x), (y,y), (y,z)
	h(2,1) = h(2,1)
	h(2,2) = m2 *s + (1-m2) *p
	h(2,3) = m*n*sp
	! (z,x), (z,y), (z,z)
	h(3,1) = n*l*sp
	h(3,2) = n*m*sp
	h(3,3) = n2 *s + (1-n2) *p

	return
	end 	subroutine slatkospp
	!.....................................................




	!.....................................................
	! Hamiltonian matrix elements between p-d orbitals using SK method
	subroutine slatkospd(r, sk, h)
	implicit none
	double precision, dimension(3), intent(in) :: r
	double precision, dimension(2), intent(in) :: sk
	double precision, dimension(3,norbtm), intent(out) :: h
	! local
	double precision :: rr
	double precision :: l,m,n ! direction cosines
	double precision :: lmn, s, p, sq3, lp,mp,np
	double precision :: 	lm, lm2, nlm2
	double precision :: l2,m2,n2

	s = sk(1); p = sk(2); ! s=sigma_pd, p=pi_pd
	sq3 = dsqrt(3.0d0);

	rr = norm2(r);
	! direction cosines; r is from first to the second atom
	l = r(1)/rr;
	m = r(2)/rr;
	n = r(3)/rr;

	l2=l**2; m2=m**2; n2=n**2;

	!..................................................	
	! with t2g
	!..................................................
	lmn = l*m*n;
	lp = (1.d0-2*l**2)*p
	mp = (1.d0-2*m**2)*p
	np = (1.d0-2*n**2)*p	

	! x, xy
	h(1,1) = sq3 *l2 *m *s + m*lp !(1.d0-2*l**2)*p
	! x, yz
	h(1,2) = (sq3 *s - 2*p)*lmn
	! x, zx
	h(1,3) = sq3 *l2 *n *s + n*lp !(1.d0-2*l**2)*p

	! y, xy
	h(2,1) = sq3 *m2 *l *s + l* mp !(1.d0-2*m**2)*p
	! y, yz
	h(2,2) = sq3 *m2 *n *s + n* mp !(1.d0-2*m**2)*p
	! y, zx
	h(2,3) = h(1,2) !(sq3 *s - 2*p)*lmn

	! z, xy
	h(3,1) = h(1,2) !(sq3 *s - 2*p)*lmn
	! z, yz
	h(3,2) = sq3 *n2 *m *s + m* np !(1.d0-2*n**2)*p
	! z, zx
	h(3,3) = sq3 *n2 *l *s + l* np ! (1.d0-2*n**2)*p
	!..................................................

	!..................................................	
	! with eg
	!..................................................
	if(norbtm==5) then
	 lm = (l2 - m2);
	 lm2 = 0.5d0*sq3*lm*s
	 nlm2 = (n2 - 0.5d0*(l2 + m2)) *s
	 ! x, x^2-y^2 
	 h(1,4) = l* (lm2 +(1.d0-lm)*p)
	 ! y, x^2-y^2 
	 h(2,4) = m* (lm2 -(1.d0+lm)*p)
	 ! z, x^2-y^2 
	 h(3,4) = n* (lm2 -lm*p)

	 ! x, 3z^2-r^2
	 h(1,5) = l*(nlm2 - sq3*n2*p)
	 ! y, 3z^2-r^2
	 h(2,5) = m*(nlm2 - sq3*n2*p)
	 ! z, 3z^2-r^2
	 h(3,5) = n*(nlm2 - sq3*(l2+m2)*p)
	endif
	!..................................................

	return
	end subroutine slatkospd
	!.....................................................




	!.....................................................
	! Hamiltonian matrix elements between d-d orbitals using SK method
	subroutine slatkosdd(r, sk, h)
	implicit none
	double precision, dimension(3), intent(in) :: r
	double precision, dimension(3), intent(in) :: sk
	double precision, dimension(norbtm,norbtm), intent(out) :: h
	! local
	double precision :: rr
	double precision :: l,m,n ! direction cosines
	double precision :: lmn, s, p, d, sq3
	double precision :: l2,m2,n2, lmp,lm,lm2	

	h(:,:) = 0.0d0;
	s = sk(1); p = sk(2); d = sk(3); ! s=sigma_dd, p=pi_dd, d = delta_dd
	sq3 = dsqrt(3.0d0);

	rr = norm2(r);
	! direction cosines; r is from first to the second atom
	l = r(1)/rr;
	m = r(2)/rr;
	n = r(3)/rr;

	!write(*,*)'sk, rr = ', sk, rr


	!..................................................	
	! (t2g,t2g)
	!..................................................
	l2 = l**2; m2 = m**2; n2 = n**2;
	! xy,xy
	h(1,1) = 3*l2*m2 *s + (l2 + m2 - 4*l2*m2)*p + (n2 + l2*m2)*d
	! xy, yz
	h(1,2) = 3*l*m2*n*s + l*n*(1.d0-4*m2)*p + l*n*(m2-1.d0)*d
	! xy, zx
	h(1,3) = 3*l2*m*n*s + m*n*(1.d0-4*l2)*p + m*n*(l2-1.d0)*d

	! perm:   x,y,z===>y,z,x    ::    l,m,n ==> m,n,l
	! yz, xy : perm of (xy, zx)
	h(2,1) = 3*m2*n*l*s + n*l*(1.d0-4*m2)*p + n*l*(m2-1.d0)*d
	! yz, yz : perm of (xy, xy)
	h(2,2) = 3*m2*n2 *s + (m2 + n2 - 4*m2*n2)*p + (l2 + m2*n2)*d
	! yz, zx : perm of (xy, yz)
	h(2,3) = 3*m*n2*l*s + m*l*(1.d0-4*n2)*p + m*l*(n2-1.d0)*d

	! zx,xy : perm of (yz,zx)
	h(3,1) = 3*n*l2*m*s + n*m*(1.d0-4*l2)*p + n*m*(l2-1.d0)*d
	! zx,yz : perm of (yz,xy)
	h(3,2) = 3*n2*l*m*s + l*m*(1.d0-4*n2)*p + l*m*(n2-1.d0)*d
	! zx,zx : perm of (yz,yz)
	h(3,3) = 3*n2*l2 *s + (n2 + l2 - 4*n2*l2)*p + (m2 + n2*l2)*d
	!..................................................

	! could enclose the part below in an if statement
	! or make two seperate routines, 
	! one for norbtm=5 (t2g only), one for norbtm=5 (full d-space)
	! ******** proceed only if norbtm=5 ***********************
	if(norbtm==3) return
	! *********************************************************
	
	!..................................................	
	! (t2g,eg)
	!..................................................
	lm = (l2 - m2);
	lm2 = 1.5d0*lm *s;
	! xy, x^2-y^2
	h(1,4) = l*m*lm2 - 2*l*m*lm*p + 0.5d0*l*m*lm*d
	! yz, x^2-y^2
	h(2,4) = m*n*lm2 - m*n*(1.d0+2*lm)*p + m*n*(1.d0+0.5d0*lm)*d
	! zx, x^2-y^2
	h(3,4) = n*l*(lm2 + (1.d0-2*lm)*p + (1.d0+0.5d0*lm)*d)
	!..................................................
	lm = l2 + m2;
	lm2 = sq3*(n**2 -0.5d0*lm)*s;
	! xy, 3z^2-r^2
	h(1,5) = l*m*lm2 - 2*sq3*l*m*n2*p + 0.5d0*sq3*l*m*(1.d0+n2)*d
	! yz, 3z^2-r^2
	h(2,5) = m*n*lm2 + sq3*m*n*(lm - n2)*p - 0.5d0*sq3*m*n*lm*d ! last term: l+m^2  or l^2+m^2 ?? check
	! zx, 3z^2-r^2
	h(3,5) = l*n*lm2 + sq3*l*n*(lm - n2)*p - 0.5d0*sq3*l*n*lm*d
	!..................................................


	!..................................................	
	! (eg,eg)
	!..................................................
	lmp = (l2 + m2);
	lm = (l2 - m2);
	lm2 =lm**2;	

	! x^2-y^2, x^2-y^2
	h(4,4) = 0.75d0*lm2*s + (lmp - lm2)*p +(n2 + 0.25*lm2)*d
	! x^2-y^2, 3z^2-r^2
	h(4,5) = 0.5d0*sq3*lm*(n2 - 0.5d0*lmp)*s -sq3*n2*lm*p + 
     .                            0.25d0*sq3*(1.d0+n2)*lm*d
	! 3z^2-r^2, 3z^2-r^2
	h(5,5) = (n2 - 0.5d0*lmp)**2*s + 3*n2*lmp*p + 0.75d0*lmp**2*d
	!..................................................	



	! the rest of the cases from swapping atom 1,2: 
	! r ===> -r: l,m,n = -l,-m,-n
	! even in l,m,n, so equal.
	h(4,1) = h(1,4);
	h(4,2) = h(2,4);
	h(4,3) = h(3,4);
	h(5,1) = h(1,5);
	h(5,2) = h(2,5);
	h(5,3) = h(3,5);
	h(5,4) = h(4,5); 


	!write(*,*)' dd: h = ', h

	return
	end 	subroutine slatkosdd
	!.....................................................



	end 	module skoster
