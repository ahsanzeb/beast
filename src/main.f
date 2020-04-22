
	program perovskite
	use modmain
	implicit none
	integer :: il, io, i
	double precision:: phi0
	
	
	write(*,'(a)') "Give Number of layers & phi and hit enter:"
	read(*,*) nlayers, phi0

	if(phi0 > 45.0) then
	 write(*,'("O atoms have to cross to make this big rotation!")')
	 stop 
	else
	 write(*,'(a)')
     .  'Octraherda size will be rescaled to keep "a" fixed.'
	endif


	!nlayers = 2;
	noctl = 2;
	noct = noctl*nlayers;
	a = 7.0d0;
	nspecies = 3;
	
	
	allocate(oct(nlayers,noctl))
	allocate(phi(nlayers))
	! phi 
	phi = 0.0d0;

	! read/set phi:
	do il=1,nlayers
	 if(mod(il,2)==0) then
	  phi(il) = phi0*pi/180.0d0;
	 else
	  phi(il) = -phi0*pi/180.0d0;
	 endif
	 oct(il,1)%phi = +phi(il)
	 oct(il,2)%phi = -phi(il) 
	end do


	! set the basic cubic structure
	! sqrt(2) x sqrt(2) x nlayers cell
	! square cell rotated by 45 wr.r.t x,y
	! lattice vectors
	avec(1,:) = (/1.0d0, -1.0d0, 0.0d0/)*a
	avec(2,:) = (/1.0d0, +1.0d0, 0.0d0/)*a
	avec(3,:) = (/0.0d0,  0.0d0, 1.0d0/)*a*nlayers

	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)

	!ainv = matmul(avec,ainv)
	!write(*,'(3f10.3)') ainv(1,:)
	!write(*,'(3f10.3)') ainv(2,:)
	!write(*,'(3f10.3)') ainv(3,:)

	! atomic positions in cartesian
	! set pos of B
	do il=1,nlayers
	 oct(il,1)%rb = (/0.5d0,-0.5d0,1.d0*(il-1)/)*a;
	 oct(il,2)%rb = (/0.5d0,+0.5d0,1.d0*(il-1)/)*a;
	enddo

	! oxygens
	do il=1,nlayers
	 do io=1, noctl
	  	oct(il,io)%xo(1,:) = (/0.5d0,0.0d0,0.0d0/)*a;
	  	oct(il,io)%xo(2,:) = (/0.0d0,0.5d0,0.0d0/)*a;
	  	oct(il,io)%xo(3,:) = (/0.0d0,0.0d0,0.5d0/)*a;
	  oct(il,io)%lo = 0.5d0 * a;
	 end do
	end do
	
	! rotate and set abs pos of oxygen
	do il=1,nlayers
	 do io=1, noctl
	  call rotate(il,io)
	 end do
	end do

	! write GEOMETRY.OUT for visualisation with VESTA
	call writegeom()



	stop
	end 	program perovskite
