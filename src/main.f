
	program perovskite
	use modmain
	implicit none
	integer :: il, io, i
	
	
	write(*,'(a)') "Give Number of layers and hit enter:"
	read(*,*) nlayers



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
	  phi(il) = 10*pi/180.0d0;
	 else
	  phi(il) = -10*pi/180.0d0;
	 endif
	 oct(il,1)%phi = +phi(il)
	 oct(il,2)%phi = -phi(il) 
	end do

	write(*,'("Octraherda rotated by an angle of ",f10.5)') phi(1)

	! set the basic cubic structure
	! lattice vectors
	avec(1,:) = (/a*noctl, 0.0d0, 0.0d0/)
	avec(2,:) = (/a*1,a*1, 0.0d0/)
	avec(3,:) = (/0.0d0, 0.0d0, a*nlayers/)

	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)

	!ainv = matmul(avec,ainv)


	write(*,'(3f10.3)') ainv(1,:)
	write(*,'(3f10.3)') ainv(2,:)
	write(*,'(3f10.3)') ainv(3,:)

	! atomic positions in cartesian
	! set pos of B
	do il=1,nlayers
	 oct(il,1)%rb = (/0.d0,0.d0,1.d0*(il-1)/)*a;
	 oct(il,2)%rb = (/1.d0,0.d0,1.d0*(il-1)/)*a;
	enddo

	! oxygens
	do il=1,nlayers
	 do io=1, noctl
	  	oct(il,io)%xo(1,:) = (/0.5d0,0.0d0,0.0d0/)*a;
	  	oct(il,io)%xo(2,:) = (/0.0d0,0.5d0,0.0d0/)*a;
	  	oct(il,io)%xo(3,:) = (/0.0d0,0.0d0,0.5d0/)*a;
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
