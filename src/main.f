
	program perovskite
	use modmain
	implicit none
	integer :: il, io, i
	
	


	nlayers = 2;
	noctl = 2;
	noct = noctl*nlayers;
	a = 2*7.0d0;
	nspecies = 3;
	
	
	allocate(oct(nlayers,noctl))
	allocate(phi(nlayers))

	phi =10.0d0*pi/180.0d0;

	! read/set phi
	do il=1,nlayers
	 oct(il,1)%phi = +phi(il)
	 oct(il,2)%phi = -phi(il) 
	end do

	! set the basic cubic structure
	! set pos of B
	oct(1,1)%rb = (/0.0,0.0,0.0/);
	oct(1,2)%rb = (/0.5,0.0,0.0/);

	oct(2,1)%rb = (/0.0,0.0,0.5/);
	oct(2,2)%rb = (/0.5,0.0,0.5/);
	! oxygens
	do il=1,nlayers
	 do io=1, noctl
	  	oct(il,io)%xo(1,:) = (/0.25,0.0,0.0/);
	  	oct(il,io)%xo(2,:) = (/-0.25,0.5,0.0/);
	  	oct(il,io)%xo(3,:) = (/0.0,0.0, 0.25/);
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
