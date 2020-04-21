
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
	 write(*,'("Octraherda size will be rescaled to keep a fixed.")')
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
	! lattice vectors
	avec(1,:) = (/a*noctl, 0.0d0, 0.0d0/)
	avec(2,:) = (/a*1,a*1, 0.0d0/)
	avec(3,:) = (/0.0d0, 0.0d0, a*nlayers/)

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
	 oct(il,1)%rb = (/0.d0,0.d0,1.d0*(il-1)/)*a;
	 oct(il,2)%rb = (/1.d0,0.d0,1.d0*(il-1)/)*a;
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





	! TM/O atoms and their neighbours
	allocate(tm(nlayers,noctl))
	allocate(ox(nlayers,noctl*3))

	! frist nns of TM B atom are 6 O atoms
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always

	  tm(il,io)%r = oct(il,io)%rb
	  allocate(tm(il,io)%nn1(6)) ! O
	 ! allocate(tm(il,io)%nn2(6)) ! B/TM

	  if(io==1) then
	  	 tm(il,io)%ia = (il-1)*noctl*8 + io
	   ! a1,a2,a3 denote the lattice vectors
	   ! O from nearby cells. 
	   tm(il,io)%nn1(4)%ia = tm(il,io)%ia + 5 ! 6, ! O_x of B_x in -a1
	   tm(il,io)%nn1(4)%r = oct(il,2)%ro(1,:) - avec(1,:)
	   tm(il,io)%nn1(5)%ia = tm(il,io)%ia + 6 ! 7, ! O_y of B_x in -a2
	   tm(il,io)%nn1(5)%r = oct(il,2)%ro(2,:) - avec(2,:)
	  	else
	   ! second TM/B atom in the layer
	  	 tm(il,io)%ia = tm(il,1)%ia + 4; 
	   ! a1,a2,a3 denote the lattice vectors
	   ! O from nearby cells. 
	   tm(il,io)%nn1(4)%ia = tm(il,1)%ia + 1 ! 2, ! O_x of B_-x
	   tm(il,io)%nn1(4)%r = oct(il,1)%ro(1,:)
	   tm(il,io)%nn1(5)%ia = tm(il,1)%ia + 2 ! 3, ! O_y of B_-x in -a2+a1
	   tm(il,io)%nn1(5)%r = oct(il,1)%ro(2,:) -avec(2,:)+avec(1,:)
	  	endif

	  do i=1,3 ! O belonging to the unit cell
	   tm(il,io)%nn1(i)%ia = tm(il,io)%ia + ioo;
	   tm(il,io)%nn1(i)%r = oct(il,io)%ro(i,:)
	  end do
	  ! 3rd O could belong to the unit cell or maybe in a cell on below.
	  if(il==1) then ! O_z from the cell below, i.e., -a3
	   tm(il,io)%nn1(6)%ia = nlayers*nactl*8 + tm(il,io)%ia + 3
	   tm(il,io)%nn1(6)%r = oct(il,io)%ro(3,:) - avec(3,:)
	  else ! belongs to the unit cell
	   tm(il,io)%nn1(6)%ia = tm(il-1,io)%ia + 3 ! lower layer O_z
	   tm(il,io)%nn1(6)%r = oct(il-1,io)%ro(3,:)
	  endif

	 end do
	end do




	! second nns of TM B atom are 6 TM B (or B'?) atoms
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always

	  !tm(il,io)%r = oct(il,io)%rb
	  !allocate(tm(il,io)%nn1(6)) ! O
	  allocate(tm(il,io)%nn2(6)) ! B/TM

	  if(io==1) then
	  	 tm(il,io)%ia = (il-1)*noctl*8 + io !  1st atom in a lyer
	   ! TM belonging to the unit cell or otherwise, all B' atoms
	   do i=1,6
	    tm(il,io)%nn2(i)%ia = tm(il,io)%ia + 4; ! 5th atom in a lyer
	   end do
		else ! io=2
	  	 tm(il,io)%ia = (il-1)*noctl*8 + 5;! 5th atom in a lyer
	   do i=1,6
	    tm(il,io)%nn2(i)%ia = tm(il,io)%ia - 4; ! 1st atom in a lyer
	   end do
	  endif
	  ! same layer
	  tm(il,io)%nn2(1)%r = oct(il,io)%rb + (1.0d0,0.0d0,0.0d0)*a
	  tm(il,io)%nn2(2)%r = oct(il,io)%rb - (1.0d0,0.0d0,0.0d0)*a
	  tm(il,io)%nn2(3)%r = oct(il,io)%rb + (0.0d0,1.0d0,0.0d0)*a
	  tm(il,io)%nn2(4)%r = oct(il,io)%rb - (0.0d0,1.0d0,0.0d0)*a
	  ! top/bottom layers
	  tm(il,io)%nn2(5)%r = oct(il,io)%rb + (0.0d0,0.0d0,1.0d0)*a
	  tm(il,io)%nn2(6)%r = oct(il,io)%rb - (0.0d0,0.0d0,1.0d0)*a

	 end do
	end do

























	stop
	end 	program perovskite
