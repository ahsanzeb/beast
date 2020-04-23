
	module modmain
	implicit none

	integer :: nlayers, natomsl, natoms
	integer :: noctl, noct
	integer :: nspecies, nsptm !
	integer :: norbtm, norbo ! number of spatial orbitals on a TM/O atom
	integer :: norbtms, norbos ! number of spin-space orbitals on a TM/O atom

	integer, allocatable, dimension(:,:) :: atom2orb
	double precision :: a ! lattice constant of a single formula unit cubic cell
	double precision, parameter :: pi = 3.141592653589793d0

	! SK parameters
	double precision, allocatable, dimension(:,:):: skbo ! spd, ppd
	double precision, allocatable, dimension(:,:,:):: skbb ! sdd, pdd, ddd
	double precision, dimension(2):: skoo ! spp, ppp
	
	type :: octahedra
		!integer :: ntot
		double precision :: phi, lo,lor,lort
		double precision, dimension(3) :: rb ! position of B atom, absolute position in a unit cell.
		double precision, dimension(3,3) :: xo ! relative to xb, position of oxygen atoms in cubic structure
		double precision, dimension(3,3) :: xor ! relative to xb, position of oxygen atoms after rotation
		double precision, dimension(3,3) :: ro ! absolute position in a unit cell, position of oxygen atoms after rotation
	end type octahedra

	type(octahedra), allocatable, dimension(:,:) :: oct

	double precision, allocatable, dimension(:) :: phi 

	! lattice vectors
	double precision, dimension(3,3) :: avec, ainv

	type :: nneighbours
	 integer :: ia ! atom index, or index of equalent atom inside the unit cell (if this atom is outside the unit cell)
	 integer :: i, j ! orbital range indices: start, end; or equiv atom's orb ranges
	 double precision, dimension(3) :: r,dr ! position, true position of this nn.
	double precision, allocatable, dimension(:,:) :: h ! hamiltonian matrix elements between orbitals of two beighbours
	end type nneighbours

	type :: atoms
	 character(len=2) :: label ! 'Ir', 'O', etc...
	 integer :: ia ! atom number, index in full list of atoms
	 integer :: i, j ! orbital range indices: start, end
	 !integer :: typ ! orbital type, 1= p, 2=d, only one type allowed.
	 double precision, dimension(3) :: r ! position
	 type(nneighbours), allocatable, dimension(:) :: nn1, nn2 ! first and second nns, inside the unit cell or outside it.
	end type atoms

	type(atoms), allocatable, dimension(:,:) :: tm ! TM
	type(atoms), allocatable, dimension(:,:,:) :: ox ! Oxygen


	contains

	!..............................................................

	subroutine writegeom()
	implicit none
	integer :: fnum, i,io,il


		
	open(fnum,file='GEOMETRY.OUT',form='FORMATTED')
	write(fnum,*)
	write(fnum,'("scale")')
	write(fnum,'(" 1.0")')
	write(fnum,*)
	write(fnum,'("avec")')
	write(fnum,'(3G18.10)') avec(1,:)
	write(fnum,'(3G18.10)') avec(2,:)
	write(fnum,'(3G18.10)') avec(3,:)
	write(fnum,*)
	write(fnum,'("atoms")')
	write(fnum,'(I4,T40," : nspecies")') nspecies

	! all oxygen
	write(fnum,'("''","O ","''",T40," : spfname")') 
	write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")') 3*noct
	do il=1, nlayers
	 do io =1, noctl
	 	do i=1,3
		 write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%ro(i,1:3),0.,0.,0.
	  end do
	 end do
	end do


	! assume different B atoms in each layer 
	if(mod(noctl,2)==0) then ! even number of octaherda in a layer
	! assume two species
	! first species:
	write(fnum,'("''",A,"''",T40," : spfname")') "Ir"
	write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")') noct/2
	do il=1, nlayers
	 do io =1, noctl,2
		write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%rb, 0.,0.,0.
	 end do
	end do
	! second species:
	write(fnum,'("''",A,"''",T40," : spfname")') "Ti"
	write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")') noct/2
	do il=1, nlayers
	 do io =2, noctl,2
		write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%rb, 0.,0.,0.
	 end do
	end do
	else ! odd
	! one species only
	write(fnum,'("''",A,"''",T40," : spfname")') "Ir"
	write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")')noct
	do il=1, nlayers
	 do io =1, noctl
	  write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%rb, 0.,0.,0.
	 end do
	end do
	endif

	close(fnum)

	return
	end 	subroutine writegeom
	!..............................................................
	
	!..............................................................
	! rotate octahedra
	subroutine rotate(il,io)
	implicit none
	integer, intent(in) :: il,io
	double precision :: phi, dl, sphi, cphi
	integer :: i,j
	double precision :: v(3)

	phi = oct(il,io)%phi;
	sphi = dsin(phi);
	cphi = dcos(phi);
	!write(*,'("phi, sphi,cphi = ",3f10.5)') phi, sphi, cphi

	do i = 1,2 ! O atoms
	 oct(il,io)%xor(i,1) =  
     .  cphi*oct(il,io)%xo(i,1) + sphi*oct(il,io)%xo(i,2)
	 oct(il,io)%xor(i,2) =  
     . -sphi*oct(il,io)%xo(i,1) + cphi*oct(il,io)%xo(i,2)
	 oct(il,io)%xor(i,3) = 	 oct(il,io)%xo(i,3); ! z comp

	! rescale distance to keep B-O-B shift in angle with boht B the same.
	oct(il,io)%lor = oct(il,io)%lo * 1.0d0/cphi;
	dl = oct(il,io)%lor - oct(il,io)%lo
	! rescale  x,y comp
	if(i==1) then ! atom along x
	 oct(il,io)%xor(i,1) = oct(il,io)%xor(i,1)  + dl * cphi;
	 oct(il,io)%xor(i,2) = oct(il,io)%xor(i,2)  + dl * sphi;
	elseif(i==2)then ! i==2, atom along y
	 oct(il,io)%xor(i,1) = oct(il,io)%xor(i,1)  - dl * sphi;
	 oct(il,io)%xor(i,2) = oct(il,io)%xor(i,2)  + dl * cphi;
	endif

	end do
	! atom on z axis: not rotated
	oct(il,io)%xor(3,:) = oct(il,io)%xo(3,:)

	! set abs value of oxygen position after rotation.
	do i = 1,3
	 oct(il,io)%ro(i,:) =	oct(il,io)%rb(:) + oct(il,io)%xor(i,:);
	end do
	
	! ro cartesian to fractional
	do i=1,3
	 !write(*,'(a)') '-----------------'		
	 !write(*,'(3f10.4)') oct(il,io)%ro(i,:)
	 call r3mv(transpose(ainv),oct(il,io)%ro(i,:),v)
	 oct(il,io)%ro(i,:) = v
	 !write(*,'(3f10.4)') oct(il,io)%ro(i,:)
	end do
	! central B atom
	call r3mv(transpose(ainv),oct(il,io)%rb,v)
	oct(il,io)%rb = v

	return
	end 	subroutine rotate
	!..............................................................
	! for orthogonal avec, we can used transpose, but have to use inverse for general cases.
	subroutine transform(s, v, v2)
	implicit none
	integer, intent(in) :: s
	double precision, dimension(3), intent(in) :: v
	double precision, dimension(3), intent(out) :: v2

	if (s==+1) then
		v2 = matmul(avec,v)
	elseif(s==-1)then
		v2 = matmul(transpose(avec),v)
	else
		stop "Error(transform): wrong input s..."
	endif

	return
	end 	subroutine transform

	!..............................................
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
	! copeid from elk-6.2.8
	subroutine r3minv(a,b)
	implicit none
	real(8), intent(in) :: a(3,3)
	real(8), intent(out) :: b(3,3)
	real(8) t1
	t1=a(1,2)*a(2,3)*a(3,1)-a(1,3)*a(2,2)*a(3,1)+a(1,3)*a(2,1)*a(3,2)
     .-a(1,1)*a(2,3)*a(3,2)+a(1,1)*a(2,2)*a(3,3)-a(1,2)*a(2,1)*a(3,3)
	if (abs(t1).lt.1.d-40) then
	 write(*,*)
	 write(*,'("Error(r3minv): singular matrix")')
	 write(*,*)
	 stop
	end if
	t1=1.d0/t1
	b(1,1)=t1*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
	b(2,1)=t1*(a(2,3)*a(3,1)-a(2,1)*a(3,3))
	b(3,1)=t1*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
	b(1,2)=t1*(a(1,3)*a(3,2)-a(1,2)*a(3,3))
	b(2,2)=t1*(a(1,1)*a(3,3)-a(1,3)*a(3,1))
	b(3,2)=t1*(a(1,2)*a(3,1)-a(1,1)*a(3,2))
	b(1,3)=t1*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
	b(2,3)=t1*(a(1,3)*a(2,1)-a(1,1)*a(2,3))
	b(3,3)=t1*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
	return
	end subroutine

!..........................................................
! For first nearest neighbours of TM, 6 O atoms:
! sets indices (true or equiv in the cell) and positions (true)
!..........................................................
	subroutine settmnn1()
	implicit none
	double precision, dimension(3) :: z, a1, a2, a3
	integer :: il, io,jo,i,j
	
	 !z = (/0.d0,0.d0,1.d0/)*a;
	 a1 = avec(1,:);
	 a2 = avec(2,:);
	 a3 = avec(3,:);

	allocate(tm(nlayers,noctl))
	!.....................................................
	! frist nns of TM B atom are 6 O atoms
	!.....................................................
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always

	  tm(il,io)%r = oct(il,io)%rb
	  allocate(tm(il,io)%nn1(6)) ! O

	  if(io==1) then
	  	 tm(il,io)%ia = (il-1)*noctl*8 + io
	   ! a1,a2,a3 denote the lattice vectors
	   ! O from nearby cells. 
	   tm(il,io)%nn1(4)%ia = tm(il,io)%ia + 5 ! 6, ! O_x of B2 in -a2
	   tm(il,io)%nn1(4)%r = oct(il,2)%ro(1,:) - a2
	   tm(il,io)%nn1(5)%ia = tm(il,io)%ia + 6 ! 7, ! O_y of B2 in -a2+a1
	   tm(il,io)%nn1(5)%r = oct(il,2)%ro(2,:) -a2+a1
	  	else ! io=2
	   ! second TM/B atom in the layer
	  	 tm(il,io)%ia = tm(il,1)%ia + 4; 
	   ! O_y of B1 in the cell
	   tm(il,io)%nn1(4)%ia = tm(il,1)%ia + 2 ! 3, O_y of B1
	   tm(il,io)%nn1(4)%r = oct(il,1)%ro(2,:)
	   ! O from nearby cell. 
	   tm(il,io)%nn1(4)%ia = tm(il,1)%ia + 1 ! 2, ! O_x of B1 in -a1
	   tm(il,io)%nn1(4)%r = oct(il,1)%ro(1,:) - a1
	  	endif

	  do i=1,3 ! O belonging to the unit cell
	   tm(il,io)%nn1(i)%ia = tm(il,io)%ia + i;
	   tm(il,io)%nn1(i)%r = oct(il,io)%ro(i,:)
	  end do
	  ! 3rd O could belong to the unit cell or maybe in a cell on below.
	  if(il==1) then ! O_z from the cell below, i.e., -a3
	   tm(il,io)%nn1(6)%ia = nlayers*noctl*8 + tm(il,io)%ia + 3
	   tm(il,io)%nn1(6)%r = oct(il,io)%ro(3,:) - a3
	  else ! belongs to the unit cell, always, even for il=nlayers
	   tm(il,io)%nn1(6)%ia = tm(il-1,io)%ia + 3 ! lower layer O_z
	   tm(il,io)%nn1(6)%r = oct(il-1,io)%ro(3,:)
	  endif

	 end do
	end do

	return
	end subroutine settmnn1
!..........................................................
! For the second nearest neighbours of TM, 6 TM atoms:
! sets indices (true or equiv in the cell) and positions (true).
!..........................................................
	subroutine settmnn2()
	implicit none
	double precision, dimension(3) :: x, y, z
	integer :: il, io,jo,i,j

	x = (/1.0d0,0.0d0,0.0d0/)*a;
	y = (/0.0d0,1.0d0,0.0d0/)*a
	z = (/0.0d0,0.0d0,1.0d0/)*a;

	!allocate(tm(nlayers,noctl))
	!.....................................................
	! second nns of TM B atom are 6 TM B (or B'?) atoms
	!.....................................................
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
	  tm(il,io)%nn2(1)%r = oct(il,io)%rb + x
	  tm(il,io)%nn2(2)%r = oct(il,io)%rb - x
	  tm(il,io)%nn2(3)%r = oct(il,io)%rb + y
	  tm(il,io)%nn2(4)%r = oct(il,io)%rb - y
	  ! top/bottom layers
	  tm(il,io)%nn2(5)%r = oct(il,io)%rb + z
	  tm(il,io)%nn2(6)%r = oct(il,io)%rb - z
	 end do
	end do

	return
	end subroutine settmnn2
!..........................................................


!..........................................................
! For the first nearest neighbours of O, 2 TM atoms:
! sets indices (true or equiv in the cell) and positions (true).
!..........................................................
	subroutine setonn1()
	implicit none
	double precision, dimension(3) :: z, a1, a2, a3
	integer :: il, io,jo,i,j
	
	 z = (/0.d0,0.d0,1.d0/)*a;
	 a1 = avec(1,:);
	 a2 = avec(2,:);
	 !a3 = avec(3,:);

	allocate(ox(nlayers,noctl,3))
	!.....................................................
	! frist nns of O atoms are 2 TM atoms
	!.....................................................
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always
	  do i=1,3
	   ox(il,io,i)%ia = tm(il,io)%ia + i ! O index
	   ox(il,io,i)%r = oct(il,io)%ro(i,:)
	   allocate(ox(il,io,i)%nn1(2)) ! 2 TM atoms
	   !...................................
	   ! set indices of the two 1st nns TM
	   !...................................
	   ! TM of the same octahedra:
	   ox(il,io,i)%nn1(1)%ia = tm(il,io)%ia 
	   ! TM of the other octahedra:
	   if(i<3) then ! same layer   
	    if(io==1) then
	     ox(il,io,i)%nn1(2)%ia = tm(il,io)%ia + 4;
	    else ! io=2
	     ox(il,io,i)%nn1(2)%ia = tm(il,io)%ia - 4;
	    endif
	   else ! i=3: TM in the top layer
	    ox(il,io,i)%nn1(1)%ia = tm(il,io)%ia ! TM of the same octahedra
	    if(il < nlayers) then ! second nn TM is inside the unit cell
	     ox(il,io,i)%nn1(2)%ia = tm(il,io)%ia + 8
	    else ! il=nlayers, second nn TM  of O_z is outside the unit cell
	     ox(il,io,i)%nn1(2)%ia = tm(1,io)%ia ! 1st layer's periodic image
	    endif
	   endif ! i<3
	   !...................................

	  end do ! i
	 end do ! io
	   !...................................
	   ! set positions of the two 1st nns TM
	   !...................................
	   ! TM of the same octahedra:
	   ! O_x of B1
	   ox(il,1,1)%nn1(1)%r = oct(il,1)%rb
	   ox(il,1,1)%nn1(2)%r = oct(il,2)%rb + a1
	   ! O_y of B1
	   ox(il,1,2)%nn1(1)%r = oct(il,1)%rb
	   ox(il,1,2)%nn1(2)%r = oct(il,2)%rb
	   ! O_z of B1
	   ox(il,1,3)%nn1(1)%r = oct(il,1)%rb
	   ox(il,1,3)%nn1(2)%r = oct(il,1)%rb + z


	   ! O_x of B2
	   ox(il,2,1)%nn1(1)%r = oct(il,2)%rb
	   ox(il,2,1)%nn1(2)%r = oct(il,1)%rb + a2
	   ! O_y of B2
	   ox(il,2,2)%nn1(1)%r = oct(il,2)%rb
	   ox(il,2,2)%nn1(2)%r = oct(il,1)%rb +a2-a1
	   ! O_z of B2
	   ox(il,2,3)%nn1(1)%r = oct(il,2)%rb
	   ox(il,2,3)%nn1(2)%r = oct(il,2)%rb + z

	end do ! il
	!.....................................................

	return
	end subroutine setonn1

!..........................................................
! For the second nearest neighbours of O, 8 O atoms:
! sets indices (true or equiv in the cell) and positions (true).
!..........................................................
	subroutine setonn2()
	implicit none
	double precision, dimension(3) :: z, a1, a2, a3
	integer :: il, io,jo,i,j
	
	 z = (/0.d0,0.d0,1.d0/)*a;
	 a1 = avec(1,:);
	 a2 = avec(2,:);
	 a3 = avec(3,:);

	!allocate(ox(nlayers,noctl,3))
	!.....................................................
	! frist nns of O atoms are 2 TM atoms
	!.....................................................

	! il=1 and il=nlayers have O that have 2nd nns in layers outside the cell
	! treat them seperately?
	do il=1,nlayers

	 ! allocate nn2(:)
	 do io=1,noctl ! noctl = 2 always
	  do i=1,3
	   allocate(ox(il,io,i)%nn2(8)) ! 8 O atoms
	  end do
	 end do ! io


	 !..................................................
	 ! nns of O with B1 octahedra
	 !..................................................

	 !..................................................
	 ! O_x of B1: ox(il,io,i)%ia = 2
	 io=1; i=1;
	 jo = 2;
	 ox(il,io,i)%nn1(1)%ia = ox(il,io,i)%ia + 1 ! 3, same octa
	 ox(il,io,i)%nn1(2)%ia = ox(il,io,i)%ia + 5 ! 7, -a1-a2
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia + 2 ! 4, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn1(4)%ia = ox(nlayers,io,i)%ia + 2 ! 4' 
	  ox(il,io,i)%nn1(8)%ia = ox(nlayers,io,i)%ia + 6 ! 8' (+a1)
	 else ! inside unit cell
	  ox(il,io,i)%nn1(4)%ia = ox(il-1,io,i)%ia + 2 ! 4'
	  ox(il,io,i)%nn1(8)%ia = ox(il-1,io,i)%ia + 6 ! 8' (+a1)
	 endif
	 ox(il,io,i)%nn1(5)%ia = ox(il,io,i)%ia + 1 ! 3, +a1
	 ox(il,io,i)%nn1(6)%ia = ox(il,io,i)%ia + 5 ! 7, +a1
	 ox(il,io,i)%nn1(7)%ia = ox(il,io,i)%ia + 6 ! 8, +a1

	 ox(il,io,i)%nn1(1)%r = oct(il,io)%ro(2,:) ! 3, same octa
	 ox(il,io,i)%nn1(2)%r = oct(il,jo)%ro(2,:) -a1-a2 ! 7, -a1-a2
	 ox(il,io,i)%nn1(3)%r = oct(il,io)%ro(3,:)  ! 4, same cell
	 ox(il,io,i)%nn1(4)%r = oct(il,io)%ro(3,:) -z ! 4', layer below
	 ox(il,io,i)%nn1(5)%r = oct(il,io)%ro(2,:) +a1 ! 3, +a1
	 ox(il,io,i)%nn1(6)%r = oct(il,jo)%ro(2,:) +a1 ! 7, +a1
	 ox(il,io,i)%nn1(7)%r = oct(il,jo)%ro(3,:) +a1 ! 8, +a1
	 ox(il,io,i)%nn1(8)%r = oct(il,jo)%ro(3,:) -z +a1 ! 8', +a1, layer below
	 !..................................................
	 ! O_y of B1: ox(il,io,i)%ia = 3
	 io=1; i=2;
	 jo = 2;
	 ox(il,io,i)%nn1(1)%ia = ox(il,io,i)%ia - 1 ! 2, same octa
	 ox(il,io,i)%nn1(2)%ia = ox(il,io,i)%ia + 3 ! 6, -a2
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia + 1 ! 4, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn1(4)%ia = ox(nlayers,io,i)%ia + 1 ! 4' 
	  ox(il,io,i)%nn1(8)%ia = ox(nlayers,io,i)%ia + 5 ! 8' 
	 else ! inside unit cell
	  ox(il,io,i)%nn1(4)%ia = ox(il-1,io,i)%ia + 1 ! 4'
	  ox(il,io,i)%nn1(8)%ia = ox(il-1,io,i)%ia + 5 ! 8'
	 endif
	 ox(il,io,i)%nn1(5)%ia = ox(il,io,i)%ia + 3 ! 6, same cell
	 ox(il,io,i)%nn1(6)%ia = ox(il,io,i)%ia - 1 ! 2, -a2
	 ox(il,io,i)%nn1(7)%ia = ox(il,io,i)%ia + 5 ! 8, same cell

	 ox(il,io,i)%nn1(1)%r = oct(il,io)%ro(1,:) ! 2, same octa
	 ox(il,io,i)%nn1(2)%r = oct(il,jo)%ro(1,:) -a2 ! 6, -a2
	 ox(il,io,i)%nn1(3)%r = oct(il,io)%ro(3,:)  ! 4, same cell
	 ox(il,io,i)%nn1(4)%r = oct(il,io)%ro(3,:) -z ! 4', layer below
	 ox(il,io,i)%nn1(5)%r = oct(il,io)%ro(1,:) ! 6, same cell
	 ox(il,io,i)%nn1(6)%r = oct(il,jo)%ro(1,:) -a2 ! 2, -a2
	 ox(il,io,i)%nn1(7)%r = oct(il,jo)%ro(3,:) ! 8, same cell
	 ox(il,io,i)%nn1(8)%r = oct(il,jo)%ro(3,:) -z ! 8',layer below
	 !..................................................
	 ! O_z of B1: ox(il,io,i)%ia = 4
	 io=1; i=3;
	 jo = 2;
	 ! same layer
	 ox(il,io,i)%nn1(1)%ia = ox(il,io,i)%ia - 2 ! 2, same octa
	 ox(il,io,i)%nn1(2)%ia = ox(il,io,i)%ia - 1 ! 3, same octa
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia + 2 ! 6, -a2
	 ox(il,io,i)%nn1(4)%ia = ox(il,io,i)%ia + 3 ! 7, -a2+a1
	 ! layer above
	 if(il<nlayers) then ! the layer above is inside the cell
	  ox(il,io,i)%nn1(5)%ia = ox(il,io,i)%nn1(1)%ia + 8! 2, +z
	  ox(il,io,i)%nn1(6)%ia = ox(il,io,i)%nn1(2)%ia + 8! 3, +z
	  ox(il,io,i)%nn1(7)%ia = ox(il,io,i)%nn1(3)%ia + 8! 6, +z -a2
	  ox(il,io,i)%nn1(8)%ia = ox(il,io,i)%nn1(4)%ia + 8! 7, +z -a2+a1
	 else ! il=nlayers: nns are in the image of the first layer in the cell
	  ox(il,io,i)%nn1(5)%ia = 2 !ox(1,io,i)%nn1(1)%ia ! 2, +z
	  ox(il,io,i)%nn1(6)%ia = 3 !ox(1,io,i)%nn1(2)%ia ! 3, +z
	  ox(il,io,i)%nn1(7)%ia = 6 !ox(1,io,i)%nn1(3)%ia ! 6, +z -a2
	  ox(il,io,i)%nn1(8)%ia = 7 !ox(1,io,i)%nn1(4)%ia ! 7, +z -a2+a1
	 endif	

	 ox(il,io,i)%nn1(1)%r = oct(il,io)%ro(1,:) ! 2, same octa
	 ox(il,io,i)%nn1(2)%r = oct(il,io)%ro(2,:) ! 3, same octa
	 ox(il,io,i)%nn1(3)%r = oct(il,jo)%ro(1,:) -a2 ! 6, -a2
	 ox(il,io,i)%nn1(4)%r = oct(il,jo)%ro(2,:) -a2+a1 ! 7, -a2+a1
	 ox(il,io,i)%nn1(5)%r = oct(il,io)%ro(1,:) +z ! 2, +z
	 ox(il,io,i)%nn1(6)%r = oct(il,io)%ro(2,:) +z !3, +z
	 ox(il,io,i)%nn1(7)%r = oct(il,jo)%ro(1,:) +z - a2! 6, -a2 +z
	 ox(il,io,i)%nn1(8)%r = oct(il,jo)%ro(2,:) +z -a2+a1 ! 7, -a2+a1 +z
	 !..................................................

	 !..................................................
	 ! nns of O with B2 octahedra
	 !..................................................

	 !..................................................
	 ! O_x of B2: ox(il,io,i)%ia = 6
	 io=2; i=1;
	 jo = 1;
	 ox(il,io,i)%nn1(1)%ia = ox(il,io,i)%ia + 1 ! 7, same octa
	 ox(il,io,i)%nn1(2)%ia = ox(il,io,i)%ia - 3 ! 3, other octa, same cell
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia + 2 ! 8, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn1(4)%ia = ox(nlayers,io,i)%ia + 2 ! 8' 
	  ox(il,io,i)%nn1(8)%ia = ox(nlayers,io,i)%ia - 2 ! 4' (+a2)
	 else ! inside unit cell
	  ox(il,io,i)%nn1(4)%ia = ox(il-1,io,i)%ia + 2 ! 8'
	  ox(il,io,i)%nn1(8)%ia = ox(il-1,io,i)%ia - 2 ! 4' (+a2)
	 endif
	 ox(il,io,i)%nn1(5)%ia = ox(il,io,i)%ia + 1 ! 7, +a1
	 ox(il,io,i)%nn1(6)%ia = ox(il,io,i)%ia - 3 ! 3, +a2
	 ox(il,io,i)%nn1(7)%ia = ox(il,io,i)%ia - 2 ! 4, +a2

	 ox(il,io,i)%nn1(1)%r = oct(il,io)%ro(2,:) ! 7, same octa
	 ox(il,io,i)%nn1(2)%r = oct(il,jo)%ro(2,:) ! 3, other octa, same cell
	 ox(il,io,i)%nn1(3)%r = oct(il,io)%ro(3,:) ! 8, same cell
	 ox(il,io,i)%nn1(4)%r = oct(il,io)%ro(3,:) -z ! 8', layer below
	 ox(il,io,i)%nn1(5)%r = oct(il,io)%ro(2,:) +a1 ! 7, +a1
	 ox(il,io,i)%nn1(6)%r = oct(il,jo)%ro(2,:) +a2 ! 3, +a2
	 ox(il,io,i)%nn1(7)%r = oct(il,jo)%ro(3,:) +a2 ! 4, +a2
	 ox(il,io,i)%nn1(8)%r = oct(il,jo)%ro(3,:) -z +a2 ! 4', +a2, layer below
	 !..................................................
	 ! O_y of B2: ox(il,io,i)%ia = 7
	 io=2; i=2;
	 jo = 1;
	 ox(il,io,i)%nn1(1)%ia = ox(il,io,i)%ia - 1 ! 6, same octa
	 ox(il,io,i)%nn1(2)%ia = ox(il,io,i)%ia - 5 ! 2, -a1
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia + 1 ! 8, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn1(4)%ia = ox(nlayers,io,i)%ia + 1 ! 8' 
	  ox(il,io,i)%nn1(8)%ia = ox(nlayers,io,i)%ia - 3 ! 4', +a2-a1
	 else ! inside unit cell
	  ox(il,io,i)%nn1(4)%ia = ox(il-1,io,i)%ia + 1 ! 8'
	  ox(il,io,i)%nn1(8)%ia = ox(il-1,io,i)%ia - 3 ! 4', +a2-a1
	 endif
	 ox(il,io,i)%nn1(5)%ia = ox(il,io,i)%ia - 1 ! 6, -a1
	 ox(il,io,i)%nn1(6)%ia = ox(il,io,i)%ia - 5 ! 2, +a2-a1
	 ox(il,io,i)%nn1(7)%ia = ox(il,io,i)%ia - 3 ! 4, +a2-a1

	 ox(il,io,i)%nn1(1)%r = oct(il,io)%ro(1,:) ! 6, same octa
	 ox(il,io,i)%nn1(2)%r = oct(il,jo)%ro(1,:) -a1 ! 2, -a1
	 ox(il,io,i)%nn1(3)%r = oct(il,io)%ro(3,:)  ! 8, same cell
	 ox(il,io,i)%nn1(4)%r = oct(il,io)%ro(3,:) -z ! 8', layer below
	 ox(il,io,i)%nn1(5)%r = oct(il,io)%ro(1,:) -a1 ! 6, -a1
	 ox(il,io,i)%nn1(6)%r = oct(il,jo)%ro(1,:) +a2-a1 ! 2, +a2-a1
	 ox(il,io,i)%nn1(7)%r = oct(il,jo)%ro(3,:) +a2-a1 ! 4, +a2-a1
	 ox(il,io,i)%nn1(8)%r = oct(il,jo)%ro(3,:) -z +a2-a1 ! 4',layer below,+a2-a1
	 !..................................................
	 ! O_z of B2: ox(il,io,i)%ia = 8
	 io=2; i=3;
	 jo = 1;
	 ! same layer
	 ox(il,io,i)%nn1(1)%ia = ox(il,io,i)%ia - 2 ! 6, same octa
	 ox(il,io,i)%nn1(2)%ia = ox(il,io,i)%ia - 1 ! 7, same octa
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia - 6 ! 2, -a1
	 ox(il,io,i)%nn1(3)%ia = ox(il,io,i)%ia - 5 ! 3, same cell
	 ! layer above
	 if(il<nlayers) then ! the layer above is inside the cell
	  ox(il,io,i)%nn1(5)%ia = ox(il,io,i)%nn1(1)%ia + 8 ! 6, same octa, +z
	  ox(il,io,i)%nn1(6)%ia = ox(il,io,i)%nn1(1)%ia + 8 ! 7, same octa, +z
	  ox(il,io,i)%nn1(7)%ia = ox(il,io,i)%nn1(1)%ia + 8 ! 2, -a1, +z
	  ox(il,io,i)%nn1(8)%ia = ox(il,io,i)%nn1(1)%ia + 8 ! 3, same cell, +z
	 else ! il=nlayers: nns are in the image of the first layer in the cell
	  ox(il,io,i)%nn1(5)%ia = 6 ! 6, +z
	  ox(il,io,i)%nn1(6)%ia = 7 ! 7, +z
	  ox(il,io,i)%nn1(7)%ia = 2 ! 2, +z -a1
	  ox(il,io,i)%nn1(8)%ia = 3 ! 3, +z
	 endif	

	 ox(il,io,i)%nn1(1)%r = oct(il,io)%ro(1,:) ! 6, same octa
	 ox(il,io,i)%nn1(2)%r = oct(il,io)%ro(2,:) ! 7, same octa
	 ox(il,io,i)%nn1(3)%r = oct(il,jo)%ro(1,:) -a1 ! 2, -a1
	 ox(il,io,i)%nn1(4)%r = oct(il,jo)%ro(2,:)  ! 3, same cell
	 ox(il,io,i)%nn1(5)%r = oct(il,io)%ro(1,:) +z ! 6, +z
	 ox(il,io,i)%nn1(6)%r = oct(il,io)%ro(2,:) +z !7, +z
	 ox(il,io,i)%nn1(7)%r = oct(il,jo)%ro(1,:) +z - a1! 2, -a1 +z
	 ox(il,io,i)%nn1(8)%r = oct(il,jo)%ro(2,:) +z ! 3, +z
	 !..................................................


	end do ! il
	!.....................................................

	return
	end subroutine setonn2
	!.....................................................
	! TM orbitals are indexed first. O later, so TM orbitals make 
	! one block that we can project onto, if we like.
	subroutine mapatom2orbs()
	implicit none
	integer :: il, io,jo,i,j

	allocate(atom2orb(2,natoms))
	! norbtm = number of d orbitals on TM atoms. 6 t2g or all 10?
	! TM atoms
	i1 = 0;
	do il=1,nlayers
	 do io=1,noctl
	  ia = tm(il,io)%ia
    atom2orb(1,ia) = i1 + 1; ! start index
    atom2orb(2,ia) = i1 + norbtms; ! end index
	  i1 = i1 + norbtms;
	 end do
	end do

	! O atoms
	do il=1,nlayers
	 do io=1,noctl
	  do i=1,3
	   ia = ox(il,io,i)%ia
     atom2orb(1,ia) = i1 + 1; ! start index
     atom2orb(2,ia) = i1 + norbos; ! end index
	   i1 = i1 + norbtms;
	  end do
	 end do
	end do

	return
	end 	subroutine mapatom2orbs
	!.....................................................

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

	s = sk(1); p = sk(2); ! s=sigma_pp, p=pi_pp
	sp = s-p;

	rr = norm2(r);
	! direction cosines; r is from first to the second atom
	l = r(1)/rr;
	m = r(2)/rr;
	n = r(3)/rr;

	! (x,x), (x,y), (x,z)
	h(1,1) = l**2 *s + (1-l**2) *p
	h(1,2) = l*m*sp
	h(1,3) = l*n*sp
	! (y,x), (y,y), (y,z)
	h(2,1) = h(2,1)
	h(2,2) = m**2 *s + (1-m**2) *p
	h(2,3) = m*n*sp
	! (z,x), (z,y), (z,z)
	h(3,1) = n*l*sp
	h(3,2) = n*m*sp
	h(3,3) = n**2 *s + (1-n**2) *p

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

	s = sk(1); p = sk(2); ! s=sigma_pd, p=pi_pd
	sq3 = dsqrt(3);

	rr = norm2(r);
	! direction cosines; r is from first to the second atom
	l = r(1)/rr;
	m = r(2)/rr;
	n = r(3)/rr;

	!..................................................	
	! with t2g
	!..................................................
	lmn = l*m*n;
	lp = (1.d0-2*l**2)*p
	mp = (1.d0-2*m**2)*p
	np = (1.d0-2*n**2)*p	

	! x, xy
	h(1,1) = sq3 *l**2 *m *s + m*lp !(1.d0-2*l**2)*p
	! x, yz
	h(1,2) = (sq3 *s - 2*p)*lmn
	! x, zx
	h(1,3) = sq3 *l**2 *n *s + n*lp !(1.d0-2*l**2)*p

	! y, xy
	h(2,1) = sq3 *m**2 *l *s + l* mp !(1.d0-2*m**2)*p
	! y, yz
	h(2,2) = sq3 *m**2 *n *s + n* mp !(1.d0-2*m**2)*p
	! y, zx
	h(2,3) = h(1,2) !(sq3 *s - 2*p)*lmn

	! z, xy
	h(3,1) = h(1,2) !(sq3 *s - 2*p)*lmn
	! z, yz
	h(3,2) = sq3 *n**2 *m *s + m* np !(1.d0-2*n**2)*p
	! z, zx
	h(3,3) = sq3 *n**2 *l *s + l* np ! (1.d0-2*n**2)*p
	!..................................................

	!..................................................	
	! with eg
	!..................................................
	if(norbtm==5) then
	 lm = (l**2 - m**2);
	 lm2 = 0.5d0*sq3*lm*s
	 nlm2 = (n**2 - 0.5d0*(l**2 + m**2)) *s
	 ! x, x^2-y^2 
	 h(1,4) = l* (lm2 +(1.d0-lm)*p)
	 ! y, x^2-y^2 
	 h(2,4) = m* (lm2 -(1.d0+lm)*p)
	 ! z, x^2-y^2 
	 h(3,4) = n* (lm2 -lm*p)

	 ! x, 3z^2-r^2
	 h(1,5) = l*(nlm2 - sq3*n**2*p)
	 ! y, 3z^2-r^2
	 h(2,5) = m*(nlm2 - sq3*n**2*p)
	 ! z, 3z^2-r^2
	 h(3,5) = n*(nlm2 - sq3*(l**2+m**2)*p)
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
	double precision :: lmn, s, p, sq3, lp,mp,np

	s = sk(1); p = sk(2); d = sk(3); ! s=sigma_dd, p=pi_dd, d = delta_dd
	sq3 = dsqrt(3);

	rr = norm2(r);
	! direction cosines; r is from first to the second atom
	l = r(1)/rr;
	m = r(2)/rr;
	n = r(3)/rr;

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


	!..................................................	
	! (t2g,eg)
	!..................................................
	lm = (l**2 - m**2);
	lm2 = 1.5d0*lm *s;
	! xy, x^2-y^2
	h(1,4) = l*m*lm2 - 2*l*m*lm*p + 0.5d0*l*m*lm*d
	! yz, x^2-y^2
	h(2,4) = m*n*lm2 - m*n*(1.d0+2*lm)*p + m*n*(1.d0+0.5d0*lm)*d
	! zx, x^2-y^2
	h(3,4) = n*l*(lm2 + (1.d0-2*lm)*p + (1.d0+0.5d0*lm)*d)
	!..................................................
	lm = l**2 + m**2;
	lm2 = sq3*(n**2 -0.5d0*lm)*s;
	n2 = n**2;
	! xy, 3z^2-r^2
	h(1,5) = l*m*lm2 - 2*sq3*l*m*n2*p + 0.5d0*sq3*l*m*(1.d0+n2)*d
	! yz, 3z^2-r^2
	h(2,5) = m*n*lm2 + sq3*m*n*(lm - n2)*p - 0.5d0*sq3*m*n*lm*d ! last term: l+m^2  or l^2+m^2 ?? check
	! zx, 3z^2-r^2
	h(3,5) = l*n*lm2 + sq3*l*n*(lm - n2)*p - 0.5d0*sq3*l*n*lm*d
	!..................................................



	


	
	

















	return
	end 	subroutine slatkosdd



	! readinput: allocate space for SK and read them...
	! 	allocate(skbo(nsptm,2))
	! allocate(skbb(nsptm,nsptm,3))







	

	! px: with eg
	! x, x^2-y^2

	! x, 3z^2-r^2

	
	! sk(1) = sigma_pd; sk(2)=pi_pd

! H(x,xy) = dsqrt(3) l**2 * m * sk(1) + m*(1-2*l**2)*sk(2)
! H(x,x^2-y^2)=0.5*dsqrt(3)*l*(l**2-m**2)*sk(1) + l*(1-l**2-m**2)*sk(2)













	end 	module modmain
