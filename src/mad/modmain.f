
	module modmain
	implicit none

	integer :: nlayers, natomsl, natoms, nspin
	integer :: noctl, noct
	integer :: nspecies, nsptm !
	integer :: norbtm, norbo ! number of spatial orbitals on a TM/O atom
	integer :: norbtms, norbos ! number of spin-space orbitals on a TM/O atom
	integer :: ntot ! hilbert space size
	logical :: tmnn2, oxnn2, singlesk, lsoc, lhu, lspin
	integer :: ntottm
	double precision :: qtot 
	double complex, parameter :: iota = dcmplx(0.0d0,1.0d0)

	integer :: maxscf ! max scf iterations
	double precision :: toldm ! = 1.0d-6
	double precision :: beta0, betamax 
	integer :: mtype

	integer :: nwplot
	logical :: lpdos, lbc, lusevmat, lgs
	double precision, parameter :: fourpi=12.566370614359172954d0

	integer, parameter, dimension(-2:2) :: m2i=(/1,2,5,3,4/) ! Mth eleme of m2i is index of the corresponding d orbital in our code.; M as in Fernandez-Seivane et al. JPCM 2006.
	double complex, dimension(10,10) :: Hsoc ! TM soc

	
	!integer, allocatable, dimension(:) :: layersp ! layer TM species
	double precision, allocatable, dimension(:) :: soc ! TM soc
	double precision, allocatable, dimension(:) :: nds ! number of elec in d orbitals in TM atom (not ion).

	double precision, allocatable, dimension(:,:,:,:,:) :: gcmat ! Gaunt coeff matrix

	double precision :: reducebf
		
	double precision:: phi0

	double complex, allocatable, dimension(:,:):: hk
	double precision, allocatable, dimension(:,:) :: eval
	double precision, allocatable, dimension(:,:) :: hii ! onsite matrix elements of h
	double complex, allocatable, dimension(:,:,:):: evec
	double complex, allocatable, dimension(:,:,:):: hksave


	!double precision, allocatable, dimension(:,:,:,:) :: kscf
	!double precision, allocatable, dimension(:,:) :: kband
	double precision, allocatable, dimension(:,:) :: vvlp1d
	integer :: nv, np
	double precision, allocatable, dimension(:,:) :: vvl, vpl
	double precision, allocatable, dimension(:)  :: dv,dp


	double precision, allocatable, dimension(:,:) :: kgrid ! grid in irreducible wedge of BZ
	double precision, allocatable, dimension(:) :: wk ! weights
	double precision, allocatable, dimension(:,:) :: wke ! weights*occupations
	double precision :: temp, efermi, wknorm

	
	integer :: ntotk, nk1, nk2, nk3

	
	integer, allocatable, dimension(:,:) :: atom2orb
	integer, allocatable, dimension(:) :: atom2species
	double precision :: a ! lattice constant of a single formula unit cubic cell
	double precision, parameter :: pi = 3.141592653589793d0

	! SK parameters
	double precision, allocatable, dimension(:,:):: skbo ! spd, ppd
	double precision, allocatable, dimension(:,:,:):: skbb ! sdd, pdd, ddd
	double precision, dimension(2):: skoo ! spp, ppp
	double precision, allocatable, dimension(:):: onsite ! onsite energies of each atom, Oxygen at position 0, TM at 1:nsptm


	type :: hubbardparam
	 double precision :: U,J
	 double precision, dimension(0:4) :: Fk ! slater integrals, although we need F0,F2,F4 only, let's keep index k standard.
	double precision, allocatable, dimension(:,:,:,:):: Vee ! matrix for all atoms of a given species, using given U&J, via Fk and gcmat.
	end type hubbardparam

	
	type(hubbardparam), allocatable, dimension(:) :: Hub !  for each TM species

	type :: octahedra
		!integer :: ntot
		double precision :: phi, lo,lor,lort
		double precision, dimension(3) :: rb, rbf ! position of B atom, absolute position in a unit cell.
		double precision, dimension(3,3) :: xo ! relative to xb, position of oxygen atoms in cubic structure
		double precision, dimension(3,3) :: xor ! relative to xb, position of oxygen atoms after rotation
		double precision, dimension(3,3) :: ro, rof ! absolute position in a unit cell, position of oxygen atoms after rotation
	 !double precision, allocatable, dimension(:,:):: dm ! dm of TM atoms 
	end type octahedra

	type(octahedra), allocatable, dimension(:,:) :: oct

	double precision, allocatable, dimension(:) :: phi 

	! lattice vectors
	double precision, dimension(3,3) :: avec, ainv

	double precision, dimension(3,3) :: bvec, cvec
	double precision :: omega,omegabz

	double complex, dimension(5,5) :: Ur, Uz ! for transformation between real and complex Ylm

	double complex, dimension(10,10) :: Ulm2j ! for transformation between real (l,m)x(1/2,spin) and (J,Jz)

	type :: nneighbours
	 integer :: ia ! atom index, or index of equalent atom inside the unit cell (if this atom is outside the unit cell)
	 integer :: is ! TM species index
	 integer :: i, j ! orbital range indices: start, end; or equiv atom's orb ranges
	 double precision, dimension(3) :: r,dr ! position, true position of this nn.
	double precision, allocatable, dimension(:,:) :: h ! hamiltonian matrix elements between orbitals of two beighbours
	end type nneighbours

	type :: oatoms
	 character(len=2) :: label ! 'Ir', 'O', etc...
	 integer :: ia ! atom index in full list of atoms.
	 integer :: is ! TM species index
	 integer :: i, j ! orbital range indices: start, end
	 !integer :: typ ! orbital type, 1= p, 2=d, only one type allowed.
	 double precision, dimension(3) :: r ! position
	 type(nneighbours), allocatable, dimension(:) :: nn1, nn2 ! first and second nns, inside the unit cell or outside it.
	 !double precision, allocatable, dimension(:,:,:,:):: dm ! dm of TM atoms
	 !double precision, allocatable, dimension(:,:):: vmat, vmatold ! dm of TM atoms
	 !double precision, dimension(3):: mag ! magnetisation 
	 !double precision, dimension(3):: mfix	! fixed spin moment
	 !double precision, dimension(3):: beff ! effective magnetic field to fix the mom		  
	end type oatoms


	type :: tmatoms
	 character(len=2) :: label ! 'Ir', 'O', etc...
	 integer :: ia ! atom index in full list of atoms.
	 integer :: is ! TM species index
	 integer :: i, j ! orbital range indices: start, end
	 !integer :: typ ! orbital type, 1= p, 2=d, only one type allowed.
	 double precision, dimension(3) :: r ! position
	 type(nneighbours), allocatable, dimension(:) :: nn1, nn2 ! first and second nns, inside the unit cell or outside it.
	 double precision, allocatable, dimension(:,:,:,:):: dm ! dm of TM atoms
	 double precision, allocatable, dimension(:,:):: vmat, vmatold ! dm of TM atoms
	 double precision, dimension(3):: mag ! magnetisation 
	 double precision, dimension(3):: mfix	! fixed spin moment
	 double precision, dimension(3):: beff ! effective, sum of all types of bfields
	 double precision, dimension(3):: bext	! external magnetic field at TM site
	end type tmatoms

	type(tmatoms), allocatable, dimension(:,:) :: tm ! TM
	type(oatoms), allocatable, dimension(:,:,:) :: ox ! Oxygen

	! fixed spin moment calculation:
	! some variables can be moved to fixmom module
	double precision, dimension(3) :: bfsmc
	double precision, dimension(3) :: momtot, momfix ! total unit cell mom for fix mom calculations
	double precision :: taufsm
	integer :: fsmtype ! type of fixed spin mom: 1,2,3 or -1,-2,- 3: elk notation.
	double precision, dimension(3) :: bfieldc ! external global Bfield
	logical :: lbfields




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
		 write(fnum,'(3F14.8,"  ",3F12.8)')oct(il,io)%rof(i,1:3),0.,0.,0.
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
		write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%rbf, 0.,0.,0.
	 end do
	end do
	! second species:
	write(fnum,'("''",A,"''",T40," : spfname")') "Ti"
	write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")') noct/2
	do il=1, nlayers
	 do io =2, noctl,2
		write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%rbf, 0.,0.,0.
	 end do
	end do
	else ! odd
	! one species only
	write(fnum,'("''",A,"''",T40," : spfname")') "Ir"
	write(fnum,'(I4,T40," : natoms; atpos, bfcmt below")')noct
	do il=1, nlayers
	 do io =1, noctl
	  write(fnum,'(3F14.8,"  ",3F12.8)') oct(il,io)%rbf, 0.,0.,0.
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
	 oct(il,io)%rof(i,:) = v

	! test... randomise O pos a bit.
!	oct(il,io)%ro(i,:)=oct(il,io)%ro(i,:)
!     .   +(/rand(0),rand(0),rand(0)/)*a/2
	!write(*,'(3f10.4)') (/rand(0),rand(0),rand(0)/)

	 !write(*,'(3f10.4)') oct(il,io)%ro(i,:)
	end do
	! central B atom
	call r3mv(transpose(ainv),oct(il,io)%rb,v)
	oct(il,io)%rbf = v


	! test... randomise O pos a bit.
	!oct(il,io)%rb = 	oct(il,io)%rb + (/rand(0),rand(0),rand(0)/)




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
	
	 z = (/0.d0,0.d0,1.d0/)*a;
	 a1 = avec(1,:);
	 a2 = avec(2,:);
	 a3 = avec(3,:);

	!allocate(tm(nlayers,noctl))
	!.....................................................
	! frist nns of TM B atom are 6 O atoms
	!.....................................................
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always

	  tm(il,io)%r = oct(il,io)%rb
	  allocate(tm(il,io)%nn1(6)) ! O

	  if(io==1) then
	  	 tm(il,io)%ia = (il-1)*8 + io
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
	   tm(il,io)%nn1(5)%ia = tm(il,1)%ia + 1 ! 2, ! O_x of B1 in -a1
	   tm(il,io)%nn1(5)%r = oct(il,1)%ro(1,:) - a1
	  	endif

	  do i=1,3 ! O belonging to the unit cell
	   tm(il,io)%nn1(i)%ia = tm(il,io)%ia + i;
	   tm(il,io)%nn1(i)%r = oct(il,io)%ro(i,:)
	  end do
	  ! 3rd O could belong to the unit cell or maybe in a cell on below.
	  if(il==1) then ! O_z from the cell below, i.e., -a3
	   tm(il,io)%nn1(6)%ia = (nlayers-1)*8 + tm(il,io)%ia + 3
	   tm(il,io)%nn1(6)%r = oct(il,io)%ro(3,:) - z
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
	  	 tm(il,io)%ia = (il-1)*8 + io !  1st atom in a lyer
	   ! TM belonging to the unit cell or otherwise, all B' atoms
	   do i=1,6
	    tm(il,io)%nn2(i)%ia = tm(il,io)%ia + 4; ! 5th atom in a lyer
	   end do
		else ! io=2
	  	 tm(il,io)%ia = (il-1)*8 + 5;! 5th atom in a lyer
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
	subroutine setoxnn1()
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
	end subroutine setoxnn1

!..........................................................
! For the second nearest neighbours of O, 8 O atoms:
! sets indices (true or equiv in the cell) and positions (true).
!..........................................................
	subroutine setoxnn2()
	implicit none
	double precision, dimension(3) :: z, a1, a2, a3
	integer :: il, io,jo,i,j
	
	 z = (/0.d0,0.d0,1.d0/)*a;
	 a1 = avec(1,:);
	 a2 = avec(2,:);
	 a3 = avec(3,:);
	 !? distances could also be written in terms of a along x,y,z directions as in setmnnn2.

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
	 ox(il,io,i)%nn2(1)%ia = ox(il,io,i)%ia + 1 ! 3, same octa
	 ox(il,io,i)%nn2(2)%ia = ox(il,io,i)%ia + 5 ! 7, -a1-a2
	 ox(il,io,i)%nn2(3)%ia = ox(il,io,i)%ia + 2 ! 4, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn2(4)%ia = ox(nlayers,io,i)%ia + 2 ! 4' 
	  ox(il,io,i)%nn2(8)%ia = ox(nlayers,io,i)%ia + 6 ! 8' (+a1)
	 else ! inside unit cell
	  ox(il,io,i)%nn2(4)%ia = ox(il-1,io,i)%ia + 2 ! 4'
	  ox(il,io,i)%nn2(8)%ia = ox(il-1,io,i)%ia + 6 ! 8' (+a1)
	 endif
	 ox(il,io,i)%nn2(5)%ia = ox(il,io,i)%ia + 1 ! 3, +a1
	 ox(il,io,i)%nn2(6)%ia = ox(il,io,i)%ia + 5 ! 7, +a1
	 ox(il,io,i)%nn2(7)%ia = ox(il,io,i)%ia + 6 ! 8, +a1

	 ox(il,io,i)%nn2(1)%r = oct(il,io)%ro(2,:) ! 3, same octa
	 ox(il,io,i)%nn2(2)%r = oct(il,jo)%ro(2,:) -a1-a2 ! 7, -a1-a2
	 ox(il,io,i)%nn2(3)%r = oct(il,io)%ro(3,:)  ! 4, same cell
	 ox(il,io,i)%nn2(4)%r = oct(il,io)%ro(3,:) -z ! 4', layer below
	 ox(il,io,i)%nn2(5)%r = oct(il,io)%ro(2,:) +a1 ! 3, +a1
	 ox(il,io,i)%nn2(6)%r = oct(il,jo)%ro(2,:) +a1 ! 7, +a1
	 ox(il,io,i)%nn2(7)%r = oct(il,jo)%ro(3,:) +a1 ! 8, +a1
	 ox(il,io,i)%nn2(8)%r = oct(il,jo)%ro(3,:) -z +a1 ! 8', +a1, layer below
	 !..................................................
	 ! O_y of B1: ox(il,io,i)%ia = 3
	 io=1; i=2;
	 jo = 2;
	 ox(il,io,i)%nn2(1)%ia = ox(il,io,i)%ia - 1 ! 2, same octa
	 ox(il,io,i)%nn2(2)%ia = ox(il,io,i)%ia + 3 ! 6, -a2
	 ox(il,io,i)%nn2(3)%ia = ox(il,io,i)%ia + 1 ! 4, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn2(4)%ia = ox(nlayers,io,i)%ia + 1 ! 4' 
	  ox(il,io,i)%nn2(8)%ia = ox(nlayers,io,i)%ia + 5 ! 8' 
	 else ! inside unit cell
	  ox(il,io,i)%nn2(4)%ia = ox(il-1,io,i)%ia + 1 ! 4'
	  ox(il,io,i)%nn2(8)%ia = ox(il-1,io,i)%ia + 5 ! 8'
	 endif
	 ox(il,io,i)%nn2(5)%ia = ox(il,io,i)%ia + 3 ! 6, same cell
	 ox(il,io,i)%nn2(6)%ia = ox(il,io,i)%ia - 1 ! 2, -a2
	 ox(il,io,i)%nn2(7)%ia = ox(il,io,i)%ia + 5 ! 8, same cell

	 ox(il,io,i)%nn2(1)%r = oct(il,io)%ro(1,:) ! 2, same octa
	 ox(il,io,i)%nn2(2)%r = oct(il,jo)%ro(1,:) -a2 ! 6, -a2
	 ox(il,io,i)%nn2(3)%r = oct(il,io)%ro(3,:)  ! 4, same cell
	 ox(il,io,i)%nn2(4)%r = oct(il,io)%ro(3,:) -z ! 4', layer below
	 ox(il,io,i)%nn2(5)%r = oct(il,io)%ro(1,:) ! 6, same cell
	 ox(il,io,i)%nn2(6)%r = oct(il,jo)%ro(1,:) -a2 ! 2, -a2
	 ox(il,io,i)%nn2(7)%r = oct(il,jo)%ro(3,:) ! 8, same cell
	 ox(il,io,i)%nn2(8)%r = oct(il,jo)%ro(3,:) -z ! 8',layer below
	 !..................................................
	 ! O_z of B1: ox(il,io,i)%ia = 4
	 io=1; i=3;
	 jo = 2;
	 ! same layer
	 ox(il,io,i)%nn2(1)%ia = ox(il,io,i)%ia - 2 ! 2, same octa
	 ox(il,io,i)%nn2(2)%ia = ox(il,io,i)%ia - 1 ! 3, same octa
	 ox(il,io,i)%nn2(3)%ia = ox(il,io,i)%ia + 2 ! 6, -a2
	 ox(il,io,i)%nn2(4)%ia = ox(il,io,i)%ia + 3 ! 7, -a2+a1
	 ! layer above
	 if(il<nlayers) then ! the layer above is inside the cell
	  ox(il,io,i)%nn2(5)%ia = ox(il,io,i)%nn2(1)%ia + 8! 2, +z
	  ox(il,io,i)%nn2(6)%ia = ox(il,io,i)%nn2(2)%ia + 8! 3, +z
	  ox(il,io,i)%nn2(7)%ia = ox(il,io,i)%nn2(3)%ia + 8! 6, +z -a2
	  ox(il,io,i)%nn2(8)%ia = ox(il,io,i)%nn2(4)%ia + 8! 7, +z -a2+a1
	 else ! il=nlayers: nns are in the image of the first layer in the cell
	  ox(il,io,i)%nn2(5)%ia = 2 !ox(1,io,i)%nn2(1)%ia ! 2, +z
	  ox(il,io,i)%nn2(6)%ia = 3 !ox(1,io,i)%nn2(2)%ia ! 3, +z
	  ox(il,io,i)%nn2(7)%ia = 6 !ox(1,io,i)%nn2(3)%ia ! 6, +z -a2
	  ox(il,io,i)%nn2(8)%ia = 7 !ox(1,io,i)%nn2(4)%ia ! 7, +z -a2+a1
	 endif	

	 ox(il,io,i)%nn2(1)%r = oct(il,io)%ro(1,:) ! 2, same octa
	 ox(il,io,i)%nn2(2)%r = oct(il,io)%ro(2,:) ! 3, same octa
	 ox(il,io,i)%nn2(3)%r = oct(il,jo)%ro(1,:) -a2 ! 6, -a2
	 ox(il,io,i)%nn2(4)%r = oct(il,jo)%ro(2,:) -a2+a1 ! 7, -a2+a1
	 ox(il,io,i)%nn2(5)%r = oct(il,io)%ro(1,:) +z ! 2, +z
	 ox(il,io,i)%nn2(6)%r = oct(il,io)%ro(2,:) +z !3, +z
	 ox(il,io,i)%nn2(7)%r = oct(il,jo)%ro(1,:) +z - a2! 6, -a2 +z
	 ox(il,io,i)%nn2(8)%r = oct(il,jo)%ro(2,:) +z -a2+a1 ! 7, -a2+a1 +z
	 !..................................................

	 !..................................................
	 ! nns of O with B2 octahedra
	 !..................................................

	 !..................................................
	 ! O_x of B2: ox(il,io,i)%ia = 6
	 io=2; i=1;
	 jo = 1;
	 ox(il,io,i)%nn2(1)%ia = ox(il,io,i)%ia + 1 ! 7, same octa
	 ox(il,io,i)%nn2(2)%ia = ox(il,io,i)%ia - 3 ! 3, other octa, same cell
	 ox(il,io,i)%nn2(3)%ia = ox(il,io,i)%ia + 2 ! 8, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn2(4)%ia = ox(nlayers,io,i)%ia + 2 ! 8' 
	  ox(il,io,i)%nn2(8)%ia = ox(nlayers,io,i)%ia - 2 ! 4' (+a2)
	 else ! inside unit cell
	  ox(il,io,i)%nn2(4)%ia = ox(il-1,io,i)%ia + 2 ! 8'
	  ox(il,io,i)%nn2(8)%ia = ox(il-1,io,i)%ia - 2 ! 4' (+a2)
	 endif
	 ox(il,io,i)%nn2(5)%ia = ox(il,io,i)%ia + 1 ! 7, +a1
	 ox(il,io,i)%nn2(6)%ia = ox(il,io,i)%ia - 3 ! 3, +a2
	 ox(il,io,i)%nn2(7)%ia = ox(il,io,i)%ia - 2 ! 4, +a2

	 ox(il,io,i)%nn2(1)%r = oct(il,io)%ro(2,:) ! 7, same octa
	 ox(il,io,i)%nn2(2)%r = oct(il,jo)%ro(2,:) ! 3, other octa, same cell
	 ox(il,io,i)%nn2(3)%r = oct(il,io)%ro(3,:) ! 8, same cell
	 ox(il,io,i)%nn2(4)%r = oct(il,io)%ro(3,:) -z ! 8', layer below
	 ox(il,io,i)%nn2(5)%r = oct(il,io)%ro(2,:) +a1 ! 7, +a1
	 ox(il,io,i)%nn2(6)%r = oct(il,jo)%ro(2,:) +a2 ! 3, +a2
	 ox(il,io,i)%nn2(7)%r = oct(il,jo)%ro(3,:) +a2 ! 4, +a2
	 ox(il,io,i)%nn2(8)%r = oct(il,jo)%ro(3,:) -z +a2 ! 4', +a2, layer below
	 !..................................................
	 ! O_y of B2: ox(il,io,i)%ia = 7
	 io=2; i=2;
	 jo = 1;
	 ox(il,io,i)%nn2(1)%ia = ox(il,io,i)%ia - 1 ! 6, same octa
	 ox(il,io,i)%nn2(2)%ia = ox(il,io,i)%ia - 5 ! 2, -a1
	 ox(il,io,i)%nn2(3)%ia = ox(il,io,i)%ia + 1 ! 8, same cell
	 ! The layer below, 4'/8'
	 if(il==1) then ! image of the top most layer in the cell
	  ox(il,io,i)%nn2(4)%ia = ox(nlayers,io,i)%ia + 1 ! 8' 
	  ox(il,io,i)%nn2(8)%ia = ox(nlayers,io,i)%ia - 3 ! 4', +a2-a1
	 else ! inside unit cell
	  ox(il,io,i)%nn2(4)%ia = ox(il-1,io,i)%ia + 1 ! 8'
	  ox(il,io,i)%nn2(8)%ia = ox(il-1,io,i)%ia - 3 ! 4', +a2-a1
	 endif
	 ox(il,io,i)%nn2(5)%ia = ox(il,io,i)%ia - 1 ! 6, -a1
	 ox(il,io,i)%nn2(6)%ia = ox(il,io,i)%ia - 5 ! 2, +a2-a1
	 ox(il,io,i)%nn2(7)%ia = ox(il,io,i)%ia - 3 ! 4, +a2-a1

	 ox(il,io,i)%nn2(1)%r = oct(il,io)%ro(1,:) ! 6, same octa
	 ox(il,io,i)%nn2(2)%r = oct(il,jo)%ro(1,:) -a1 ! 2, -a1
	 ox(il,io,i)%nn2(3)%r = oct(il,io)%ro(3,:)  ! 8, same cell
	 ox(il,io,i)%nn2(4)%r = oct(il,io)%ro(3,:) -z ! 8', layer below
	 ox(il,io,i)%nn2(5)%r = oct(il,io)%ro(1,:) -a1 ! 6, -a1
	 ox(il,io,i)%nn2(6)%r = oct(il,jo)%ro(1,:) +a2-a1 ! 2, +a2-a1
	 ox(il,io,i)%nn2(7)%r = oct(il,jo)%ro(3,:) +a2-a1 ! 4, +a2-a1
	 ox(il,io,i)%nn2(8)%r = oct(il,jo)%ro(3,:) -z +a2-a1 ! 4',layer below,+a2-a1
	 !..................................................
	 ! O_z of B2: ox(il,io,i)%ia = 8
	 io=2; i=3;
	 jo = 1;
	 ! same layer
	 ox(il,io,i)%nn2(1)%ia = ox(il,io,i)%ia - 2 ! 6, same octa
	 ox(il,io,i)%nn2(2)%ia = ox(il,io,i)%ia - 1 ! 7, same octa
	 ox(il,io,i)%nn2(3)%ia = ox(il,io,i)%ia - 6 ! 2, -a1
	 ox(il,io,i)%nn2(4)%ia = ox(il,io,i)%ia - 5 ! 3, same cell
	 ! layer above
	 if(il<nlayers) then ! the layer above is inside the cell
	  ox(il,io,i)%nn2(5)%ia = ox(il,io,i)%nn2(1)%ia + 8 ! 6, same octa, +z
	  ox(il,io,i)%nn2(6)%ia = ox(il,io,i)%nn2(1)%ia + 8 ! 7, same octa, +z
	  ox(il,io,i)%nn2(7)%ia = ox(il,io,i)%nn2(1)%ia + 8 ! 2, -a1, +z
	  ox(il,io,i)%nn2(8)%ia = ox(il,io,i)%nn2(1)%ia + 8 ! 3, same cell, +z
	 else ! il=nlayers: nns are in the image of the first layer in the cell
	  ox(il,io,i)%nn2(5)%ia = 6 ! 6, +z
	  ox(il,io,i)%nn2(6)%ia = 7 ! 7, +z
	  ox(il,io,i)%nn2(7)%ia = 2 ! 2, +z -a1
	  ox(il,io,i)%nn2(8)%ia = 3 ! 3, +z
	 endif	

	 ox(il,io,i)%nn2(1)%r = oct(il,io)%ro(1,:) ! 6, same octa
	 ox(il,io,i)%nn2(2)%r = oct(il,io)%ro(2,:) ! 7, same octa
	 ox(il,io,i)%nn2(3)%r = oct(il,jo)%ro(1,:) -a1 ! 2, -a1
	 ox(il,io,i)%nn2(4)%r = oct(il,jo)%ro(2,:)  ! 3, same cell
	 ox(il,io,i)%nn2(5)%r = oct(il,io)%ro(1,:) +z ! 6, +z
	 ox(il,io,i)%nn2(6)%r = oct(il,io)%ro(2,:) +z !7, +z
	 ox(il,io,i)%nn2(7)%r = oct(il,jo)%ro(1,:) +z - a1! 2, -a1 +z
	 ox(il,io,i)%nn2(8)%r = oct(il,jo)%ro(2,:) +z ! 3, +z
	 !..................................................


	end do ! il
	!.....................................................

	return
	end subroutine setoxnn2
	!.....................................................
	! TM orbitals are indexed first. O later, so TM orbitals make 
	! one block that we can project onto, if we like.
	subroutine mapatom2orbs()
	implicit none
	integer :: il, io,i, i1, ia

	allocate(atom2orb(2*nspin,natoms))
	! norbtm = number of d orbitals on TM atoms. 6 t2g or all 10?
	! TM atoms
	i1 = 0;
	do il=1,nlayers
	 do io=1,noctl
	  ia = tm(il,io)%ia
	  atom2orb(1,ia) = i1 + 1; ! start index, up spin
	  i1 = i1 + norbtm;
	  atom2orb(2,ia) = i1; ! end index, up spin
		if(nspin==2) then
	   atom2orb(3,ia) = i1 + 1; ! start index, down spin
	   i1 = i1 + norbtm;
	   atom2orb(4,ia) = i1; ! end index, down spin
		endif
	  !write(*,*)'ia, i2 = ',ia, i1
	 end do
	end do

	! O atoms
	do il=1,nlayers
	 do io=1,noctl
	  do i=1,3
	   ia = ox(il,io,i)%ia
	   atom2orb(1,ia) = i1 + 1; ! start index
	   i1 = i1 + norbo;
	   atom2orb(2,ia) = i1; ! end index
		 if(nspin==2) then
	    atom2orb(3,ia) = i1 + 1; ! start index, down spin
	    i1 = i1 + norbo;
	    atom2orb(4,ia) = i1; ! end index, down spin
		 endif
	  end do
	 end do
	end do

	return
	end 	subroutine mapatom2orbs
!.....................................................
! atom index to species index.... used for TM atoms to identify their SK param
! if SK param are set equal for all TM, we can actually set js=1 where atom2species() map is used. but, it cost nothing so does not matter.
	subroutine mapatom2species()
	implicit none
	integer :: il, io, i, ia, is

	allocate(atom2species(natoms))
	atom2species(:) = -1;

	! TM ia to is:
	do il=1, nlayers
	do io=1, noctl
	 if(io==1)then
	  ia = (il-1)*8 + 1 ! TM index of this octahedra
	 else
	  ia = (il-1)*8 + 5 ! TM index of this octahedra
	 endif
	 atom2species(ia) = tm(il,io)%is !layersp(il)
	end do
	end do
	!write(*,*)'Warning(mapatom2species): setting TM atom2species(:)=1'

	!write(*,'(1000i5)')atom2species
	return
	end 	subroutine mapatom2species
!.....................................................

	! readinput: allocate space for SK and read them...
	! 	allocate(skbo(nsptm,2))
	! allocate(skbb(nsptm,nsptm,3))



	



!..................................................................
! copied from ELK 6.2 version:
!..................................................................
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
	subroutine reciplat(avec,bvec,omega,omegabz)
! !INPUT/OUTPUT PARAMETERS:
!   avec    : lattice vectors (in,real(3,3))
!   bvec    : reciprocal lattice vectors (out,real(3,3))
!   omega   : unit cell volume (out,real)
!   omegabz : Brillouin zone volume (out,real)
! !DESCRIPTION:
!   Generates the reciprocal lattice vectors from the real-space lattice vectors
!   \begin{align*}
!     {\bf b}_1&=\frac{2\pi}{s}({\bf a}_2\times{\bf a}_3)\\
!     {\bf b}_2&=\frac{2\pi}{s}({\bf a}_3\times{\bf a}_1)\\
!     {\bf b}_3&=\frac{2\pi}{s}({\bf a}_1\times{\bf a}_2)
!   \end{align*}
!   and finds the unit cell volume $\Omega=|s|$, where
!   $s={\bf a}_1\cdot({\bf a}_2\times{\bf a}_3)$, and the Brillouin zone volume
!   $\Omega_{\rm BZ}=(2\pi)^3/\Omega$.
!
	implicit none
! arguments
	real(8), intent(in) :: avec(3,3)
	real(8), intent(out) :: bvec(3,3)
	real(8), intent(out) :: omega,omegabz
! local variables
	real(8), parameter :: twopi=6.2831853071795864769d0
	real(8) t1
	call r3cross(avec(:,2),avec(:,3),bvec(:,1))
	call r3cross(avec(:,3),avec(:,1),bvec(:,2))
	call r3cross(avec(:,1),avec(:,2),bvec(:,3))
	t1=avec(1,1)*bvec(1,1)+avec(2,1)*bvec(2,1)+avec(3,1)*bvec(3,1)
! unit cell volume
	omega=abs(t1)
	if (omega.lt.1.d-6) then
	  write(*,*)
	  write(*,'("Error(reciplat) omega too small : ",G18.10)') omega
	  write(*,'(" Lattice vectors may be collinear")')
	  write(*,*)
	  stop
	end if
	bvec(:,:)=(twopi/t1)*bvec(:,:)
! Brillouin zone volume
	omegabz=(twopi**3)/omega
	return
	end subroutine
!..................................................................
	subroutine plotpt1d(cvec,nv,np,vvl,vpl,dv,dp)
! !INPUT/OUTPUT PARAMETERS:
!   cvec : matrix of (reciprocal) lattice vectors stored column-wise
!         (in,real(3,3))
!   nv   : number of vertices (in,integer)
!   np   : number of connecting points (in,integer)
!   vvl  : vertex vectors in lattice coordinates (in,real(3,nv))
!   vpl  : connecting point vectors in lattice coordinates (out,real(3,np))
!   dv   : cummulative distance to each vertex (out,real(nv))
!   dp   : cummulative distance to each connecting point (out,real(np))
! !DESCRIPTION:
!   Generates a set of points which interpolate between a given set of vertices.
!   Vertex points are supplied in lattice coordinates in the array {\tt vvl} and
!   converted to Cartesian coordinates with the matrix {\tt cvec}. Interpolating
!   points are stored in the array {\tt vpl}. The cummulative distances to the
!   vertices and points along the path are stored in arrays {\tt dv} and
!   {\tt dp}, respectively.
!
	implicit none
! arguments
	real(8), intent(in) :: cvec(3,3)
	integer, intent(in) :: nv,np
	real(8), intent(in) :: vvl(3,nv)
	real(8), intent(out) :: vpl(3,np),dv(nv),dp(np)
! local variables
	integer i,j,k,m,n
	real(8) vl(3),vc(3)
	real(8) dt,f,t1
! alloctable arrays
	real(8), allocatable :: seg(:)
	if (nv.lt.1) then
	  write(*,*)
	  write(*,'("Error(plotpt1d): nv < 1 : ",I8)') nv
	  write(*,*)
	  stop
	end if
	if (np.lt.nv) then
	  write(*,*)
	  write(*,'("Error(plotpt1d): np < nv : ",2I8)') np,nv
	  write(*,*)
	  stop
	end if
	! special case of 1 vertex
	if (nv.eq.1) then
	  dv(1)=0.d0
	  dp(:)=0.d0
	  do i=1,np
	    vpl(:,i)=vvl(:,1)
	  end do
	  return
	end if
	allocate(seg(nv))
! find the length of each segment and total distance
	dt=0.d0
	do i=1,nv-1
	  dv(i)=dt
	  vl(:)=vvl(:,i+1)-vvl(:,i)
	  call r3mv(cvec,vl,vc)
	  seg(i)=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
	  dt=dt+seg(i)
	end do
	dv(nv)=dt
	! add small amount to total distance to avoid 0/0 condition
	dt=dt+1.d-8
! number of points to use between vertices
	n=np-nv
! construct the interpolating path
	k=0
	do i=1,nv-1
	  t1=dble(n)*seg(i)/dt
	  m=nint(t1)
	  if ((m.gt.n).or.(i.eq.(nv-1))) m=n
	  do j=1,m+1
	    k=k+1
	    f=dble(j-1)/dble(m+1)
	    dp(k)=dv(i)+f*seg(i)
	    vpl(:,k)=vvl(:,i)*(1.d0-f)+vvl(:,i+1)*f
	  end do
	  dt=dt-seg(i)
	  n=n-m
	end do
	dp(np)=dv(nv)
	vpl(:,np)=vvl(:,nv)
	deallocate(seg)
	return
	end subroutine
!..................................................................
	subroutine r3cross(x,y,z)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector 1 (in,real(3))
!   y : input vector 2 (in,real(3))
!   z : output cross-product (out,real(3))
! !DESCRIPTION:
!   Returns the cross product of two real 3-vectors.
	implicit none
! arguments
	real(8), intent(in) :: x(3),y(3)
	real(8), intent(out) :: z(3)
	z(1)=x(2)*y(3)-x(3)*y(2)
	z(2)=x(3)*y(1)-x(1)*y(3)
	z(3)=x(1)*y(2)-x(2)*y(1)
	return
	end subroutine
!..................................................................




!-----------------------------------
	subroutine ddiag(ntot,H,eval,ijob,nev)
	implicit none
	integer, intent(in) :: ntot,ijob,nev
	!double precision, dimension(ntot,ntot), intent(inout):: H	
	!double precision, dimension(ntot), intent(inout):: W
	double precision, dimension(ntot,ntot), intent(inout):: H	
	double precision, dimension(ntot), intent(out):: eval
	! local
	double precision, dimension(ntot) :: W
	integer :: lwork != 3*ntot ! LWORK >= max(1,3*N-1)
	double precision, dimension(3*ntot):: WORK
	integer :: info,j
	double precision ::starttime, endtime,starttime1, endtime1
	EXTERNAL	 DSYEV ! LAPACK/BLAS

	!write(*,'(/,a)') " diagonal: diagonalising H ... "
	!starttime1 = clock()

	lwork = 3*ntot ! LWORK >= max(1,3*N-1)

	if(1==0) then
	write(*,*)'************************************'
	do j=1,ntot
		write(*,'(10000E8.1)') H(j,:)
	end do
	write(*,*)'************************************'
	endif
	
	
	CALL DSYEV('Vectors','Upper',ntot,H,ntot,W,WORK,LWORK,INFO)
	! Check for convergence.
	IF( INFO.GT.0 ) THEN
		write(*,'(/,a)')
     .   'diagonal: The algorithm failed to compute eigenvalues.'
		STOP
	END IF


	!eig(ijob)%evec = H(1:ntot,1:nev)
	eval = W(1:nev)



	return
	END subroutine ddiag
!---------------------------------------	


!-----------------------------------
	subroutine zdiag(ntot,H,eval,ijob,nev)
	implicit none
	integer, intent(in) :: ntot,ijob,nev
	!double precision, dimension(ntot,ntot), intent(inout):: H	
	!double precision, dimension(ntot), intent(inout):: W
	double complex, dimension(ntot,ntot), intent(inout):: H	
	double precision, dimension(ntot), intent(out):: eval
	!double complex, dimension(ntot,ntot), intent(out) :: evec
	! local
	double precision, dimension(ntot) :: W, RWORK(3*ntot-2)
	integer :: lwork, LWMAX  != 3*ntot ! LWORK >= max(1,3*N-1)
	double complex, dimension(3*ntot):: WORK
	integer :: info,j
	double precision ::starttime, endtime,starttime1, endtime1
	EXTERNAL ZHEEV ! LAPACK/BLAS
	INTRINSIC INT, MIN

	LWMAX = 3*ntot;
	
	!write(*,'(/,a)') " diagonal: diagonalising H ... "
	!starttime1 = clock()

	!lwork = 3*ntot ! LWORK >= max(1,3*N-1)

	if(1==0) then
	write(*,*)'************************************'
	do j=1,ntot
		write(*,'(10000E8.1)') H(j,:)
	end do
	write(*,*)'************************************'
	endif

	! Query the optimal workspace.
	LWORK = -1
	CALL ZHEEV('Vectors','Upper',ntot,H,ntot,W,WORK,LWORK,RWORK,INFO)
	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
	
	CALL ZHEEV('Vectors','Upper',ntot,H,ntot,W,WORK,LWORK,RWORK,INFO)
	! Check for convergence.
	IF( INFO.GT.0 ) THEN
		write(*,'(/,a)')
     .   'diagonal: The algorithm failed to compute eigenvalues.'
		STOP
	END IF


	!eig(ijob)%evec = H(1:ntot,1:nev)
	eval = W(1:nev)



	return
	END subroutine zdiag
!---------------------------------------	
! for details, see Fernandez-Seivane et al. J. Phys.: Condens. Matter 18 (2006) 7999–8013
! their expression for real spherical harmonics has a different sign for M=+1,-1. than the standard, wiki
! https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
! http://www.theochem.ru.nl/~pwormer/Knowino/knowino.org/wiki/Spherical_harmonics.html expression near the end of the page.
! we will see later if their results for Vso need any correction or if their choice is cnsistent and we have to modify their results for our choice of real Ylm.
C *********************************************************************
C SIESTA:
C SPIN-ORBIT INTERACTION ---  ON-SITE APPROXIMATION 
C
C The spin-orbit hamiltonian has the form:
C 
C                        | Lz             Lx - i Ly |
C  <i|HSO|j> = <i| V(r)  |                          | |j>
C                        | Lx + i Ly      - Lz      |
C
C                        | i*L1(li,mi,mj)         L2(li,mi,mj)-i*L3(li,mi,mj) |
C  <i|HSO|j> = 0.5 M(li) |                                                    |
C                        | -L2(li,mi,mj)-i*L3(li,mi,mj)       -i*L1(li,mi,mj) |
C
C where M(li) is the radial part, and Li(li,mi,mj) are the radial bits
C 
C Hence,  <i|Lz|j> =  i L1
C         <i|Lx|j> = -i L3
C         <i|Ly|j> =  i L2
C

	subroutine mkHsoc()
	implicit none
	! local
	double complex, dimension(5,5) :: Lz, Ldn, Lup
	double precision, dimension(3) :: L
	!double complex, parameter :: iota = dcmplx(0.0d0,1.0d0);
	integer :: m,n

	! blocks in Hsoc
	Lz = dcmplx(0.0d0,0.0d0);
	Lup = dcmplx(0.0d0,0.0d0);
	Ldn = dcmplx(0.0d0,0.0d0);

	! calc blocks in terms of L1,L2,L3 real variables as in siesta
	do m=-2,2,1
	 do n=-2,2,1
	  call int_so_ang(m,n,L)
	  Lz(m2i(m),m2i(n)) =  iota * L(1)
	  Ldn(m2i(m),m2i(n)) = -iota* L(3) + L(2)
	  Lup(m2i(m),m2i(n)) = -iota* L(3) - L(2)
	 end do
	end do

	! set blocks to their respective positons:
	Hsoc(1:5,1:5) = Lz;
	Hsoc(6:10,6:10) = -Lz;
	Hsoc(1:5,6:10) = Ldn;
	Hsoc(6:10,1:5) = Lup;

	if(1==0) then
	write(*,*)"Lz in real Ylm basis (our M-order):"
	do m=1,5
	 write(*,'(10f10.6)') Lz(m,1:5)
	end do
	write(*,*)"Lup in real Ylm basis (our M-order):"
	do m=1,5
	 write(*,'(10f10.6)') Lup(m,:)
	end do
	endif
	

	return
	end 	subroutine mkHsoc


C *********************************************************************
C
C Subroutine to calculate the spin-orbit angular integral
C Calculates L1(li,mi,mj), L2(li,mi,mj) and L3(li,mi,,mj) 
C
C *********************************************************************
	subroutine int_so_ang(mi, mj, L)
	implicit none
	integer, intent(in)  :: mi, mj
	double precision,intent(out) :: L(3)
	double precision::  La, Lb
	double precision, parameter :: one = 1.d0, two = 2.d0
	double precision, parameter :: six = 6.d0
	integer :: li

	li = 2; ! always d-orbitals in our case
	L(1:3)= 0.0d0

	La = sqrt(li*(li+1.d0)/2.d0)
	Lb = sqrt(li*(li+1.d0)-2.d0)/2.d0

	if((mi+mj).EQ.0) L(1) = mj*1.0d0

	select case ( mi )
      case ( 0 )
         select case ( mj )
         case ( -1 )
            L(3) = La
         case ( 1 )
            L(2) = La
         end select
      case ( 1 )
         select case ( mj )
         case ( -2 )
            L(3) =  Lb
         case ( 0 )
            L(2) = -La
         case ( 2 )
            L(2) =  Lb
         end select
      case ( -1 )
         select case ( mj )
         case ( -2 )
            L(2) =  Lb
         case ( 0 )
            L(3) = -La
         case ( 2 )
            L(3) = -Lb
         end select
      case ( 2 )
         select case ( mj )
         case ( -1 )
            L(3) =  Lb
         case ( 1 )
            L(2) = -Lb
         end select
      case ( -2 )
         select case ( mj )
          case ( -1 )
            L(2) = -Lb
         case ( 1 )
            L(3) = -Lb
         end select
      end select

	! Ahsan
	! taking care of the difference sign convention by Fernandez-Seivane; 
	! our Y_{2,1/-1}^{our/standard} = - Y_{2,1/-1}^{Fernandez}
	! this might also be needed in siesta
	! I dont know if siesta uses Fernandez convention or the standard one.
	L = L * Ylmsgns(mi,mj) 	

	return
	end subroutine int_so_ang
C *********************************************************************
! not the most efficient way, but does not matter much.
	double precision function Ylmsgns(m,n)
	implicit none
	integer, intent(in) :: m,n

	Ylmsgns = 1.0d0
	if(m == 1 .or. m==-1 ) Ylmsgns = -1.0d0 * Ylmsgns
	if(n == 1 .or. n==-1 ) Ylmsgns = -1.0d0 * Ylmsgns

	return
	end function Ylmsgns
C *********************************************************************
! Hubbard U for d electrons:
! copied from elk-6.2.8/src/genfdu.f90 
! !INPUT/OUTPUT PARAMETERS:
!   i : DFT+U entry (in,integer)
!   u : parameter U (inout,real)
!   j : parameter J (inout,real)
!   f : Slater parameters (inout,real)
! !DESCRIPTION:
!   Calculate the Slater parameters for DFT+$U$ calculation with different
!   approaches, see  {\it Phys. Rev. B} {\bf 80}, 035121 (2009). The relations
!   among Slater and Racah parameters are from E.U. Condon and G.H. Shortley,
!   {\it The Theory of Atomic Spectra},  The University Press, Cambridge (1935).
C *********************************************************************
	subroutine setUrUz()
	implicit none
	double precision, parameter :: sqrt2inv=1.0d0/dsqrt(2.0d0)
	!double complex, parameter :: iota = dcmplx(0.0d0,1.0d0)
	! m2i=(/1,2,5,3,4/) ! Mth eleme of m2i is index of the corresponding d orbital in our code.
	double complex :: ab(5,5)
	integer :: i
	
	Ur(1:5,1) = (/1.0d0,0.0d0,0.0d0, 0.0d0,-1.0d0/)*iota*sqrt2inv ! r_{-2}
	Ur(1:5,2) = (/0.0d0,1.0d0,0.0d0, 1.0d0, 0.0d0/)*iota*sqrt2inv ! r_{-1}
	Ur(1:5,5) = (/0.0d0,0.0d0,1.0d0, 0.0d0, 0.0d0/) ! r_{0}
	Ur(1:5,3) = (/0.0d0,1.0d0,0.0d0,-1.0d0, 0.0d0/)*sqrt2inv ! r_{+1}
	Ur(1:5,4) = (/1.0d0,0.0d0,0.0d0, 0.0d0, 1.0d0/)*sqrt2inv ! r_{+2}
	
	Uz = conjg(transpose(Ur));

	if(1==0) then
	! r_{0} in z basis
	ab = 0.0d0; ab(3,3)=1.0d0; !matmul(Uz,Ur)

	ab = matmul(Uz,ab);
	ab = matmul(ab,Ur);
	! r_{0} in our real basis. ===> should become ab(5,5)=1.0, ab(else,else)=0.
	
	write(*,*)'------- real(Uz.Ur) --------- '
	do i=1,5
	 write(*,'(10f5.2)') real(ab(i,:))
	end do
	write(*,*)'------- Im(Uz.Ur) --------- '
	do i=1,5
	 write(*,'(10f5.2)') aimag(ab(i,:))
	end do
	write(*,*)'----------------------------'	
	endif
	
	return
	end 	subroutine setUrUz
C *********************************************************************
	
	subroutine setUlm2j()
	implicit none
	double precision, parameter :: sqrt2inv=1.0d0/dsqrt(2.0d0)
	double complex, dimension(10,10) :: U
	integer :: i,j

	!---------------------------------------------------------------------
	! U: our real Ylm to complex Ylm
	U=0.0d0
	U(1:5,1:5) = Ur; ! for up spin
	U(6:10,6:10) = Ur; ! for down spin
	! spins are not coupled in our representation so the two off-diagonal blocks of U are zero. 
	!---------------------------------------------------------------------
	! set complex Ylm,spin to J obtained by ClebschGordan[] in Mathematica
	! with spin up block of (complex) Ylm first: copied bwlow (as mathematica array rules after conversion to sparse):
	!{{1,6}->1,{2,1}->1/Sqrt[5],{2,7}->2/Sqrt[5],{3,2}->Sqrt[2/5],{3,8}->Sqrt[3/5],{4,3}->Sqrt[3/5],{4,9}->Sqrt[2/5],{5,4}->2/Sqrt[5],{5,10}->1/Sqrt[5],{6,5}->1,{7,1}->-(2/Sqrt[5]),{7,7}->1/Sqrt[5],{8,2}->-Sqrt[(3/5)],{8,8}->Sqrt[2/5],{9,3}->-Sqrt[(2/5)],{9,9}->Sqrt[3/5],{10,4}->-(1/Sqrt[5]),{10,10}->2/Sqrt[5],{_,_}->0}
	Ulm2j(:,:)=0.0d0;	
	Ulm2j(1,6)=1.0d0; 
	Ulm2j(2,1)=1.0d0/dSqrt(5.0d0); 
	Ulm2j(2,7)=2.0d0/dSqrt(5.0d0);
	Ulm2j(3,2)=dSqrt(2.0d0/5.0d0); 
	Ulm2j(3,8)=dSqrt(3.0d0/5.0d0);
	Ulm2j(4,3)=dSqrt(3.0d0/5.0d0); 
	Ulm2j(4,9)=dSqrt(2.0d0/5.0d0);
	Ulm2j(5,4)=2.0d0/dSqrt(5.0d0); 
	Ulm2j(5,10)=1.0d0/dSqrt(5.0d0);
	Ulm2j(6,5)=1.0d0; 
	Ulm2j(7,1)=-(2.0d0/dSqrt(5.0d0));
	Ulm2j(7,7)=1.0d0/dSqrt(5.0d0); 
	Ulm2j(8,2)=-dSqrt((3.0d0/5.0d0));
	Ulm2j(8,8)=dSqrt(2.0d0/5.0d0);
	Ulm2j(9,3)=-dSqrt((2.0d0/5.0d0));
	Ulm2j(9,9)=dSqrt(3.0d0/5.0d0);
	Ulm2j(10,4)=-(1.0d0/dSqrt(5.0d0));
	Ulm2j(10,10)=2.0d0/dSqrt(5.0d0)
	!---------------------------------------------------------------------
	! operating from left on a state in real Ylm and up/dn basis,
	! U changes the basis to complex Ylm and up/dn spin
	! Ulm2j changes this complex Ylm,spin to J.
	! combine the two transformation and saving as Ulm2j below:
	! so now Ulm2j is from real Ylm,spin to J:
	Ulm2j = matmul(Ulm2j,U);

	return
	end 	subroutine setUlm2j
!======================================================================
!	given atomic occupations:
!	subroutine setupdm()
!	implicit none
	! do not assume that TM atom gives its 2 d electrons to O
	! occupations on diagonal of dm of TM:
!	do il=1,nlayers
!	 do io=1,noctl
!	  allocate(tm(il,io)%vmat(norbtms,norbtms))
!	  allocate(tm(il,io)%vmatold(norbtms,norbtms))
!	  allocate(tm(il,io)%dm(norbtm, nspin,norbtm, nspin))
!	  tm(il,io)%dm = 0.0d0;
!	  is = layersp(il);
!	  do i=1,norbtm
!	   tm(il,io)%dm(i,i) = 1.0d-1 * (nds(is)+tm(il,io)%mz) ! average occupation of up
!	   tm(il,io)%dm(5+i,5+i) = 1.0d-1 * (nds(is)-tm(il,io)%mz) ! down
!	  end do
!	 end do
!	end do
!	return
!	end 	subroutine setupdm
!======================================================================


	end 	module 
