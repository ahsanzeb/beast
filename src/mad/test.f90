program test
use modmain
use readinput

use esvar
use estatic, only: setmadvar, initmadelung, getatomic
use mtbmpol


	integer :: il, io, i, ia
	integer :: ik,ib
	double precision, dimension(3) :: kvec








	! read input file 'input.in'
	call input()

	! set nspin
	if(lspin) then
	 nspin=2
	elseif(lsoc .or. lhu) then
	 nspin = 2;
	 write(*,*)"main: setting nspin=2 as lsoc/lhu = T"
	else
	 nspin=1
	endif



		
	!write(*,'(a)') "Give Number of layers & phi and hit enter:"
	!read(*,*) nlayers, phi0
	!nlayers=5; phi0=10.0d0

	if(any(phi(:) > 45.0)) then
	 write(*,'("O atoms have to cross to make this big rotation!")')
	 stop 
	else
	 write(*,'(a)') 'Octraherda size will be rescaled to keep "a" fixed.'
	endif

	!tmnn2 = .true.;
	!oxnn2 = .true.;


	!nlayers = 2;
	noctl = 2;
	noct = noctl*nlayers;
	natoms = noct*4;
	a = 7.0d0;


	!nsptm =1;
	nspecies = nsptm;
	
	norbtm = 5; norbtms = norbtm*nspin;
	norbo = 3; norbos = norbo*nspin;

	ntot = noct*norbtms + noct*3*norbos;
	ntottm = noct*norbtms;

	write(*,*)'ntot, ntottm, qtot = ', ntot, ntottm, qtot
	!i = noct*1 + noct*6*norbos
	!write(*,*)'TM d^1: nelec = ', i
	!write(*,*)'if only one spin, filled bands ~ ',0.5*i



	
	allocate(oct(nlayers,noctl))
	!allocate(phi(nlayers))
	! phi 
	!phi = 0.0d0;

	! read/set phi:
	do il=1,nlayers
	 !if(mod(il,2)==0) then
	 ! phi(il) = phi0*pi/180.0d0;
	 !else
	 ! phi(il) = -phi0*pi/180.0d0;
	 !endif
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
	  	oct(il,io)%xo(1,:) = (/0.5d0,0.0d0,0.0d0/)*a
 !    .   +(/rand(0),rand(0),rand(0)/)*a
	  	oct(il,io)%xo(2,:) = (/0.0d0,0.5d0,0.0d0/)*a
!     .   +(/rand(0),rand(0),rand(0)/)*a
	  	oct(il,io)%xo(3,:) = (/0.0d0,0.0d0,0.5d0/)*a
!     .   +(/rand(0),rand(0),rand(0)/)*a
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

	! set nearest neighbours:
	call settmnn1()
	call setoxnn1()
	if(tmnn2) call settmnn2()
	if(oxnn2) call setoxnn2()

	! atom to orbitals map
	call mapatom2orbs()
	! atom to species (TM) map
	call mapatom2species()

 call reciplat(avec,bvec,omega,omegabz)

 ldip = 2;
 !--------------------------------------------------------
 ! once before SCF cycle starts:
 call setmadvar()
 call initmadelung()
 !--------------------------------------------------------

 !--------------------------------------------------------
 ! dueing each SCF iteration:
 ! calculate rhoc etc to be used for Qmpol
 call getatomic()
 ! calculate Qmpol
 call tbmpol()
 ! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
 ! call tbeseld()
 !--------------------------------------------------------


 


 write(*,*)'test program ended...'
end program test
