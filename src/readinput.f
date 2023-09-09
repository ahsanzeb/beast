
! copied a lot from ELK readinput
! given energies in eV, other variables in a.u. 
! energy: U, J, qpol, soc, skparam,...
! length in bohr, 
! Note: a.u. of magnetic field is a really big unit.



! ylm rotatuions codes:
! https://github.com/polarch/Spherical-Harmonic-Transform/tree/master
!https://github.com/tak0kada/ivanic

	module readinput
	use modmain

	implicit none

	logical :: lsys

	contains

	subroutine input()
	implicit none
	! local
	character(40) :: block
	integer :: iostat
	logical :: file_exists
	integer :: il,i,j

	integer :: ios, is, maxsp, ion_d
	double precision :: r1, delB, DelO

	!From monopoles ionic crystal, with d^0 ionic config of TM atoms,
	! and p^6 Oxygen config. [i.e., -4 and +2 charges at TM and O]
	delO = -0.99161304 ! Hartree. 
	DelB = 1.8376130 ! Hartree

	lcage = .false.

	Ham4EDRIXS = .false.
	
	ewalda = 1.d-7; !50.0d0 ! at the moment, k-space sum is discarded so welda has to be infinitesimal to make real space sum equal to the actual sum. of course, this also needs a larger n_r for the convergence! ewaldnr > 5 is okay, converged to 2 decimal palces for NaCl structure, calculated Madelung ==> 1.7457834017323322; converged  is 1.747565 which requires ewaldnr~30.
	
	ewaldnr = 10;
	ewaldnk=10;
	nround = 8;
	noctl = 2;
	lsys = .false.

	lscf = .true.
	ltmgs= .false.
	
	tmnn2= .false.;
	oxnn2= .false.;
	lsoc = .false.;
	lhu = .false.;
	nlayers = 1;
	!phi0 = 0.0d0;
	lspin = .false.;
	lbc = .false.;
	lbands = .false.;
	lusevmat = .false.;
	lgs = .true.;

	temp = 0.025d0 * eV2Har;

	maxscf = 100

	beta0 = 0.3d0;
	betamax = 0.5d0;
	mtype = 2
	toldm = 1.0d-5;

	! pdos
	nwplot = 300;
	mine = -1.0d0; ! -1.0 Hartree
	maxe = +1.0d0; ! +1.0 Hartree
	lpdos = .false.;
	
	reducebf = 0.5;
	bfieldc = 0.0d0;
	fsmtype = 0
	momfix = 0.0d0;
	taufsm = 0.01d0;

	lbfields = .true.
	lesH = .true.
	
	Dcf = 0.0d0;
	Del112 = 0.0d0;
	Del222 = 0.0d0;
	Del224 = 0.0d0;

	nv =2; ! vertices
	np =10; ! points
	allocate(vvlp1d(3,nv))
	vvlp1d(:,1) = (/0.0d0,0.0d0,0.0d0/)
	vvlp1d(:,2) = (/0.5d0,0.0d0,0.0d0/)

	! k-grid
	nk1=1; nk2=1; nk3=1;

	a0 = 7.0; !6.68119
	xsf = .false.
	!tolnns = 0.1 * a0;
	
!-------------------------------------------
! readinput:
	open(50,file='input.in', action='read', iostat=ios)
	if (ios.ne.0) then
	write(*,*)
	write(*,'("Error(readinput): error opening input.in")')
	write(*,*)
	stop
	end if
10	 continue
	read(50,*,end=30) block
	! check for a comment
	if ((scan(trim(block),'!').eq.1).or.
     .           (scan(trim(block),'#').eq.1)) goto 10
	select case(trim(block))

!	case('nlayers')
!	  read(50,*,err=20) nlayers
!	case('phi')
!	  read(50,*,err=20) phi0

	case('system')
	 lsys = .true.
	 read(50,*,err=20) nlayers, nsptm, xsf
	 natoms = nlayers * 8; ! 2 octahedra per layer
	 if(nsptm > nlayers) then
	  write(*,*)"Warning: only first nlayers species will be used."
	  write(*,*)"a layer has one type of atoms, just for simpliciy."
	 endif

	allocate(tm(nlayers,noctl))
	! Hardness and onsite:
	allocate(hardU(0:nsptm))
	allocate(onsite(0:nsptm))
	 
	 ! species index in each layer
	 !allocate(layersp(nlayers))
	 !read(50,*,err=20) (layersp(il), il=1,nlayers)
	 read(50,*,err=20) ((tm(il,i)%is,i=1,2), il=1,nlayers)

	 write(*,*) 'readinp: tm species, layer 1:', (tm(1,i)%is,i=1,2)
	 
	 maxsp = 0
	 do il=1,nlayers
	  maxsp = max(maxsp, tm(il,1)%is)
	  !tm(il,2)%is = tm(il,1)%is ! both atoms belong to the same species
	 end do
	 
	 if(maxsp > nsptm) then
		write(*,'("nsptm smaller than used in layers!")')
		stop
	 endif
	 ! octahedra rotation angle for each layers
	 !allocate(phi(nlayers))
	 !read(50,*,err=20) (phi(il),il=1,nlayers) ! degrees


	case('LatticeConstant','a0')
	   read(50,*,err=20) a0

	case('TiltRotation')
	   read(50,*,err=20) theta !, phii
			phii = theta

	case('RunBeastParam') ! qA, lam, d222, phi, theta
	! assuming that this flag will be at the end of the file so that its 
	! var are already allocated (soc) and never overwritten.
	   read(50,*,err=20) Hub(1)%U, Hub(1)%J, ion_d, 
     .                   qA, soc(1), Del222, phii, theta
	   soc(1) = soc(1)*eV2Har
	   Dcf(5,2) = Del222
	   nds(1) = (6-qa) + ion_d
	  ! set Hub%Fk and allocate Hub%Vee
	  call 	setFk(1, eV2Har) !Hub(1)%U,Hub(1)%J,

	case('Ham4EDRIXS') ! T/F?, cages?, qA, phi, theta
	! assuming that this flag will be at the end of the file so that its 
	! var are never overwritten.
	   read(50,*,err=20) Ham4EDRIXS, lcage, qA, phii, theta
			
	case('Bfield')
		read(50,*,err=20) bfieldc(1:3)
		
	case('OnsiteBfield')
	  call sysfirst()
	  do il=1,nlayers
	   read(50,*,err=20) tm(il,1)%bext(1:3), tm(il,2)%bext(1:3)
	  end do
	case('EwaldNr')
	   read(50,*,err=20) ewaldnr ! ewalda, ewaldnr, ewaldnk

	case('RoundStrux')
	   read(50,*,err=20) nround

	case('FixedSpinMomentType')
	   read(50,*,err=20) fsmtype

	case('FixedMoment')
	   read(50,*,err=20) momfix(1:3)		

	case('OnsiteFixedMoment')
	  call sysfirst()
	  do il=1,nlayers
	   read(50,*,err=20) tm(il,1)%mfix(1:3), tm(il,2)%mfix(1:3)
	  end do

	case('taufsm')
	 read(50,*,err=20) taufsm
	  if (taufsm.lt.0.d0) then
	   write(*,*)
	   write(*,'("Error(readinput): taufsm < 0 : ",G18.10)') taufsm
	   write(*,*)
	   stop
	  end if

	case('UseSTATE','usestate', 'usedata')
	 read(50,*,err=20) lusevmat
	 if(lusevmat) then
	  lgs = .false.
	  write(*,'(a)')"UseVMAT = T, assuming we dont want GS...?"
	 endif
	case('SpinPolarised')
	 read(50,*,err=20) lspin

	case('reducebf')
	 read(50,*,err=20) reducebf

	case('CalculationType','Hes')
	 read(50,*,err=20) lesH, lscf, ltmgs 
	 !if(ltmgs) lscf = .false.

	case('CrystalField')
	 read(50,*,err=20) lcage, Del112, Del222, Del224
	 Dcf(2,1) = Del112 ! Oxygen
	 Dcf(5,2) = Del222 ! TM
	 Dcf(7,2) = Del224 ! TM
	 !Dcf = Dcf * eV2Har 

	! Delta = 1.d5 gives 7349.89 splitting in bands
	! to make 1:2, rescale the Delta by: 
	!Dcf = Dcf * 2*1.d5/7349.89d0
	!write(*,*)"readinp: CrystalField and SOC values:"
	!write(*,*)"readinp: 2 level sys: V produces 2V splitting"
	 
	 !read(50,*,err=20) r1
	 !read(50,*,err=20) (Dcf(is,1), is=1,7) ! Oxygen
	 !read(50,*,err=20) (Dcf(is,2), is=1,7) ! all TM
	 !Dcf = Dcf*r1

	!......................................................................
	case('TMSpecies', 'Species')
	 call sysfirst()
	 read(50,*,err=20) lhu, lsoc
	 if(lhu) then
	  allocate(Hub(nsptm))
	  ! species hubbard U and J: given in eV
	  read(50,*,err=20) (Hub(is)%U, is=1,nsptm)
	  read(50,*,err=20) (Hub(is)%J, is=1,nsptm)
	  ! set Hub%Fk and allocate Hub%Vee
	  	do is=1,nsptm
	   call 	setFk(is, eV2Har) !Hub(is)%U,Hub(is)%J,
	   allocate(Hub(is)%Vee(5,5,5,5))
	   !write(*,'(a,i5,5f10.5)')"is, Hub(is)%fk = ",is,Hub(is)%fk
	 end do
	else
	 read(50,*,err=20) 
	 read(50,*,err=20) 
	endif	 
	 ! species spin-orbit
	 allocate(soc(nsptm))
	 soc = 0.0d0
	 if(lsoc) then
	  read(50,*,err=20) (soc(il), il=1,nsptm) ! SOC given in eV
	  	! Convert to Hartree from eV
	  soc = soc * eV2Har
	else
	 read(50,*,err=20) 
	endif
	 allocate(nds(nsptm))
	 nds = 0.0d0
	 read(50,*,err=20) (nds(il), il=1,nsptm)
	 read(50,*,err=20) qa
	 
	 qtot = 0.0d0
	 ! total electrons in the unit cell
	 do il=1,nlayers ! 2.0* for two octaherda per layer
	  do i=1,2
	  qtot=qtot + ( (nds(tm(il,i)%is) + 0.0d0) + ! +0 for TM s-orbital
     .                          3.0d0*4.0d0 + qa) !2.0d0) ! +2 for Ca/Sr A-site
		enddo
	 end do
	!......................................................................
	
	case('BandCharacter','bc')
	 read(50,*,err=20) lbc

	case('pdos')
	 read(50,*,err=20) lpdos, mine, maxe, nwplot
	!case('nwplot')
	! read(50,*,err=20) nwplot

	case('temp')
	  read(50,*,err=20) temp  
	  temp = temp * eV2Har
	case('maxscf')
	  read(50,*,err=20) maxscf 
	case('beta0') ! scf mixing weight
	  read(50,*,err=20) beta0
	case('betamax') ! scf mixing weight max
	  read(50,*,err=20) betamax
	case('MixType','mtype','mixtype') ! scf mixing type
	  read(50,*,err=20) mtype ! 0: Linear mixing, 1: adoptive
	case('scftol','SCFtol') ! scf mixing weight max
	  read(50,*,err=20) toldm

	case('skparam') ! SK parameters given in eV
	 call sysfirst()
	! single set of SK param for all TM species? 2nd nns of TM? 2nd nns of O?
	 read(50,*,err=20) singlesk, tmnn2, oxnn2
	 allocate(skbo(nsptm,2)) 	! 2: sigma_pd, pi_pd
	 allocate(skbb(nsptm,nsptm,3))	! 3: sigma_dd, pi_dd, delta_dd
	 skbo= 0.0d0; 	 skbb= 0.0d0; 	 skoo= 0.0d0;
	 if(singlesk) then
	   read(50,*,err=20) skbo(1,1:2) ! TM-O / O-TM 1st nns
	   read(50,*,err=20) skbb(1,1,1:3) ! TM-TM 2nd nns
	   read(50,*,err=20) skoo(1:2) ! O-O 2nd nns	  

	  ! eV to Hartree:
	  skbo = skbo * eV2Har;
	  skbb = skbb * eV2Har;
	  skoo = skoo * eV2Har;

	 ! set full arrays:
	 do i=1,nsptm
	 	skbo(i,1:2) = skbo(1,1:2)
	  do j=1,nsptm
	   skbb(i,j,1:3)=skbb(1,1,1:3)
	  end do
	 end do

	 else
	 	write(6,'("Error(readinp): singlesk=F not implemented yet.")')	
	 	stop  
	 endif	 
!	 do i=1,ntmsp1
!	  read(50,*,err=20) skbo(i,1:2) ! TM-O / O-TM 1st nns
!		do j=1,ntmsp1
!	   read(50,*,err=20) skbb(i,j,1:3) ! TM-TM 2nd nns
!		end do
!	 end do 
!	 read(50,*,err=20) skbo(i,1:2) ! O-O 2nd nns


	case('HardnessAndOnsite') ! Onsite energies given in eV
	call sysfirst()
	hardU = 0.0d0
	read(50,*,err=20) (hardU(i),i=1,nsptm) ! eta_TM1, eta_TM2, ...	
	onsite= 0.0d0; 
	!onsite(0) = +100.0 
	!write(*,*)'readinp: test... setting O onsite to +100.0'
	read(50,*,err=20) (onsite(i),i=1,nsptm)

	! fake splitting in B1:
	!read(50,*,err=20) (delta(i),i=1,5)
	!delta = delta * eV2Har;
	
	! eV to Ryd
	onsite = onsite * eV2Har;
	hardU = hardU * eV2Har;

	! renormalised onsite energies due to electrostatic interaction 
	! and charge transfers that are actually used in the code:
	!onsite(0	) = -DelO; !to mk  effective onsite of O ==> 0.0d
	!do i=1,nsptm ! Oxygen hardness=0
	!	onsite(i) = onsite(i) - delB + 4.d0*hardU(i);
	!	write(*,*) 'i, onsite(i) = ', i, onsite(i)
	!end do
	

!	case('hardness','Hardness', 'eta') ! hardness energies given in eV
!	read(50,*,err=20) (hardU(i),i=0,nsptm) ! eta_O, eta_TM1, eta_TM2, ...
!	! eV to Ryd


	case('USOCLoops') ! Onsite energies given in eV
	call sysfirst()
	read(50,*,err=20) (isploop(i),i=1,2)
	if(maxval(isploop) > nsptm ) then
	 write(*,*)"Warning(readinput): USOCLoops: max(isploop) > nsptm"
	 if(isploop(1) .ne. 1) then
	  stop "Error(readinput): isploop(1) .ne. 1 "
	 endif
	endif
	write(*,*) "reading: uloop(1:6). only two tm species?"
	read(50,*,err=20) (uloop(i), i=1,6)
	read(50,*,err=20) (sloop(i), i=1,6)
	uloop = dble(uloop)
	sloop = dble(sloop)
	! to stop the loops in case we have increments = 0 in the input file.
	 if(dabs(uloop(3))<1.0d-10) uloop(3) = 10.0d0
	 if(dabs(uloop(6))<1.0d-10) uloop(6) = 10.0d0
	 if(dabs(sloop(3))<1.0d-10) sloop(3) = 10.0d0
	 if(dabs(sloop(6))<1.0d-10) sloop(6) = 10.0d0
	

	case('kgrid')
	read(50,*,err=20) nk1,nk3
	


	case('plot1d') ! bandlines
	  read(50,*,err=20) lbands, nv,np ! number of bandline segments, total k-points
	  if (nv.lt.1) then
	    write(*,*)
	    write(*,'("Error: nv < 1 : ",I8)') nv
	    write(*,*)
	    stop
	  end if
	  if (np.lt.nv) then
	    write(*,*)
	    write(*,'("Error: np < nv : ",2I8)') np,nv
	    write(*,*)
	    stop
	  end if
	  if (allocated(vvl)) deallocate(vvl)
	  allocate(vvl(3,nv))
	 do i=1,nv
	    read(50,*,err=20) vvl(:,i)
	 end do
!-------------------------------------------
	case('')
	goto 10
	case default
	write(*,*)
	write(*,'("Error: invalid block name : ",A)') trim(block)
	write(*,*)
	stop
	end select
	goto 10
20	 continue
	write(*,*)
	write(*,'("Error: error reading from input.in")')
	write(*,'("Problem in ''",A,"'' block")') trim(block)
	write(*,'("Check input convention in manual!! :D")')
	write(*,*)
	stop
30	 continue
	close(50)


	if(.not. allocated(onsite)) then
	 allocate(onsite(0:nsptm))
	 onsite = 0.0d0
	 write(*,*)
     . "Warning: onsite not provided, using 0.0 for all atoms!"
	endif


	! set nspin
	if(lspin) then
	 nspin=2
	elseif(lsoc .or. lhu) then
	 nspin = 2;
	 write(*,*)"main: setting nspin=2 as lsoc/lhu = T"
	else
	 nspin=1
	endif

	if(xsf) then
	 write(*,*) "========>>>>> reading XSF file..... "
	 ! read structure info from a file and assign atom index 
	 ! according to the location, similar to our octahedra convention.
	 call readstruc()
	endif
	
	return
	end 	subroutine input
!-------------------------------------------

	subroutine sysfirst()
	implicit none
	if(.not. lsys) then
	 write(6,*)"Please give system details first in input.in"
	 write(6,*)"Sorry, I need to set nsptm before many others!"
	 stop
	endif
	end 	subroutine sysfirst
!-------------------------------------------



!=====================================================================
	subroutine readstruc()
	implicit none
	! local
	character(40) :: block
	integer :: iostat
	logical :: file_exists
	integer :: il, io

	integer :: nasio, i,j,k, ios
	character(2), allocatable, dimension(:) :: splabel, splabel2
	double precision, allocatable, dimension(:,:) :: xasio,xasiof
	double precision, dimension(3) :: r
	double precision :: d2min, tol, ang2au, kz, kzhalf
	double precision, allocatable, dimension(:) :: d2j
	ang2au = 1.8897259886d0;
	
	open(50,file='crystal.xsf', action='read', iostat=ios)
	if (ios.ne.0) then
	write(*,*)
	write(*,'("Error(readinput): error opening crystal.xsf")')
	write(*,*)
	stop
	end if
10	 continue
	read(50,*,end=30) block
	! check for a comment
	if ((scan(trim(block),'!').eq.1).or.
     .           (scan(trim(block),'#').eq.1)) goto 10
	select case(trim(block))

!-------------------------------------------
	case('CRYSTAL') ! XSF format, lattice vectors
		continue
	case('PRIMVEC') ! XSF format, lattice vectors
		do i=1,3
	   read(50,*,err=20) avec(:,i)
	  end do

	  avec = avec * ang2au;

		a0 = dsqrt(0.5d0*(avec(1,1)**2 + avec(2,1)**2+avec(3,1)**2))
		write(*,*) "=====>>> lattice scale a0 = a1/sqrt(2) = ", a0
	case('PRIMCOORD') ! XSF format, atomic coordinates
		read(50,*,err=20) nasio ! Sr Ir O: structure na_SIO
		if(nasio /= nlayers*10) then
	   write(*,*)
     . "Error: XSF inconsistent with nlayers in input.in"
	   write(*,*) 'Natoms in XSF = ', nasio
	   write(*,*) 'nlayers*10 = ', nlayers*10
	   stop
	  endif
		allocate(splabel(nasio), xasio(nasio,3))! ,xasiof(nasio,3))
		!allocate(splabel2(nasio))

		do i=1,nasio
	   read(50,*,err=20) splabel(i), xasio(i,1:3) ! assuming Sr, Ir and O. will be reassigned using info in system block
	   !write(*,*) splabel(i), xasio(i,1:3)
	  end do

	  xasio = xasio * ang2au;
!-------------------------------------------
	case('')
	goto 10
	case default
	write(*,*)
	write(*,'("Error: invalid block name : ",A)') trim(block)
	write(*,*)
	stop
	end select
	goto 10
20	 continue
	write(*,*)
	write(*,'("Error: error reading from input.in")')
	write(*,'("Problem in ''",A,"'' block")') trim(block)
	write(*,'("Check input convention in manual!! :D")')
	write(*,*)
	stop
30	 continue
	close(50)




! ---------------------------------------------------------------------
! some variables that were previously set in getstructure()
!
	noctl = 2;
	noct = noctl*nlayers;
	natoms = noct*4;
	a = a0; !7.0d0;
	nspecies = nsptm;
	
	norbtm = 5; norbtms = norbtm*nspin;
	norbo = 3; norbos = norbo*nspin;

	ntot = noct*norbtms + noct*3*norbos;
	ntottm = noct*norbtms;

	write(*,*)'ntot, ntottm, qtot = ', ntot, ntottm, qtot

	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)
! ---------------------------------------------------------------------

	allocate(posA(nlayers*2,3))
	allocate(pos(nlayers*8,3))

	write(*,*) "natoms, nlayers = ",natoms, nlayers 
	! set positions of Sr/A atoms:
	j=0
	do i=1,nasio
		if(splabel(i)=='Sr') then
			j = j+1
			posA(j,:) = xasio(i,:);
		endif
	end do
	if(j < nlayers*2) stop 'Number of Sr/A atoms found < nlayers*2 '

	! set positions of Ir/B atoms:
	j=1; k=0
	do i=1,nasio
		if(splabel(i)=='Ir') then
			k = k + 1;
			pos(j,:) = xasio(i,:);
			j = j+4
		endif
	end do
	if(k < nlayers*2) stop 'Number of Ir/B atoms found < nlayers*2 '


	! set positions of O atoms:
	k=0; j=0;
	do i=1,nasio
		if(splabel(i)=='O') then
			k = k + 1;
			j = j + 1;
			if(mod(j-1,4)==0) j = j + 1; ! skip TM atom indexes: 1,5,9,13,.... 
			pos(j,:) = xasio(i,:);
		endif
	end do
	if(k < nlayers*2*3) stop 'Number of O atoms found < nlayers*6 '

	!write(*,*)  ' ===========>>>> before return.... '
	!do i=1,natoms
	! write(*,'(i5,3f10.5)') i, pos(i,:)
	!end do


	return 
	! lines below are disabled
	! not assigning indexes above according to octahedra or layers. 






	! start from:
	! position of Ir/O atoms in a cubic system: pos in fractional coordinates
	! all atoms inside the unit cell, similar to what spacegroup generates
	kz = 1.d0/dble(nlayers);
	kzhalf = kz/2.0d0;
	do il=1,nlayers
		i = (il-1)*8
		! TM 1
	  pos(i+1,:) = (/ 0.5d0, 0.0d0, (il-1)*kz /);
		! Oxygen
	  pos(i+2,:) = (/ 0.75d0, 0.75d0, (il-1)*kz  /)
	  pos(i+3,:) = (/ 0.75d0, 0.25d0, (il-1)*kz  /)
	  pos(i+4,:) = (/ 0.5d0, 0.0d0, (il-1)*kz + kzhalf  /)

		! TM 2
	  pos(i+5,:) = (/ 0.0d0, 0.5d0, (il-1)*kz /);
		! Oxygen
	  pos(i+6,:) = (/ 0.25d0, 0.25d0, (il-1)*kz  /)
	  pos(i+7,:) = (/ 0.25d0, 0.75d0, (il-1)*kz  /)
	  pos(i+8,:) = (/ 0.0d0, 0.5d0, (il-1)*kz + kzhalf  /)
	enddo



	if(1==0)then
	! fractional to cartesian:
		write(*,*)'CUBIC structure pos:'
	do i=1,natoms,4
	  call r3mv(avec,pos(i,:)+(/0.5,0.5,0.5/), r)
		write(*,'(a,3f10.5)') 'Sr', r
	end do
		do i=1,natoms,4
	    call r3mv(avec,pos(i,:), r)
			write(*,'(a,3f10.5)') 'Ir',r
		end do
		do i=1,natoms
			if(mod(i-1,4) /= 0) then
	     call r3mv(avec,pos(i,:), r)
			 write(*,'(a,3f10.5)') 'O', r
			endif
		end do

	 stop 'readinp: testing... '
	 return
	endif





	!write(*,*) 'pos: '
	!do i=1,natoms
	!	write(*,'(3f10.4)') pos(i,:)
	!end do
	!write(*,*)'-------------------'




	! cartesian to fractional
	call r3minv(avec,ainv) ! inverse matrix
	do i=1, nasio
		call r3mv(ainv,xasio(i,:),xasiof(i,:))
		write(*,'(a,3f10.4)') splabel(i), xasiof(i,:)
	end do



	allocate(d2j(nasio))
	! find actual positions (of TM/O) but with the reference structure indexing.
	! i.e., in every layaer, 1-8 atoms: TM-O-O-O-TM-O-O-O
	do i=1, natoms ! reference structure
		d2j(:) = 1.0d5; ! large number
		do j=1,nasio ! input structure
			if(splabel(j) /= 'Sr') then ! Ir or O
			 r = pos(i,:) - xasiof(j,:)
			 d2j(j) = r(1)**2 +	r(2)**2 + r(3)**2;
			endif
		end do
		j = minloc(d2j,1)
		write(*,'(a,2i4,1f8.3)')'i,j, d2j(j) = ',i,j, d2j(j)
		pos(i,:) = xasiof(j,:);
	end do

	return
	end subroutine readstruc
!=====================================================================

	end module readinput
