
! copied a lot from ELK readinput
! given energies in eV, other variables in a.u. 
! energy: U, J, qpol, soc, skparam,...
! length in bohr, 
! Note: a.u. of magnetic field is a really big unit.

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

	integer :: ios, is, maxsp
	double precision :: r1


	noctl = 2;
	lsys = .false.

	tmnn2= .false.;
	oxnn2= .false.;
	lsoc = .false.;
	lhu = .false.;
	nlayers = 1;
	!phi0 = 0.0d0;
	lspin = .false.;
	lbc = .false.;
	lusevmat = .false.;
	lgs = .true.;

	temp = 0.025d0;

	maxscf = 100

	beta0 = 0.3d0;
	betamax = 0.5d0;
	mtype = 2
	toldm = 1.0d-5;

	nwplot = 100;

	reducebf = 0.5;
	bfieldc = 0.0d0;
	fsmtype = 0
	momfix = 0.0d0;
	taufsm = 0.01d0;

	lbfields = .true.

	nv =2; ! vertices
	np =10; ! points
	allocate(vvlp1d(3,nv))
	vvlp1d(:,1) = (/0.0d0,0.0d0,0.0d0/)
	vvlp1d(:,2) = (/0.5d0,0.0d0,0.0d0/)

	! k-grid
	nk1=1; nk2=1; nk3=1;

!-------------------------------------------
! readinput:
	open(50,file='input.in', action='read')
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
	 read(50,*,err=20) nlayers, nsptm
	 if(nsptm > nlayers) then
	  write(*,*)"Warning: only first nlayers species will be used."
	  write(*,*)"a layer has one type of atoms, just for simpliciy."
	 endif

	allocate(tm(nlayers,noctl))
	 
	 ! species index in each layer
	 !allocate(layersp(nlayers))
	 !read(50,*,err=20) (layersp(il), il=1,nlayers)
	 read(50,*,err=20) (tm(il,1)%is, il=1,nlayers)
	 maxsp = 0
	 do il=1,nlayers
	  maxsp = max(maxsp, tm(il,1)%is)
	  tm(il,2)%is = tm(il,1)%is ! both atoms belong to the same species
	 end do
	 
	 if(maxsp > nsptm) then
		write(*,'("nsptm smaller than used in layers!")')
		stop
	 endif
	 ! octahedra rotation angle for each layers
	 allocate(phi(nlayers))
	 read(50,*,err=20) (phi(il),il=1,nlayers) ! degrees

	case('Bfield')
		read(50,*,err=20) bfieldc(1:3)
		
	case('OnsiteBfield')
	  call sysfirst()
	  do il=1,nlayers
	   read(50,*,err=20) tm(il,1)%bext(1:3), tm(il,2)%bext(1:3)
	  end do

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
	  	 ! Convert to Rydberg from eV
	  	 Hub(is)%U = Hub(is)%U * eV2Ryd
	  	 Hub(is)%J = Hub(is)%J * eV2Ryd

	  	 Hub(is)%fk(:) = 0.0d0
	   Hub(is)%fk(0)=Hub(is)%U
	   ! r1 = F(4)/F(2), see PRB 52, R5467 (1995)
	   r1=0.625d0
	   Hub(is)%fk(2)=(14.d0*Hub(is)%J)/(1.d0+r1)
	   Hub(is)%fk(4)=Hub(is)%fk(2)*r1
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
	  	! Convert to Rydberg from eV
	  soc = soc * eV2Ryd
	else
	 read(50,*,err=20) 
	endif
	 allocate(nds(nsptm))
	 nds = 0.0d0
	 read(50,*,err=20) (nds(il), il=1,nsptm)
	 qtot = 0.0d0
	 ! total electrons in the unit cell
	 do il=1,nlayers ! 2.0* for two octaherda per layer
	  qtot=qtot + 2.0d0*( (nds(tm(il,1)%is) + 2.0d0) + ! +2 for TM s
     .                          3.0d0*4.0d0 + 2.0d0) ! +2 for Ca/Sr A-site
	 end do
	!......................................................................
	
	case('BandCharacter','bc')
	 read(50,*,err=20) lbc

	case('pdos')
	 read(50,*,err=20) lpdos
	case('nwplot')
	 read(50,*,err=20) nwplot

	case('temp')
	  read(50,*,err=20) temp  
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
	  if(tmnn2)then
	  ! TM-O / O-TM 1st nns; ! TM-TM 2nd nns
	   read(50,*,err=20) skbo(1,1:2), skbb(1,1,1:3) 
		else
	   read(50,*,err=20) skbo(1,1:2)
	  endif
	  if(oxnn2) then
	 	 read(50,*,err=20) skoo(1:2) ! O-O 2nd nns
	  else
	 	 read(50,*,err=20) ! avoid error if data present but switch off using oxnn2=T/F
	  endif

	  ! eV to Rydberg:
	  skbo = skbo * eV2Ryd;
	  skbb = skbb * eV2Ryd;
	  skoo = skoo * eV2Ryd;
	  
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


	case('onsite') ! Onsite energies given in eV
	call sysfirst()
	allocate(onsite(0:nsptm))
	onsite= 0.0d0
	read(50,*,err=20) (onsite(i),i=0,nsptm)
	! eV to Ryd
	onsite = onsite * eV2Ryd;

	case('kgrid')
	read(50,*,err=20) nk1,nk3
	


	case('plot1d') ! bandlines
	  read(50,*,err=20) nv,np ! number of bandline segments, total k-points
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




	end module readinput
