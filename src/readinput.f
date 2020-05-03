
! copied a lot from ELK readinput

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

	integer :: ios, is
	double precision :: r1

	lsys = .false.

	tmnn2= .false.;
	oxnn2= .false.;
	lsoc = .false.;
	lhu = .false.;
	nlayers = 1;
	!phi0 = 0.0d0;
	lspin = .false.

	temp = 0.025d0;

	maxscf = 100
	beta = 0.5d0; 

	reducebf = 0.5;
	bfield = 0.0d0;

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
	 
	 ! species index in each layer
	 allocate(layersp(nlayers))
	 read(50,*,err=20) (layersp(il), il=1,nlayers)
	 if(maxval(layersp) > nsptm) then
		write(*,'("nsptm smaller than used in layers!")')
		stop
	 endif
	 ! octahedra rotation angle for each layers
	 allocate(phi(nlayers))
	 read(50,*,err=20) (phi(il),il=1,nlayers) ! degrees
	 allocate(nds(nsptm))
	 nds = 0.0d0
	 read(50,*,err=20) (nds(il), il=1,nsptm)
	 qtot = 0.0d0
	 ! total electrons in the unit cell
	 do il=1,nlayers ! 2.0* for two octaherda per layer
	  qtot=qtot + 2.0d0*( (nds(layersp(il)) + 2.0d0) +  ! +2 for TM s
     .                          3.0d0*4.0d0 + 2.0d0) ! +2 for Ca/Sr site
	 end do

	case('SpinPolarised')
	 read(50,*,err=20) lspin

	case('bf')
	 read(50,*,err=20) bfield

	case('reducebf')
	 read(50,*,err=20) reducebf

	case('SpinOrbit','SOC','soc')
	 call sysfirst()
	 read(50,*,err=20) lsoc
	 ! species spin-orbit
	 allocate(soc(nsptm))
	 soc = 0.0d0
	 if(lsoc) then
	  read(50,*,err=20) (soc(il), il=1,nsptm)
	else
	 read(50,*,err=20) 
	endif

	case('HubbardUJ', 'UJ')
	 call sysfirst()
	 read(50,*,err=20) lhu
	 if(lhu) then
	  allocate(Hub(nsptm))
	  ! species hubbard U and J
	  read(50,*,err=20) (Hub(is)%U, is=1,nsptm)
	  read(50,*,err=20) (Hub(is)%J, is=1,nsptm)
	  ! set Hub%Fk and allocate Hub%Vee
	  	do is=1,nsptm
	  	 Hub(is)%fk(:) = 0.0d0
	   Hub(is)%fk(0)=Hub(is)%U
	   ! r1 = F(4)/F(2), see PRB 52, R5467 (1995)
	   r1=0.625d0
	   Hub(is)%fk(2)=(14.d0*Hub(is)%J)/(1.d0+r1)
	   Hub(is)%fk(4)=Hub(is)%fk(2)*r1
	   allocate(Hub(is)%Vee(5,5,5,5))
	   !write(*,'(a,i5,5f10.5)')"is, Hub(is)%fk = ",is,Hub(is)%fk
	 end do
	 call setUrUz() ! sets Uz and Ur matrices
	else
	 read(50,*,err=20) 
	 read(50,*,err=20) 
	endif	 

	case('temp')
	  read(50,*,err=20) temp  
	case('maxscf')
	  read(50,*,err=20) maxscf 
	case('beta') ! scf mixing weight
	  read(50,*,err=20) beta

	case('skparam')
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


	case('onsite')
	call sysfirst()
	allocate(onsite(0:nsptm))
	onsite= 0.0d0
	read(50,*,err=20) (onsite(i),i=0,nsptm)

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
