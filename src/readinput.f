
! copied a lot from ELK readinput

	module readinput
	use modmain

	implicit none



	contains

	subroutine input()
	implicit none
	! local
	character(40) :: block
	integer :: iostat
	logical :: file_exists
	logical :: parameters
	double precision :: memory
	integer :: i, nevdm, nexdm, anev
	logical :: fixrhoexdm

	integer :: ios

	nlayers = 1;
	phi0 = 0.0d0;
	
	nv =2; ! vertices
	np =10; ! points
	allocate(vvlp1d(3,nv))
	vvlp1d(:,1) = (/0.0d0,0.0d0,0.0d0/)
	vvlp1d(:,2) = (/0.5d0,0.0d0,0.0d0/)
	
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

	case('nlayers')
	  read(50,*,err=20) nlayers
	case('phi')
	  read(50,*,err=20) phi0
	
	

	case('bandlines')
	  read(50,*,err=20) nv,np
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

	return
	end 	subroutine input
!-------------------------------------------





	end module readinput
