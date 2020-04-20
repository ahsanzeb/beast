
	module modmain
	implicit none

	integer :: nlayers, natomsl, natoms
	integer :: noctl, noct
	integer :: nspecies
	double precision :: a ! lattice constant of a single formula unit cubic cell
	double precision, parameter :: pi = 3.141592653589793d0

	type :: octahedra
		!integer :: ntot
		double precision :: phi
		double precision, dimension(3) :: rb ! position of B atom, absolute position in a unit cell.
		double precision, dimension(3,3) :: xo ! relative to xb, position of oxygen atoms in cubic structure
		double precision, dimension(3,3) :: xor ! relative to xb, position of oxygen atoms after rotation
		double precision, dimension(3,3) :: ro ! absolute position in a unit cell, position of oxygen atoms after rotation
	end type octahedra

	type(octahedra), allocatable, dimension(:,:) :: oct

	double precision, allocatable, dimension(:) :: phi 

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
	write(fnum,'(3G18.10)') a, 0.0, 0.0
	write(fnum,'(3G18.10)') a/2,a/2, 0.0
	write(fnum,'(3G18.10)') 0.0, 0.0, a
	write(fnum,*)
	write(fnum,'("atoms")')
	write(fnum,'(I4,T40," : nspecies")') nspecies

	! all oxygen
	write(fnum,'("''",A,"''",T40," : spfname")') 'O'
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
	!double precision, dimension(2,2) :: rot ! rotation matrix, rot about z axis
	double precision :: phi
	integer :: i,j

	phi = oct(il,io)%phi;
	!rot(1,1) = dcos(phi); rot(1,2) = dsin(phi); 
	!rot(2,1) =-dsin(phi); rot(2,2) = dcos(phi); 

	! convert atom coordinates to cartesian, rotate, and convert back to fractional.

	! avec(:,:) 

	
	do i = 1,2
	 oct(il,io)%xor(i,1) =  
     .  dcos(phi)*oct(il,io)%xo(i,1) + dsin(phi)*oct(il,io)%xo(i,2)
	 oct(il,io)%xor(i,2) =  
     . -dsin(phi)*oct(il,io)%xo(i,1) + dcos(phi)*oct(il,io)%xo(i,2)
	end do

	! not rotated
	oct(il,io)%xor(3,:) = oct(il,io)%xo(3,:)

	! set abs valie of oxygen position after rotation.
	do i = 1,3
	 oct(il,io)%ro(i,:) =	oct(il,io)%rb(:) + oct(il,io)%xor(i,:);
	end do

			 
	return
	end 	subroutine rotate
	!..............................................................








	end 	module modmain
