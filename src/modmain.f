
	module modmain
	implicit none

	integer :: nlayers, natomsl, natoms
	integer :: noctl, noct
	integer :: nspecies
	double precision :: a ! lattice constant of a single formula unit cubic cell
	double precision, parameter :: pi = 3.141592653589793d0

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
	else ! i==2, atom along y
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
	!..............................................






	end 	module modmain
