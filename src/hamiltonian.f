
	module hamiltonian
	use modmain

	implicit none





	contains

!----------------------------------------------------------------
	subroutine getHk(kvec,hk)
	implicit none
	!integer, intent(in):: ik
	double precision, dimension(3), intent(in) :: kvec
	double complex, dimension(ntot,ntot), intent(out) :: hk
	double complex :: eikr
	double precision :: kr
	integer :: i1,i2, j1,j2,is,js,ia,ja
	integer :: il,io,ii,i,k

	hk(:,:) = 0.0d0
	!...................................................................
	! TM-O  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl

	  !is = tm(il,io)%is; ! species index of TM atom
	  ! atom index, and orbital range
	  ia = tm(il,io)%ia;
	  i1 = atom2orb(1,ia); i1 = atom2orb(2,ia); ! orbital ranges
	  !i1 = tm(il,io)%i; i2 = tm(il,io)%j;
	  !write(*,*) 'il,io, ia = ',il,io, ia

	  do k = 1,6 ! 1st nns O
	   !write(*,*) 'k, tm(il,io)%nn1(k)%ia = ',k,tm(il,io)%nn1(k)%ia
	   ja = tm(il,io)%nn1(k)%ia;	   
	   j1 = atom2orb(1,ja); j2 = atom2orb(2,ja); ! orbital ranges

	   !write(*,*) 'ja, j1,j2 = ',ja,j1,j2 

	   kr = DOT_PRODUCT(kvec,tm(il,io)%nn1(k)%dr)
	   eikr = dcmplx(dcos(kr),dsin(kr))
	   hk(i1:i2,j1:j2) = hk(i1:i2,j1:j2) !+ tm(il,io)%nn1(k)%h * eikr
	  end do ! k
	 end do ! io
	end do ! il

	write(*,*)'-------------------- a'

	!...................................................................
	! 	TM-TM (2nd neighbours) 
	!...................................................................
	if(tmnn2) then
	do il=1,nlayers
	 do io=1, noctl
	  	ia = tm(il,io)%ia;
	  i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
	  !write(*,*) 'il,io, ia = ',il,io, ia
	  !write(*,*) 'i1,i2 = ',i1,i2

	  !i1 = tm(il,io)%i; i2 = tm(il,io)%j;
	  !	write(*,*)'is,norbtm = ',is, norbtm
	  do k = 1,6 ! 2nd nns TM
	   ja = tm(il,io)%nn2(k)%ia;	   

	   !write(*,*) 'k, tm(il,io)%nn2(k)%ia = ',k,tm(il,io)%nn2(k)%ia

	   j1 = atom2orb(1,ja); j2 = atom2orb(2,ja); ! orbital ranges

	   !write(*,*) 'ja, j1,j2 = ',ja,j1,j2 

	  ! write(*,*) tm(il,io)%nn2(k)%h

	   kr = DOT_PRODUCT(kvec,tm(il,io)%nn2(k)%dr)
	   eikr = dcmplx(dcos(kr),dsin(kr))
	   hk(i1:i2,j1:j2) = hk(i1:i2,j1:j2) + tm(il,io)%nn2(k)%h *eikr
	  end do ! k
	 end do ! io
	end do ! il
	endif
	!...................................................................

	write(*,*)'-------------------- b'


	!...................................................................
	! O-TM  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  ia = ox(il,io,ii)%ia;
	  	i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
	  !i1 = ox(il,io,ii)%i; i2 = ox(il,io,ii)%j; 
	  do k = 1,2 ! 1st nns TM
	   ja = ox(il,io,ii)%nn1(k)%ia;	   
	   j1 = atom2orb(1,ja); j2 = atom2orb(2,ja); ! orbital ranges
	   kr = DOT_PRODUCT(kvec,ox(il,io,ii)%nn1(k)%dr)
	   eikr = dcmplx(dcos(kr),dsin(kr))
	   hk(i1:i2,j1:j2) = hk(i1:i2,j1:j2) + ox(il,io,ii)%nn1(k)%h *eikr
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	!...................................................................
	write(*,*)'-------------------- c'

	!...................................................................
	! O-O  (2nd neighbours)
	!...................................................................
	if(oxnn2) then
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  ia = ox(il,io,ii)%ia;
	  	i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
	  !i1 = ox(il,io,ii)%i; i2 = ox(il,io,ii)%j; 
	  !write(*,*) 'il,io,ii, ia = ',il,io,ii, ia
	  !write(*,*) 'i1,i2 = ',i1,i2
	  do k = 1,8 ! 2nd nns O
	   ja = ox(il,io,ii)%nn2(k)%ia;	   
	   !write(*,*) 'k, ja, j1,j2 = ',k, ja,j1,j2 
	   j1 = atom2orb(1,ja); j2 = atom2orb(2,ja); ! orbital ranges
	   kr = DOT_PRODUCT(kvec,ox(il,io,ii)%nn2(k)%dr)
	   eikr = dcmplx(dcos(kr),dsin(kr))
	   hk(i1:i2,j1:j2) = hk(i1:i2,j1:j2) + ox(il,io,ii)%nn2(k)%h *eikr
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	endif
	!...................................................................

	write(*,*)'-------------------- d'




	
	return
	end 	subroutine getHk
!----------------------------------------------------------------






	end module hamiltonian
	
