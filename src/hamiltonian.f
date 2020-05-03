
	module hamiltonian
	use modmain

	implicit none


	contains

!----------------------------------------------------------------
	subroutine getHk(ik,kvec,hk,iscf)
	implicit none
	integer, intent(in):: ik, iscf
	double precision, dimension(3), intent(in) :: kvec
	double complex, dimension(ntot,ntot), intent(out) :: hk
	double complex :: eikr, eikr2
	double precision :: kr
	integer :: i1,i2, j1,j2,is,js,ia,ja
	integer :: il,io,ii,i,k,k2
	integer :: i3,i4,j3,j4
	! local
	double precision :: t0
	
	! hksave: to keep a copy of hk without vmat part; to build hk during scf loop
	! vmatold: to keep a copy of vmat of previous iteration
	!        for mixing during scf loop

	!===================================================================
	! if iscf > 1:
	! H(k) calculated from H(k)_{save} 
	! Only vmat part updated by mixing old with new.
	!===================================================================
	if(iscf > 1) then ! use hksave and vmat to construct hk, & return
	 hk = hksave(:,:,ik); ! important to reset, because hk has eigenvectors on return from diag()
	if(lhu) then
	 t0 = 1.0d0 - beta
	 do il=1,nlayers
	  do io=1, noctl
	   ia = tm(il,io)%ia;
	   !is = atom2species(ia);
	   i1 = atom2orb(1,ia); !i2 = atom2orb(2,ia); ! orbital ranges
	   i4 = atom2orb(4,ia);
	   ! Hubbard e-e interaction matrix
	   !hk(i1:i4,i1:i4) = tm(il,io)%vmat(:,:);
	   if(1==0 .and. il==1 .and. io==2 .and. ik==1) then
	    write(*,*) "io=2: i1,i4 = ",i1,i4
	    	write(*,*) "tm(il,io)%vmat: "
	    do i2=1,10
	    write(*,'(10000e10.2)') tm(il,io)%vmat(i2,:)
	    end do
	   endif
	   ! simple linear mixing, only vmat part in Hk is updated
	   if(iscf==2) then ! we only have vmat
	    hk(i1:i4,i1:i4) = hk(i1:i4,i1:i4) + tm(il,io)%vmat
	   else ! we have vmatold to mix!
	    hk(i1:i4,i1:i4) = hk(i1:i4,i1:i4) + beta*tm(il,io)%vmat
     .                                 + t0*tm(il,io)%vmatold
	   endif
	   ! update vmatold
	   tm(il,io)%vmatold = tm(il,io)%vmat !


	   if(1==0 .and. il==1 .and. io==2 .and. ik==1) then
	    write(*,*) "hk: "
	    do i2=1,10
	    write(*,'(10000e10.2)') hk(i1-1+i2,i1:i4)
	    end do
	   endif
	   
	  end do ! io
	 end do ! il
	!elseif ! no Hubbard U, just copy Hk from Hsave
	 !hk = hsave ; already reset above outside lhu if block
	 !write(*,'(a)')"Warning: H(k) is nothing to update."
	 !write(*,'(a)')"                Why trying scf .... ?"
	endif ! lhu
	
	 return
	endif ! iscf>1


	!write(*,*) '**************************************** iscf=',iscf
	!===================================================================
	! H(k) calculated (element by element) only if iscf==1 or 0
	! also sets H(k)_{save}
	!===================================================================
	hk(:,:) = 0.0d0

	! onsite matrix elements
	! can be done once and reused for every k point
	! but, will cost memory and adding full space hii to hij 
	! might cost more than setting only relevant matrix elements as done below.
	do il=1,nlayers
	 do io=1, noctl
	  ia = tm(il,io)%ia;
	  i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
	  is = atom2species(ia);
	  do i=i1,i2
	   hk(i,i) = hk(i,i) + dcmplx(onsite(is),0.0d0)
	  end do
	  if(nspin==2) then ! down spin block
	   do i=i1+norbtm,i2+norbtm
	    hk(i,i) = hk(i,i) + dcmplx(onsite(is),0.0d0)
	   end do
	  endif	  

		do ii=1,3
	   ia = ox(il,io,ii)%ia;
	   i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
	   do i=i1,i2
	    hk(i,i) = hk(i,i) + dcmplx(onsite(0),0.0d0)
	   end do
	   if(nspin==2) then ! down spin block
	    do i=i1+norbo,i2+norbo
	     hk(i,i) = hk(i,i) + dcmplx(onsite(0),0.0d0)
	    end do
	   endif
		end do
	 end do ! io
	end do ! il

	!write(*,*)'onsite(0) = ',onsite(0)

	!...................................................................
	! TM-O  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl

	  !is = tm(il,io)%is; ! species index of TM atom
	  ! atom index, and orbital range
	  ia = tm(il,io)%ia;
	  i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
	  !i1 = tm(il,io)%i; i2 = tm(il,io)%j;
	  !write(*,'(a,5i5)') 'il,io, ia , i1,i2= ',il,io, ia,i1,i2

	  do k = 1,6 ! 1st nns O
	   !write(*,*) 'k, tm(il,io)%nn1(k)%ia = ',k,tm(il,io)%nn1(k)%ia
	   ja = tm(il,io)%nn1(k)%ia;	   
	   j1 = atom2orb(1,ja); j2 = atom2orb(2,ja); ! orbital ranges

	   !write(*,*) 'ja, j1,j2 = ',ja,j1,j2 

	   kr = DOT_PRODUCT(kvec,tm(il,io)%nn1(k)%dr)  
	   eikr = dcmplx(dcos(kr),dsin(kr))

	   hk(i1:i2,j1:j2) =  hk(i1:i2,j1:j2) + tm(il,io)%nn1(k)%h * eikr
!	  if(ia ==1 .and. k==3 ) then
!	  	 k2=6
!	   kr = DOT_PRODUCT(kvec,tm(il,io)%nn1(k2)%dr)  
!	   eikr2 = dcmplx(dcos(kr),dsin(kr))
!	   write(*,'(a,1000G10.3)')' e* : ',aimag(eikr),aimag(eikr2)
!	   write(*,'(a,1000G10.3)')'k=3:h ',tm(il,io)%nn1(k)%h
!	   write(*,'(a,1000G10.3)')'k=6:h ', tm(il,io)%nn1(k2)%h
!	   write(*,'(a,1000G10.3)')'k=3:h*e ',tm(il,io)%nn1(k)%h*eikr
!	   write(*,'(a,1000G10.3)')'k=6:h*e ', tm(il,io)%nn1(k2)%h*eikr2
!	  endif




!	if(1==0) then
!	do k2=k,6
!	   kr = DOT_PRODUCT(kvec,tm(il,io)%nn1(k2)%dr)  
!	   eikr2 = dcmplx(dcos(kr),dsin(kr))
!	 if(k2 .ne. k .and. tm(il,io)%nn1(k)%ia==tm(il,io)%nn1(k2)%ia !.and.
!!     . norm2(abs(tm(il,io)%nn1(k)%h*eikr+tm(il,io)%nn1(k2)%h*eikr2))
!!     .     < 1.0d-10 
!     .    ) then
!	    write(*,'(a,10i5)') 'k,k2, il,io, ia, ja/ja2 =',
!     . k,k2,il,io, ia,ja !,tm(il,io)%nn1(k2)%ia
!	 endif
!	end do
!	endif




!	   if (norm2(abs(hk(i1:i2,j1:j2))) < 1.0d-10 .and. 1==0)then
!	   	write(*,'(a,5i5)') 'il,io, ia , i1,i2= ',il,io, ia,i1,i2
!	    write(*,'(a,5i5)')'k, ja, j1,j2 = ',k,ja,j1,j2
!	    !write(*,'(1000f6.2)')(norm2(tm(il,io)%nn1(k)%h(:,js)),js=1,3)
!	    !write(*,'(1000f6.2)')(norm2(abs(hk(i1:i2,j1-1+js))),js=1,3)
!	    !write(*,'(a,1000f6.2)')':h=', tm(il,io)%nn1(k)%h*eikr
!	   endif

	   !if (tm(il,io)%nn1(k)%ia==4)then
	   !write(*,'(a,1000f6.2)')'eikr =', real(eikr), aimag(eikr)
	   !endif


!	   if (tm(il,io)%nn1(k)%ia==4 .and. 1==0)then
!	    write(*,'(a,3f10.3)') 'kvec = ', kvec
!	    write(*,'(a,3f10.3)') 'Oz with B1: dr=', tm(il,io)%nn1(k)%dr
!	   	!write(*,'(a,5i5)') 'il,io, ia , i1,i2= ',il,io, ia,i1,i2
!	    !write(*,'(a,5i5)')'k, ja, j1,j2 = ',k,ja,j1,j2
!	    !write(*,'(1000f6.2)')(norm2(tm(il,io)%nn1(k)%h(:,js)),js=1,3)
!	    !write(*,'(1000f6.2)')(norm2(abs(hk(i1:i2,j1-1+js))),js=1,3)
!	    write(*,'(a,1000f6.2)')'h = ', tm(il,io)%nn1(k)%h
!	    write(*,'(a,1000f6.2)')'eikr =', real(eikr), aimag(eikr)
!	    write(*,'(a,1000f6.2)')'hk = ', tm(il,io)%nn1(k)%h*eikr
!	   endif

!	   if (tm(il,io)%nn1(k)%ia==8 .and. 1==0)then
!	    write(*,'(a,3f10.3)') 'Oz with B2: dr=', tm(il,io)%nn1(k)%dr
!	   	!write(*,'(a,5i5)') 'il,io, ia , i1,i2= ',il,io, ia,i1,i2
!	    !write(*,'(a,5i5)')'k, ja, j1,j2 = ',k,ja,j1,j2
!	    !write(*,'(1000f6.2)')(norm2(tm(il,io)%nn1(k)%h(:,js)),js=1,3)
!	    !write(*,'(1000f6.2)')(norm2(abs(hk(i1:i2,j1-1+js))),js=1,3)
!	    write(*,'(a,1000f6.2)')':h=', tm(il,io)%nn1(k)%h*eikr
!	   endif

	  end do ! k
	 end do ! io
	end do ! il

	!write(*,*)'-------------------- a: tmnn2=',tmnn2

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
	   !write(*,*) 'kr = ',kr
	   eikr = dcmplx(dcos(kr),dsin(kr))
	   hk(i1:i2,j1:j2) = hk(i1:i2,j1:j2) + tm(il,io)%nn2(k)%h *eikr
	  end do ! k
	 end do ! io
	end do ! il
	endif
	!...................................................................

	!write(*,*)'-------------------- b'

	!write(*,*)'===> O-TM disabled ....'
	!if(1==0) then
	!...................................................................
	! O-TM  (1st neighbours)
	!...................................................................
	!if(1==0) then
	!write(*,*)'O-TM disabled.... '
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
	!endif
	!...................................................................
	!write(*,*)'-------------------- c'
	!endif
	
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

	!write(*,*)'-------------------- d'


!	if(norm2(kvec) < 1.0d-6 .and. 1==0) then
!	do il=1,nlayers
!	 do io=1, noctl
!		do ii=1,3
!	   ia = ox(il,io,ii)%ia;
!	   i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges
!	    write(*,*)'O space: H(k=0): ia = ',ia 
!	    write(*,'(10000f10.5)') (norm2(abs(hk(i,1:ntottm))),i=i1,i2)
!	  write(*,'(10000f10.5)')(norm2(abs(hk(i,ntottm+1:ntot))),i=i1,i2)
!		end do
!	 end do ! io
!	end do ! il
!	endif

	!write(*,*) 'nspin = ',nspin
	if(nspin==2) then
	 ! make a copy of above Hk at the spin-down indices...
	 do ia=1,natoms
	  i1 = atom2orb(1,ia); ! orbital ranges, spin up
	  i2 = atom2orb(2,ia); 
		i3 = atom2orb(3,ia); ! spin down
		i4 = atom2orb(4,ia); 
	  do ja=1,natoms
	   j1 = atom2orb(1,ja); ! orbital ranges, spin up
	   j2 = atom2orb(2,ja); 
		 j3 = atom2orb(3,ja); ! spin down
		 j4 = atom2orb(4,ja); 
		 ! copy spin up blocks to spin down blocks
	   hk(i3:i4,j3:j4) = hk(i1:i2,j1:j2)
	   ! add local atomic spin orbit coupling blocks of TM atoms
	   if(ia==ja .and. lsoc) then
	    is = atom2species(ia);
	    if (is > 0) then ! TM atoms, for O atom2species gives -1.
	     ! up-up block
	     hk(i1:i2,j1:j2) = hk(i1:i2,j1:j2) + 
     .          soc(is)*Hsoc(1:norbtm,1:norbtm)
	     ! dn-dn block
	     hk(i3:i4,j3:j4) = hk(i3:i4,j3:j4) + 
     .          soc(is)*Hsoc(6:5+norbtm,6:5+norbtm)
	     ! up-dn block; this block is 0 at this stage
	     hk(i1:i2,j3:j4) = soc(is)*Hsoc(1:norbtm,6:5+norbtm)
	     ! dn-up block; this block is 0 at this stage; 
	     ! ? not needed because we use Upper Triangular in diag routines
	     hk(i3:i4,j1:j2) = soc(is)*Hsoc(6:5+norbtm,1:norbtm)
	    endif ! is > 0
	   endif !ia==ja
	  end do ! ja
	 end do ! ia
	endif ! nspin==2

	!===================================================================
  ! save hk for scf loop, hamiltonian mixing will only update vmat part.... 
	if(iscf ==1) then ! setting isscf=0 for band calculations avoids saving these arrays
	 hksave(:,:,ik) = hk;
	endif
	!===================================================================


	! band structure: assumes SCF has been performaed...
	! or maybe, vmat read from some file.... written after scf in a previous run...
	
	if(iscf==0 .and. lhu) then
	 do il=1,nlayers
	  do io=1, noctl
	   ia = tm(il,io)%ia;
	   i1 = atom2orb(1,ia);
	   i4 = atom2orb(4,ia);
	   ! Converged Hubbard e-e interaction matrix
	   hk(i1:i4,i1:i4) = hk(i1:i4,i1:i4) + tm(il,io)%vmat(:,:);	   
	  end do ! io
	 end do ! il
	endif

	return
	end 	subroutine getHk
!----------------------------------------------------------------


	end module hamiltonian
	
