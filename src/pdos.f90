
module pdos
use modmain

implicit none




contains

!===========================================================================
subroutine getpdos(ef) !mine, maxe, nwplot)
implicit none
!integer, intent(in) :: nwplot
!double precision, intent(in) :: mine, maxe
real(8), intent(in) :: ef
real(4), allocatable :: bc(:,:,:,:,:)
real(4), allocatable :: bcj(:,:,:,:),bcl(:,:,:,:)

real(8), allocatable :: w(:), es(:,:), dos(:)
real(8), allocatable :: gc(:,:,:), gcj(:,:),gcl(:,:)

integer :: natomtm, iw, itm
double precision :: nbymaxe, dw
integer :: ik, ist, il, io, ia,i1,i4
double complex, dimension(10) :: wf,wf2

natomtm = natoms/4;

allocate(es(ntotk,ntot)) 
allocate(w(nwplot))
allocate(dos(nwplot))
allocate(gc(nwplot,norbtm,nspin))
allocate(gcj(nwplot,norbtm*nspin))
allocate(gcl(nwplot,norbtm*nspin))
allocate(bc(norbtm,nspin,natomtm,ntot,ntotk))
allocate(bcj(norbtm*nspin,natomtm,ntot,ntotk))
allocate(bcl(norbtm*nspin,natomtm,ntot,ntotk))

! to get index easily of a state in the energy bins array
!mine = minval(eval);
!maxe = maxval(eval);

nbymaxe = (nwplot-1)/(maxe-mine);
es = 1.0d0 +  nbymaxe*(eval-ef - mine); ! lies in range 1 to ne.

dw = (maxe - mine)/dble(nwplot-1);
do iw=1,nwplot
 w(iw) = mine + (iw-1)*dw
end do

!--------------------------------------------------------------
! total dos
dos = 0.0d0
do ist=1,ntot
 do ik=1,ntotk
   iw = int(es(ik,ist))
   if(iw > 0 .and. iw <= nwplot) then
    dos(iw) = dos(iw) + wk(ik) ! total dos
   end if
 end do
end do
open(52,file='DOS.OUT',form='FORMATTED',position='append')
do iw=1,nwplot
 write(52,'(20G18.10)') w(iw), dos(iw)
end do
write(52,'("     ")')
write(52,'("     ")') ! two empty lines for gnu indexing
close(52)
!--------------------------------------------------------------

if(lpdos) then

 ! band character for real Ylm/spin and for J,Jz basis
 itm=0;
 do il=1,nlayers
 do io=1,noctl
  itm = itm + 1;
  ia= tm(il,io)%ia;
  i1=atom2orb(1,ia)
  i4=atom2orb(4,ia)
  do ist=1,ntot
   do ik=1,ntotk
    wf = evec(ik,i1:i4,ist)
    ! band character in real Ylm,spin basis (our basis)
	  bc(:,1,itm,ist,ik) = zabs(wf(1:5))**2  ! up
	  bc(:,2,itm,ist,ik) = zabs(wf(6:10))**2 ! down
    ! change the basis to J,Jz
	  wf2 = matmul(Ulm2j,wf) 
	  bcj(:,itm,ist,ik) = zabs(wf2)**2; ! in J,Jz basis 
	  ! first 6: J=5/2, Jz=-5/2,-3/2,-1/2,1/2,3/2,5/2
	  ! later 4: J=3/2, Jz=-3/2,-1/2,1/2,3/2

		! local Ham eigenstates: local ham given by getHtmsUJ(), eig calc in scf.f90
		wf2 = matmul(tm(il,io)%UT, wf)
	  bcl(:,itm,ist,ik) = zabs(wf2)**2; ! in local eigenstate bases


   end do
  end do
 end do
 end do

 ! write bc, bcj to a file??? use the above code to compute bc for bands, make a routine or just copy the code 
 !---------------------------------------------------------------
 ! calc and write pdos 
 dos = 0.0d0; ! for all Oxygen atoms
 open(50,file='PDOS.OUT',form='FORMATTED',position='append')
 open(51,file='PDOSJ.OUT',form='FORMATTED',position='append')
 open(52,file='PDOSL.OUT',form='FORMATTED',position='append')

 do itm=1,natomtm
 gc = 0.0d0
 gcj = 0.0d0;
 gcl = 0.0d0;
 do ist=1,ntot
  do ik=1,ntotk
   iw = int(es(ik,ist))
   if(iw > 0 .and. iw <= nwplot) then
    gc(iw,:,1) = gc(iw,:,1) + bc(:,1,itm,ist,ik)*wk(ik)
    gc(iw,:,2) = gc(iw,:,2) + bc(:,2,itm,ist,ik)*wk(ik)
    gcj(iw,:) = gcj(iw,:) + bcj(:,itm,ist,ik)*wk(ik)
    gcl(iw,:) = gcl(iw,:) + bcl(:,itm,ist,ik)*wk(ik)

   endif
  end do
 end do

 ! write output file
 do iw=1,nwplot
  dos(iw) = (1.0d0 - sum(gcj(iw,:)))/dble(3*natomtm)  ! Oxygen atoms average 
  write(50,'(21G18.10)') w(iw), dos(iw), gc(iw,:,1), gc(iw,:,2)
  write(51,'(21G18.10)') w(iw), dos(iw), gcj(iw,:)
  write(52,'(21G18.10)') w(iw), dos(iw), gcl(iw,:)

 end do
 write(50,'("     ")')
 write(50,'("     ")') ! two empty lines for gnu indexing
 write(51,'("     ")')
 write(51,'("     ")') ! two empty lines for gnu indexing
 write(52,'("     ")')
 write(52,'("     ")') ! two empty lines for gnu indexing

 end do ! atoms
 close(50)
 close(51)
 close(52)


endif ! lpdos
 deallocate(bc,bcj,es,w,gc,gcj,dos,gcl,bcl)

return
end subroutine getpdos
!===========================================================================


!===========================================================================
subroutine getbc()
implicit none
real(4), allocatable :: bc(:,:,:,:,:)
real(4), allocatable :: bcj(:,:,:,:), bcl(:,:,:,:)
integer :: natomtm, itm, ib
integer :: ik, ist, il, io, ia,i1,i4, ispin
double complex, dimension(10) :: wf, wf2

natomtm = natoms/4;

allocate(bc(norbtm,nspin,natomtm,ntot,np))
allocate(bcj(norbtm*nspin,natomtm,ntot,np))
allocate(bcl(norbtm*nspin,natomtm,ntot,np))

 ! band character for real Ylm/spin and for J,Jz basis
 itm=0;
 do il=1,nlayers
 do io=1,noctl
  itm = itm + 1;
  ia= tm(il,io)%ia;
  i1=atom2orb(1,ia)
  i4=atom2orb(4,ia)
  do ist=1,ntot
   do ik=1,np
    wf = evec(ik,i1:i4,ist)
    ! band character in real Ylm,spin basis (our basis)
	  bc(:,1,itm,ist,ik) = zabs(wf(1:5))**2  ! up
	  bc(:,2,itm,ist,ik) = zabs(wf(6:10))**2 ! down
    ! change the basis to J,Jz
	  wf2 = matmul(Ulm2j,wf) 
	  bcj(:,itm,ist,ik) = zabs(wf2)**2; ! in J,Jz basis 
	  ! first 6: J=5/2, Jz=-5/2,-3/2,-1/2,1/2,3/2,5/2
	  ! later 4: J=3/2, Jz=-3/2,-1/2,1/2,3/2

		! local Ham eigenstates: local ham given by getHtmsUJ(), eig calc in scf.f90
		wf2 = matmul(tm(il,io)%UT, wf)
	  bcl(:,itm,ist,ik) = zabs(wf2)**2; ! in local eigenstate bases

   end do
  end do
 end do
 end do

	! format: dp, e, bc_{all O orbitals}, bc_{TM 1 up}, bc_{TM 1 dn}, bc_{TM 2 up},...
	!-------------------------------------------
	open(10,file='BANDC.OUT',form='FORMATTED',action='write')
	do ib=1,ntot
	 do ik=1,np
	  write(10,'(2G20.8,3x,100000e10.2)') dp(ik), eval(ik,ib)-efermi, &
	    (1.0d0 - sum(bc(:,:,:,ib,ik)))/dble(3*natomtm) , &
      ((bc(:,ispin,itm,ib,ik),ispin=1,nspin), itm=1,natomtm)	  
	 end do
	 write(10,'(2f15.8)')
	 write(10,'(2f15.8)') 
	end do
	close(10)
	!-------------------------------------------
	! format: dp, e, bc_{all O orbitals}, bc_{TM 1 up}, bc_{TM 1 dn}, bc_{TM 2 up},...
	open(10,file='BANDCJ.OUT',form='FORMATTED',action='write')
	do ib=1,ntot
	 do ik=1,np
	  write(10,'(2G20.8,3x,100000e10.3)') dp(ik), eval(ik,ib)-efermi, &
      (1.0d0 - sum(bc(:,:,:,ib,ik)))/dble(3*natomtm), &
      (bcj(:,itm,ib,ik), itm=1,natomtm) 
	 end do
	 write(10,'(2f15.8)')
	 write(10,'(2f15.8)') 
	end do
	close(10)
	!-------------------------------------------

	open(10,file='BANDCL.OUT',form='FORMATTED',action='write')
	do ib=1,ntot
	 do ik=1,np
	  write(10,'(2G20.8,3x,100000e10.3)') dp(ik), eval(ik,ib)-efermi, &
      (1.0d0 - sum(bc(:,:,:,ib,ik)))/dble(3*natomtm), &
      (bcl(:,itm,ib,ik), itm=1,natomtm) 
	 end do
	 write(10,'(2f15.8)')
	 write(10,'(2f15.8)') 
	end do
	close(10)
	!-------------------------------------------
 deallocate(bc,bcj)
 deallocate(bcl)

return
end subroutine getbc
!===========================================================================


end module pdos
