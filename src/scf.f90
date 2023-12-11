
module scf
use modmain
use hamiltonian
use fermi
use hubbard
use mixing
use fixmom
use mtbmpol
use estatic
use mtbeseld

implicit none



contains

!=====================================================================
subroutine tmgroundstate()
implicit none
integer :: il,io,i, ia,is
double precision, dimension(10) :: eval
double precision :: engyes
double complex, dimension(10,10) :: ham


 write(6,'("Calculating TM atoms GS .... ")')

	! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
 call tbeseld(lesH,engyes)

open(804,file='evals-tm.dat',action='write',position='append')
open(805,file='evecs-tm.dat',action='write',position='append')

 !-------- -------- -------- -------- -------- -------------- 
 ! Hamiltonian, eigenstates and eigenvalues: on kgrid in BZ
 !-------- -------- -------- -------- -------- -------------- 

 call getHtms()
 !write(*,*) '100 * h_B1: --------------------'
 !write(*, '(20f15.6)') 100*tm(1,1)%ham(:,:) ! evec

 call writeHamOctaFrame()

 do il=1,1 !nlayers
  do io=1,1 ! noctl
   ia = tm(il,io)%ia;
	 is = atom2species(ia);

	 ham = tm(il,io)%ham
	 !write(*,*) '------------'
	 !write(*,'(10f8.3)') dble(ham)

   call zdiag(10,ham,eval,0,10)

	 write(804,'(10f18.12)') eval* 27.21138505d0 ! *1/eV2har
	 write(805,*) ham ! evec

  end do ! io
 end do !il
	
 close(804)
 close(805)

return
end subroutine tmgroundstate
!=====================================================================





!====================================================================
	subroutine writeHamOctaFrame()
	implicit none
	integer ::i,ia,il,io,is
	character*50 :: fname ='beast-tm-ham-octa-frame.dat'
	character*50 :: fname2 ='beast-tm-ham.dat'

	double precision, dimension(5,5) :: h, R, RT

	open(10,file=trim(fname),action='write',position='append')
	open(11,file=trim(fname2),action='write',position='append')

	do il=1, 2
	do io=1, 2
   ia = tm(il,io)%ia;
	 is = atom2species(ia);

	 h = tm(il,io)%ham(1:5,1:5) ! just up spin block
	 write(11,'(5f18.12)') h * 27.21138505d0 ! write 5x5 matrices as an array of 25 elements

	 i = (il-1)*2 + io
	 R = Rlmax(5:9,5:9,i) ! Rlmax is rotation matrix for Ylm of TM atom i
	 RT = transpose(R);
	  ! transform to rotated frame of the octagedron i 
	  ! R^T.h.R
	  h = matmul(RT,h)
	  h = matmul(h,R) ! 
	
	 write(10,'(5f18.12)') h * 27.21138505d0 ! write 5x5 matrices as an array of 25 elements

  end do !io
  end do ! il

	close(10)
	close(11)	
	return
	end 	subroutine WriteHamOctaFrame
!----------------------------------------------------------------






!=====================================================================
subroutine tmscfground()
implicit none
double precision, dimension(10) :: eval, r
double precision :: ddm, engyadu, ebands, engyes, energyb, edu
double complex, dimension(10,10) :: ham
double precision :: qtm
integer :: il,io,i, ia,is, iscf, k, ib, i1,j1,j, m1,m2

double precision, dimension(1,10) ::wke, eval1
double complex, dimension(2,2) :: x






 write(6,'("Calculating TM atoms SCF GS .... ")')

	! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
 call tbeseld(lesH,engyes)


! initialise the mixer
 call mixerifc(0, ddm) ! allocate nu, mu, f, beta arrays

do iscf = 1, maxscf
  if(mod(iscf,100)==0) write(6,'(a,i5)') "  SCF iteration ", iscf
 !--------------------------------------------------------

 call getHtmsUJ(iscf)

 do il=1,nlayers
  do io=1, noctl
   ia = tm(il,io)%ia;
	 is = atom2species(ia);
	 ham = tm(il,io)%ham
   call zdiag(10,ham,eval,0,10)



 !-------- -------- -------- -------- -------- -------------- 
 ! find the Fermi level and wke = occupation weighted by wk 
 !-------- -------- -------- -------- -------- -------------- 
 eval1(1,:) = eval ! do we need this new array?
 qtm = q0(is) - (6-qa)
 call fermid(1, (/1.0d0/), 1.0d0, 10, eval1, temp, &
                      qtm, wke, efermi, nspin, ebands)
 !write(*,'(a,f15.6)') "N_electron = ", qtot
 !write(*,'(a, 3f15.10)') "Fermi energy = ", efermi

	! Hubbard U & Hund's J: tm%vmat
	!write(*,*) 'scf: qa, qaa, int(qa), ne = ',qa, qaa, int(qa), ne
	call mkvmattm(iscf, il, io, is, ham, wke(1,:), edu)


  ! calculate rhoc etc to be used for Qmpol: 
  ! [call getatomic() in full system case]
  atm(ia)%qs = 0.0d0
  do i=1,5
   do k=1,nspin
   atm(ia)%qs(k) = atm(ia)%qs(k) + dble(tm(il,io)%dm(i,k,i,k))
   end do
  end do
!  write(*,*)'scf: ia, atm(ia)%qs = ',atm(ia)%qs 
	! set rhoc = dm; reshape first.
	call mkdm2rhoc(tm(il,io)%dm, atm(ia)%rhoc) 

  end do ! io
 end do !il


	! Oxygen: set atm%rhoc and atm%qs
	! assume point charges, all 6 orbitals filled
	do ib=1, nbas
	 if(mod(ib-1,4)==0) cycle ! it's a tm atoms
	 atm(ib)%rhoc = 0.0d0
	 do i=2,4 ! ilm range for l=1
	  atm(ib)%rhoc(i,i) = 1.0d0
	 end do
	 atm(ib)%qs(:) = 3.0d0 ! 3 up and 3 down spin electrons
	end do ! ib

!  write(6,'(a,i5)') "1  ============ iscf=", iscf

  ! calculate Qmpol
  call tbmpol()
  ! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
	call tbeseld(lesH, engyes)
  ! now we can add atm%dh to the hamiltonian....
!  write(6,'(a,i5)') "2  ============ iscf=", iscf

 !-------- -------- -------- -------- -------- -------------- 	
  ! Hubbard U potential matrices for TM atoms
  ! already done in il,io loops
  !-------- -------- -------- -------- -------- -------------- 	
  ! mix vmat & Vm (multipoles)
  !-------- -------- -------- -------- -------- -------------- 	
  call mixpack(iscf, .true.) ! vectorise tm%vmat for mixing
  call mixerifc(iscf, ddm) ! does mixing
  call mixpack(iscf, .false.) ! tm%vmat from vectorised form.
  !-------- -------- -------- -------- -------- -------------- 
  ! check convergence of scf:
  !-------- -------- -------- -------- -------- -------------- 
  if(ddm < toldm) then ! iscf .ne. 1 .and. 
	 write(6,'("SCF coverged in ",i5," iterations!")') iscf
	 write(6,'("tolerance, change in Vmat & dH = ", 2e20.6)') &
	              toldm, ddm
	 exit ! exit scf loop
	elseif(mod(iscf,100)==0) then
   write(6,'("Absolute change in Vmat & dH = ", 2e20.6)') ddm
  endif
 !-------- -------- -------- -------- -------- -------------- 
 ! calculate magnetic fields for next iteration:
 call fsmbfield(iscf, energyb)

!write(6,'("Ebs, Euj, Eq, Eb = ", 4e20.6)') ebands, engyadu, engyes, energyb
!write(6,'("Etot = ", 1e20.6)') ebands + engyadu + engyes + energyb


end do ! iscf







 do il=1,nlayers
 do io=1,2

	write(*,'(a)') 'tm%vmat: real n imag parts:'
	write(*,'(10f7.4)') dble(tm(il,io)%vmat)
	write(*,'(a)')'---------------------------'	
	write(*,'(10f7.4)') dimag(tm(il,io)%vmat)

 
  ia = tm(il,io)%ia;
  atm(ia)%qs = 0.0d0
  x = 0.0d0
  do i=1,2
   do j=1,2
    do i1=1,5
    x(i,j) = x(i,j) + tm(il,io)%dm(i1,i,i1,j) ! diagonal in orbital index
    end do
   end do
   atm(ia)%qs(i) = atm(ia)%qs(i) + dble(x(i,i))
  end do
  write(*,'(a)') 'il,io, dm_spin:'
  write(*,'(2f8.5,10x,2f8.5)') x
  write(*,'(a,f10.5)') '======>   m_z: ', atm(ia)%qs(1)-atm(ia)%qs(2)
 end do
 end do


! do m1=1,5
 !do m2=1,5
!  write(*,'(5f7.4)') (Hub(1)%Vee(m1,m2,m1,m2), m2=1,5)
 !end do
! end do




!  write(*,*)'scf: U: ',Hub(1)%U
!  write(*,*)'scf: J: ',Hub(1)%J


! diagonalise the last tm%ham to get eigenvalues and eigenvectors [Self-consistent]
 open(804,file='evals-tm-scf.dat',action='write',position='append')
 open(805,file='evecs-tm-scf.dat',action='write',position='append')

 do il=1,nlayers
  do io=1, noctl
	 ham = tm(il,io)%ham
   call zdiag(10,ham,eval,0,10)
	 write(804,'(10f18.12)') eval
	 write(805,*) ham ! evec
  end do ! io
 end do !il
 	
 close(804)
 close(805)

return
end subroutine tmscfground
!=====================================================================

subroutine mkdm2rhoc(v,u)
implicit none
double complex, dimension(5,2,5,2), intent(in) :: v
double complex, dimension(5,5), intent(out) :: u
integer :: i,i1,ispin,j,jspin,j1
	! unfold vmat to directly add in hk
	u = 0.0d0
	do i=1,norbtm
	do j=1,norbtm
	 do ispin=1,nspin ! diagonal spin
		u(j,i)  = u(j,i) + v(j,ispin,i,ispin)
   end do
  end do
  end do

return
end subroutine mkdm2rhoc

!=====================================================================
subroutine nscfgroundstate()
implicit none
integer :: ik
double precision, dimension(3) :: kvec
double precision :: ddm, engyadu, ebands, engyes, energyb
integer :: il,io,i


!-------- BZ integration -------- -------- -------- -------- 
 if (.not. allocated(kgrid)) call mkkgrid(nk1,nk3) ! makes kgrid & wk for BZ integration, sets ntotk
 if (.not. allocated(hk)) allocate(hk(ntot,ntot))
 ! if larger sys or too large kgrid, then we can save mem by using sparse format for hksave, (and even for hk...)
 if (.not. allocated(hksave)) allocate(hksave(ntot,ntot,ntotk)) ! for mixing.... in hamiltonain module
 if (.not. allocated(eval)) allocate(eval(ntotk,ntot))
 if (.not. allocated(evec)) allocate(evec(ntotk,ntot,ntot))
 if (.not. allocated(wke)) allocate(wke(ntotk,ntot))


 write(6,'("Calculating non-selfconsistent GS .... ")')

! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
 call tbeseld(lesH,engyes)



do i=1,2 !natoms,4
	write(*,*)"scf.f90: BEFORE SCF STARTS: MONOPOLES: atom = ",i 
	write(*,'(5f10.5)') atm(i)%dh
enddo

!open(804,file='evals.dat',action='write',position='append')
 !-------- -------- -------- -------- -------- -------------- 
 ! Hamiltonian, eigenstates and eigenvalues: on kgrid in BZ
 !-------- -------- -------- -------- -------- -------------- 
 do ik= 1,ntotk
  kvec = kgrid(1:3,ik) ! already in cartesian coordinates
  ! to disable Hubbard term in getHk() but calc all other terms
  call getHk(ik,kvec, hk, -1) ! iscf=-1

	if(ik==1) then
	 open(123,file='ham.OUT',action='write')
	 write(123,'(112f6.2/)') hk
	 close(123)
	endif
  
  call zdiag(ntot,hk,eval(ik,:),ik,ntot)
  evec(ik,:,:) = Hk ! eigenvector are columns of Hk
  !call useeigdata(ntot,hk) ! calc whatever we like at this k-point
  ! will take a weighted average after the k-loop completes.

!	write(804,*) eval(ik,:)
  
 end do ! ik
! close(804)
 
 !-------- -------- -------- -------- -------- -------------- 
 ! find the Fermi level and wke = occupation weighted by wk 
 !-------- -------- -------- -------- -------- -------------- 
 call fermid(ntotk, wk, wknorm, ntot, eval, temp, &
                      qtot, wke, efermi, nspin, ebands)
 !write(*,'(a,f15.6)') "N_electron = ", qtot
 write(*,'(a, 3f15.10)') "Fermi energy = ", efermi

	open(456,file='ebands.dat',action='write',position='append')
	write(456,*) ebands
	close(456)

 !--------------------------------------------------------
  ! calculate rhoc etc to be used for Qmpol
  call getatomic()
  ! calculate Qmpol
  call tbmpol()

  ! TM%dm; iscf=1 to trick to allocate and set tm%dm 
  call mkvmat(1, engyadu) ! uses global evec


	open(114,file='dh_B1.dat', action='write',position='append')
	write(114,'(5f15.10)') atm(1)%dh
	close(114)


	open(114,file='mom-B1.dat', action='write',position='append')
	write(114,'(5f18.10)') tm(1,1)%mag
	close(114)

	open(114,file='pop-B1.dat', action='write',position='append')
	write(114,'(10f8.4)') (dble(tm(1,1)%dm(i,1,i,1)),i=1,5), (dble(tm(1,1)%dm(i,2,i,2)),i=1,5)
	write(114,'(10f8.4)') (dble(tm(1,2)%dm(i,1,i,1)),i=1,5), (dble(tm(1,2)%dm(i,2,i,2)),i=1,5)
	write(114,'(10f8.4)') 
	close(114)


do i=1,natoms
	if( mod(i-1,4) == 0 ) then
	 write(*,*)"scf.f90: TM atom: ",i 
	 write(*,'(a,i5,25f7.3)')"scf.f90: Qm: atom, i = ",i,qmpol(:,i)
	 !write(*,'(5f10.5)') atm(i)%dh
	else
		write(*,*)"O atom: ",i 
		!write(*,'(3f10.5)') atm(i)%dh
		write(*,'(a,i5,25f7.3)')"scf.f90: Qm: atom, i = ",i,qmpol(:,i)
	endif
enddo

return
end subroutine nscfgroundstate
!=====================================================================




!=====================================================================
subroutine groundstate()
implicit none
integer :: iscf, ik
double precision, dimension(3) :: kvec
double precision :: ddmold, ddm, engyadu, ebands, engyes, energyb
integer :: il,io,i
double precision, dimension(10) :: evall
integer:: jspn, ispn
! some large number
ddmold = 1.0d8;

write(*,*)'===============================      ntot = ',ntot

!-------- BZ integration -------- -------- -------- -------- 
 if (.not. allocated(kgrid)) call mkkgrid(nk1,nk3) ! makes kgrid & wk for BZ integration, sets ntotk

 !write(*,*) 'scf: kgrid = ', kgrid

 if (.not. allocated(hk)) allocate(hk(ntot,ntot))
 ! if larger sys or too large kgrid, then we can save mem by using sparse format for hksave, (and even for hk...)
 if (.not. allocated(hksave)) allocate(hksave(ntot,ntot,ntotk)) ! for mixing.... in hamiltonain module
 if (.not. allocated(eval)) allocate(eval(ntotk,ntot))
 if (.not. allocated(evec)) allocate(evec(ntotk,ntot,ntot))
 if (.not. allocated(wke)) allocate(wke(ntotk,ntot))
! if(lhu)then
!	allocate(tm(il,io)%vmat(norbtms, nspin,norbtms, nspin))
!	allocate(tm(il,io)%dm(norbtms, nspin,norbtms, nspin))
! endif

if(.not.lhu) then
 maxscf =1
 write(6,'("HubbardU = F ==> only 1 SCF cycle")')
else
 write(6,'("Starting SCF loop .... ")')
! if(lusevmat) then
!  call readvmat()
! else
! ! given atomic densities/occupations:
! call setupdm()
! endif
endif

! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
 call tbeseld(lesH,engyes)



do i=1,2 !natoms,4
	write(*,*)"scf.f90: BEFORE SCF STARTS: MONOPOLES: atom = ",i 
	write(*,'(5f10.5)') atm(i)%dh
enddo





! initialise the mixer
 call mixerifc(0, ddm) ! allocate nu, mu, f, beta arrays

do iscf = 1, maxscf
 write(6,'(a,i5)') "===================>>> SCF iteration ", iscf


	!write(*,'(5f18.13)') (atm(1)%dh(i,i),i=5,9)






 !-------- -------- -------- -------- -------- -------------- 
 ! Hamiltonian, eigenstates and eigenvalues: on kgrid in BZ
 !-------- -------- -------- -------- -------- -------------- 
 do ik= 1,ntotk
  kvec = kgrid(1:3,ik) ! already in cartesian coordinates
  call getHk(ik,kvec, hk, iscf)

	if(ik==1) then
	 open(123,file='ham.OUT',action='write')
	 write(123,'(112f6.2/)') hk
	 close(123)
	endif
  
  call zdiag(ntot,hk,eval(ik,:),ik,ntot)
  evec(ik,:,:) = Hk ! eigenvector are columns of Hk
  !call useeigdata(ntot,hk) ! calc whatever we like at this k-point
  ! will take a weighted average after the k-loop completes.
 end do ! ik

 
 !-------- -------- -------- -------- -------- -------------- 
 ! find the Fermi level and wke = occupation weighted by wk 
 !-------- -------- -------- -------- -------- -------------- 
 call fermid(ntotk, wk, wknorm, ntot, eval, temp, &
                      qtot, wke, efermi, nspin, ebands)
 !write(*,'(a,f15.6)') "N_electron = ", qtot
 write(*,'(a, 3f15.10)') "Fermi energy = ", efermi

 !if(iscf==20) write(*,'(a, 100000f6.1)')'eval=',eval

 !-------- -------- -------- -------- -------- -------------- 	
 ! something observables to be calculated? weighted averages here.
 ! call averages(efermi)
 !-------- -------- -------- -------- -------- -------------- 	

 !--------------------------------------------------------
 if (lhu) then
  ! dueing each SCF iteration:
  ! calculate rhoc etc to be used for Qmpol
  call getatomic()
  ! calculate Qmpol
  call tbmpol()

  !write(*,*)'scf: qmpol_0 = ',qmpol(1,:)
  
  ! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
	call tbeseld(lesH, engyes)
  ! now we can add atm%dh to the hamiltonian....
 !-------- -------- -------- -------- -------- -------------- 	
  ! Hubbard U potential matrices for TM atoms
  call mkvmat(iscf, engyadu) ! uses global evec

	!write(*,'(a)')"scf: B1: vmat:"
	!write(*,'(10f8.4)') dble(tm(1,1)%vmat)


  !-------- -------- -------- -------- -------- -------------- 	
  ! mix vmat & Vm (multipoles)
  !-------- -------- -------- -------- -------- -------------- 	
  call mixpack(iscf, .true.) ! vectorise tm%vmat for mixing
  call mixerifc(iscf, ddm) ! does mixing
  call mixpack(iscf, .false.) ! tm%vmat from vectorised form.

	!write(6,*) 'iscf, ddm  = ', iscf, ddm
  !-------- -------- -------- -------- -------- -------------- 
  ! check convergence of scf:
  !-------- -------- -------- -------- -------- -------------- 
  if(ddm < toldm) then ! iscf .ne. 1 .and. 
	 write(6,'("SCF coverged in ",i5," iterations!")') iscf
	 write(6,'("tolerance, change in Vmat & dH = ", 2e20.6)') &
	              toldm, ddm
	 !write(6,'("Ebs + Euj = ", 2e20.6)') ebands, engyadu

	 !write(*,*)"CFM(2,2,4,:): -----------" 
	 !write(*,'(10f7.2)') CFM(2,2,4,:)

	write(*,*)"scf: hubbard.f, occ of individual orb n_m in FLL DCC"

	 !write(*,*)"scf.f90: real of tm%vmat"
	 !write(*,'(10f8.4)') dble(tm(1,1)%vmat) ! 
	 !write(*,*)"scf.f90: imag of tm%vmat"
	 !write(*,'(10f8.4)') dimag(tm(1,1)%vmat) ! 

	 !write(*,'(a)') ' - - - - - - - '
	 !write(*,'(10f10.5)') tm(2,1)%vmat
 
	 exit ! exit scf loop
	else
   write(6,'("Absolute change in Vmat & dH = ", 2e20.6)') ddm
  endif
  !-------- -------- -------- -------- -------- -------------- 	
  !ddmold = ddm
 endif

 
 !-------- -------- -------- -------- -------- -------------- 
 ! calculate magnetic fields for next iteration:
 call fsmbfield(iscf, energyb)

write(6,'("Ebs, Euj, Eq, Eb = ", 4e20.6)') ebands, engyadu, engyes, energyb
write(6,'("Etot = ", 1e20.6)') ebands + engyadu + engyes + energyb

end do! iscf



!-----------------------------------------------------------------
! calc local ham and its eigenvectors for pdoc abd band characters
 call getHtmsUJ(10) ! iscf>1 will include tm%vmat
 
 open(115,file='localeig.dat', action='write',position='append')
 do il=1,nlayers
  do io=1, noctl
   call zdiag(10,tm(il,io)%ham,evall,0,10)
   tm(il,io)%UT = transpose(dconjg(tm(il,io)%ham))

   ! write only 5 eval, up spin 
   write(115,'(a,10f15.10)') 'Eigenvalues, spin up: ',evall(1:10:2)
   ! in bases aligned with lattice/global axes
   ! write(115,'(20f15.10)') tm(il,io)%ham

   ! assuming spin degen, separate diff spin comp
   call spinseparate(tm(il,io)%ham)
   ! local eigenvectors in local rotated bases 
   i = (il-1)*2 + io;

   !write(115,'(a)')'------ global axes bases ---------'
   !write(115,'(5f15.10)') dble(tm(il,io)%ham(1:5,1:5))

   tm(il,io)%ham = matmul(Rl2T(:,:,i),tm(il,io)%ham)
   !write(115,'(20f15.10)') tm(il,io)%ham  
   ! write only up spin vec
   !tm(il,io)%ham(1:5,1:5) = matmul(Rl2(1:5,1:5,i),tm(il,io)%ham(1:5,1:5))
   !write(115,'(a)')'------ local axes bases ---------'
   write(115,'(5f15.10)') dble(tm(il,io)%ham(1:5,1:5))

	 !write(*,'(a,i10)') ' Rl2T.Rl2:    ===========> i = ',i
   !write(*,'(10f6.2)')  matmul(Rl2T(:,:,i),Rl2(:,:,i))

   
 enddo
 enddo
 close(115)
!-----------------------------------------------------------------


! save tm%dm data for orbital ordering... [and something ike wigner function??]
! tm(il,io)%dm(norbtm, nspin,norbtm, nspin)
 open(115,file='tmdms.dat', action='write',position='append')
 do il=1,nlayers
  do io=1, noctl
   ! per tm atom: 3 blocks each dim 10x10
   do ispn=1,2
   do jspn=ispn,2
    write(115,'(10f15.10)') tm(il,io)%dmc(:,ispn,:,jspn) ! complex: 10x10 dble values blocks
   end do
   end do 
 enddo
 enddo
 close(115)















	open(114,file='dh_B1.dat', action='write',position='append')
	write(114,'(5f15.10)') atm(1)%dh
	close(114)


	open(114,file='mom-B1.dat', action='write',position='append')
	write(114,'(5f18.10)') tm(1,1)%mag
	close(114)

	open(114,file='pop-B1.dat', action='write',position='append')
	write(114,'(10f8.4)') (dble(tm(1,1)%dm(i,1,i,1)),i=1,5), (dble(tm(1,1)%dm(i,2,i,2)),i=1,5)
	write(114,'(10f8.4)') (dble(tm(1,2)%dm(i,1,i,1)),i=1,5), (dble(tm(1,2)%dm(i,2,i,2)),i=1,5)
	write(114,'(10f8.4)') 
	close(114)


! 
!do i=1,natoms,4
!	write(*,*)"scf.f90: TM atom: ",i 
!	write(*,'(5f10.5)') atm(i)%dh
!enddo


do i=1,natoms
	if( mod(i-1,4) == 0 ) then
	 write(*,*)"scf.f90: TM atom: ",i 
	 write(*,'(a,i5,25f7.3)')"scf.f90: Qm: atom, i = ",i,qmpol(:,i)
	 !write(*,'(5f10.5)') atm(i)%dh
	else
		write(*,*)"O atom: ",i 
		!write(*,'(3f10.5)') atm(i)%dh
		write(*,'(a,i5,25f7.3)')"scf.f90: Qm: atom, i = ",i,qmpol(:,i)
	endif
enddo


if (1==0) then
	ddm = 0.0d0
	do i=5,9
		ddm = ddm + atm(1)%dh(i,i)
	end do
	write(*,'(a,f16.12)')'===> sum of diag dh = ',ddm
endif


 call mixerifc(-1, ddm) ! deallocate nu, mu, f, beta arrays

!write(6,'("Ebs, Euj, Eq, Eb = ", 4e20.6)') ebands, engyadu, engyes, energyb


!-------- -------- -------- -------- -------- -------------- 

! write engyadu to output, use it alongwith occupied states energy to get total energy....


! if evec, eval not needed anymore:
! for band structure calculations, new kpoint list:
!deallocate(hk, eval,evec,hksave, wke)
 
! save converged vmat that can be read in a future calculation, 
! avoids scf cycle in that calculation, even if it has a different k-grid.
!-------------------------------------------
if(lhu) then
  write(*,'(a)')"Writing e-e interaction matrices & Beff in STATE.OUT"
	open(10,file='STATE.OUT',form='FORMATTED',action='write')
	do il=1,nlayers
	 do io=1, noctl
    write(10,'(10i5)') il, io, tm(il,io)%ia, &
                       tm(il,io)%is, atom2orb(1:4,tm(il,io)%ia)
    do i=1,10
	   write(10,'(13G20.8)')tm(il,io)%vmat(i,1:10),tm(il,io)%beff(1:3)
	  end do 
	 end do ! io
	end do ! il
  close(10)
else
  write(*,'(a)')"lhu = F : Not writing STATE.OUT file"
endif
!-------------------------------------------








return
end subroutine groundstate
!=====================================================================




!=====================================================================
! makes kgrid and weight of k points in the irreducible wedge of BZ of tetragonal lattice.
! n1 point along x/y; n3 points along z.
	subroutine mkkgrid(n1,n3)
	implicit none
	integer, intent(in) :: n1,n3
	! local
	integer :: i,j,k, ind, nkslice, i1,i2
	double precision, dimension(3) :: k1,k2,k3,b1,b2,b3

	ntotk = (n1*(n1+1))/2 *n3; ! in irreducible wedge of the BZ of tetragonal lattice
	
	allocate(kgrid(3,ntotk))
	allocate(wk(ntotk))

	! gamma point calculation: 
	if(ntotk==1) then
		kgrid(:,1) = (/0.0d0,0.0d0,0.0d0 /);
		wk(1) = 1.0d0
		wknorm = 1.d0;

	write(*,*) 'scf: gamma point calculation!'
	return

	endif

	
	b1 = 0.5d0*bvec(:,1)/n1
	b2 = 0.5d0*bvec(:,2)/n1
	b3 = 0.5d0*bvec(:,3)/n3

! make a general triangular slice
	ind = 0;
	do i=0,n1-1
	 k1 = i*b1
	 do j=i,n1-1
	  k2 = j*b2
	   ind = ind + 1;
	   kgrid(:,ind) = k1 + k2; ! cartesian comp. 
	   if(i == 0) then
	    if(j == 0) then
	     wk(ind) = 0.125d0 ! 1/8
	    elseif(j == n1-1) then
	     wk(ind) = 0.25d0 ! 1/4 ! (1,1) corner 
	    else
	     wk(ind) = 0.5d0 ! 1/2 ! on x=0 line
	    endif
	   elseif(i == j .or. j== n1-1) then ! on diagonal (x=y) or right side (y=1); 
	   ! i==j==n1-1 reset to correct value 1/8 below after i,j loops
	    wk(ind) = 0.5d0 ! 1/2
	   else ! deep inside the wedge
	    wk(ind) = 1.0d0 ! 1
	   endif
	 end do
	end do
	! last point is at i=j= n1-1: change weight to 1/8
	wk(ind) = 0.125d0 ! 1/8

	nkslice = (n1*(n1+1))/2;

	! z> 0, z<1 slices:
	! use the slice to make slices at z>0
	i1=1;i2=1;
	
	do k=1,n3-1
	  k3 = k*b3
	  i1 = k*nkslice + 1;
	  i2 = (k+1)*nkslice;
	  do i=1,3
	   kgrid(i,i1:i2) = kgrid(i,1:nkslice) + k3(i);
	 end do
	 wk(i1:i2) = wk(1:nkslice)
	end do

	! z=0 slice; ! kgrid same; wk half of the above
	wk(1:nkslice) = 0.5d0 * wk(1:nkslice)
	! z=1 slice; 
	!write(*,*) 'scf: i1,i2: ',i1,i2

	wk(i1:i2) = 0.5d0 * wk(1:nkslice) ! i1,i2 set to correct value in the last iter of k loop above

	!write(*,'(10000f10.3)') wk

	wknorm = 1.0d0/sum(wk(:));

	!write(*,'(a,10000f10.5)') 'scf: wk = ',wknorm, wk

	
	return
	end subroutine mkkgrid
!=====================================================================


	subroutine spinseparate(vs)
	implicit none
	double complex, dimension(10,10), intent(inout) :: vs
	double complex, dimension(10,5) :: us
	double complex, dimension(5) :: v1, v2
	double precision:: norm1, norm2
	integer:: i, j1,j2

	j1=0; j2=0;
	do i=1,10
		v1 = vs(1:5,i) ! ith column
		v2 = vs(6:10,i)
		norm1 = dsqrt(sum(zabs(v1)**2))
		norm2 = dsqrt(sum(zabs(v2)**2))

		! set the state up or down depending on which spin has larger comp
		if(norm1 >= norm2) then
			vs(1:5,i) = v1/norm1;
			vs(6:10,i) = 0.0d0;
			j1=j1+1;
			us(j1,:) = vs(1:5,i); 
		else
			vs(6:10,i) = v2/norm2;
			vs(1:5,i) = 0.0d0
			j2=j2+1;
			us(5+j2,:) = vs(6:10,i);
		endif
	end do

	vs=0.0d0;
	vs(1:5,1:5) = us(1:5,:)
	vs(6:10,6:10) = us(6:10,:)

	return	
	end 	subroutine spinseparate
	
	



end module scf
