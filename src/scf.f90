
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
subroutine groundstate()
implicit none
integer :: iscf, ik
double precision, dimension(3) :: kvec
double precision :: ddmold, ddm, engyadu, ebands, engyes, energyb
integer :: il,io,i

engyes = 0.0d0
! some large number
ddmold = 1.0d8;

!-------- BZ integration -------- -------- -------- -------- 
 call mkkgrid(nk1,nk3) ! makes kgrid & wk for BZ integration, sets ntotk

 !write(*,*) 'scf: kgrid = ', kgrid

 allocate(hk(ntot,ntot))
 ! if larger sys or too large kgrid, then we can save mem by using sparse format for hksave, (and even for hk...)
 allocate(hksave(ntot,ntot,ntotk)) ! for mixing.... in hamiltonain module
 allocate(eval(ntotk,ntot))
 allocate(evec(ntotk,ntot,ntot))
 allocate(wke(ntotk,ntot))
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
 call tbeseld(engyes)


! initialise the mixer
 call mixerifc(0, ddm) ! allocate nu, mu, f, beta arrays


do iscf = 1, maxscf
 write(6,'("SCF iteration ",i5)') iscf
 !-------- -------- -------- -------- -------- -------------- 
 ! Hamiltonian, eigenstates and eigenvalues: on kgrid in BZ
 !-------- -------- -------- -------- -------- -------------- 
 do ik= 1,ntotk
  kvec = kgrid(1:3,ik) ! already in cartesian coordinates
  call getHk(ik,kvec, hk, iscf)
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
  ! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
  if(lesH) call tbeseld(engyes)
  ! now we can add atm%dh to the hamiltonian....
 !-------- -------- -------- -------- -------- -------------- 	
  ! Hubbard U potential matrices for TM atoms
  call mkvmat(iscf, engyadu) ! uses global evec
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
  if(ddm < toldm) then
	 write(6,'("SCF coverged in ",i5," iterations!")') iscf
	 write(6,'("tolerance, change in Vmat & dH = ", 2e20.6)') &
	              toldm, ddm
	 !write(6,'("Ebs + Euj = ", 2e20.6)') ebands, engyadu
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
	wk(i1:i2) = 0.5d0 * wk(1:nkslice) ! i1,i2 set to correct value in the last iter of k loop above

	!write(*,'(10000f10.3)') wk

	wknorm = 1.0d0/sum(wk(:));

	return
	end subroutine mkkgrid
!=====================================================================





end module scf
