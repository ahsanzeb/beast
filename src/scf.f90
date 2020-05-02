
module scf
use modmain
use hamiltonian
use fermi
use hubbard

implicit none





contains

!=====================================================================
subroutine groundstate()
implicit none
integer :: iscf, ik
double precision, dimension(3) :: kvec
double precision :: ddmold, ddm

! some large number
ddmold = 100.0d0;

!-------- BZ integration -------- -------- -------- -------- 
 call mkkgrid(nk1,nk3) ! makes kgrid & wk for BZ integration, sets ntotk
	
 allocate(hk(ntot,ntot))
 allocate(eval(ntotk,ntot))
 allocate(evec(ntotk,ntot,ntot))
 allocate(wke(ntotk,ntot))
! if(lhu)then
!	allocate(tm(il,io)%vmat(norbtms, nspin,norbtms, nspin))
!	allocate(tm(il,io)%dm(norbtms, nspin,norbtms, nspin))
! endif

write(6,'("Starting SCF loop .... ")')
do iscf = 1, maxscf
 write(6,'("SCF iteration ",i5)') iscf
 !-------- -------- -------- -------- -------- -------------- 
 ! Hamiltonian, eigenstates and eigenvalues: on kgrid in BZ
 !-------- -------- -------- -------- -------- -------------- 
 do ik= 1,ntotk
  kvec = kgrid(1:3,ik) ! already in cartesian coordinates
  call getHk(kvec, hk, iscf)
  call zdiag(ntot,hk,eval(ik,:),ik,ntot)
  evec(ik,:,:) = Hk ! eigenvector are columns of Hk
  !call useeigdata(ntot,hk) ! calc whatever we like at this k-point
  ! will take a weighted average after the k-loop completes.
 end do ! ik
 !-------- -------- -------- -------- -------- -------------- 
 ! find the Fermi level and wke = occupation weighted by wk 
 !-------- -------- -------- -------- -------- -------------- 
 call fermid(ntotk, wk, ntot, eval, temp, qtot, wke, efermi, nspin)
 write(*,'(a,f15.6)') "N_electron = ", qtot
 write(*,'(a,f25.10)') "Fermi energy = ", efermi
 !-------- -------- -------- -------- -------- -------------- 
 ! check convergence of scf:
 !-------- -------- -------- -------- -------- -------------- 

	


 !-------- -------- -------- -------- -------- -------------- 	
 ! calc occupations of all states, and then weighted averages now.
 ! call averages(efermi)
 ! Hubbard U potential matrices for TM atoms
 if(lhu) then
  call mkvmat(iscf, ddm)
  	write(6,*) "SCF mixing.... not implemented yet... do it."
  if(dabs(ddmold - ddm) < toldm) then
	 write(6,'("SCF coverged in ",i5," iterations!")') iscf
	 write(6,'("norm(dm), delta(norm(dm)) = ", 2f20.15)') &
	              ddm, dabs(ddmold - ddm)
	 exit ! exit scf loop
  endif
  ddmold = ddm
 endif
 !-------- -------- -------- -------- -------- -------------- 


 !
end do! iscf
!-------- -------- -------- -------- -------- -------------- 

! if evec, eval not needed anymore:
! for band structure calculations, new kpoint list:
if(allocated(eval)) deallocate(eval)
if(allocated(hk)) deallocate(hk)
 









return
end subroutine groundstate
!=====================================================================

end module scf
