
module scf
use modmain
implicit none





contains

!=====================================================================
subroutine groundstate()
implicit none
integer :: iscf, ik

!-------- BZ integration -------- -------- -------- -------- 
 call mkkgrid(nk1,nk3) ! makes kgrid & wk for BZ integration, sets ntotk
	
 allocate(hk(ntot,ntot))
 allocate(eval(ntotk,ntot))
 allocate(evec(ntotk,ntot,ntot))
 allocate(wke(ntotk,ntot))

write(6,'("Starting SCF loop .... ")')
do iscf = 1, maxscf
 write(6,'("SCF iteration ",i5)'),iscf
 !-------- -------- -------- -------- -------- -------------- 
 ! Hamiltonian, eigenstates and eigenvalues: on kgrid in BZ
 !-------- -------- -------- -------- -------- -------------- 
 do ik= 1,ntotk
  kvec = kgrid(1:3,ik) ! already in cartesian coordinates
  call getHk(kvec, hk)
  call zdiag(ntot,hk,eval(ik,:),evec(ik,:,:),ik,ntot)
  !call useeigdata(ntot,hk) ! calc whatever we like at this k-point
  ! will take a weighted average after the k-loop completes.
 end do ! ik
 !-------- -------- -------- -------- -------- -------------- 
 ! find the Fermi level and wke = occupation weighted by wk 
 !-------- -------- -------- -------- -------- -------------- 
 call fermid(ntotk, wk, ntot, eval, temp, qtot, wke, efermi, nspin)
 write(*,'(a,f15.6)') "N_electron = ", qtot
 write(*,'(a,f15.6)') "Fermi energy = ", efermi
 !-------- -------- -------- -------- -------- -------------- 
 ! check convergence of scf:
 !-------- -------- -------- -------- -------- -------------- 

	


 !-------- -------- -------- -------- -------- -------------- 	
 ! calc occupations of all states, and then weighted averages now.
 ! call averages(efermi)
 ! Hubbard U potential matrices for TM atoms
 if(lhu) then
  call mkvmat()
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
