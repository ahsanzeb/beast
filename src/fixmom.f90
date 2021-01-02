
module fixmom
! use modmain, only: 

implicit none
real(8), parameter :: cb = 1.d0 ! g_e/4c 
real(8), allocatable, dimension(:,:,:) :: bfsmcmt


contains
!======================================================================
subroutine fsmbfield(iscf, energyb)
! !USES:
use modmain
implicit none
integer, intent(in) :: iscf
double precision, intent(out) :: energyb
! local variables
real(8) v1(3),v2(3),t1, b
integer :: il, io

!write(*,*) 'called: fsmbfield(iscf)..... 1'

! energy of physical global field




!.............................................................
! energy due to magnetic field: to be subtracted from the total energy.
!.............................................................
! E = - mu.B
energyb = 0.0d0;
! (Both mag and field from the previour iteration)
if (iscf > 1) then
 do il=1,nlayers
  do io=1,noctl
   energyb = energyb + sum(tm(il,io)%mag * tm(il,io)%beff)
	end do
 end do
end if


!.............................................................
! external fields:
!.............................................................
! initialise tm(il,io)%beff with external global and site fields
 do il=1,nlayers
  do io=1,noctl
   tm(il,io)%beff = cb*(tm(il,io)%bext + bfieldc)
	end do
end do
! reduce B's for the next iteration
Bfieldc = Bfieldc * reducebf
do il=1,nlayers
 do io=1,noctl
  tm(il,io)%bext = tm(il,io)%bext * reducebf
 end do
end do
!.............................................................
! fixed spin moment, constraint fields:
!.............................................................
if(iscf==1) then
  if (.not. allocated(bfsmcmt)) allocate(bfsmcmt(3,nlayers,noctl))
endif

!write(*,*) 'called: fsmbfield(iscf)..... 2 '

bfsmc = 0.0d0;

if (nspin == 1 .or. fsmtype == 0) return

! determine the global effective field
if ((abs(fsmtype).eq.1).or.(abs(fsmtype).eq.3)) then
  bfsmc = taufsm* (momtot - momfix)
 ! make sure that the constraining field is perpendicular to the fixed moment
 ! for fixed direction calculations (Y. Kvashnin and LN)
 write(*,'(a, 3f15.10)')'bfsmc = ', bfsmc
 if (fsmtype.lt.0) call r3vo(momfix,bfsmc)
 ! add to effective:
 write(*,'(a, 3f15.10)')'bfsmc = ', bfsmc
 do il=1,nlayers
  do io=1,noctl
   tm(il,io)%beff = tm(il,io)%beff + bfsmc
	end do
end do
end if

 !write(*,*) 'called: fsmbfield(iscf)..... 3 '

!. . . . . . . . . . . . . . . . . . . . . . . . .
if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
  ! determine the muffin-tin fields for fixed local moments
  do il=1,nlayers
    do io=1,noctl
      ! if any component is >= 1000 then do not fix the moment
      t1=sum(abs(tm(il,io)%mfix))
      if (t1 .ge. 1000.d0) cycle
      v2 = tm(il,io)%mag - tm(il,io)%mfix
      bfsmcmt(:,il,io) =  taufsm*v2
      ! fixed spin direction
      if (fsmtype.lt.0) call r3vo(tm(il,io)%mfix,bfsmcmt(:,il,io))
      tm(il,io)%beff = tm(il,io)%beff + bfsmcmt(:,il,io)
                         
    !write(*,'(a, 3f15.10)')'tm(il,io)%beff = ', tm(il,io)%beff
    end do
  end do
end if

!write(*,*) 'called: fsmbfield(iscf)..... 4 '
return
end subroutine
!EOC
!======================================================================
subroutine r3vo(x,y)
implicit none
! arguments
real(8), intent(in) :: x(3)
real(8), intent(inout) :: y(3)
! local variables
real(8) t1,t2
t1=x(1)**2+x(2)**2+x(3)**2
if (t1.lt.1.d-8) return
t2=(x(1)*y(1)+x(2)*y(2)+x(3)*y(3))/t1
y(:)=y(:)-t2*x(:)
return
end subroutine
!======================================================================
!======================================================================


end module fixmom
