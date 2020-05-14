
module mixing
implicit none

! variable with save attribute that are only required by routines in this module 
! can be defined here without creating any dependence of this module on modmain
real(8), allocatable, dimension(:) :: nu, mu, beta, f
!real(8), allocatable, dimension(:) :: v
integer :: n

contains
!======================================================================
! style and some code adopted from elk-6.2.8
! elk copyright: J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
!======================================================================
subroutine mixpack(iscl,tpack)
use modmain, only: nlayers, noctl, tm, natoms
use esvar, only: atm, types2norb
implicit none
! arguments
integer, intent(in) :: iscl
logical, intent(in) :: tpack
! local variables
complex(8), dimension(100) :: vz
integer :: i1,i2,i3, il,io, ia, it, norb

! reshape function order: rows first, usual fortran order.
! caution: 10x10 to 100, and 100 to 10x10 should use the same order.

if(tpack) then ! complex matrix to double vector 
 nu = 0.0d0
 i1 = 0;
 do il=1,nlayers
  do io=1,noctl
   i2 = i1 + 100;
   vz = reshape(tm(il,io)%vmat, (/100/));
   nu(i1+1:i2) = dble(vz);
   i1 = i2; i2 = i1 + 100; 
   nu(i1+1:i2) = aimag(vz);
   i1 = i2;
  end do
 end do

 ! append madelung term of the hamiltonian:
 do ia=1,natoms
  it = atm(ia)%it;
  norb = types2norb(it) * types2norb(it);
  i2 = i1 + norb;
  nu(i1+1:i2) = reshape(atm(ia)%dh, (/ norb /));
  i1 = i2;
 end do
 
else ! double vector to complex matrix
 i1 = 0;
 do il=1,nlayers
  do io=1,noctl
   !tm(il,io)%vmat = 0.0d0
   i2 = i1 + 100;
   i3 = i2 + 100;
   tm(il,io)%vmat = dcmplx(reshape(nu(i1+1:i2),  (/10,10/)), & 
                           reshape(nu(i2+1:i3),  (/10,10/))   )
   i1=i3;
  end do
 end do

 ! madelung term
 do ia=1,natoms
  it = atm(ia)%it;
  norb = types2norb(it);
  i2 = i1 + norb*norb;
  atm(ia)%dh = reshape(nu(i1+1:i2), (/ norb, norb /));
  i1 = i2;
 end do

endif

if (iscl == 1) then
mu = 0.0d0; ! nu;
endif 

return
end subroutine
!======================================================================
subroutine mixerifc(iscl,d)
use modmain, only: beta0, betamax, mtype, natoms
implicit none
! arguments
integer, intent(in) :: iscl
real(8), intent(out) :: d
real(8) :: t1, t0
integer :: i

if(iscl==0) then ! initialise mixer
 if (mtype==0) then
  n = 50*natoms + (natoms/4)*25 + (natoms/4)*3*9 ! TM vmat, TM dh, O dh
  allocate(nu(n))
  allocate(mu(n))
  d=1.d0
  return
 elseif(mtype==1) then
  n = 50*natoms + (natoms/4)*25 + (natoms/4)*3*9 ! TM vmat, TM dh, O dh
  allocate(nu(n))
  allocate(mu(n))
  allocate(beta(n))
  allocate(f(n))
  f(:)=0.d0
  beta(:)=beta0
  d=1.d0
  return
 else	
  write(*,*)'Error(mixerifc): wrong mtype...!'
 endif
elseif(iscl==-1) then ! deallocate work arrays
 if (mtype==0) then
  deallocate(nu)
  deallocate(mu)
  return
 elseif(mtype==1) then
  deallocate(nu)
  deallocate(mu)
  deallocate(beta)
  deallocate(f)
  return
 endif
endif

! call routines to do mixing....
if(mtype==0) then
 ! linear mixing
 t0=1.d0-beta0
 d=0.d0
 do i=1,n
  t1=nu(i)-mu(i)
  nu(i)=beta0*nu(i)+t0*mu(i)
  d=d+t1**2
  mu(i)=nu(i)
 end do
 d=sqrt(d/dble(n))
elseif(mtype==1)then
 ! adaptive linear mixing
 d=0.d0
 do i=1,n
  t1=nu(i)-mu(i)
  d=d+t1**2
  if (t1*f(i).ge.0.d0) then
    beta(i)=beta(i)+beta0
    if (beta(i).gt.betamax) beta(i)=betamax
  else
    beta(i)=0.5d0*(beta(i)+beta0)
  end if
  f(i)=t1
  nu(i)=beta(i)*nu(i)+(1.d0-beta(i))*mu(i)
  mu(i)=nu(i)
 end do
 d=sqrt(d/dble(n))
else
  write(*,*)
  write(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
  write(*,*)
  stop
endif

!write(*,*) 'mixerifc: d = ', d


return
end subroutine
!======================================================================

end module mixing
