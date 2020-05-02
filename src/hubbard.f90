
module Hubbard
implicit none

!======================================================================
! once in the program for a given set of U&J
!======================================================================
!1. calc Slater-integrals Fk for each TM species using given U & J
!2. calc Vee matrices for each TM species using Fk and gcmat
! note: Vee will be in complex spherical harmonics
!======================================================================

!======================================================================
! at each SCF iteration, for each TM atom
!======================================================================
!1. calc dmr in d-space of each individual TM atoms
!note: dmr will be in our basis of real spherical harmonics
! 2. change the basis of dmr to the complex Ylm, let's call it dmc
! 3. calc vmatc using Vee and dmc
! 4. change the basis of vmatc to our real Ylm, let's call it vmatr
! we can now use vmatr in construction of hamiltonian, hk.
!======================================================================

contains
!======================================================================
! Make Vee from Fk and gaunt matrix gcmat; for all TM species
!======================================================================
! gcmat calculated by a call in main.f:
! make gaunt coefficients matrix that will be used to calc Vee 
!allocate(gcmat(0:2,-2:2,-2:2,-2:2,-2:2)) ! for l=2, d orbitals
!call mkdgaunt(2,gcmat)

! Hub%fk set in readinput, Hub%Vee allocated there as well.
! allocate(Hub(nsptm))
! allocate(Hub(is)%Vee(-2:2,-2:2,-2:2,-2:2))
subroutine mkvee(l)
implicit none
integer, intent(in) :: l
! local
integer :: m1,m2,m3,m4,k
double precision :: sum1
double precision, parameter :: fourpi=12.566370614359172954d0

do m1=-l,l
  do m2=-l,l
    do m3=-l,l
      do m4=-l,l
        ! all species in one go
        do is=1,nsptm
         sum1=0.d0
         do k=0,2*l,2
          sum1=sum1+Hub(is)%Fk(k)*gcmat(k/2,m4,m3,m2,m1) 
         end do
         Hub(is)%Vee(m1,m3,m2,m4)=fourpi*sum1
        end do ! is

      end do
    end do
  end do
end do

return
end subroutine mkvee
!======================================================================
! at each k, during k-loop
subroutine mktmdmk()
implicit none
integer :: is,il,io,ik
integer :: i,i1,i2, ist
double complex, dimension(norbtms, nspin) :: wf
double complex, dimension(norbtms, nspin,norbtms, nspin) :: dm, dmc

do il=1,nlayers
 do io=1,noctl
  ia = tm(il,io)%ia;
  is = atom2species(ia);
  i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges, up spin
  i3 = atom2orb(3,ia); i4 = atom2orb(4,ia); ! orbital ranges, down spin
	!...............................................................
	! calculate dm
	!...............................................................
	dm = 0.0d0;
  do ik=1,ntotk
   do ist=1,ntot
    ! assuming eigenvectors of Hk are arranged as columns
    wf(1:norbtm,1) = evec(ik,i1:i2,ist) 
    wf(1:norbtm,2) = evec(ik,i3:i4,ist)
    do j=1,norbtm
     do jspin=1,nspin
      do i=1,norbtm
       do ispin=1,nspin
       ! oct(il,ia)%
        dm(j,jpsin,i,ispin) = dm(j,jpsin,i,ispin) &
         + conjg(wf(j,jspin))*wf(i,ispin) * wke(ik,ist)
       end do
      end do
     end do
    end do
   end do ! ist
  end do ! ik
	!...............................................................
	! convert dm to complex spherical harmonics
	! mind our order of real harmoncis: m2i list; i2m list?
	!...............................................................
	do ispin=1,nspin
	 do jspin=1,nspin
		call rtozflm(n,dm(:,jspin,:,isin),dmc(:,jspin,:,isin))
	 end do
	end do





	!...............................................................
	! vmat in complex spherical harmonics using dm and vee
	!...............................................................
	! 


	!...............................................................
	! convert vmat to real spherical harmonics, our basis
	!...............................................................




 end do ! io
end do !il






return
end subroutine mktmdmk

!======================================================================

subroutine rtozflm(n,fr,fz)
implicit none
double complex, dimension(n,n), intent(in) :: fr
double complex, dimension(n,n), intent(out) :: fz
 fz = matmul(Ur,fr);
 fz = matmul(fz,Uz);
return
end subroutine ztorflm
!======================================================================
subroutine ztorflm(n,fz,fr)
implicit none
double complex, dimension(n,n), intent(in) :: fz
double complex, dimension(n,n), intent(out) :: fr
 fr = matmul(Uz,fz);
 fr = matmul(fr,Ur);
return
end subroutine ztorflm
!======================================================================








subroutine genfdu(i,u,j,f)
implicit none

u=ujdu(1,i)
j=ujdu(2,i)
f(:)=fdu(:,i)
lambda=lambdadu(i)

f(0)=u
! r1 = F(4)/F(2), see PRB 52, R5467 (1995)
r1=0.625d0
f(2)=(14.d0*j)/(1.d0+r1)
f(4)=f(2)*r1

ujdu(1,i)=u
ujdu(2,i)=j
fdu(:,i)=f(:)
lambdadu(i)=lambda





















end module Hubbard
