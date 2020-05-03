
module Hubbard
use modmain
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
integer :: i1,i2,i3,i4, is
double precision :: sum1
double precision, parameter :: fourpi=12.566370614359172954d0

do m1=-l,l
  i1 = l+m1+1
  do m2=-l,l
   i2 = l+m2+1
    do m3=-l,l
      i3 = l+m3+1
      do m4=-l,l
        i4 = l+m4+1
        ! all species in one go
        do is=1,nsptm
         sum1=0.d0
         do k=0,2*l,2
          sum1=sum1+Hub(is)%Fk(k)*gcmat(k/2,m4,m3,m2,m1) 
         end do
         Hub(is)%Vee(i1,i3,i2,i4)=fourpi*sum1
        end do ! is

      end do
    end do
  end do
end do

return
end subroutine mkvee
!======================================================================
! calculates density matrices of tm atoms in their d-orbital-spin space
! calculates electron-electron interaction potential matrices for all atoms
! also calculates magnetisation of tm atoms
! IN: global evec & wke ( & global Uz, Ur, and a lot of other indexing arrays, and sizes, etc.. )
subroutine mkvmat(iscf,ddm)
implicit none
integer, intent(in) :: iscf
double precision, intent(out) :: ddm
integer :: is,il,io,ik
integer :: i,j,i1,i2, j1,ist, i3, i4, ia, ispin,jspin
double complex, dimension(norbtm, nspin) :: wf
double complex, dimension(norbtm, nspin,norbtm, nspin) :: dm, dmc !dm in real, compelx Ylm
double complex, dimension(norbtm, nspin,norbtm, nspin) :: vmat, vmatc
double precision, dimension(3) :: mag
!logical, save :: first = .true.

ddm = 0.0d0

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
    wf = 0.0d0
    wf(1:norbtm,1) = evec(ik,i1:i2,ist) 
    wf(1:norbtm,2) = evec(ik,i3:i4,ist)
    do j=1,norbtm
     do jspin=1,nspin
      do i=1,norbtm
       do ispin=1,nspin
        dm(j,jspin,i,ispin) = dm(j,jspin,i,ispin) &
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
		call rtozflm(dm(:,jspin,:,ispin),dmc(:,jspin,:,ispin))
	 end do
	end do
	!...............................................................
	! vmat in complex spherical harmonics using dm and vee
	! FLL double counting correction, also calc mag
	!...............................................................
	call genvmat(is,norbtm,nspin,dmc,vmatc,mag)
	!...............................................................
	! convert vmat to real spherical harmonics, our basis
	!...............................................................
	do ispin=1,nspin
	 do jspin=1,nspin
		call ztorflm(vmatc(:,jspin,:,ispin),vmat(:,jspin,:,ispin))
	 end do
	end do
	! we can now perhaps add this vmat to our hamiltonian...
	!...............................................................
	! assign the results to global variables of TM atom concerned.

	if(iscf==1)then
   allocate(tm(il,io)%vmat(norbtms,norbtms))
   allocate(tm(il,io)%vmatold(norbtms,norbtms))
   allocate(tm(il,io)%dm(norbtm, nspin,norbtm, nspin))
  endif
	tm(il,io)%mag = mag
	tm(il,io)%dm = dm
	!tm(il,io)%vmat = vmat

	ddm = ddm + norm2(dble(dm))
	
	! unfold vmat to directly add in hk
	do i=1,norbtm ! 5
	do ispin=1,nspin
   i1 = (ispin-1)*norbtm + i
	 do j=1,norbtm
	 do jspin=1,nspin
    j1 = (jspin-1)*norbtm + j
		tm(il,io)%vmat(j1,i1)  = vmat(j,jspin,i,ispin)
   end do
   end do
  end do
  end do

 end do ! io
end do !il

return
end subroutine mkvmat

!======================================================================

subroutine rtozflm(fr,fz)
implicit none
double complex, dimension(5,5), intent(in) :: fr
double complex, dimension(5,5), intent(out) :: fz
 fz = matmul(Ur,fr);
 fz = matmul(fz,Uz);
return
end subroutine rtozflm
!======================================================================
subroutine ztorflm(fz,fr)
implicit none
double complex, dimension(5,5), intent(in) :: fz
double complex, dimension(5,5), intent(out) :: fr
 fr = matmul(Uz,fz);
 fr = matmul(fr,Ur);
return
end subroutine ztorflm
!======================================================================
! calculates electron-electron interaction potential matrices for all atoms
! also calculates magnetisation of tm atoms
! dm & vmat both in complex Ylm basis; adapted from ELK code
subroutine genvmat(is,norbtm,nspin,dm,vmat,mg)
implicit none
integer, intent(in) :: is,norbtm,nspin
double complex, dimension(norbtm,nspin,norbtm,nspin), intent(in) ::dm
double complex, dimension(norbtm,nspin,norbtm,nspin), intent(out) ::vmat
double precision, dimension(3), intent(out) :: mg
! local
integer ispn,jspn
integer m1,m2,m3,m4,nm
complex(8) z1,z2
double precision :: U, J, n
!complex(8), parameter :: iota = dcmplx(0.0d0,1.0d0)
double complex, dimension(nspin,nspin) :: dms

!-----------------------------------------------------
! spin density matrix
!-----------------------------------------------------
dms(:,:)=0.d0
do ispn=1,nspin
do jspn=1,nspin
 do m1=1,5
  dms(ispn,jspn)=dms(ispn,jspn)+dm(m1,ispn,m1,jspn)
 end do
end do
end do
! trace over spin
n=dble(dms(1,1))
if (nspin==2) n=n+dble(dms(2,2))
! magnetisation
if (nspin==2) then
 mg(:)=0.d0
 mg(3)=dble(dms(1,1)-dms(2,2))
! non-collinear terms
 mg(1)=dble(dms(1,2)+dms(2,1))
 mg(2)=dble(iota*(dms(1,2)-dms(2,1)))
 write(*,'(a,10f15.6)') 'mom: ', mg
end if
!-----------------------------------------------------
! vmat
!-----------------------------------------------------
vmat = 0.0d0
! calculation of DFT+U potential and energy
! begin loops over m1 and m2
do m1=1,5 ! for magnetic qunatum number ms = -2:2
do m2=1,5
! begin loops over m3 and m4
 do m3=1,5
 do m4=1,5
  do ispn=1,nspin
  do jspn=1,nspin
   z1=dm(m2,ispn,m1,ispn)*dm(m4,jspn,m3,jspn)
   z2=dm(m4,jspn,m1,ispn)*dm(m2,ispn,m3,jspn)
   !engyadu(ia,i)=engyadu(ia,i)+dble(z1-z2)*Hub(is)%Vee(m1,m3,m2,m4)
   vmat(m1,ispn,m2,ispn)=vmat(m1,ispn,m2,ispn) &
                 + dm(m4,jspn,m3,jspn)*Hub(is)%Vee(m1,m3,m2,m4)
   vmat(m1,ispn,m4,jspn)=vmat(m1,ispn,m4,jspn) &
                 - dm(m2,ispn,m3,jspn)*Hub(is)%Vee(m1,m3,m2,m4)
  end do
  end do
! end loops over m3 and m4
 end do
 end do
! end loops over m1 and m2
end do
end do
!-----------------------------------------------------
! double counting correction: FLL
!-----------------------------------------------------
! non-collinear case
! correction to the potential
! U, J,  dms, n
U = Hub(is)%U
J = Hub(is)%J
do m1=1,5
 vmat(m1,1,m1,1)=vmat(m1,1,m1,1) &
             -u*(n-0.5d0)+j*(dms(1,1)-0.5d0)
 vmat(m1,2,m1,2)=vmat(m1,2,m1,2) &
             -u*(n-0.5d0)+j*(dms(2,2)-0.5d0)
 vmat(m1,1,m1,2)=vmat(m1,1,m1,2)+j*dms(1,2)
 vmat(m1,2,m1,1)=vmat(m1,2,m1,1)+j*dms(2,1)
end do
!-----------------------------------------------------

return
end subroutine genvmat
!======================================================================


end module Hubbard
