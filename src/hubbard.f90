
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


!write(*,*)'Hubb: F(k): ', Hub(1)%Fk(:)

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

!write(*,*)"========>>>> hubbard: testing, setting Vee=0"
!do is=1,nsptm
!Hub(is)%Vee = 0.0d0
!end do

return
end subroutine mkvee
!======================================================================
! calculates density matrices of tm atoms in their d-orbital-spin space
! calculates electron-electron interaction potential matrices for all atoms
! also calculates magnetisation of tm atoms
! IN: global evec & wke ( & global Uz, Ur, and a lot of other indexing arrays, and sizes, etc.. )
subroutine mkvmat(iscf,edu) ! ddm,
implicit none
integer, intent(in) :: iscf
double precision, intent(out) :: edu ! ddm
integer :: is,il,io,ik
integer :: i,j,i1,i2, j1,ist, i3, i4, ia, ispin,jspin
double complex, dimension(norbtm, nspin) :: wf
double complex, dimension(norbtm, nspin,norbtm, nspin) :: dm, dmc !dm in real, compelx Ylm
double complex, dimension(norbtm, nspin,norbtm, nspin) :: vmat, vmatc
double precision, dimension(3) :: mag
double precision :: edua
!logical, save :: first = .true.

!ddm = 0.0d0
edu = 0.0d0;
edua = 0.0d0;
momtot = 0.0d0

do il=1,nlayers
 do io=1,noctl
  ia = tm(il,io)%ia;
  is = atom2species(ia);
  i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges, up spin
  i3 = atom2orb(3,ia); i4 = atom2orb(4,ia); ! orbital ranges, down spin
	!...............................................................
	! calculate dm
	!...............................................................
	dm = 0.0d0; dmc = 0.0d0;
  do ik=1,ntotk
   do ist=1,ntot
    ! assuming eigenvectors of Hk are arranged as columns
    wf = 0.0d0
    wf(1:norbtm,1) = evec(ik,i1:i2,ist) 
    wf(1:norbtm,2) = evec(ik,i3:i4,ist)
    dmc = 0.0d0;
    do j=1,norbtm
     do jspin=1,nspin
      do i=1,norbtm
       do ispin=1,nspin
        dmc(j,jspin,i,ispin) = conjg(wf(j,jspin))*wf(i,ispin) ! dmc used as a dummy here.
       end do
      end do
     end do
    end do
    !if(ik==1 .and. ist>72 .and. il==1 .and. io==1) then
    	!	write(*,*)'hubbard: dm_{ik,ist} = ', il,io
    	!	write(*,'(10f8.3)') cdabs(dmc)
		!endif
    dm = dm + dmc * wke(ik,ist);
   end do ! ist
  end do ! ik
  dmc = 0.0d0; ! reset.

	!write(*,*)'hubbard: dm: il,io = ', il,io
	!write(*,'(10f8.3)') dble(dm)

	!...............................................................
	! convert dm to complex spherical harmonics
	! mind our order of real harmoncis: m2i list; i2m list?
	!...............................................................
	do ispin=1,nspin
	 do jspin=1,nspin
		call rtozflm(dm(:,jspin,:,ispin),dmc(:,jspin,:,ispin))
	 end do
	end do
	! save dmc
	tm(il,io)%dmc = dmc
	!...............................................................
	! vmat in complex spherical harmonics using dm and vee
	! FLL double counting correction, also calc mag
	!...............................................................
	call genvmat(is,norbtm,nspin,dmc,vmatc,mag,edua)
	edu = edu + edua;
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

	if(iscf==1 .and. .not. allocated(tm(il,io)%vmat))then
   allocate(tm(il,io)%vmat(norbtms,norbtms))
   allocate(tm(il,io)%vmatold(norbtms,norbtms))
   allocate(tm(il,io)%dm(norbtm, nspin,norbtm, nspin))
   tm(il,io)%vmatold = 0.0d0
  endif
	tm(il,io)%mag = mag
	tm(il,io)%dm = dm
	momtot = momtot + mag;

	
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

	!ddm = ddm + norm2(dble(tm(il,io)%vmat-tm(il,io)%vmatold))

	!write(*,'(a,2i5)')'il,io, tm(il,io)%vmat', il,io
	!write(*,'(10f8.3)')tm(il,io)%vmat

 end do ! io
end do !il

!ddm = ddm/dble(norbtms*0.25d0*natoms)

return
end subroutine mkvmat

!======================================================================

!======================================================================
! Hubbard interaction: for TM-only case: when we solve tm atoms separately:
!======================================================================
! calculates density matrices of tm atoms in their d-orbital-spin space
! calculates electron-electron interaction potential matrices for all atoms
! also calculates magnetisation of tm atoms
! IN: global evec & wke ( & global Uz, Ur, and a lot of other indexing arrays, and sizes, etc.. )
subroutine mkvmattm(iscf, il, io, is, evec, wke, edu) ! ddm,
implicit none
integer, intent(in) :: iscf, il, io, is
double complex, dimension(10, 10), intent(in) :: evec
double precision, dimension(10), intent(in) :: wke
double precision, intent(out) :: edu ! ddm
integer :: i,j,i1,i2, j1,ist, i3, i4, ia, ispin,jspin
double complex, dimension(norbtm, nspin) :: wf
double complex, dimension(norbtm, nspin,norbtm, nspin) :: dm, dmc !dm in real, compelx Ylm
double complex, dimension(norbtm, nspin,norbtm, nspin) :: vmat, vmatc
double precision, dimension(3) :: mag
double precision :: edua
!logical, save :: first = .true.

!ddm = 0.0d0
edu = 0.0d0;
edua = 0.0d0;
momtot = 0.0d0

	!...............................................................
	! calculate dm
	!...............................................................
	dm = 0.0d0;
   do ist=1, 10
    wf = 0.0d0
    wf(1:norbtm,1) = evec(1:5,ist) ! assuming eigenvec in columns
    wf(1:norbtm,2) = evec(6:10,ist)
    dmc = 0.0d0;
    do j=1,norbtm
     do jspin=1,nspin
      do i=1,norbtm
       do ispin=1,nspin
        dmc(j,jspin,i,ispin) = conjg(wf(j,jspin))*wf(i,ispin) ! dmc used as a dummy here.
       end do
      end do
     end do
    end do
    dm = dm + dmc* wke(ist);
   end do ! ist
  dmc = 0.0d0; ! reset.

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
	call genvmat(is,norbtm,nspin,dmc,vmatc,mag,edua)
	edu = edu + edua;
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

	if(iscf==1 .and. .not. allocated(tm(il,io)%vmat))then
   allocate(tm(il,io)%vmat(norbtms,norbtms))
   allocate(tm(il,io)%vmatold(norbtms,norbtms))
   allocate(tm(il,io)%dm(norbtm, nspin,norbtm, nspin))
   tm(il,io)%vmatold = 0.0d0
  endif
	tm(il,io)%mag = mag
	tm(il,io)%dm = dm
	momtot = momtot + mag;

	
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

return
end subroutine mkvmattm

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
subroutine genvmat(is,norbtm,nspin,dm,vmat,mg,engyadu)
implicit none
integer, intent(in) :: is,norbtm,nspin
double complex, dimension(norbtm,nspin,norbtm,nspin), intent(in) ::dm
double complex, dimension(norbtm,nspin,norbtm,nspin), intent(out) ::vmat
double precision, dimension(3), intent(out) :: mg
double precision, intent(out) :: engyadu
! local
integer ispn,jspn
integer m1,m2,m3,m4,nm
complex(8) z1,z2
double precision :: U, J, n, sum1, edc, v
!complex(8), parameter :: iota = dcmplx(0.0d0,1.0d0)
double complex, dimension(nspin,nspin) :: dms
double precision, dimension(norbtm) :: n0
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
!write(*,*)"hubbard: testing: using average occupation of the orbital n0"
!n = n/dble(nspin*norbtm) ! average occupation of the orbital


!write(*,*)"hubbard: testing: should it be occ of individual orb n_m???"
!n0 = 0.0d0
!do m1=1,norbtm
!	do ispn=1,nspin
!		n0(m1)= n0(m1) + dble(dm(m1,ispn,m1,ispn)) ! diag of dm should be real
!	end do
!end do

!write(*,'(a,5f10.4)')'n0 : ', n0

! magnetisation
if (nspin==2) then
 mg(:)=0.d0
 mg(3)=dble(dms(1,1)-dms(2,2))
! non-collinear terms
 mg(1)=dble(dms(1,2)+dms(2,1))
 mg(2)=dble(iota*(dms(1,2)-dms(2,1)))
 !write(*,'(a,10f15.6)') 'mom: ', mg
end if
!-----------------------------------------------------
! vmat
!-----------------------------------------------------
vmat = 0.0d0;
engyadu=0.0d0;
! calculation of DFT+U potential and energy
! begin loops over m1 and m2
do m1=1,5 ! for magnetic qunatum number ms = -2:2
do m2=1,5
! begin loops over m3 and m4
 do m3=1,5
 do m4=1,5
  v = Hub(is)%Vee(m1,m3,m2,m4)
  do ispn=1,nspin
  do jspn=1,nspin
   z1=dm(m2,ispn,m1,ispn)*dm(m4,jspn,m3,jspn)
   z2=dm(m4,jspn,m1,ispn)*dm(m2,ispn,m3,jspn)
   engyadu = engyadu + dble(z1-z2)*v
   vmat(m1,ispn,m2,ispn)=vmat(m1,ispn,m2,ispn) &
                 + dm(m4,jspn,m3,jspn)*v
   vmat(m1,ispn,m4,jspn)=vmat(m1,ispn,m4,jspn) &
                 - dm(m2,ispn,m3,jspn)*v
  end do
  end do
! end loops over m3 and m4
 end do
 end do
! end loops over m1 and m2
end do
end do
! multiply energy by factor 1/2
engyadu=0.5d0*engyadu
!-----------------------------------------------------
! double counting correction: FLL
!-----------------------------------------------------
! non-collinear case
! correction to the potential
! U, J,  dms, n
!write(*,*)'hubbard: testing: setting u,j=0 for DC correction... '

U = Hub(is)%U*ev2har
J = Hub(is)%J*ev2har
! U,J given in eV in input, and kept in eV.

! correction to the energy
edc=0.5d0*u*n*(n-1.d0)
edc=edc-0.5d0*j*dble(dms(1,1)*(dms(1,1)-1.d0))
edc=edc-0.5d0*j*dble(dms(2,2)*(dms(2,2)-1.d0))
edc=edc-0.5d0*j*dble(dms(1,2)*dms(2,1))
edc=edc-0.5d0*j*dble(dms(2,1)*dms(1,2))
engyadu=engyadu-edc

do m1=1,5
 vmat(m1,1,m1,1)=vmat(m1,1,m1,1) &
             -u*(n-0.5d0)+j*(dms(1,1)-0.5d0)
 vmat(m1,2,m1,2)=vmat(m1,2,m1,2) &
             -u*(n-0.5d0)+j*(dms(2,2)-0.5d0)
 vmat(m1,1,m1,2)=vmat(m1,1,m1,2)+j*dms(1,2)
 vmat(m1,2,m1,1)=vmat(m1,2,m1,1)+j*dms(2,1)
end do

!-----------------------------------------------------
! trace of dmat times vmat
sum1=0.d0
do ispn=1,nspin
 do m1=1,5
  do jspn=1,nspin
   do m2=1,5
    sum1=sum1+dble(dm(m1,ispn,m2,jspn)*vmat(m2,jspn,m1,ispn))
   end do
  end do
 end do
end do
! subtract contribution to the energy of vmat potential
engyadu=engyadu-sum1


return
end subroutine genvmat
!======================================================================








! To compare with Kugel-Khomskii's model, we need to have their Hubbard model
! that has simple form (usual, rotationaly non-invariant type).

subroutine genvmatkk()
implicit none

return
end subroutine genvmatkk
!======================================================================




end module Hubbard
