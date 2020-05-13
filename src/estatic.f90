module estatic
use esvar
!use tbesel

implicit none

contains

!======================================================================
! We want to keep this module and its related routines seperate from 
! the rest of out code.
! This routines sets the global variables of this module
! It can be combined with initmadelung()
!======================================================================
subroutine setmadvar()
use modmain, only: oct, natoms, noctl, nlayers, avec, bvec, a, &
                  twopi, omega, nsptm, atom2species, nds
implicit none
integer :: ilm, i,l, m, io, il, ib, itm, is, ia, it, i1, i2

ldip = 2; ! dipole corrections, yes!

! ilm ranges for electronic orbitals in the basis set:
! allocate rhoc for these species/classes accordingly, 
! e.g., for Oxygen: allocate(rhoc(ilm12(1,0):ilm12(2,0))) etc.
! Oxygen: only p-orbitals: index ilm: 2-4
ilm12(1,1) = 2; ilm12(2,1) = 4;
! TM: only d-orbitals: index ilm: 5-9
ilm12(1,2) = 5; ilm12(2,2) = 9;
! two types: O & TM from itypes(ic/ib) ===> 1,2; 1 for O, 2 for TM
! set types from ic map
types2norb(1) = 3; ! O
types2norb(2) = 5; ! TM

nclass = nsptm + 1;
allocate(species2type(0:nsptm))
species2type(0) = 1; ! O: species number 0 has type number 1
species2type(1:nsptm) = 2; ! TM: species number >0 have type number 2

!integer, intent(in) :: ntot,ntms,nox, maxlv
nbas = natoms;
lmxl = 2  ! use > 2, using dummy here.
nclass = nsptm + 1;

! some sizes:
nlmi = (lmxl+1)*(lmxl+1);
lmxst = 2* lmxl;
nlm = (lmxst+1)*(lmxst+1);

allocate(atm(natoms))
do ib=1,nbas
 is = atom2species(ib);
 it = species2type(is); 
 atm(ib)%is = is
 atm(ib)%it = it; 
 i1 = ilm12(1,it); i2 = ilm12(2,it)
 allocate(atm(ib)%rhoc(i1:i2, i1:i2))
end do

allocate(qmpol(nlmi,nbas))

allocate(qpol(7,0:nsptm)) ! Crystal field constats for various species

qpol = 1.0d0; ! dummy. put it in modmain so that readinp can set it directly from input.in/default?

allocate(q0(0:nsptm)) ! neutral atom number of electrons
q0(0) = 4; ! Sr/A atom in perovskite gives 2 electrons; how to include them?
do is=1,nsptm
q0(1:nsptm) = nds
end do


! direct lattive
nxd = 10; nzd = 10; ! find a suitable number, check tbe code's method to find it, or its default.
! reciprocal lattive
nxg=10; nzg=10;

! lattice vectors:
s_lat%plat = avec/a
s_lat%alat = a
s_lat%vol = omega; ! nlayers * 2.0d0 * a**3;
s_lat%awald = 2.0d0 ! using dummy here.
s_lat%qlat = bvec*a/twopi;

vol = omega

! set positions of all atoms
allocate(s_lat%pos(3,nbas))
do il=1,nlayers
 do io=1,noctl ! 2
  ib = (il-1)*noctl + io; ! octrahedron number
  itm = (ib-1)*4 + 1 ! 4 atoms per octahedron
  s_lat%pos(1:3,itm) = oct(il,1)%rb/a
  do i=1,3
   s_lat%pos(1:3,itm+i) = oct(il,io)%ro(1,i)/a
  end do
 end do
end do


! glat, dlat : k-space and direct space lattice translational vectors for ewald sums
!       awald = s_lat%awald
!      nkd = s_lat%nkd
!      nkq = s_lat%nkq
!      gam = s_lat%gam
!      as = s_lat%as
!      ewtol = s_lat%tol
!     rpad = s_lat%rpad
!      nkdmx = s_lat%nkdmx
!      nkqmx = s_lat%nkqmx


! set up: glat, dlat; lattice vectors for ewald sums
 nkd = nxd* nxd* nzd
 allocate(dlat(3,nkd))
 call setEwaldvecs(nxd, nxd, nzd, avec, dlat)

 nkg = nxg* nxg* nzg
 allocate(glat(3,nkg))
 call setEwaldvecs(nxg, nxg, nzg, bvec, glat)




! set ll array to be used in hstrd() in mkstrxd.f90: 
allocate(ll(nlm))
ll = 0;
ilm = 0;
do l=0,lmxst
 do m = -l,l
  ilm = ilm + 1
  ll(ilm) = l
 end do
end do



! allocate and set crystal field. getM() can now just read this CFM.
 call setCFM()


return
end subroutine setmadvar

!======================================================================
! only one, before the SCF loop
! calculates the structure matrix B_{RL,R'L'}
! Ewald summation in real and reciprocal spaces
!======================================================================
subroutine initmadelung()
use esvar
use mscg
use mmakcg9
use mmkstrxd
use msylm, only: sylmnc
implicit none
	
 ! cy, gaunt coeff, struxidx and struxd

 call sylmnc(cy,16)
 call scg(9,cg,indxcg,jcg)
 call makcg9(indxcg,jcg,cg,gaunt,ak)

 !write(*,*)'Gaunt calculated.... '
 !write(*,'(a,100f10.3)') 'cg(1:100) = ', cg(1:100)

 allocate(struxd(nlmi,nlmi,nbas,nbas))
 !C   --- Structure constants for Ewald sums ---
 !nlmq = nlm; nlmq1 = nlm;
 call mkstrxd()

 !write(*,'(a)')'Structure matrix:'
 !write(*,'(a,1000e12.3)') 'struxd(:,:,1,1) = ',struxd(:,:,1,1)
 !write(*,'(a,1000e12.3)') 'struxd(:,:,2,1) = ',struxd(:,:,2,1)
 !write(*,'(a,1000e12.3)') 'struxd(:,:,8,1) = ',struxd(:,:,8,1)

 !allocate(rhoc(nl**4*nbas))
 !rhoc = 0.0_8

return
end 	subroutine initmadelung

!===============================================
subroutine setEwaldvecs(nx, ny, nz, c, v)
implicit none
integer, intent(in):: nx, ny, nz
double precision, dimension(3,3), intent(in) :: c
double precision, dimension(3,nx*ny*nz), intent(out) :: v
integer :: i,j,k, ind
double precision, dimension(3) :: px, py

ind =0;
do i=1, nx
 px = i*c(:,1)
 do j=1, ny
  py = j*c(:,2)
  do k=1, nz
   ind = ind + 1
   v(:,ind) = px + py + k*c(:,3)
  end do
 end do
end do

return
end subroutine setEwaldvecs
!===============================================
! getatomic() calculates atm%rhoc and atm%qup/qdn etc
subroutine getatomic()
use modmain, only: evec, wke, atom2species, atom2orb, norbtm, nspin, &
                   ntotk, ntot
implicit none
! IN: global evec & wke ( & a lot of other indexing arrays, and sizes, etc.. )
integer :: ik, ist
integer :: i,j,i1,i2, j1, i3, i4, ia,is, ispin, it, norb
double complex, dimension(norbtm, nspin) :: wf
double complex, allocatable, dimension(:,:) :: dmx
double precision, dimension(3) :: mag
double precision, dimension(2) :: qs
double precision :: dmij, qij
integer :: norbold

	write(*,*)'getatomic: testing... '
	
 norbold = 5;
 norb = norbold;
 allocate(dmx(norb,norb))

 do ia=1,nbas

  is = atom2species(ia); 
  i1 = atom2orb(1,ia); i2 = atom2orb(2,ia); ! orbital ranges, up spin
  i3 = atom2orb(3,ia); i4 = atom2orb(4,ia); ! orbital ranges, down spin
	it = species2type(is); 
	norb = types2norb(it);

	if(norb /= norbold) then
	 deallocate(dmx)
	 allocate(dmx(norb,norb))
	 norbold = norb
	endif

	! calculate rhoc
	atm(ia)%rhoc = 0.0d0;
  do ik=1,ntotk
   do ist=1,ntot
    ! assuming eigenvectors of Hk are arranged as columns
    wf = 0.0d0;
    wf(1:norb,1) = evec(ik,i1:i2,ist) 
    wf(1:norb,2) = evec(ik,i3:i4,ist)
    dmx = 0.0d0;
    qs = 0.0d0;

    do j=1,norb
     do i=1,norb
      dmij = 0.0d0;
      do ispin=1,nspin
       qij = conjg(wf(j,ispin)) * wf(i,ispin);
       dmij = dmij + qij
       qs(ispin) = qs(ispin) + qij
      end do
      dmx(j,i) = dmx(j,i) +  dmij
     end do ! i
    end do ! j
    
    atm(ia)%rhoc = atm(ia)%rhoc + dmx * wke(ik,ist) ! rhoc is allocated in a diff way, but i think it should work... 
    do ispin=1,nspin
     atm(ia)%qs(ispin) = atm(ia)%qs(ispin) + qs(ispin)* wke(ik,ist)
    end do
    
   end do ! ist
  end do ! ik


 end do ! ib

 deallocate(dmx)


return
end subroutine getatomic
!===============================================

end module estatic
