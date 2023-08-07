module estatic
use esvar

implicit none

contains

!======================================================================
! We want to keep this module and its related routines seperate from 
! the rest of our code.
! This routines sets the global variables of this module
! It can be combined with initmadelung()
!======================================================================
subroutine setmadvar()
use modmain, only: oct, natoms, noctl, nlayers, avec, bvec, a, &
                  twopi, omega, nsptm, atom2species, nds, nspin, Dcf, &
                  ainv, ewalda, ewaldnr, ewaldnk, hardU, nround, &
                  pos, posA, noct, qa
implicit none
integer :: ilm, i,l, m, io, il, ib, itm, is, ia, it, i1, i2
!integer :: nvevEwalsR, nvevEwalsK
double precision :: rcut

! round struxd and strucdA to nround-th decimal.
decimal = dble(10**nround);

ldip = 0; ! dipole corrections, yes!

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

!nclass = nsptm + 1;
allocate(species2type(0:nsptm));
species2type(0) = 1; ! O: species number 0 has type number 1
species2type(1:nsptm) = 2; ! TM: species number >0 have type number 2

noctl=2; ! =1 for fake cubic
nbas = natoms;
lmxl = 4  ! max l for potential... 
nclass = nsptm + 1;
nbasA = nbas/4; ! A-sites = number of octahedra


! some sizes:
nlmi = (lmxl+1)*(lmxl+1);
lmxst = 2* lmxl; !lmx_structure_matrix
nlm = (lmxst+1)*(lmxst+1); 
nsp = nspin;

write(*,*) 'nlmi,nlm = ',nlmi,nlm
!write(*,*)'natoms = ',natoms
allocate(qmpol(nlmi,nbas))
qmpol = 0.0d0;

allocate(qref(0:nsptm))
allocate(atm(natoms))
do ib=1,nbas
 is = atom2species(ib);
 it = species2type(is); 
 
 !write(*,*) 'ib, is, it : ',ib, is, it
 
 atm(ib)%is = is
 atm(ib)%it = it; 
 i1 = ilm12(1,it); i2 = ilm12(2,it)
 allocate(atm(ib)%rhoc(i1:i2, i1:i2))
 allocate(atm(ib)%dh(i1:i2, i1:i2))
 ! set initial guess: O 2e added; TM 4e lost.
 if(is==0) then ! set as monopoles
  qmpol(1,ib) =  2.0d0; !2.0d0; !2.0d0! 0.0d0; ! +2.0d0
 else
! FCC NaCl sturc if only B are +- charged
! madelung const: 1.746941
! if(ib==1 .or. ib==13) then 
!	qmpol(1,ib) =  1.0d0
! else
! 	qmpol(1,ib) =  1.0d0
! endif
  qmpol(1,ib) = -(6.0d0-qa) ! +1.0d0 !1.0d0 ! -4.0d0
 endif
 qref(is) = qmpol(1,ib) ! reference charge state for hardness term
end do

! extra normal atoms
!do io=1,noct
!	qmpol(1,noct*4+io) = -2.0d0
!end do

allocate(hard(0:nsptm))
hard = hardU ! hardU [modmain] read from the input file [readinput]


allocate(qpol(7,0:nsptm)) ! Crystal field constats for various species

qpol(:,0) = Dcf(:,1) ! O
do is=1,nsptm
 qpol(:,is) = Dcf(:,2) ! TM
end do

qmpolA = -qa; !-4.0d0 ! -2.0d0 !-1.0d0 ! -2.0d0;

allocate(q0(0:nsptm)) ! neutral atom number of electrons
q0(0) =  4.0; ! Oxygen q0 in p orbitals ! Sr/A atom in perovskite gives 2 electrons; how to include them?
do is=1,nsptm
q0(1:nsptm) = nds  !+ 2.0d0; ! +2 for TM s electrons;
                            ! for qmpol, q0 is assumed shperical symmetric, so consistent with s orbit.
end do

!nvevEwalsR = ewaldnr;
!nvevEwalsK = ewaldnk;

! direct lattice:
nxd = ewaldnr; 
nzd = nxd; !max(1,int(nxd*dsqrt(2.0d0)/dble(nlayers)));
! r_cut:
rcut = dble(nxd*dsqrt(2.0d0)*a);

! reciprocal lattive
nxg=ewaldnk; 
nzg= nxg !max(1,int(dble(nxg*nlayers)/dsqrt(2.0d0)));

! recommended value for alpha, ewald convergence parameter (named s_lat%awald below). 
! https://wanglab.hosted.uark.edu/DLPOLY2/node114.html
! s_lat%awald = 0.32/r_cut
!Also see: http://ambermd.org/Questions/ewald.html


!write(*,*) 'a = ',a,  ' avec  below'
!write(*,'(3f10.5)') avec
!write(*,'(3f10.5)') ainv
!write(*,'(3f10.5)') matmul(avec,ainv)

! lattice vectors:
s_lat%plat = avec/a ! *1/a to make plat dimensionless
s_lat%alat = a 
s_lat%vol = omega; ! nlayers * 2.0d0 * a**3;
! s_lat%awald has dimensions of [L^-1], i.e., 1/alat
s_lat%awald = ewalda*0.32d0/rcut; ! r_cut ~ nvevEwals*a;
s_lat%qlat = ainv*a !bvec*(a/twopi); ! plat.qlat = I and not 2pi for directshortn()



write(*,*) 'nlayers = ',nlayers
write(*,*) 'nxd, nzd, nxg, nzg, alpha*a = ',nxd, nzd, nxg, nzg, s_lat%awald*a
!write(*,*) 'nvevEwalsR*dsqrt(2.0d0)*a, 1/s_lat%awald = '
!write(*,*) nvevEwalsR*dsqrt(2.0d0)*a, ewalda * 0.32d0

vol = omega

! set positions of all atoms
! dimensionless positions (rescaled by lattice constant)
allocate(s_lat%pos(3,nbas))
allocate(s_lat%posA(3,nbasA))

s_lat%pos = transpose(pos)/a;
s_lat%posA = transpose(posA)/a;

if (1==0) then
do il=1,nlayers
 do io=1,noctl ! 2
  ib = (il-1)*noctl + io; ! octrahedron number
  itm = (ib-1)*4 + 1 ! 4 atoms per octahedron
  ! set locations of TM:
  s_lat%pos(1:3,itm) = oct(il,io)%rb/a
  ! set locations of A-atoms:
  s_lat%posA(1:3,ib) = s_lat%pos(1:3,itm) + (/0.5d0, 0.5d0, 0.5d0 /);
  !write(*,*) 'itm, r_tm = ',itm, s_lat%pos(1:3,itm)
  !write(*,*) 'ib, r_A = ',ib,s_lat%posA(1:3,ib)

  ! set locations of O: 
  do i=1,3
   s_lat%pos(1:3,itm+i) = oct(il,io)%ro(1:3,i)/a
   !write(*,*) 'i,ro = ',i,oct(il,io)%ro(1:3,i)/a
  end do
 end do
end do
end if


!do ib=1,nbas
! write(*,'(a,i5,3f10.5)') 'atom, r = ', ib, s_lat%pos(:,ib)
!end do

!do ib=1,nbasA
! write(*,'(a,i5,3f10.5)') 'A: atom, r = ', ib, s_lat%posA(:,ib)
!end do


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
 nkd = (2*nxd+1)*(2*nxd+1)*(2*nzd+1);
 allocate(dlat(3,nkd))
 call setEwaldvecs(nxd, nzd, nkd, avec, dlat)

 nkg = (2*nxg+1)*(2*nxg+1)*(2*nzg+1);
 allocate(glat(3,nkg))
 call setEwaldvecs(nxg, nzg,nkg,bvec, glat)








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

 !CFM = 0.0d0
 !write(*,*)'testing: setting CFM=M_{l1,l2,l}=0'


return
end subroutine setmadvar

!======================================================================
! only once, before the SCF loop
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
integer :: i,j
	
 ! cy, gaunt coeff, struxidx and struxd

 call sylmnc(cy,16)
 call scg(9,cg,indxcg,jcg)
 call makcg9(indxcg,jcg,cg,gaunt,ak)

 open(190,file='gaunt.dat',action='write')
	do i=1,9
	 do j=1,9
		write(190,*) gaunt(i,j,1:25)
	 enddo
	enddo
 close(190)

 
 !write(*,*)'Gaunt calculated.... '
 !write(*,'(a,100f10.3)') 'cg(1:100) = ', cg(1:100)

 allocate(struxd(nlmi,nlmi,nbas,nbas))
 allocate(struxdA(nlmi,nbas,nbasA))

 !C   --- Structure constants for Ewald sums ---
 !nlmq = nlm; nlmq1 = nlm;
 call mkstrxd()
 call mkstrxdA() ! makes struxdA for A-sites monopoles

 !write(*,'(a)')'Structure matrix:'
 !write(*,'(a,1000e12.3)') 'struxd(:,:,1,1) = ',struxd(:,:,1,1)
 !write(*,'(a,1000e12.3)') 'struxd(:,:,2,1) = ',struxd(:,:,2,1)
 !write(*,'(a,1000e12.3)') 'struxd(:,:,8,1) = ',struxd(:,:,8,1)

 !allocate(rhoc(nl**4*nbas))
 !rhoc = 0.0_8

return
end 	subroutine initmadelung

!===============================================
subroutine setEwaldvecs(nx, nz, nn, c, v)
implicit none
integer, intent(in):: nx, nz, nn
double precision, dimension(3,3), intent(in) :: c
double precision, dimension(3,nn), intent(out) :: v
integer :: i,j,k, ind
double precision, dimension(3) :: px, py

ind =0;
do i=-nx, nx
 px = i*c(:,1)
 do j=-nx, nx
  py = j*c(:,2)
  do k=-nz, nz
    ind = ind + 1
    v(:,ind) = px + py + k*c(:,3)
  end do
 end do
end do

! swap places of first vector with the one with zero length (at ind given below).
px = v(:,1)
! index of zero vector:
ind = nx*(2*nx+1)*(2*nz+1) +  nx*(2*nz+1) + (nz + 1);
v(:,1) = v(:,ind);
v(:,ind) = px;

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
double complex :: dmij, qij
!double precision 
integer :: norbold

 !write(*,*)'getatomic: testing... '
	
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
 	dmx = 0.0d0;

	atm(ia)%qs = 0.0d0; ! setting TM q=0 means no s electrons on TM
!	if(.false.) then
! 	if(is==0) then ! Oxygen atoms
!	 atm(ia)%qs = 0.0d0;
!	else ! TM atoms: +2 monopoles for s-orbitals electrons that can go to Oxygen or remain on the parent TM atom.
!	 if(nspin==1) then
!	  atm(ia)%qs(1) = -2.0d0; ! spinless case, for two s electrons
!	 else
!	  atm(ia)%qs(:) = -1.0d0; ! spin polarised case, for one e of either spin belonging to s orbital
!	 endif
!	endif
!	endif
	
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
     !qs(ispin) = qs(ispin) + dble(qij)
     
     do i=1,norb
      dmij = 0.0d0;
      do ispin=1,nspin
       qij = conjg(wf(j,ispin)) * wf(i,ispin);
       dmij = dmij + qij
       if(i==j) qs(ispin) = qs(ispin) + dble(qij)
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

  !write(*,'(a,3x,i5,3x,10000f8.2)')'ia, qs', ia, atm(ia)%qs	   

 end do ! ia

 deallocate(dmx)


return
end subroutine getatomic
!===============================================

end module estatic
