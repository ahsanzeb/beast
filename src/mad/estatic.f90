module estatic
use esvar
implicit none

contains

!======================================================================
! We want to keep this module and its related routines seperate from 
! the rest of out code.
! This routines sets the global variables of this module
! It can be combined with initmadelung()
!======================================================================
subroutine setmadvar()
use modmain, only: oct, ntot, noctl, nlayers
implicit none
integer :: ilm, i,l, m, io, il, ib


!integer, intent(in) :: ntot,ntms,nox, maxlv
nbas = ntot;
!lmxl = maxlv;


! lattice vectors:
!s_lat%plat = avec ! check if its columns of rows? ! lattc.f: plat(*,k) holds lattice vector k
!s_lat%alat = a
!s_lat%vol = vol !?
!s_lat%awald = awald ! ?

!s_lat%qlat = bvec ! check if rows or columns? + 2pi factors?
!s_lat%det = det ! det of avec?

! set positions of all atoms
do il=1,nlayers
 do io=1,noctl ! 2
  ib = (il-1)*noctl + io; ! octrahedron number
  s_lat%pos(1:3,ib) = oct(il,1)%rb
  do i=1,3
   s_lat%pos(1:3,ib+i) = oct(il,io)%ro(1,i)
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




! some sizes:
nlmi = (lmxl+1)*(lmxl+1);
lmxst = 2* lmxl;
nlm = (lmxst+1)*(lmxst+1);

! set ll array to be used in hstrd() in mkstrxd.f90: 
allocate(ll(nlmi))
ll = 0;
ilm = 0;
do l=0,lmxl
 do m = -l,l
  ilm = ilm + 1
  ll(ilm) = l
 end do
end do

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

 allocate(struxd(nlmi,nlmi,nbas,nbas))
 !C   --- Structure constants for Ewald sums ---
 !nlmq = nlm; nlmq1 = nlm;
 call mkstrxd()

 !allocate(rhoc(nl**4*nbas))
 !rhoc = 0.0_8

return
end 	subroutine initmadelung



end module estatic
