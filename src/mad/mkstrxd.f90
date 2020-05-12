module mmkstrxd
use esvar
use mrcnsl0
use mhstra
implicit none
contains
	

subroutine mkstrxd() !s_ctrl, ipc, s_lat, tbc, nlmq, nlmq1, lmxl, ldip, dlat, nkd, glat, nkg, &
                  !& indxcg, jcg, cg, cy, struxd, struxidx,  struxsize)
!- Make lookup table of strux and radial derivatives
! ----------------------------------------------------------------------
!i Inputs:
!i  ldip  : 3 include 'spherical' dipole correction to Ewald
!i        : 2 include 'slab' dipole correction to Ewald
!i        : any other number - skip the dipole correction
!i   nbas,bas,awld,alat,vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy
!i   nlmq1: leading dimension of strx, nlmq1=(ll(nlmq)+1)**2
!i   nlmq : max L-cutoff for multipoles, second dimension of strx
!i   lmxl : l-cutoff for multipoles for each class of species
!o Outputs:
!o  struxd : coefficients B_LL'(R'-R) (structure constants)
!r Remarks
!r   Instead of making B on each sc iteration, the program prepares the whole B
!r   in the beginning of the self-consistency cycle to be used as a lookup table.
!r   This results in faster calculation at an expence of memory.
!r
!r   Calling mkstrx is set as the default option. To call instead hstr on each
!r   iteration (and therefore to reduce memory but increase computational time),
!r   use switch --sfly
!r
!r   Efficiency issues: there are two symmetry properties of zero energy
!r   structure constants that are used in mkstrx:
!r     B_LL'(R'-R) = B_L'L(R-R')                (1)
!r     B_LL'(R'-R) = (-1)^(l+l') B_L'L(R'-R)    (2)
!r   same properties hold for the radial derivative of B.
!r
!r This symmetry however is not employed because the complexity that a
!r  parallel implementation will require and the diminishing gains + increased overhead with.
!r Curently the matrices and the indices are distributed in contiguous blocks across the long
!r side and no communication whatsoever is needed. For a process with id pid, the block starts
!r at atom tbc%esamap(pid)+1 and ends at tbc%esamap(pid+1) inclusive. For atom ib from this range
!r there is a column of neighbours in struxidx(:,ib-tbc%esamap(pid)). For each of pair (jb,ib)
!r a block of size (nlmj,nlmi1) is preallocated starting at struxd(struxidx(jb,ib-tbc%esamap(pid)))
!r The block size for dstrx is (nlmj,nlmi)
!r
!b Bugs
!b   mkstrx is not written in the most efficient way and could be further refined wrt to
!b   amount of calculated elements of strux. At the moment mkstrx is a result of certain
!b   trade off between performance and clarity of the program.
!b
!u Updates
!u           2013 (DP)  strx moved to indexed linear array with auxiliary index for better compatibility with C
!u           2013 (DP)  Fortran pointer based strx for memory efficiency
!u           2012 (DP)  Parallell distributed strx & dstrx
!u    05 Mar 2010 (SL)  optimization (make strux up to Lmax and for triangle
!u                      jb<=ib only)
!u    19 Jan 2010 (ATP) first written
! ----------------------------------------------------------------------

   !use tbprl
   !use structures
   implicit none

   integer :: pib,jb,lmax,lmxf,i1mach,iprint,ib
   integer :: li,li1,lj, ilm,jlm, i, cmol, cpv, nbas1, u, pi, sz
   real(dp) :: tau(3),taua(3), hl(100),bl(100), plat(3,3), qlat(3,3), alat, &
                              & awald, det   
      
   plat  = s_lat%plat
   alat  = s_lat%alat
   vol   = s_lat%vol
   awald = s_lat%awald

   !call dinv33(s_lat%plat,0,qlat,det)
   qlat = s_lat%qlat ! maz
   det = s_lat%det


   do  ib = 1, nbas
      do  jb = 1, nbas ! can we restrict jb to ib:nbas ? 
         tau = s_lat%pos(1:3,jb)-s_lat%pos(1:3,ib)
         call directshortn(tau,plat,qlat)
         call rcnsl0(tau, awald, alat, hl) ! lmxst,nlm,alat,glat,nkg,dlat,nkd,vol,cy(1:nlm)
         call hstra(struxd(:,:,jb,ib), hl) ! hstra is efficient: uses symmetry of B
      enddo
   end do

return
end subroutine mkstrxd



subroutine directshortn(p,plat,qlat)
!       Shorten vector 'p' to cell 'plat' using precomputed ilat = plat^-1
!       'plat' vectors shall be stored in the columns.
         implicit none
         real(dp), intent(in) :: plat(3,3), qlat(3,3)
         real(dp), intent(inout) :: p(3)
         real(dp) :: r(3)

         r = matmul(qlat, p)
         ! Shorten r to a vector within a cell described by the identity matrix I.
         ! Preserve direction on boundary cases (inclusive boundary).
         r = r + sign(1.0d0, r)*floor(0.5d0 - abs(r))
         p = matmul(plat, r)
return
end subroutine directshortn 

end 	module mmkstrxd   
