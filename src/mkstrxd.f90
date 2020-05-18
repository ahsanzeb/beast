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
   real(dp) :: tau(3), hl(nlm), awald, taux(3)
   real(dp) ::  plat(3,3), qlat(3,3), alat

   plat  = s_lat%plat
   qlat = s_lat%qlat
   alat  = s_lat%alat
   awald = s_lat%awald

	write(*,'(a,10000f10.4)') 'plat, qlat, vol, awald= ',plat, qlat, vol, awald
  write(*,'(a,1000i5)') 'nkg, nkd = ', nkg, nkd
  
   do  ib = 1, nbas
      do  jb = 1, nbas ! can we restrict jb to ib:nbas ? 
         !write(*,'(a,2i5)') 'ib, jb = ',ib,jb
         tau = s_lat%pos(1:3,jb)-s_lat%pos(1:3,ib);
         !taux = tau;
         !write(*,'(a,10000f10.2)') 'tau = ',tau 
         call directshortn(tau,plat,qlat)
         !write(*,'(a,3f7.3,3x,3f7.3)')'tau_in, tau_out = ',taux, tau
         call rcnsl0(alat*tau, awald, hl) ! lmxst,nlm,alat,glat,nkg,dlat,nkd,vol,cy(1:nlm)

          !write(*,'(a,2i5,1000e10.1)') 'ib,jb, hl = ',ib,jb, hl

         !if(ib==1 .and. jb == 2) then
           !write(*,'(a,100000f10.5)') 'hl = ',hl
           !write(*,'(a,100000f10.5)') 'cy = ',cy         
         !endif

         ! struxd(:,:,jb,ib): ~ Ylm of ib at jb (second dim) & Ylm of jb at ib (first dim); 
         call hstra(struxd(:,:,jb,ib), hl) ! hstra is efficient: uses symmetry of B  
               
         !do ilm=1,nlmi
         ! write(*,'(a,i5,1000e10.1)') 'ilm, struxd(:,ilm,jb,ib) = ',ilm, struxd(:,ilm,jb,ib)
         !end do

      enddo
   end do

! write(*,'(a)')'..... .... .... ....... ..... .... .... '
! write(*,'(a)')'mkstrucd.f90: '
! do ib=1,nbas
!  write(*,'(a,100f10.5)') (norm2(struxd(:,:,ib,jb) - &
!                       transpose(struxd(:,:,jb,ib))), jb=1,nbas)
! end do
! write(*,'(a)')'..... .... .... ....... ..... .... .... '

 if(1==0) then
 write(*,'(a)')'..... .... .... ....... ..... .... .... '
 write(*,'(a)')'mkstrucd.f90: '
 do ib=1,nbas
  write(*,'(100f15.5)') (norm2(struxd(:,:,ib,jb)), jb=1,nbas)
 end do
 write(*,'(a)')'..... .... .... ....... ..... .... .... '
 endif
 

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



subroutine directshortnx(p,plat,qlat)
!       Shorten vector 'p' to cell 'plat' using precomputed ilat = plat^-1
!       'plat' vectors shall be stored in the columns.
         implicit none
         real(dp), intent(in) :: plat(3,3), qlat(3,3)
         real(dp), intent(inout) :: p(3)
         real(dp) :: r(3)
         write(*,*)'-------------------------'
	       !write(*,*) 'p= ',p
	       write(*,*) p
         r = matmul(qlat, p)
	       !write(*,*) 'r=matmul(qlat, p):'
	       write(*,*) r
         ! Shorten r to a vector within a cell described by the identity matrix I.
         ! Preserve direction on boundary cases (inclusive boundary).
         r = r + sign(1.0d0, r)*floor(0.5d0 - abs(r))
	       !write(*,*) 'r=r + sign(1.0d0, r)*floor(0.5d0 - abs(r)):'
	       write(*,*) r

         p = matmul(plat, r)
	       !write(*,*) 'r=matmul(plat, p):'
	       write(*,*) p

return
end subroutine directshortnx 






!==============================================================================
! strux involving A sites: asymmetric: A only as jb atom, other TM/O only as ib atom
!==============================================================================
subroutine mkstrxdA() 
implicit none

   integer :: pib,jb,lmax,lmxf,i1mach,iprint,ib
   integer :: li,li1,lj, ilm,jlm, i, cmol, cpv, nbas1, u, pi, sz
   real(dp) :: tau(3), hl(nlm), awald ,taux(3)
   real(dp) ::  plat(3,3), qlat(3,3), alat
      
   plat  = s_lat%plat
   qlat = s_lat%qlat
   alat  = s_lat%alat
   awald = s_lat%awald

	write(*,'(a,100f10.4)')'plat = ', plat
		write(*,'(a,100f10.4)')'qlat = ', qlat
	write(*,'(a,100f10.4)')'plat.qlat = ', matmul(plat,qlat)
	
	write(*,'(a,10000f10.4)') 'plat, qlat, vol, awald= ',plat, qlat, vol, awald
  write(*,'(a,1000i5)') 'nkg, nkd = ', nkg, nkd

   ! ib= normal TM/O atoms; jb=A-site atoms; 
   ! ib,jb switched below: see tau & hstraA call.
   do  ib = 1, nbas
      do  jb = 1, nbasA ! 
         tau = s_lat%pos(1:3,ib)-s_lat%posA(1:3,jb); ! ib, jb switched: 
                                                     ! Ylm(tau): centred around A-site, at r_{ib}.
         call directshortn(tau,plat,qlat)
         !taux = s_lat%pos(1:3,ib);                                                    
         !if(jb==1)call directshortnx(taux,plat,qlat)
         if(jb==1)write(*,'(a,i5,f7.3)')'ia, |tau| = ',ib, norm2(tau)

         call rcnsl0(alat*tau, awald, hl)
         call hstraA(struxdA(:, ib,jb), hl) ! ib, jb switched because we need 1:nlmi for ib; to preserve the structure of struxdA:struxd.
      enddo
   end do

 write(*,'(a)')'..... .... .... ....... ..... .... .... '
 write(*,'(a)')'mkstrucd.f90 A: '
 do ib=1,nbas
  write(*,*) (struxdA(:,ib,jb), jb=1,nbasA)
 end do
 write(*,'(a)')'..... .... .... ....... ..... .... .... '



return
end subroutine mkstrxdA

end 	module mmkstrxd   
