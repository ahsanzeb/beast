	module mhstra
	use mtbshfl
	use esvar
	implicit none
	contains
	

!==============================================================================
! maz: copies from hstra.f:
!==============================================================================
	subroutine hstra(strx,hl)
!C- Make structure constants from reduced strux at energy zero taking
!C- advantage of symmetry properties of zero energy strux
!C ----------------------------------------------------------------------
!Ci Inputs:
!Ci  ldip  : 3 include 'spherical' dipole correction to Ewald
!Ci        : 2 include 'slab' dipole correction to Ewald
!Ci        : any other number - skip the dipole correction
!Ci  nlmf,nlm :make B_LL' for L = 1, nlmf, L'= 1, nlm
!Ci  nlmq1 : leading dimensions of B, nlmq1=(ll(nlmq)+1)**2
!Ci  nlmq  : max L-cutoff for multipoles, leading dimension of B'
!Ci  hl    : solid Hankels for given point tau
!Ci  cg,indxcg,jcg : Clebsch Gordan coefficients in condensed form
!Ci  vol   : unit cell volume
!Co Outputs:
!Co  strx  : coefficients B_LL'(tau) (structure constants)
!Co  dstrx : derivatives of structure constants (x dB_LL'/dx) at x = tau
!Co          if pv = F dstrx is not touched
!Cr Remarks:
!Cr  strx are zero energy structure constants B_LL' are calculated as:
!Cr    B_LL'(tau) = 4\pi \sum_L" (-1)^l' C_LL'L" H_l"(tau)
!Cr  where C_LL'L" are Gaunt coefficients (C_LL'L" = \int Y_L * Y_L' * Y_L")
!Cr  and H_L are solid Hankels at zero energy
!Cr    H_L(tau) = h_l"(|tau|) Y_L"(tau)
!Cr  prepared by either rcnsl0.f (periodic branch, MOL = .F.) or
!Cr  soldh.f (non-periodic branch, MOL = .T.).
!Cr
!Cr  If pv is set, returns tau*d/dtau B in drstrx
!Cr
!Cr  The program is an accelerated version of hstr.f which takes
!Cr  advantage of the symmetry of B_LL':
!Cr          B_L'L(tau) = B_LL'(tau) * (-1)^(l+l')   (*)
!Cr  (*), in turn, is a consequence of the property
!Cr  Y_L(-tau) = (-1)^l Y_L(tau)
!Cr
!Cr  hstra first calculates the lower left triangle of B_LL'(L = 1, nlmf,
!Cr  L' = 1, L) and then fills up the upper right triangle using (*).
!Cr  For the sake of clarity, the program is restricted to nlmf >= nlm.
!Cr  No such restriction is imposed in hstr.f though.
!Cr
!Cu Updates
!Cu    05 Mar 2010 (SL) created from hstr.f
!C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      !integer, intent(in) :: ldip
      double precision, intent(in) :: hl(nlm)
      double precision, intent(out) :: strx(nlmi,nlmi)
C Local Parameters
      integer :: icg,icg1,icg2,ii,ilm,indx,klm,l,lk,llm,lm,lp,mlm
      integer, parameter :: lmxx=12
      double precision sig(0:lmxst),fourpi,fpibv,sumx

!       call tcn('hstra: make strux')

      fourpi = 16.0d0*datan(1.0d0)

C --- (-1)^l ---
      sig(0) = 1d0
      if (lmxst > 0) then
        do  l = 1, lmxst
          sig(l) = -sig(l-1)
        enddo
      endif
C --- add together Gaunt-sums ---
      do  mlm = 1, nlmi
        lm = ll(mlm)
        do  klm = 1, mlm
          lk = ll(klm)
          sumx = 0.0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2+min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  icg = icg1, icg2
            llm = jcg(icg)
            lp = ll(llm)
C...        Only lp = lm + lk contribute
            if (lp /= lm+lk) cycle
            sumx = sumx + cg(icg)*hl(llm)
          enddo
          strx(mlm,klm) = sumx*fourpi*sig(lk)
        enddo
      enddo
C --- Add the remaining off-diagonal terms using B_LL' = -1^(l+l')*B_L'L
      if (nlmi > 1) then
        do  mlm = 1, nlmi-1
          lm = ll(mlm)
          do  klm = mlm+1, nlmi
            lk = ll(klm)
            strx(mlm,klm) =  (sig(lm)*sig(lk))*strx(klm,mlm)
          enddo
        enddo
      endif

C --- the following includes extra p terms 'implicitly' ---
      if (ldip == 2 .or. ldip == 3) then
        if (nlmi > 1) then
          fpibv = fourpi/vol
          do  ilm = 2, 4
            strx(ilm,ilm) = strx(ilm,ilm) - fpibv
          enddo
        endif
      endif

      call tbshfl(0,nlmi,nlmi,nlmi,strx)
      call strfac(0,nlmi,nlmi,nlmi,strx)

!       call tcx('hstra: make strux')
      end subroutine hstra
	end 	module mhstra
