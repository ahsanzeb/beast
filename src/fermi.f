! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
	module fermi
	use modmain
	
	implicit none

	public :: fermid, stepf

	CONTAINS

	subroutine fermid( NK, WK, NE, E, temp, qtot, WKE, EF, nspin)

C *********************************************************************
C Finds the Fermi energy and the occupation weights of states.
C Written by J.M.Soler. August'96.
C Simple single excitation introduced by E. Artacho August 2002.
C Alternative occupation functions introduced by P. Ordejon, May'03.
C ********** INPUT ****************************************************
C INTEGER nspin    : Number of different spin polarizations (1 or 2)
C INTEGER maxspn   : Maximum number of different spin polarizations (1 or 2)
C                    for E and WKE matrices dimensions
C INTEGER NK       : Number of K-points
C REAL*8  WK(NK)   : Sampling weights of k-points (must sum 1)
C INTEGER maxe     : First dimension of E and WKE
C INTEGER NE       : Number of bands
C REAL*8  E(maxe,maxspn,NK) : State eigenvalues
C REAL*8  temp     : Temperature (in the same units of E)
C REAL*8  qtot     : Total valence charge (number of electrons)
C ********** OUTPUT ***************************************************
C REAL*8  WKE(maxe,maxspn,NK) : Occupations multiplied by k-point weights
C                               (sum qtot)
C REAL*8  EF                 : Fermi energy
C *********************************************************************

C
C  Modules
C

C Passed variables
	integer, intent(in)  :: nk, ne, nspin
	double precision, intent(in) :: wk(nk), e(nk,ne), temp, qtot
	double precision, intent(out) :: wke(nk,ne), ef

C Local variables
	integer        :: ie
	integer        :: ief
	integer        :: ik
	integer        :: ispin
	integer        :: iter
	integer,  save :: nitmax = 150
	logical,  save :: blread = .false.
	double precision, save :: tol = 1.0d-10
	double precision:: sumq, emin, emax, t, drange, wkebuf, w, eik
	double precision :: tbnspn, wknorm


C Zero occupancies, including those not explicitly
C calculated here if ne < maxe
	wke(1:nk,1:ne) = 0.0d0

	tbnspn = 2.0d0/dble(nspin);
	wknorm = 1.0d0/sum(wk(:));
	tbnspn = tbnspn * wknorm;
	
C Determine Fermi level
	sumq = 0.0d0
	emin = e(1,1)
	emax = e(1,1)
	do ik = 1,nk
	 do ie = 1,ne
	  wke(ik,ie) = wk(ik)*tbnspn
	  sumq = sumq + wke(ik,ie)
	  emin = min(emin,e(ik,ie))
	  emax = max(emax,e(ik,ie))
	 enddo
	enddo

	ef = emax

	if (sumq.lt.qtot) then
	  write(6,*) 'Fermid: Not enough states'
	  write(6,*) 'Fermid: qtot,sumq=',qtot,sumq
	endif

	T = max(temp,1.d-6)
	drange = T*dsqrt(-dlog(tol*0.01d0))
	emin = emin - drange
	emax = emax + drange
	
	do iter = 1,nitmax
	
	 ef = 0.5d0*(emin + emax)
	 sumq = 0.0d0
	 do ik = 1,nk
	  do ie = 1,ne
	   wke(ik,ie) = wk(ik)*stepf((e(ik,ie)-ef)/T)*tbnspn
	   sumq = sumq + wke(ik,ie)
	  enddo
	 enddo

C If the Fermi level was found..................... 
	 if (dabs(sumq-qtot).lt.tol) then
	  return
	 endif

	 if (sumq.le.qtot) emin = ef
	 if (sumq.ge.qtot) emax = ef
	 
	enddo

	write(6,*) 'Fermid: Iteration has not converged.'
	write(6,*) 'Fermid: qtot,sumq=',qtot,sumq
	stop

	end subroutine fermid

!---------------------------------------------------------
	double precision function stepf(X)
	double precision, intent(in) :: x
C Fermi-Dirac distribution
	if (x.gt.100.D0) then
	 stepf = 0.D0
	elseif (x.lt.-100.d0) then
	 stepf = 1.d0
	else
	 stepf = 1.d0 / ( 1.d0 + dexp(x) )
	endif

	end function stepf
!---------------------------------------------------------

! makes kgrid and weight of k points in the irreducible wedge of BZ of tetragonal lattice.
! n1 point along x/y; n3 points along z.
	subroutine mkkgrid(n1,n3)
	implicit none
	integer, intent(in) :: n1,n3
	! local
	integer :: i,j,k, ind, nkslice, i1,i2
	double precision, dimension(3) :: k1,k2,k3,b1,b2,b3

	ntotk = (n1*(n1+1))/2 *n3; ! in irreducible wedge of the BZ of tetragonal lattice
	
	allocate(kgrid(3,ntotk))
	allocate(wk(ntotk))
	b1 = 0.5d0*bvec(:,1)/n1
	b2 = 0.5d0*bvec(:,2)/n1
	b3 = 0.5d0*bvec(:,3)/n3

! make a general triangular slice
	ind = 0;
	do i=0,n1-1
	 k1 = i*b1
	 do j=i,n1-1
	  k2 = j*b2
	   ind = ind + 1;
	   kgrid(:,ind) = k1 + k2; ! cartesian comp. 
	   if(i == 0) then
	    if(j == 0) then
	     wk(ind) = 0.125d0 ! 1/8
	    elseif(j == n1-1) then
	     wk(ind) = 0.25d0 ! 1/4 ! (1,1) corner 
	    else
	     wk(ind) = 0.5d0 ! 1/2 ! on x=0 line
	    endif
	   elseif(i == j .or. j== n1-1) then ! on diagonal (x=y) or right side (y=1); 
	   ! i==j==n1-1 reset to correct value 1/8 below after i,j loops
	    wk(ind) = 0.5d0 ! 1/2
	   else ! deep inside the wedge
	    wk(ind) = 1.0d0 ! 1
	   endif
	 end do
	end do
	! last point is at i=j= n1-1: change weight to 1/8
	wk(ind) = 0.125d0 ! 1/8

	nkslice = (n1*(n1+1))/2;

	! z> 0, z<1 slices:
	! use the slice to make slices at z>0
	do k=1,n1-1
	  k3 = k*b3
	  i1 = k*nkslice + 1;
	  i2 = (k+1)*nkslice;
	  do i=1,3
	   kgrid(i,i1:i2) = kgrid(i,1:nkslice) + k3(i);
	 end do
	 wk(i1:i2) = wk(1:nkslice)
	end do

	! z=0 slice; ! kgrid same; wk half of the above
	wk(1:nkslice) = 0.5d0 * wk(1:nkslice)
	! z=1 slice; 
	wk(i1:i2) = 0.5d0 * wk(1:nkslice) ! i1,i2 set to correct value in the last iter of k loop above

	!write(*,'(10000f10.3)') wk

	return
	end subroutine mkkgrid
C *********************************************************************

	end module fermi
