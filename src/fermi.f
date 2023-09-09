! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
	module fermi
	!use modmain
	
	implicit none

	public :: fermid, stepf

	CONTAINS

	subroutine fermid(NK,WK,wknorm,NE,E,temp, 
     .                    qtot, WKE, EF, nspin, ebands)

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
	double precision, intent(out) :: wke(nk,ne), ef, ebands
	double precision, intent(in) :: wknorm
C Local variables
	integer        :: ie
	integer        :: ief
	integer        :: ik
	integer        :: ispin
	integer        :: iter
	integer,  save :: nitmax = 200
	logical,  save :: blread = .false.
	double precision, save :: tol = 1.0d-8
	double precision:: sumq, emin, emax, t, drange, wkebuf, w, eik
	double precision :: tbnspn

C Zero occupancies, including those not explicitly
C calculated here if ne < maxe
	wke(1:nk,1:ne) = 0.0d0
	ebands = 0.0d0

	tbnspn = 2.0d0/dble(nspin);
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
	  !calc energy of occupied states
	  do ik=1,nk
	   do ie=1,ne
	    ebands = ebands + wke(ik,ie)*e(ik,ie)
	   end do
	  end do
	  
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

	end module fermi
