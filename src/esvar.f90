module esvar
implicit none

integer, parameter :: dp = 8

integer :: indxcg(7400),jcg(62200)
double precision :: cg(62200)
double precision :: cy(289),gaunt(9,9,25),ak(9,9,9,9,3)
!integer, allocatable, dimension(:,:) :: struxidx
double precision, allocatable, dimension(:,:,:,:) :: struxd
double precision, allocatable, dimension(:,:,:) :: struxdA ! for A-sites atoms
double precision :: decimal ! to round the struxd, struxdA

double precision, dimension(1:25,1:25,4) :: Rlmax ! up to l=4; rotation matrix for spherical harmonics; for all 4 B atoms
logical :: cage
double precision :: qaa, fVoctO, fVoctAB


!nbas,nsp,nl,nlmq,s_pot%qnu,s_ctrl%ipc,lmxl,gaunt,qpol,rho,rhoc,qmpol,mmom
integer :: nbas, nsp, nclass
integer :: nbasA ! A-sites atoms = total number of octahedra.

integer :: nxd, nzd
integer :: nxg, nzg

integer :: nkg, nkd
real(8), allocatable, dimension(:,:) :: glat, dlat

integer :: lmxl, nlmi ! lmxl = max lm for potential, same for all species.
integer :: lmxst, nlm ! max lm for structure matrix; lmxst = 2*lmxl
integer, allocatable, dimension(:) :: ll

!logical :: ldip ! dipole correction?
integer :: ldip

real(8) :: vol

! atmomic species related:
real(8), allocatable, dimension(:,:,:,:) :: CFM ! species: crystal field matrix: only a few elements are nonzero... 
real(8), allocatable, dimension(:) :: q0  ! species: number of electrons when neutral
integer, dimension(2,2) :: ilm12
! Oxygen: only p-orbitals: index ilm: 2-4
!ilm12(1,0) = 2; ilm12(2,0) = 4;
! TM: only d-orbitals: index ilm: 5-9
!ilm12(1,1) = 5; ilm12(2,1) = 9;
integer, dimension(2) :: types2norb
integer, allocatable, dimension(:) :: species2type

! qmpol does not have to be a part of gatoms object, as it has the same dimension for all atoms.
double precision, allocatable, dimension(:,:) :: qmpol ! size=(nlmi,nbas)
double precision, allocatable, dimension(:,:) :: qpol ! Crstal field constant Delta_{l',l'',l} in Paxton's notes
double precision :: qmpolA ! monopole charge on A-site.
double precision, allocatable, dimension(:) :: qref ! size=(nbas)

double precision, allocatable, dimension(:):: hard ! hardness of each atom, p orb of Oxygen at position 0, d orb of TM at 1:nsptm


! do we need gatoms here or would it suffice to have it in modmain??
	type :: gatoms
	 integer :: ia ! atom index in full list of atoms.
	 integer :: is ! TM species index
	 integer :: it ! type: 1 for O, 2 for TM
	 !double precision :: q0 ! number of electron in neutral configuration
	 double precision, dimension(2) ::	 qs ! spin resolved q
	 double precision, dimension(3) :: r ! position
	 double complex, allocatable, dimension(:,:):: rhoc
	 double precision, allocatable, dimension(:,:):: dh  	 
	 !allocate( rhoc( ilm12(1,it) : ilm12(2,it) ) ) 
	 double precision, dimension(3):: mag ! magnetisation 
	end type gatoms
	type(gatoms), allocatable, dimension(:) :: atm ! atoms, general.



!r --- str_lat ---
!- Lattice and parameters related to structure
! ----------------------------------------------------------------
!r  Element  Purpose
!r   alat    lattice parameter, in a.u.
!r   as      dimensionless Ewald smoothing parameter
!r   awald   Ewald smoothing parameter (lattic)
!r   plat    lattice vectors, units of alat (lattic)
!r   pos     pointer to site positions (susite)
!r   qlat    reciprocal lattice vectors, units 2pi/a (lattic)
!r   vol     cell volume
! ----------------------------------------------------------------
type str_lat
 real(8)    ::  alat, as, awald, vol
 real(8), dimension(3,3)::  plat, qlat
 real(8), allocatable ::  pos(:,:), posA(:,:)
end type str_lat

type(str_lat)::   s_lat


 contains

 !======================================================================
subroutine getM0(ilm,ilmp,ilmpp,qpol,M)
!C-
!C ----------------------------------------------------------------------
!Ci Inputs:
!Ci
!Co Outputs:
!Co
!Cr Remarks
!Cr  The tight binding parameters in qpol are as follows
!Cr  qpol(1) = M(011) = M(101)
!Cr  qpol(2) = M(112)
!Cr  qpol(3) = M(022) = M(202)
!Cr  qpol(4) = M(121) = M(211)
!Cr  qpol(5) = M(222)
!Cr  qpol(6) = M(123)
!Cr  qpol(7) = M(224)
!Cr  These are converted into Stone's definitions by multiplying the
!Cr  values from ctrl by \sqrt{4\pi/(2l+1)}
!C ----------------------------------------------------------------------
      implicit none
!C Passed Parameters
      integer ilm,ilmp,ilmpp
      double precision qpol(10),M
!C Local Variables
      integer l,lp,lpp
      double precision fourpi,sqr4pi,fac,dsqrt,datan
      fourpi = 16d0*datan(1d0)
      sqr4pi = dsqrt(fourpi)

      l = ll(ilm)
      lp = ll(ilmp)
      lpp = ll(ilmpp)
      fac = dsqrt(fourpi/(2*l + 1))
      M = 0d0
      if (lp == 1 .and. lpp == 0 .and. l == 1) M = qpol(1)
      if (lp == 0 .and. lpp == 1 .and. l == 1) M = qpol(1)
      if (lp == 1 .and. lpp == 1 .and. l == 2) M = qpol(2)
      if (lp == 2 .and. lpp == 0 .and. l == 2) M = qpol(3)
      if (lp == 0 .and. lpp == 2 .and. l == 2) M = qpol(3)
      if (lp == 2 .and. lpp == 1 .and. l == 1) M = qpol(4)
      if (lp == 1 .and. lpp == 2 .and. l == 1) M = qpol(4)
      if (lp == 2 .and. lpp == 2 .and. l == 2) M = qpol(5)
      if (lp == 1 .and. lpp == 2 .and. l == 3) M = qpol(6)
      if (lp == 2 .and. lpp == 1 .and. l == 3) M = qpol(6)
      if (lp == 2 .and. lpp == 2 .and. l == 4) M = qpol(7)
      M = fac*M
      if (l == 0) then
        if (lp == lpp) then
          M = sqr4pi
        endif
      endif
      end subroutine getM0
!======================================================================
      subroutine getM(ilm,ilmp,ilmpp,ispec,qpol,M)
!C-
!C ----------------------------------------------------------------------
!Ci Inputs:
!Ci
!Co Outputs:
!Co
!Cr Remarks
!Cr  The tight binding parameters in qpol are as follows
!Cr  qpol(1) = M(011) = M(101)
!Cr  qpol(2) = M(112)
!Cr  qpol(3) = M(022) = M(202)
!Cr  qpol(4) = M(121) = M(211)
!Cr  qpol(5) = M(222)
!Cr  qpol(6) = M(123)
!Cr  qpol(7) = M(224)
!Cr  These are converted into Stone's definitions by multiplying the
!Cr  values from ctrl by \sqrt{4\pi/(2l+1)}
!C ----------------------------------------------------------------------
      implicit none
!C Passed Parameters
      integer, intent(in) :: ilm,ilmp,ilmpp, ispec
      double precision, intent(in) ::  qpol(10)
      double precision, intent(out) ::  M
!C Local Variables
      integer l,lp,lpp

      l = ll(ilm)
      lp = ll(ilmp)
      lpp = ll(ilmpp)

      M = CFM(lpp,lp,l,ispec)
      
      return
      end subroutine getM
!======================================================================
! set crystal field for all species: TMs & O
!======================================================================
      subroutine setCFM()

! Ahsan Jun 2023: sqrt{4pi/(2l+1)} factors should be divided, not multiplied.
! calibrated by calculating madelung const of NaCl and CsCl structures
! by setting appropriate charges on TM/A atoms and comparing with converged values from litrature. 
! bug reported to TBE team:
! https://bitbucket.org/lmto/lm/issues/139/bug-in-multipole-term-of-tb-hamiltonian

!C
!C ----------------------------------------------------------------------
!Ci Inputs:
!Ci
!Co Outputs:
!Co
!Cr Remarks
!Cr  The tight binding parameters in qpol are as follows
!Cr  qpol(1) = M(011) = M(101)
!Cr  qpol(2) = M(112)
!Cr  qpol(3) = M(022) = M(202)
!Cr  qpol(4) = M(121) = M(211)
!Cr  qpol(5) = M(222)
!Cr  qpol(6) = M(123)
!Cr  qpol(7) = M(224)
!Cr  These are converted into Stone's definitions by multiplying the
!Cr  values from ctrl by \sqrt{4\pi/(2l+1)}
!C ----------------------------------------------------------------------
 implicit none
 ! local
 integer ll,l,lp,lpp, i
 double precision fourpi,sqr4pi
 double precision, dimension(4) :: fac

 integer, parameter :: lmxl=2

 fourpi = 16d0*datan(1d0)

 sqr4pi = dsqrt(fourpi)
 do l=1,4
  fac(l)= dsqrt(fourpi/dble(2*l + 1))
 end do
  
 allocate(CFM(0:lmxl,0:lmxl,0:2*lmxl,0:nclass-1)) ! O has species index of 0
 CFM(:,:,:,:) = 0.0d0;

 do i=0, nclass-1     
  do lp=0, lmxl
   CFM(lp, lp, 0,i) = sqr4pi
  end do
      CFM(1,0,1,i) = qpol(1,i)* fac(1)
      CFM(0,1,1,i) = qpol(1,i)* fac(1)
      CFM(1,1,2,i) = qpol(2,i)* fac(2) ! 
      CFM(2,0,2,i) = qpol(3,i)* fac(2)
      CFM(0,2,2,i) = qpol(3,i)* fac(2)
      CFM(2,1,1,i) = qpol(4,i)* fac(1)
      CFM(1,2,1,i) = qpol(4,i)* fac(1)
      CFM(2,2,2,i) = qpol(5,i)* fac(2)
      CFM(1,2,3,i) = qpol(6,i)* fac(3)
      CFM(2,1,3,i) = qpol(6,i)* fac(3)
      CFM(2,2,4,i) = qpol(7,i)* fac(4)
      !write(*,*)'i, CFM(2,2,2,i) = ',i, CFM(2,2,2,i)
      !write(*,*)'i, CFM(2,2,4,i) = ',i, CFM(2,2,4,i)

 end do ! i

 !write(*,*) 'esvar.f90: setCFM(),  fac : ',fac

 
 return
 end subroutine setCFM
!======================================================================

 double precision function dsum(n,dx,incx)
!c
!c     takes the sum of the values.  Adapted from:
!c     jack dongarra, linpack, 3/11/78.
!c
 double precision dx(1)
 integer i,incx,n,nincx

 dsum = 0d0
 if (n <= 0) return

 nincx = n*incx
 do i = 1, nincx, incx
  dsum = dsum + dx(i)
 end do
      
 end function dsum
!======================================================================


end module esvar
