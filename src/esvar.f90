module esvar
implicit none

integer, parameter :: dp = 8

integer :: indxcg(7400),jcg(62200)
double precision :: cg(62200)
double precision :: cy(289),gaunt(9,9,25),ak(9,9,9,9,3)
!integer, allocatable, dimension(:,:) :: struxidx
double precision, allocatable, dimension(:,:,:,:) :: struxd

!nbas,nsp,nl,nlmq,s_pot%qnu,s_ctrl%ipc,lmxl,gaunt,qpol,rho,rhoc,qmpol,mmom
integer :: nbas, nsp, nl, struxsize, nclass


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


! do we need gatoms here or would it suffice to have it in modmain??
	type :: gatoms
	 integer :: ia ! atom index in full list of atoms.
	 integer :: is ! TM species index
	 integer :: it ! type: 1 for O, 2 for TM
	 !double precision :: q0 ! number of electron in neutral configuration
	 double precision, dimension(2) ::	 qs ! spin resolved q
	 double precision, dimension(3) :: r ! position
	 double precision, allocatable, dimension(:,:):: rhoc !  dm, summed over spin.
	 !allocate( rhoc( ilm12(1,it) : ilm12(2,it) ) ) 
	 double precision, dimension(3):: mag ! magnetisation 
	end type gatoms
	type(gatoms), allocatable, dimension(:) :: atm ! atoms, general.



!r --- str_lat ---
!- Lattice and parameters related to structure
! ----------------------------------------------------------------
!r  Element  Purpose
!r   ag      pointer to symmetry group translations
!r   alat    lattice parameter, in a.u.
!r   as      dimensionless Ewald smoothing parameter
!r   avw     average MT radius, in a.u.
!r   awald   Ewald smoothing parameter (lattic)
!r   bgv     phase factor sum for symmetrization of mesh rho (supot)
!r   cg      pointer to Clebsch Gordan coeffs (setcg)
!r   cy      pointer to Ylm normalization constants (setcg)
!r   dist    deformation parameters (lattic)
!r   dlv     pointer to direct lattice vector (lattic)
!r   gam     lattice shear parms: gam, gx,gy,gz
!r   gmax    cutoff gmax for Fourier transforms
!r   gv      pointer to list of F.T. G vectors (supot)
!r   gvq     pointer to q-dependent G-vectors (sugvec)
!r   indxcg  pointer to Clebsch Gordan indxcg (setcg)
!r   igv     pointer to q-dependent G-vectors, in units of qlat (sugvec)
!r   igv2    analog of igv, but transpose, used for APWs (sugvec)
!r   ips0    pointer to first vec in star (for symm mesh rho) (supot)
!r   istab   pointer to site permutations table for group ops
!r   jcg     pointer to Clebsch Gordan jcg (setcg)
!r   kv      pointer to indices in list of F.T. G vectors (supot->gvlst)
!r   kv2     Analog of kv, for APWs (sugvec)
!r   ldist   switch specifying what kind of dist (lattdf)
!r   lsym    0 no special symmetry considerations
!r           bits 0,1:
!r           1 Make symops including SOC with spin quantized along z
!r           2 Allow operations that take z to -z
!r           bit 2
!r           4 When symmetrizing over k in the noncollinear case,
!r             do not include spinor rotation
!r   lmxst   L=cutoff when setting tolerance in Ewald sums.  If 0 => use internal default
!r   nabc    no. divisions for F.T. mesh
!r   ng      no. G vectors
!r   napw    no. G vectors for APW
!r   ngq     -
!r   nkd     no. direct latt. vecs. for Ewald sum (lattic)
!r   nkdmx   dimensioning for arrays holding latt. vecs
!r   nkq     no. reciprocal latt. vecs. for Ewald sum (lattic)
!r   nkqmx   dimensioning for arrays holding latt. vecs
!r   npgrp   Number of point symmetry group operations
!r   nsgrp   Number of space symmetry group operations
!r   nsafm   Column in (symgr, ag, istab) for special AFM symmetry, if any
!r            0 => no special AFM symmetry
!r           >0 => special AFM symmetry, but both spins calculated
!r           <0 => special AFM symmetry, spin 2 derived from spin 1
!r   plat    lattice vectors, units of alat (lattic)
!r   plat0   lattice vectors before distortion (lattic)
!r   plat2   secondary lattice vecs used in various contexts
!r   plate   order-N
!r   platl   pgf (lattic)
!r   platr   pgf (lattic)
!r   pos     pointer to site positions (susite)
!r   qlat    reciprocal lattice vectors, units 2pi/a (lattic)
!r   qlv     pointer to Ewald reciprocal lattice vectors (lattic)
!r   rpad    truncate Ewald to rpad*rmax when lattice vector
!r           list has to be padded in order to include at
!r           least one lattice vector
!r   slat    superlattice vectors
!r   symgr   pointer to symmetry group rotation matrices
!r   tol     Ewald tolerance
!r   tolft   FT mesh tolerance
!r   vol     cell volume
! ----------------------------------------------------------------
      type str_lat
      sequence

      integer    ::  ldist
      integer    ::  lmxst = 0
      integer    ::  lsym
      integer    ::  nabc  (3)
      integer    ::  ng
      integer    ::  napw
      integer    ::  ngq
      integer    ::  ngvcc ! number of lattice vectors igvcc
      integer    ::  nkd
      integer    ::  nkdmx
      integer    ::  nkq
      integer    ::  nkqmx
      integer    ::  npgrp
      integer    ::  nsgrp
      integer    ::  nsafm
      integer    ::  inull ! ensure 8-byte word boundary

!     real(8)    ::  afmt  (3)
      real(8)    ::  alat
      real(8)    ::  as
      real(8)    ::  avw
      real(8)    ::  awald
      real(8)    ::  dist  (3,3)
      real(8)    ::  gam   (4)
      real(8)    ::  gmax

      real(8)    ::  plat    (3,3)
      real(8)    ::  plat0   (3,3)
      real(8)    ::  plat2   (3,3)
      real(8)    ::  plate   (3,3)
      real(8)    ::  platl   (3,3)
      real(8)    ::  platr   (3,3)
      real(8)    ::  qlat    (3,3)
      real(8)    ::  rpad
      real(8)    ::  slat    (3,3)
      real(8)    ::  tol
      real(8)    ::  tolft
      real(8)    ::  vol

! ... Integer pointer arrays
      integer, pointer ::  indxcg(:)
      integer, pointer ::  ips0(:)
      integer, pointer ::  istab(:,:)
      integer, pointer ::  jcg(:)
      integer, pointer ::  kv(:,:)
      integer, pointer ::  igv(:,:)
      integer, pointer ::  igv2(:,:)
      integer, pointer ::  igvcc(:,:)
      integer, pointer ::  kv2(:,:)

! ... Real pointer arrays
      real(8), pointer::   ag(:,:)
      complex(8),pointer:: bgv(:)
      real(8), pointer ::  cg(:)
      real(8), pointer ::  cy(:)
      real(8), pointer ::  dlv(:,:)
      real(8), pointer ::  gv(:,:)
      real(8), pointer ::  pos(:,:)
      real(8), pointer ::  gvq(:,:)
      real(8), pointer ::  qlv(:,:)
      real(8), pointer ::  symgr(:,:)

      !type(str_symops), pointer ::  s_sym(:)

      end type str_lat







type(str_lat)::   s_lat
!type(str_ctrl)::  s_ctrl

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
      CFM(1,1,2,i) = qpol(2,i)* fac(2)
      CFM(2,0,2,i) = qpol(3,i)* fac(2)
      CFM(0,2,2,i) = qpol(3,i)* fac(2)
      CFM(2,1,1,i) = qpol(4,i)* fac(1)
      CFM(1,2,1,i) = qpol(4,i)* fac(1)
      CFM(2,2,2,i) = qpol(5,i)* fac(2)
      CFM(1,2,3,i) = qpol(6,i)* fac(3)
      CFM(2,1,3,i) = qpol(6,i)* fac(3)
      CFM(2,2,4,i) = qpol(7,i)* fac(4)
 end do ! i
 
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
