module esvar
implicit none

integer, parameter :: dp = 8

integer :: indxcg(7400),jcg(62200)
double precision :: cg(62200)
double precision :: cy(289),gaunt(9,9,25),ak(9,9,9,9,3)
!integer, allocatable, dimension(:,:) :: struxidx
double precision, allocatable, dimension(:,:,:,:) :: struxd
! rhoc complex or real? check
double complex, allocatable, dimension(:) :: rhoc 

!nbas,nsp,nl,nlmq,s_pot%qnu,s_ctrl%ipc,lmxl,gaunt,qpol,rho,rhoc,qmpol,mmom
integer :: nbas, nsp, nl, struxsize


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



end module esvar
