 module mtbeseld
 use esvar, only: nsp, nbas, nlmi, ll, ilm12, atm, struxd, & 
                  qmpol, CFM, gaunt, nbasA, struxdA

 implicit none

 contains

 
!============================================================================
 subroutine tbeseld(ecorr)
                       
!- Make L-expanded electrostatic potential, multipole moments; force
! ----------------------------------------------------------------------
!i Inputs:   rhoc=c*_RL c_RL' + cc;
!i           qpol, 10 polarisation parameters;
!i           rhon=c*_RL S c_RL' + cc (TB+U only);
!i           drhosl=c*_RL dS_RLR'L'/dR c_R'L' + cc (ovlp only)
!i           drhos: work array to make drhosl summed over all LL'
!i           gaunt: Gaunt coefficients
!i    force: if F, skip calculation of forces (see tbzint.f)
!i    pv   : if F, skip calculation of pressure (see tbzint.f)
!i    ldip : 3 include 'spherical' dipole correction to Ewald
!i         : 2 include 'slab' dipole correction to Ewald
!i         : any other number - skip the dipole correction
!i    nlmq1: leading dimension of strx and vm (see Remarks)
!i    nlmq : global L-cutoff for multipoles, leading dimension of dstrx
!i           and qmpol, second dimension of strx
!i    lmxl : species-resolved l-cutoffs for multipoles
!i    stni : Stoner parameter (l=2)
!i    idu  : for TB+U, index to which channels get additional potential
!i    uh,uj: Hubbard U's and J's for TB+U and spin pol TB-L
!i           uh is the screened U for TB+U. Unscreened U is taken from
!i           start parameters (qnu)
!i    rho  : s, p, d Mulliken charges, dimension nl,nsp,nbas (UL only)
!i    qmpol: multipole moments on each site (from tbmpole)
!i    ak   : product of Gaunt integrals for TB+U (see makcg9)
!i    struxd, dstrxd: structure constants and their radial derivatives (input) (data arrays)
!i    struxidx, dstrxidx: pointer arrays for struxd and dstrxd,
!i                        struxd(struxidx(jb,ib-tbc%esamap(pid))) points to a (nlmj,nlmi1) dimensioned block and
!i                        dstrxd(dstrxidx(jb,ib-tbc%esamap(pid))) points to a (nlmj,nlmi) dimensioned block.
!io Input/Outputs:
!io   fstrx, dfstrx: local copy of strx, dstrx. With the --sfly switch
!io                  both arrays are calculated in the program (see remarks)
!i    mmom : magnetic moments on each site (from tbmpole)
!o Outputs:
!o    rho  : s, p, d Mulliken charges, dimension nl,nsp,nbas
!o    dh   : increment to on-site hamiltonian, \Delta H_{RLRL'},
!o    pot0 : monopole potential at site R (eq. 7.81, Finnis)
!o           ---all these in Rydberg atomic units.
!o    f    : force from e'static terms
!o    fnou : force from overlap dependence of monopoles (Hubbard part)
!o    fnom : force from overlap dependence of monopoles (Madelung part)
!o    ecorr: second order correction to total energy
!o           E_2 of Finnis' book, eq (7.75); or E^U in TB+U
!o    ppdip: accumulated dipole moment, RELAX=0 in ctrl is misused
!o           to permit the dipole to be accumulated over a subset
!o           of the sites (see tbtote also)
!o    tpvq : contribution of multipoles to 3pV (pressure)
!l Local variables:
!o    vm   : Madelung potential due to charge transfer
!o    vu   : Hartree U part of the potential
!o    vj   : Stoner J part of the potential
!r Remarks
!r    The multipoles Q_L and components of electrostatic potential V_L
!r    are defined here in terms of the conventions of Stone:
!r    "Theory of intermolecular forces." If the multipole moments of the
!r    charge at sites i are Q_L then the components of the potential
!r    at sites j are defined from
!r    V(r) = \Sum_L sqrt(4pi/(2l+1)) V_L r^l Y_L(r)
!r         = \Sum_L sqrt(4pi/(2l+1))(2l+1)!! V_L J_L(r) where J_L are
!r    Bessel functions as defined by Michael (see subroutine BESSL).
!r    V_L = \Sum_L' sqrt((2l+1)*(2l'+1))/((2l+1)!!(2l'+1)!!) B_L'L Q_L'
!r    where B_L'L are Michael's structure constants. Here the Q are
!r    Stone's, not Jackson's multipole moments. They are related by
!r    Q_L=sqrt(4pi/(2l+1)) Q_L(Jackson)
!r    Jackson's are the ones used in the FP-LMTO and NFP programs.
!r    That means that the conventional multipole moments (Jackson) are
!r    these moments multiplied by sqrt{(2l+1)/4pi}. The V_L here are
!r    defined such that the potential V(r) is expanded as
!r    V(r)=\sum_L \sqrt(4pi/(2l+1)) V_L r^l Y_L(r),
!r    therefore they must be multiplied by \sqrt{4pi/(2l+1)}
!r    to get the V_L used in FP-LMTO, in which
!r    V(r)=\sum_L V_L r^l Y_L(r).
!r    The extra force from the electrostatics is \sum_L Q_L grad V_L
!r    Note that in the Finnis 2nd order theory, multipoles are those
!r    of the difference in charge with respect to free atoms. The old
!r    MRS theory referred to total electronic and core charge is no
!r    longer implemented. These theories are essentially equivalent
!r    anyway.
!r
!r    If NOUAVG=F is set in CTRL, then U is averaged over all the
!r    open channels on that site. Hence the Hubbard potential is
!r    UdN where dN is the total charge transfer on that site. Otherwise
!r    the Hubbard potential is V_l = \sum_l U_{ll'} dn_l' and dn_l is
!r    the charge tansfer in the l-channel. This latter approach is
!r    closer to the TB+U formula for E^U. In that case non diagonal
!r    U_{ll'} must be defined. For now, we use U_{ll'}=min(U_l, U_l').
!r    Other possibilities would be (U_l+U_l')/2 or \sqrt(U_l*U_l')
!r    The Hubbard energy is then (1/2)\sum_{ll'}U_{ll'}dn_l dn_l'
!r
!r    In TB+U we determine on-site Coulomb integrals indirectly using
!r    I=(U+2lJ)/(2l+1) in the d channels (l=2) and U is taken from
!r    qnu and I is taken from the I= token in SPEC. In the s and p
!r    channels we take U from qnu and set J=0. We then use F^0=U
!r    and J=(5I-U)/4=(F^2+F^4)/14 if l=2. This furnishes us with the
!r    Slater parameters F^0=U (all l); F^2=14J/2.6, F^4=1.6F^2 if l=2
!r    F^2=F^4=0 if l<2.
!r
!r    The structure constants strx and dstrx may be made within tbesel
!r    on the fly for each connecting vector (switch --sfly) or they
!r    may be passed in as a lookup table (default). The latter increases
!r    speed at the expense of memory.
!r
!r    L-dimensions of structure constants strx, nlmq1 and nlmq,
!r    correspond to (l+2)^2 and (l+1)^2, respectively, where L is the
!r    l-cutoff for multipoles, or equivalently the largest l for which
!r    Delta_l'l"l\=0. For Hamiltonian or total energy calculations
!r    nlmq suffices, but for forces We need to go one angular momentum up,
!r    hence nlmq1. No such complication arises for dstrx though.
!r
!r    Meaning of input switches force and pv:
!r    parameters force and pv do not necesserily coincide with those
!r    stored in variable ltb (in 2^4 and 2^7 registers, respectively).
!r    The latter define whether the force or pressure calculations are
!r    required in principle (and hence affect dimension nlmq1) whereas
!r    the former tell whether these should be done during current call
!r    of tbesel. Both force and pv are set to .False. in tbzint.f if
!r    self-consistency has not been reached yet.
!r
!u Upgrades
!u    As of version 8.0, tbe uses (\delta q)*U rather than q*U as the
!u    'Hubbard' potential. \delta q = q - q_0 and q_0 is computed in
!u    each channel from qnu. The Hubbard U may now be the average
!u    over all channels. The older theory (as in the MRS paper) is
!u    no longer implemented.
!u    TB+U retains all terms except the Hubbard potential which is
!u    modified to be spin (s) and orbital dependent:
!u    V^{s}_{LL'}=\sum_{L''L'''} V_{LL''L'''L'} dn^{-s}_{L''L'''} +
!u                    (V_{LL''L'''L'}-V_{LL''L'L'''}) dn^{s}_{L''L'''}
!u    n is the spin polarised density matrix (see tbfrce) and dn is
!u    the difference given by subtracting the diagonal elements q_0,
!u    the V are given in terms of the Slater radial integrals as
!u    V_{LL'L''L'''}=\sum_{k=0}^{2lmax} R^k(ll'l''l''') A^k(LL'L''L''')
!u    (k an even number less that 2l) with
!u    A^k(LL'L''L''')=4pi/2k+1 \sum_{p=-k}^{k} C_{LL'''K} C_{L'L''K}
!u    where K={kp} as L={lm} and C are the Gaunt (CG) coefficients.
!u    As a first approximation to the Slater integrals we will have
!u    R^0(llll) = F^0 = U, in all channels
!u    R^2(2222) = F^2, R^4(2222) = F^4 (ie l=2, d-channel)
!u    all others zero.
!u Updates
!u        2013 (DP)  use compact dynamically sized structure constants
!u                   and improve the parallel implementation
!u   10 Nov 11 Begin migration to f90 structures
!u    3 Mar 10 (SL)  L-summation limits, parameters force and pv
!u   19 Jan 10 (ATP) Option for strux from lookup table
!u   04 Jun 08 (ATP) Cleanup, multipole contribution to pressure
!u   17 Mar 06 (ATP) Cleanup, strip out old MRS theory, start on
!u                   non orthogonal TB-L
!u    6 Jun 05 (ATP) Added TB+U
!u   15 Feb 02 (ATP) Added MPI parallelization
!u   03 Dec 01 (ATP) bug fix
! ----------------------------------------------------------------------

 implicit none
 real(8), intent(out) :: ecorr
 double precision, allocatable :: vm(:,:)
 real(8), parameter :: pi = 4d0*datan(1d0)
 integer :: ib,jb, ic, jc, it, ilm, ilmp, ilmpp, isp
 real(8) :: M, sumV

 !---------------------------------------------------------------------
 ! Madelung potential
 !---------------------------------------------------------------------
 allocate(vm(nlmi, nbas))
 vm = 0.0d0
 do  ib = 1, nbas
  do  jb = 1, nbas
   do  ilm = 1, nlmi
    vm(ilm,ib) = vm(ilm,ib) + 2d0*sum(struxd(ilm,1:nlmi,ib,jb)*qmpol(1:nlmi,jb))
   end do
  enddo
 enddo ! ib loop
 ! Sr/Ca etc: A-sites +2e monopoles go here:
 do  ib = 1, nbas
  do  jb = 1, nbasA
   do  ilm = 1, nlmi
    ! all A-site monopoles = +2e: qmpolA(1,1:nbasA) = +2.0
    vm(ilm,ib) = vm(ilm,ib) + struxdA(ilm,ib,jb)*4.0d0
   end do
  enddo
 enddo
 !---------------------------------------------------------------------


 !---------------------------------------------------------------------
 ! hamiltonian matrix elements
 ! atm(ib)%dh: spin can be dropped...
 !---------------------------------------------------------------------
 do ib = 1, nbas
  atm(ib)%dh = 0.0d0;
  ic = atm(ib)%is ! atom2species(ib) ! 
	it = atm(ib)%it ! species2type(ic) !two types: O & TM
   do  ilmp = ilm12(1,it), ilm12(2,it) ! Hilbert space
    do  ilmpp = ilm12(1,it), ilm12(2,it)! Hilbert space
     do  ilm = 1, nlmi ! potential components
       M = CFM(ll(ilmpp),ll(ilmp),ll(ilm),ic)
       atm(ib)%dh(ilmp,ilmpp) = atm(ib)%dh(ilmp,ilmpp) + &
                      vm(ilm,ib) * M * gaunt(ilmp,ilmpp,ilm)
     enddo
    enddo ! ilmpp
   enddo ! ilmp
 enddo ! ib
 !---------------------------------------------------------------------

 !C --- electrostatic energy ---
 !C ... dQ * V :
 sumV = 0.0d0
 do  ib = 1, nbas
  do  ilm = 1, nlmi
   sumV = sumV + qmpol(ilm,ib) * vm(ilm,ib)
  enddo
 enddo
 ecorr = 0.5d0*sumV
 write(*,425) ecorr 


425	format ('   (1/2) dQ dV             : ',f12.6)

 deallocate(vm)
 
 return
 end subroutine tbeseld

end module mtbeseld
