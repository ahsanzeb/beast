 module mtbeseld
 use esvar, only: nsp, nbas, nlmi, ll, ilm12, atm, struxd, & 
                  qmpol, CFM, gaunt, nbasA, struxdA, qmpolA, s_lat, &
                  hard, qref, Rlmax, cage, qaa, fVoctO, fVoctAB
 use rotylm, only: getVoct
 implicit none

 contains




!============================================================================
 subroutine tbeseld(lesH, ecorr)
                       
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
 logical, intent(in) :: lesH
 real(8), intent(out) :: ecorr
 double precision, allocatable :: vm(:,:),vm1(:,:), vmA(:), Uq
 real(8), parameter :: pi = 4d0*datan(1d0)
 integer :: ib,jb, ic, jc, it, ilm, ilmp, ilmpp, isp,l
 real(8) :: M, sumV, average, vv
 logical, save:: first=.true.


! integer, parameter, dimension(6) :: ioct = (/2,3,4,6,7,12/)

 if(.not. lesH) then  ! only set variables to 0.0 and return
	write(*,*)"tbeseld: Hes=T ==> CrysField/multipoles dH=0"
  ecorr = 0.0d0
  do ib = 1, nbas
   atm(ib)%dh = 0.0d0;
  enddo
  return
 endif


 call getvm(cage)

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

    ! skipping monopole term: ilm=1
     do  ilm = 2,nlmi ! 21,25,4 ! potential components: for l=0,1,2,... ilm starts from 1,2,5,10,17
       M = CFM(ll(ilmpp),ll(ilmp),ll(ilm),ic);
      ! if(ib==2 .and. ilmp==ilmpp) then
       ! write(*,*)'ib=2, O atom: ilm, M = ', ilm, M
			!	if(abs(M)>1.d-5) then
			!	 write(*,*)ll(ilmpp),ll(ilmp),ll(ilm),ic
			!	endif
       ! endif
       
       atm(ib)%dh(ilmp,ilmpp) = atm(ib)%dh(ilmp,ilmpp) + &
                       vm(ilm,ib) * M * gaunt(ilmp,ilmpp,ilm)
!      if(ilmp==ilmpp .and. ilm==1) then
      ! write(*,*) 'l,lp,lpp, M, vm, gaunt = ', &
      ! ilmp,ilmpp,ilm, M,vm(ilm,ib), gaunt(ilmp,ilmpp,ilm)
!       write(*,'(a,5f10.5)') 'dh = ',atm(ib)%dh(ilmp,ilmpp)    
!			write(*,'(i5,f8.3,4f10.4)') ilmp, M, vm(ilm,ib), gaunt(ilmp,ilmpp,ilm),&
! & vm(ilm,ib) * M * gaunt(ilmp,ilmpp,ilm), vm(ilm,ib) * M *0.866
!     endif
   	
!	if(ib==1 .and. ll(ilm)==4) then
!		if(dabs(gaunt(ilmp,ilmpp,ilm))> 1.d-5) then
!	write(*,'(3i5,4f16.10)')ilmp,ilmpp,ilm, &
!   vm(ilm,ib), gaunt(ilmp,ilmpp,ilm), &
!   vm(ilm,ib)*gaunt(ilmp,ilmpp,ilm)
!	!write(*,'(10f16.10)') vm(17:25,ib)
!		endif
!	end if
			
     enddo
    enddo ! ilmpp
   enddo ! ilmp
   
!   stop "tbeself.f90: stooping....."

   
 	!write(6,'(a)') 'B1: Hmp before subtracting average :'
	!write(6,'(5f12.8/)') atm(1)%dh

	! find the average of diagonal of dh and reset it to zero.
	 average = 0.0d0
	 do ilmp=ilm12(1,it), ilm12(2,it)
    average = average + atm(ib)%dh(ilmp,ilmp)
   end do
   average = average/(ilm12(2,it)-ilm12(1,it) + 1);
	 !do ilmp=ilm12(1,it), ilm12(2,it)
   ! atm(ib)%dh(ilmp,ilmp) = atm(ib)%dh(ilmp,ilmp) - average
   !end do
   !write(*,*)"tbseld: ib, aver = ",ib, average
   
   ! add Hardness term: diagonal.; also [-average] moved here from just above.
   Uq = qmpol(1,ib)*hard(ic) !- average; ! - sign for e cahrge is taken positive here; [raising a level will decrease electron occupation]
   do ilmp=ilm12(1,it), ilm12(2,it)
    atm(ib)%dh(ilmp,ilmp) = atm(ib)%dh(ilmp,ilmp) + Uq
   end do
   !write(6,*) 'ib, q*U = ', ib, Uq

   !write(*,*) "tbeself: setting atm(ib)%dh = 0 to check H=H(U&J)"
	 !atm(ib)%dh = 0.0d0
 enddo ! ib

 	!write(*,'(a,f15.10)') 'tbseld: atm(1)%dh =', atm(1)%dh(5,5)

 !---------------------------------------------------------------------
 	!write(6,'(a)') 'B1: Hmp:'
	!write(6,'(5f12.8/)') atm(1)%dh
	!write(6,'(a,50f9.5)') 'O: Hmp = ',atm(2)%dh

	!write(6,'(a)') ' Vmp:'
	!write(6,'(9f12.8/)') vm(:,:)

 deallocate(vm, vm1)

 return
 end subroutine tbeseld

!--------------------------------
	subroutine getvm(cage)
	implicit none
	logical, intent(in) :: cage

	if(cage) then
		call getvmcage()
	else
		call getvmlattice()
	endif

	return
	end subroutine getvm
!--------------------------------
	subroutine getvmcage() ! only l=4 potential, in unrotated frame.
	implicit none
	double precision, parameter :: dsqrt57= dsqrt(5.0d0/7.0d0)
	! vm for o; vm1 for A+B

	! Oxygen octahedral potential rotated back to the original unrotated frame
	vm =0.0d0;
	do ib=1,4 ! only 4 octahedra of orthorhombic; can be used for nlayers = 2*
		vm(17:25,ib) = Rlmax(17:25,21,ib) + dsqrt57 * Rlmax(17:25,25,ib) ! vm in rotated frame; 
		vm(17:25,ib) = vm(17:25,ib) * fVoctO! fVoctO=(35.0d0/4.0d0)*(2.0d0/dc**5) ! prefactor: Pavarini, orbital order notes
	end do


	! A+B potential in original frame (it is independent of rotation)
	vm1 = 0.0d0
	do ib=1,4 ! only 4 octahedra of orthorhombic; can be used for nlayers = 2*
		vm(21,ib) = 1.0d0; vm(25,ib) = dsqrt57;
		vm(:,ib) = vm(:,ib) * fVoctAB;
	end do
	

	return
	end subroutine 
!--------------------------------
	subroutine getvmlattice()
	implicit none
 double precision, allocatable :: vm(:,:),vm1(:,:), vmA(:), Uq
 real(8), parameter :: pi = 4d0*datan(1d0)
 integer :: ib,jb, ic, jc, it, ilm, ilmp, ilmpp, isp,l
 real(8) :: M, sumV, average, vv
 logical, save:: first=.true.


 !---------------------------------------------------------------------
 ! Madelung potential
 !---------------------------------------------------------------------

 !qmpol(:,:)=0.0d0;
 !qmpolA = 0.0d0;
 
 !qmpol(:,6:8)=0.0d0;
 !write(*,*)'tbeseld.f90: testing: setting qA & qmpol = 0'

 ! struxd(ilm,jlm,ib,jb) has Ylm of jb atom at location of ib atom [ilm comp? expanded in Ylm of ib?]
 ! due to TM/O sites

 !write(*,*)"test: tbeseld: exclude O atoms: if(atm(ib)%it==1 ..."


!M = 0.0d0
!do ib=1,nbas
!	if(mod(ib-1,4) /=0 ) M = M + struxd(1,1,1,ib)
!end do

!write(*,*) 'B_{0,0} = ', M
!write(*,'(8f10.5)') struxd(1,1,1,:)
!write(*,'(8f10.5)') struxd(1,1,:,1)

! testing: set B atoms charges zero
!do jb=1,nbas
!	if(mod(jb-1,4) /= 0) qmpol(:,jb) = 0.0d0
!enddo

!write(*,*)'tbseld: qmpol(:,jb) = ',qmpol(:,1)


!m=struxd(1,1,1,1)+struxd(1,1,1,13)
!sumv=struxd(1,1,1,5)+struxd(1,1,1,9)
!write(*,*)' +, -, tot = ', m, sumv, (sumv-m)*7


 allocate(vm(nlmi, nbas))
 vm = 0.0d0

if(1==1 .and. first) then

write(*,*) "tbseld: writing B1 potential; A,B,O resolved" 
open(101,file='vm-B1.dat',action='write',position='append')

do  ib = 1,1! nbas
	! testing: exclude Ir to consider only O atoms.
 	!if(mod(jb-1,4) ==0) cycle
  
!  if(ib==1 .and. mod(jb-1,4)==0) then !
!   write(*,*)'i,j, v_ij = ',ib,jb, struxd(1,1,ib,jb)*qmpol(1,jb),struxd(1,1,jb,ib)*qmpol(1,jb)
!  endif

!	vm(:,ib) = 0.0d0
!	do jb=17,20
!    do  ilm = 17, nlmi
!     vm(ilm,ib) = vm(ilm,ib) + &
!     sum(struxd(ilm,1:nlmi,ib,jb)*qmpol(1:nlmi,jb))
!     !write(*,*) 'ilm,strx:',ilm,sum(struxd(ilm,1:nlmi,ib,jb)*qmpol(1:nlmi,jb))
!    end do
!	end do
!  write(*,*) 'A: Y40, Y44 : ',vm(21,ib),vm(25,ib)

	vm(:,ib) = 0.0d0
	do jb=1,16,4
    do  ilm = 1, nlmi
     vm(ilm,ib) = vm(ilm,ib) + &
     sum(struxd(ilm,1:nlmi,ib,jb)*qmpol(1:nlmi,jb))
    end do
	end do
  write(*,*) 'B: Y40, Y44 : ',vm(21,ib),vm(25,ib)

! write B1 potential, A,B,O resolved. 
write(101,*) vm(:,ib)
!write(101,'(a,9f12.6)')'...................................'
!write(101,'(a,9f12.6)')'vm1 ', vm(1:9,ib)
!write(101,'(a,9f12.6)')'...................................'
!write(101,'(a,9f12.6)')'B1: ',struxd(1:9,1,1,1)
!write(101,'(a,9f12.6)')'B3: ',struxd(1:9,1,9,1)
!write(101,'(a,9f12.6)')'B2: ',struxd(1:9,1,5,1)
!write(101,'(a,9f12.6)')'B4: ',struxd(1:9,1,13,1)

	vm(:,ib) = 0.0d0
	do jb=1,16
		if(mod(jb-1,4)==0) cycle
    do  ilm = 1, nlmi
     vm(ilm,ib) = vm(ilm,ib) + &
     sum(struxd(ilm,1:nlmi,ib,jb)*qmpol(1:nlmi,jb))
    end do
	end do
  write(*,*) 'O: Y40, Y44 : ',vm(21,ib),vm(25,ib)

! write B1 potential, A,B,O resolved. 
write(101,*) vm(:,ib)

end do

endif ! 1==0 and first










vm = 0.0d0
do ib=1,nbas
 	do jb=1,nbas
   do  ilm = 1, nlmi
    vm(ilm,ib) = vm(ilm,ib) + &
    sum(struxd(ilm,1:nlmi,ib,jb)*qmpol(1:nlmi,jb)) !
   end do
	end do
enddo ! ib loop

 !write(*,'(100f10.4)') vm(:,1)


! potentail due to A-site monopoles
 allocate(vm1(nlmi, nbas))
 vm1 = 0.0d0

 ! Sr/Ca etc: A-sites +2e monopoles go here:
 do  ib = 1, nbas
  do  jb = 1, nbasA
   do  ilm = 1, nlmi
    ! all A-site monopoles = +2e: qmpolA(1,1:nbasA) = -2.0 ! elec charge taken positive.
    vm1(ilm,ib) = vm1(ilm,ib) + struxdA(ilm,ib,jb)*qmpolA 
   end do
  enddo
 enddo

if(first)then
! write B1 potential, A,B,O resolved. 
	write(101,*) vm1(:,1)
 	close(101)
	first = .false.
endif

 !C --- electrostatic energy ---
 !C ... dQ * V :
 ! A atoms full term, other 1/2 due to double counting in the sum.
 sumV = 0.0d0
 do  ib = 1, nbas
  do  ilm = 1, nlmi
   sumV = sumV + qmpol(ilm,ib) * (0.5d0*vm(ilm,ib) + vm1(ilm,ib) )
  enddo
 enddo
 ecorr = sumV
 write(*,425) ecorr
 
425	format ('   (1/2) dQ dV             : ',f12.6)

!write(*,'(a,20f15.8)') 'VmB: ', vm(1,1) !/3.544907701811032d0 ! sqrt(4pi)
!write(*,'(a,20f15.8)') 'VmA: ', vm1(1,1)

!write(*,*) 'vm1(1,1) = ', vm1(1,1)
!write(*,*) 'vm(1,1) = ', vm(1,1)

!	write(*,*)'qmpolA = ',qmpolA

	write(*,'(a,100f15.10)')'B1: vm = ',vm(17:,1)
	write(*,'(a,100f15.10)')'B1: vm1 = ',vm1(17:,1)
	
 ! now combine the two terms to get full potential for Hij:
 vm = vm + vm1;

 write(*,'(a,100f10.4)') 'Total V_0 at B1: ',vm(1,1)
 
 !write(*,'(a,100f10.4)') 'Total V_0 at B2: ',vm(1,5)
 !write(*,'(a,100f10.4)') 'Total V_0 at B3: ',vm(1,9)
! write(*,'(a,100f10.4)') 'Total V_0 at B4: ',vm(1,13)

! CsCl: 1.763

!	stop 'tbseld:  stop '

 !write(*,'(a)')'tbeseld: !vm = vm + vm1; enable it. & -vm used'

!write(*,*) 'vm = vm1+vm2:' 
!write(*,*) 'vm1(1,1) = ', vm1(1,1)
!write(*,*) 'vm(1,1), alpha_CsCl, s_lat%alat ',s_lat%alat
!write(*,'(2f16.12)')  vm(1,1), vm(1,1)*dsqrt(3.0d0)/2.0d0*s_lat%alat


!write(*,*) 'tbseld: vm(1,ib=1) = ',vm(1,1)

!write(*,*) 'vm(1,1), alpha_NaCl, s_lat%alat ',s_lat%alat
!write(*,'(2f16.12)')  vm(1,1), vm(1,1)*s_lat%alat



!write(*,'(2f16.12)')  vm(1,1), vm(1,1)*s_lat%alat


!write(*,'(a,20f15.8)') 'VmA+VmB: ', vm(1,1)
!write(*,*)'Vm: l=0,1,2,3,4:'
!write(*,'(10f16.10)') vm(1:4,1)
!write(*,'(10f16.10)') vm(5:9,1)
!write(*,'(10f16.10)') vm(10:16,1)
!write(*,'(10f16.10)') vm(17:25,1)
!write(*,*) 'tbeseld: setting V40=1 and V44=sqrt(5/7) for atom 1: '
!vm(21,1) = 1.0d0;
!vm(25,1) = dsqrt(5.0d0/7.0d0);
!write(*,*) 'tbeseld: setting V44=-V44 for atom 2: '
!vm(25,2) = -vm(25,2)
 !---------------------------------------------------------------------
 !write(*,'(a)')'..... .... .... qmpol ..... .... .... '
 !do ib=1,4
  !write(*,'(i5, 100f10.5)') ib, qmpol(1,ib)
 !end do
 write(*,'(a,50f15.8)')"Layer Q0: ", (sum(qmpol(1,(ib-1)*8+1:ib*8))-4,ib=1,nbas/8)

 write(*,'(8f6.2)') qmpol(1,:)
 !write(*,*) 'TM Q0:', qmpol(1,1),qmpol(1,5)
 !write(*,*) 'TM V0:', vm(1,1), vm(1,5)
 !write(*,*) 'V0:', vm(1,:)


	! testing: Y40/Y44 ratio diff than sqrt[5/7]
 	!vm(21,1) = 1.0d0; vm(25,1) = 0.8 !dsqrt(5.0d0/7.0d0)
	return
	end subroutine 






end module mtbeseld
