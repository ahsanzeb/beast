
	module mtbmpol
	use esvar
	implicit none

	contains

      subroutine tbmpol() !nbas,nsp,nl,nlmq,qnu,ipc,lmxl,gaunt,qpol,rho,rhoc,qmpol,mmom)
C- Make multipole (and magnetic) moments from TB "density matrix"
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas,nl,nsp,qnu,ipc,gaunt
Ci   nlmq : L-cutoff for multipoles, also leading dimension of qmpol
Ci   lmxl : species-resolved l-cutoffs for multipoles
Ci   qpol : polarisation parameters (delta_ll'l'')
Ci   rho, rhoc (see tbfrce); if mixrho then these are from rhomix.
Co Outputs:
Co   qmpol: electrostatic point multipoles Q_L [Finnis, Eq.(7.76)]
Co   mmom : magnetic moments
Cr Remarks
Cr   Our multipoles are defined such that the monopole moment is the
Cr   total Mulliken charge, while the higher multipoles are approximated
Cr   as an onsite sum without reference to overlap (Finnis, p. 216)
C ----------------------------------------------------------------------
	implicit none
C Passed Parameters
!	integer, intent(in) :: nbas,nl,nlmq,nsp,ipc(nbas),lmxl(*)
!	double precision, intent(in)  :: qpol(10,*),qnu(3,0:nl-1,nsp,*),
!     .                                 gaunt(9,9,25) !,rho(nl,nsp,nbas),rhoc(nl**2,nl**2,nbas)
!	double precision, intent(out) :: qmpol(nlmi,nbas), mmom(nbas)
C Local Variables
	integer :: ib,ic,ilm,ilmp,ilmpp, it
	double precision :: M
c...deb
	double precision :: qmpt(nlmi)

c... Checks
	if (nlmi > 25)write(*,*)' tbmpol: nlmq is too big, nlmq = ',nlmi

	qmpol = 0.0d0;
	qmpt = 0.0d0;

C --- get multipole moments ---
	do  ib = 1, nbas
	 ic = atm(ib)%is ! atom2species(ib) ! class/species index
	 it = atm(ib)%it ! species2type(ic) ! =atm(ib)%it ! two types: O & TM
	 qmpol(1,ib) = atm(ib)%qs(1)
	 if (nsp == 2) then
	  atm(ib)%mag = atm(ib)%qs(1) - atm(ib)%qs(2) 
	  qmpol(1,ib) = atm(ib)%qs(1) + atm(ib)%qs(2) 
	 end if
	 qmpol(1,ib) = qmpol(1,ib) - q0(ic) ! ? atm(ib)%q0: so many O atoms, duplicate data!

	!write(*,*) 'ib, atm(ib)%qs = ',ib, atm(ib)%qs
	!write(*,'(2x,100000f6.2)') q0(ic), qmpol(:,ib)
	!write(*,*) 'atm(ib)%rhoc is complex: how can we make qmpol real?'
	!write(*,*) 'testing. using dble(atm(ib)%rhoc) ... fix it later..'

	!write (*,100)
        
	if (nlmi > 1) then
	 do ilm = 2, nlmi
	  do ilmp = ilm12(1,it), ilm12(2,it)
			! use the exchange symmetry of l' & l'' in eq 41 of Paxton's notes
			! sum_{l',l''} X_{l',l''} <====> sum_{l',l''=l',lmxl} 2*Re[X_{l',l''}]
			! as the complex rhoc terms will become 2Re[*] and M*gaunt terms are 
			! already symmetric under this exchange.
			! this also allows to use a real qmpol. above commented section had to use a complex qmpol.
			!........................................
			! diagonal in ll, l'': already real.
	   M = CFM(ll(ilmp),ll(ilmp),ll(ilm),ic)
	   qmpol(ilm,ib) = qmpol(ilm,ib) + 
     .    dble(atm(ib)%rhoc(ilmp,ilmp)) * 
     .    M * gaunt(ilmp,ilmp,ilm)
			!........................................
			! off-diagonal: combined upper & lower becomes real.
	   do  ilmpp = ilmp+1, ilm12(2,it)
	    M = CFM(ll(ilmpp),ll(ilmp),ll(ilm),ic)
	    qmpol(ilm,ib) = qmpol(ilm,ib) + 
     .    2.0d0*dble(atm(ib)%rhoc(ilmp,ilmpp)) *
     .    M * gaunt(ilmp,ilmpp,ilm)
	   enddo

	  enddo
	 enddo
	endif

	enddo ! ib

100	format ('  L''   L''''  L    l''   l''''  l      M         CG
     .      rho_L''L''''')
200	format (6(i3,2x),2(2x,f6.2,2x),2f10.6)

c...deb
C total MP moments of the cell

      if (1==0) then
        qmpt = 0d0
        do ib = 1, nbas
          qmpt(1:nlmi) = qmpt(1:nlmi) + qmpol(1:nlmi,ib)
          write(*,'(i5,2x,100000f6.2)') ib, qmpol(:,ib)
        enddo

        print *,' tbmpol: total multipole moments of the cell: '
        print *,' (maz: ill defined! diff centres summed!)'

        write(*,'(1x,10000f8.4)') qmpt
        if (nlmi >= 4) print '(1x,25f8.4)',qmpt(2:4)
        if (nlmi >= 9) print '(1x,25f8.4)',qmpt(5:9)
        if (nlmi >= 16) print '(1x,25f8.4)',qmpt(10:16)
        if (nlmi >= 25) print '(1x,25f8.4)',qmpt(17:25)
      endif
c...deb

	end subroutine tbmpol

	end module mtbmpol
