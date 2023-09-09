
	program perovskite
	use modmain
	use skoster
	use hamiltonian
	use readinput
	use fermi
	use modgaunt, only: mkdgaunt
	use Hubbard, only: mkvee, mkvmat, rtozflm
	use scf, only: groundstate
	use scf, only: nscfgroundstate, tmgroundstate, tmscfground
	use pdos

	use estatic, only: setmadvar, initmadelung
	use mtbeseld
	
	implicit none
	integer :: il, io, i, ia, iu, is, is2
	integer :: ik,ib
	double precision, dimension(3) :: kvec
	double precision :: soc1, soc2	
	double precision, dimension(10,10) :: Rl2s
	double complex, dimension(10,10) :: Hsocr
	

	write(*,'(a)')
	write(*,'(a)')
	write(*,'(a)')
	write(*,'(a)')'=================================================='
	write(*,'(a)')"*************** BEAST STARTING *******************"
	write(*,'(a)')'=================================================='

	write(*,'(a)') 'fix s_lat%awald ~ 0.32 and remove ewalds from inp'
	write(*,'(a)')'rm:struxd = dble (INT(struxd * decimal)) / decimal'
	!write(*,'(a)')'Attention: check if we are consistent?'
	!write(*,'(a)')' TB: Ry, while ELK: Hartree! ==> '
	!write(*,'(a)')'Madelung in Ry but Hubburd U/J in Haretree'
	!write(*,'(a)')'=================================================='

	write(*,*) floor(-2.6),floor(2.6)
	
	! read input file 'input.in'
	call input()


	if(lhu) then ! HubbardU
	 ! make gaunt coefficients matrix that will be used to calc Vee 
	 allocate(gcmat(0:2,-2:2,-2:2,-2:2,-2:2)) ! for l=2, d orbitals
	 call mkdgaunt(2,gcmat)
	 call mkvee(2); !l=2 for d orbitals
	endif
	!write(*,'(a,10000f10.4)') 'Gaunt = ', gcmat(:,0,0,0,0)
	!stop

	!write(*,'(a)') "Give Number of layers & phi and hit enter:"
	!read(*,*) nlayers, phi0
	!nlayers=5; phi0=10.0d0

	! set structure/geometry
	if(.not. xsf) call getstructure2()
	!write(*,*) "struc done... "

	! write GEOMETRY.OUT for visualisation with VESTA
	call writegeom()

	if(1==0) then
	open(111,file='avecs.dat',action='write',position='append')
	write(111,'(9G18.10)') avec(:,1),avec(:,2),avec(:,3)
	close(111)
	endif
	!call writegeomxsf()
	!write(*,*)'-------------------- 1'
	
	! set nearest neighbours:
	!call getneighbours()

	call getnns(natoms, pos) !,tolnns)
	!write(*,*) "getnns done... "
	
	! set maps:
	call getmaps()
	!write(*,*) "getmaps done... "
	
	! reciprocal lattice vectors
	call reciplat(avec,bvec,omega,omegabz)

	!write(*,*) 'main: test... omega,omegabz'
	!omega=1.0d0; omegabz=1.0d0;

	!write(*,'(3f10.5)') avec
	!write(*,'(3f10.5)') bvec
	!write(*,'(3f10.5)') matmul(avec,bvec)	
 !--------------------------------------------------------
 ! once before SCF cycle starts:
	if (lhu) then ! bind the calculation of multipoles with the Hubbard e-e, lhu.
	 call setmadvar()
	 call initmadelung()
	endif
 !--------------------------------------------------------
	
	! calculate Vmpol and corresping H_{i,j} due to Vmpol & Qmpol
	!call tbeseldx(omegabz)
	
	! dummy data set:
	! nsptm = number of species of TM atoms
	!allocate(skbo(nsptm,2)) 	! 2: sigma_pd, pi_pd
	!allocate(skbb(nsptm,nsptm,3))	! 3: sigma_dd, pi_dd, delta_dd
	!skbo = 0.0d0; skbb=0.0d0;
	!skbo(1,:) = (/1.0d0,0.5d0/);
		
	!skbb(1,1,:)	 = (/0.0d0,0.0d0,0.00d0/);
	!skoo = 0.0d0 ! some prob with o-o, calc... if finite skoo.

	! sets some basis transformation matrices:
	! some condition? inside an if block?
	call setUrUz() ! sets Uz and Ur matrices

	if (ham4edrixs) then ! assuming SOC is not incuded.

	 !call mkhsocl2()
	
	 if(lcage) then
	  call getvmdhcage4edrixs()
	 else
	  call getvmdhlat4edrixs()
	 endif

	 call writeB1ham4edrixs()
	 stop
	endif

	call setUlm2j() ! sets Ulm2j, uses Ur that is set in setUrUz() call above.


	! real hamiltonian matrix elements using SK method
	call realHij()

	! spin-orbit coupling universal hamiltonian in full d-orbital space of a TM atom
	call mkhsocl2()

	if(1==0) then
	open(100,file='Hsoc.dat', action='write')
	 do i=1,10
	  write(100,'(10f6.2)') dble(Hsoc(i,:))
	 end do
	 do i=1,10
	  write(100,'(10f6.2)') dimag(Hsoc(i,:))
	 end do
	close(100)
	endif

	! testing 
	if (1==0) then
	
	Rl2s = 0.0d0;
	Rl2s(1:5,1:5) = Rlmax(5:9,5:9,1)
	Rl2s(6:10,6:10) = Rlmax(5:9,5:9,1)

	Hsocr = matmul(Rl2s,Hsoc)
	Hsocr = matmul(Hsocr,transpose(Rl2s))

	open(100,file='Hsoc-rot.dat', action='write')
	 do i=1,10
	  write(100,*) dble(Hsocr(i,:))
	 end do
	 do i=1,10
	  write(100,*) dimag(Hsocr(i,:))
	 end do
	close(100)
	
	endif
	



	!call realHii() ! sets onsite hamiltonian matrix elements

	!write(*,*)'------------real Hij done -------- '

	if(1==0) then
	write(*,'("Reciprocal lattive vectors:")')
	write(*,'(3f10.3)') bvec(:,1)
	write(*,'(3f10.3)') bvec(:,2)
	write(*,'(3f10.3)') bvec(:,3)
	endif

	!write(*,*) 'lhu, lgs, lscf, ltmgs = ',lhu, lgs, lscf, ltmgs
	if(lhu) then

	 if(lgs) then
		 if (lscf) then
		  if (ltmgs) then
		   call tmscfground()
		 	 write(*,*)'main: TM SCF ground state calculated!'
		 	 stop
		  else
	     call groundstate()
	    endif
		 elseif (ltmgs) then
		 	call tmgroundstate() ! only tm atoms ham diagonalisation
		 	write(*,*)'main: TM ground state calculated!'
		 	stop
		 else
		  call nscfgroundstate() ! full system diagonalisation
		 endif
	 endif


	if(1==0) then
	 open(10,file='LayerQ0.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	 write(10,'(50f15.8)') Hub(1)%U, Hub(2)%U, 
     .                  soc,
     .    (sum( qmpol(1,(ib-1)*8+1:ib*8)-4.0d0 ),ib=1,nlayers) ! -4.0 to set neutra at 0.0 ; otherwise 2e per A sites ==> 4.0 per layer
	 close(10)
	 open(10,file='Q0.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	 write(10,'(100000f15.8)') Hub(1)%U, Hub(2)%U, 
     .                  soc, qmpol(1,:)
	 close(10)
	end if ! 1==0
	
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	elseif (nsptm > 1 ) then
	!---------------------------------------
	! initial values
	Hub(isploop(1))%U = uloop(1)	
	do while ( Hub(isploop(1))%U <=uloop(2) )
	 is = isploop(1)
	 !call 	setFk(is,Hub(is)%U,Hub(is)%J, eV2Har) ! sets Fk using U&J
	 Hub(isploop(2))%U = uloop(4)	
	 do while (Hub(isploop(2))%U <= uloop(5))
	  is2 = isploop(2)
	  !call 	setFk(is2,Hub(is2)%U,Hub(is2)%J, eV2Har) ! sets Fk using U&J
	  soc1 = sloop(1)	
	  do while (soc1 <= sloop(2))
	   soc(isploop(1)) = soc1*eV2Har

	   soc2 = sloop(4)	
	   do while (soc2 <= sloop(5))
	   soc(isploop(2)) = soc2*eV2Har
	 !---------------------------------------




	write(*,*)"**************************************"
	!write(*,'(a,50f15.8)')"U  : ", uloop
	!write(*,'(a,50f15.8)')"SOC:", sloop
	write(*,*) 'isploop:',isploop
	write(*,'(a,50f15.8)')'Us :', Hub(isploop(1))%U, Hub(isploop(2))%U
	write(*,'(a,50f15.8)')'SOCs: ', soc1, soc2


	if(lgs) call groundstate()
	open(10,file='LayerQ0.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	write(10,'(50f15.8)') Hub(isploop(1))%U, Hub(isploop(2))%U, 
     .                  soc1, soc2, 
     .    (sum( qmpol(1,(ib-1)*8+1:ib*8) ),ib=1,nlayers)
	close(10)
	open(10,file='Q0.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	write(10,'(100000f15.8)') Hub(isploop(1))%U, Hub(isploop(2))%U, 
     .                  soc1, soc2, qmpol(1,:)
	close(10)

	!.......................................................
	open(22,file='EFERMI.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	write(22,'(G20.12)') efermi
	close(22)
	!.......................................................

	if(lpdos) call getpdos(efermi) ! mine, maxe, nwplot

	 !---------------------------------------
	    soc2 = soc2 + sloop(6)
	   end do ! sp2
	   soc1 = soc1 + sloop(3)
	  end do ! sp1
	  Hub(isploop(2))%U = Hub(isploop(2))%U + uloop(6) ! u + du2
	 end do
	 Hub(isploop(1))%U = Hub(isploop(1))%U + uloop(3) ! u + du1
	end do
	!---------------------------------------
	else ! nsptm = 1 !+++++++++++++++++++++++++++++++++++++++++++
	!---------------------------------------
	! initial values
	Hub(isploop(1))%U = uloop(1)	
	do while (Hub(isploop(1))%U .le. uloop(2))
	 is = isploop(1)
	 call 	setFk(is, eV2Har) ! Hub(is)%U,Hub(is)%J, !sets Fk using U&J
	  soc1 = sloop(1)	
	  do while (soc1 .le. sloop(2))
	   soc(isploop(1)) = soc1*eV2Har
	 !---------------------------------------
	if(lgs) call groundstate()
	open(10,file='LayerQ0.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	write(10,'(50f15.8)') Hub(isploop(1))%U, Hub(isploop(2))%U, 
     .                  soc1, soc2, 
     .    (sum( qmpol(1,(ib-1)*8+1:ib*8)-4.0d0 ),ib=1,nlayers) ! -4.0 to set neutra at 0.0 ; otherwise 2e per A sites ==> 4.0 per layer
	close(10)
	open(10,file='Q0.OUT',form='FORMATTED',action='write', 
     .                              position='append')
	write(10,'(100000f15.8)') Hub(isploop(1))%U, Hub(isploop(2))%U, 
     .                  soc1, soc2, qmpol(1,:)
	close(10)
	 !---------------------------------------
	   soc1 = soc1 + sloop(3)
	  end do ! sp1
	 Hub(isploop(1))%U = Hub(isploop(1))%U + uloop(3) ! u + du1
	end do
	!---------------------------------------
	endif


	if(1==0) then
	write(*,'(a)')'TM:-------------------'
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always
	  ia = tm(il,io)%ia
	  write(*,'(a,i5,a,i5)') 'ia=',ia, '; is=',atom2species(ia)
	 end do
	end do
	end if

	
	!if(lpdos) call getpdos(efermi) ! mine, maxe, nwplot

	if(lbands) then
	 call getbands()
	 call writebands()
	 ! band character:
	 if(lbc) then
	  call getbc()
	 endif
	endif
	
	
	write(*,'(a)')'--------- complete -----------'


	contains
!======================================================================
! internal routines:
!======================================================================

	subroutine getbands()
	implicit none
	
	allocate(vpl(3,np))
	allocate(dv(nv))
	allocate(dp(np))

	call plotpt1d(bvec,nv,np,vvl,vpl,dv,dp)

	!write(*,*) vpl
	
	! vpl has kpoints along the bandlines....
	if(allocated(hk)) deallocate(hk)
	if(allocated(eval)) deallocate(eval)
	if(allocated(evec)) deallocate(evec)
	allocate(hk(ntot,ntot))
	allocate(eval(np,ntot))
	allocate(evec(np,ntot,ntot))

	open(1,file='KPOINTS.OUT',form='FORMATTED',action='write')

	do ik= 1,np
	 !call getHk((/0.5d0,0.5d0,0.5d0/),hk)
	 !write(*,'(a,i10,3f10.4)')' ik, k = ',ik,vpl(:,ik)
	 kvec = bvec(:,1)*vpl(1,ik) + 
     .    bvec(:,2)*vpl(2,ik) +
     .    bvec(:,3)*vpl(3,ik)
	 call getHk(ik,kvec, hk, 0)
	 !write(*,*)' ham done... '
	!write(1,'(3f15.6)') kvec


	if(1==0)then
	il=0
	do i=ntottm+1,ntot ! Oxygen orbitals
	if(norm2(abs(hk(i,1:ntottm))) < 1.0d-8) then
	il = il + 1
	endif
	end do
	write(*,*) 'ik, i_uncoupled = ',il
	endif
	
	if(1==0) then
	il = 0;
	do i=ntottm+1,ntot ! Oxygen orbitals
	if(norm2(abs(hk(i,1:ntot))) < 1.0d-8) then
	il = il + 1
	!write(*,'(a,i3,a,1000f6.1)') 'i=',i,' H(i,1:ntottm)=0'
	!write(*,*)'O: io = ', (i-ntottm)
	endif
	end do
	write(*,*) 'a total of ',il,' O orbitals uncoupled from TMs!'
	endif

	if(1==0)then
	open(101,file='hk.dat',action='write')
	write(101,*) dble(hk)
	write(101,*) dimag(hk)
	close(101)
	endif

	
	 call zdiag(ntot,hk,eval(ik,:),ik,ntot)
	 evec(ik,:,:) = Hk ! eigenvector are columns of Hk

	 
	! zdiag: hk contains eigenvectors (columns) on return 
	if(ik==1 .AND. 1==0) then
	if (1==1) then
	write(*,'(10000f10.2)') eval(ik,20)
	write(*,'(a,10000f10.2)') 'TM:', real(hk(1:20,21))
	write(*,'(a,10000f10.2)') 'O :',real(hk(21:ntot,21))
	endif
	!if(norm2(hk(1:36,1:20)) > 1.0d-6) 
	!write(*,*)'norm2(hk(...))=',norm2(abs(hk(1:20,21:36)))
	!write(*,*)'norm2(hk(37:,21:))=',norm2(hk(1:20,1:36))
	endif
	 
	end do ! ik
	close(1)
	
	return
	end 	subroutine getbands
!----------------------------------------------------------------------

	subroutine writebands()
	implicit none	
	open(10,file='BAND.OUT',form='FORMATTED',action='write')
	do ib=1,ntot
	 do ik=1,np
	  write(10,'(2G20.8)') dp(ik), eval(ik,ib) !-efermi
	 end do
	 write(10,'(2f15.8)')
	 write(10,'(2f15.8)') 
	end do
	close(10)

	! output the vertex location lines
	open(50,file='BANDLINES.OUT',form='FORMATTED',action='write')
	do i=1,nv
	write(50,'(2G18.10)') dv(i),0.0 !emin
	write(50,'(2G18.10)') dv(i),0.0 !emax
	write(50,'("     ")')
	end do
	close(50)

	open(10,file='eval-t2g.dat',action='write',position='append')
	 write(10,'(100G20.8)') (eval(1,ib), ib=73,84,1)
	close(10)


	return
	end 	subroutine writebands

!----------------------------------------------------------------------

	subroutine getstructure()
	implicit none
	integer :: i

	if(any(phi(:) > 0.0)) then
	if(any(phi(:) > 45.0)) then
	 write(*,'("O atoms have to cross to make this big rotation!")')
	 stop 
	else
	 write(*,'(a)')
     .  'Octraherda size will be rescaled to keep "a" fixed.'
	endif
	endif
	
	!tmnn2 = .true.;
	!oxnn2 = .true.;


	!nlayers = 2;
	noctl = 2;
	noct = noctl*nlayers;
	natoms = noct*4;
	a = a0; !7.0d0;


	!nsptm =1;
	nspecies = nsptm;
	
	norbtm = 5; norbtms = norbtm*nspin;
	norbo = 3; norbos = norbo*nspin;

	ntot = noct*norbtms + noct*3*norbos;
	ntottm = noct*norbtms;

	write(*,*)'ntot, ntottm, qtot,a = ', ntot, ntottm, qtot, a
	!i = noct*1 + noct*6*norbos
	!write(*,*)'TM d^1: nelec = ', i
	!write(*,*)'if only one spin, filled bands ~ ',0.5*i



	
	allocate(oct(nlayers,noctl))
	!allocate(phi(nlayers))
	! phi 
	!phi = 0.0d0;

	! read/set phi:
	do il=1,nlayers
	 !if(mod(il,2)==0) then
	 ! phi(il) = phi0*pi/180.0d0;
	 !else
	 ! phi(il) = -phi0*pi/180.0d0;
	 !endif
	 oct(il,1)%phi = +phi(il)
	 oct(il,2)%phi = -phi(il) 
	end do

	! set the basic cubic structure
	! sqrt(2) x sqrt(2) x nlayers cell
	! square cell rotated by 45 wr.r.t x,y
	! lattice vectors
	!
	avec(:,1) = (/1.0d0, -1.0d0, 0.0d0/)*a
	avec(:,2) = (/1.0d0, +1.0d0, 0.0d0/)*a
	avec(:,3) = (/0.0d0,  0.0d0, 1.0d0/)*a*nlayers


	!ainv = matmul(avec,ainv)
	!write(*,'(3f10.3)') ainv(1,:)
	!write(*,'(3f10.3)') ainv(2,:)
	!write(*,'(3f10.3)') ainv(3,:)

	! atomic positions in cartesian
	! set pos of B
	do il=1,nlayers
	 oct(il,1)%rb = (/0.5d0,-0.5d0,1.d0*(il-1)/)*a;
	 oct(il,2)%rb = (/0.5d0,+0.5d0,1.d0*(il-1)/)*a;
	enddo

	! oxygens
	do il=1,nlayers
	 do io=1, noctl
	  	oct(il,io)%xo(:,1) = (/0.5d0,0.0d0,0.0d0/)*a
 !    .   +(/rand(0),rand(0),rand(0)/)*a
	  	oct(il,io)%xo(:,2) = (/0.0d0,0.5d0,0.0d0/)*a
!     .   +(/rand(0),rand(0),rand(0)/)*a
	  	oct(il,io)%xo(:,3) = (/0.0d0,0.0d0,0.5d0/)*a
!     .   +(/rand(0),rand(0),rand(0)/)*a
	  oct(il,io)%lo = 0.5d0 * a;
	 end do
	 end do




	if(1==0) then
	if(nlayers>1) stop 'testing isolated octahedron: use nlayers=1'

	! set pos and posA for getnns()
	allocate(posA(nlayers*2,3))
	allocate(pos(nlayers*8,3))
	
	do il=1,nlayers
	  io=1 
	 	i = (il-1)*8 + (io-1)*4;
	  pos(i+1,:) = oct(il,io)%rb
	  pos(i+2,:) = oct(il,io)%rb + oct(il,io)%xo(:,1)
	  pos(i+3,:) = oct(il,io)%rb + oct(il,io)%xo(:,2)
	  pos(i+4,:) = oct(il,io)%rb + oct(il,io)%xo(:,3)
		posA((il-1)*2 + io,:) = oct(il,io)%rb + (/0.5d0,0.5d0,0.5d0/)*a

	  io=2 
	 	i = (il-1)*8 + (io-1)*4;
		! associate all 6 O atoms to first octahedron:
	  pos(i+2,:) = oct(il,1)%rb - oct(il,1)%xo(:,1)
	  pos(i+3,:) = oct(il,1)%rb - oct(il,1)%xo(:,2)
	  pos(i+4,:) = oct(il,1)%rb - oct(il,1)%xo(:,3)

	  pos(i+1,:) = oct(il,io)%rb
		posA((il-1)*2 + io,:) = oct(il,io)%rb + (/0.5d0,0.5d0,0.5d0/)*a

	end do

	! testing: redefine avec: change unnit cell to isolate the favourite octahedron:
	!avec(:,1) = 100* avec(:,1) 
	!avec(:,2) = 100* avec(:,2) 
	!avec(:,3) = 100* avec(:,3) 


	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)

	return
	
	endif ! 1==0







	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)

	! rotate and set abs pos of oxygen
	do il=1,nlayers
	 do io=1, noctl
	  call rotate(il,io)
	 end do
	end do


	! set pos and posA for getnns()
	allocate(posA(nlayers*2,3))
	allocate(pos(nlayers*8,3))

	do il=1,nlayers
	 do io=1, noctl
	 	i = (il-1)*8 + (io-1)*4;
	  pos(i+1,:) = oct(il,io)%rb
	  pos(i+2,:) = oct(il,io)%rb + oct(il,io)%xo(:,1)
	  pos(i+3,:) = oct(il,io)%rb + oct(il,io)%xo(:,2)
	  pos(i+4,:) = oct(il,io)%rb + oct(il,io)%xo(:,3)
		posA((il-1)*2 + io,:) = oct(il,io)%rb + (/0.5d0,0.5d0,0.5d0/)*a
	 end do
	end do


	if(1==0) then
	
	write(*,*)' avec: '
	do i=1,3
		write(*,*) avec(i,:)
	end do
	
	write(*,*)' Position of atoms: '
	do i=1,2*nlayers
		write(*,*) posA(i,:)
	end do
	do i=1,natoms
		write(*,*) pos(i,:)
	end do

	endif ! 1==0
!	 write(*,*) 'Pos:'
!	do i=1,natoms
!	 write(*,'(3f10.5)') pos(i,:)
!	end do
	
	return
	end 	subroutine getstructure
!----------------------------------------------------------------------











	subroutine getstructure2()
	implicit none
	integer :: i
	double precision :: v(3,2), u(3)


	!nlayers = 2;
	noctl = 2;
	noct = noctl*nlayers;
	natoms = noct*4 !(4+1); ! testing: add Sr to normal atoms to mk Strucx with the smae routine
	a = a0; !7.0d0;

	!nsptm =1;
	nspecies = nsptm;
	
	norbtm = 5; norbtms = norbtm*nspin;
	norbo = 3; norbos = norbo*nspin;

	ntot = noct*norbtms + noct*3*norbos;
	ntottm = noct*norbtms;

	write(*,*)'ntot, ntottm, qtot,a = ', ntot, ntottm, qtot, a
	!i = noct*1 + noct*6*norbos
	!write(*,*)'TM d^1: nelec = ', i
	!write(*,*)'if only one spin, filled bands ~ ',0.5*i



	
	allocate(oct(nlayers,noctl))

	! set the basic cubic structure
	! sqrt(2) x sqrt(2) x nlayers cell
	! square cell rotated by 45 wr.r.t x,y


	! after tilt/rotation of octahedra, unit cell/lattice vectors  rescale
	! so avec and oct(:,:)%rb will be calculated AFTER the tilt/rotation, disabling these below:
	if(1==0) then
	! lattice vectors
	!
	avec(:,1) = (/1.0d0, -1.0d0, 0.0d0/)*a
	avec(:,2) = (/1.0d0, +1.0d0, 0.0d0/)*a
	avec(:,3) = (/0.0d0,  0.0d0, 1.0d0/)*a*nlayers
	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)

	! atomic positions in cartesian
	! set pos of B
	do il=1,nlayers
	 oct(il,1)%rb = (/0.5d0,-0.5d0,1.d0*(il-1)/)*a;
	 oct(il,2)%rb = (/0.5d0,+0.5d0,1.d0*(il-1)/)*a;
	enddo

	endif


	! oxygens
	do il=1,nlayers
	 do io=1, noctl
	  	oct(il,io)%xo(:,1) = (/0.5d0,0.0d0,0.0d0/)*a
	  	oct(il,io)%xo(:,2) = (/0.0d0,0.5d0,0.0d0/)*a
	  	oct(il,io)%xo(:,3) = (/0.0d0,0.0d0,0.5d0/)*a
	  oct(il,io)%lo = 0.5d0 * a;

		! fake atoms to get the structure with rotations/tilts
	  	oct(il,io)%xo(:,4) = -oct(il,io)%xo(:,1)
	  	oct(il,io)%xo(:,5) = -oct(il,io)%xo(:,2)
	  	oct(il,io)%xo(:,6) = -oct(il,io)%xo(:,3)

	 end do
	 end do


	! rotate and set abs pos of oxygen
	call rotoctall(theta,phii,a1,a2,a3)

	! set pos and posA for getnns()
	allocate(posA(nlayers*2,3))
	allocate(pos(nlayers*(8),3)) !*(8+2)

	! Sr atoms fractional positions in layer1
	v(:,1) = (/0.5d0,0.0d0, 0.25d0/);
	v(:,2) = (/0.0d0,0.5d0, 0.25d0/);

	!write(*,*)'main: testing... oct(il,io)%rb:'
	
	do il=1,nlayers
	 do io=1, noctl
	 	i = (il-1)*8 + (io-1)*4;
	  pos(i+1,:) = oct(il,io)%rb

		!write(*,'(3f10.5)') oct(il,io)%rb

	  !pos(i+2,:) = oct(il,io)%rb + oct(il,io)%xo(:,1)
	  !pos(i+3,:) = oct(il,io)%rb + oct(il,io)%xo(:,2)
	  !pos(i+4,:) = oct(il,io)%rb + oct(il,io)%xo(:,3)
		!posA((il-1)*2 + io,:) = oct(il,io)%rb + (/a1,a2,a3/)*0.5d0

	  pos(i+2,:) = oct(il,io)%ro(:,1)
	  pos(i+3,:) = oct(il,io)%ro(:,2)
	  pos(i+4,:) = oct(il,io)%ro(:,3)

		! frac coor of Sr
		u = v(:,io) + (/0.0d0,0.0d0,(il-1)*0.5d0/);
		! cartesian coor of St atoms:
		call r3mv(avec, u , posA((il-1)*2 + io,:))

		!pos(noct*4+(il-1)*2 + io,:) = posA((il-1)*2 + io,:) ! extra normal atoms

	 end do
	end do

	!do i=1,natoms
	!	write(*,*) 'pos =', pos(i,:)
	!end do

	return
	end 	subroutine getstructure2
!----------------------------------------------------------------------








	!..............................................................

	subroutine rotoctall(th,phi,a1,a2,a3)
	use esvar, only: Rlmax
	use rotylm, only: getRlmax
	implicit none
	double precision, intent(in) ::	th, phi
	double precision, intent(out) :: a1,a2,a3
	double precision :: v(3), tt
	double precision, dimension(3) :: s,s1,s2,s5,r2,r4,r5,t,r3,t6
	integer :: i,il,io
	double precision, dimension(4) :: ths,phis
	
	if(mod(nlayers,2) /= 0) then
		write(*,*) "Error: even number of layers req for tilting!"
	endif

	! Carter/Kee/Zeb PRB 2012; (theta,phi) signs:
	! il=1: Blue = ++ , Red = --
	! il=2: Yellow= +-, Green =-+

	do il=1,nlayers,2
	 call rotoct(il  ,1, th, phi) ! Blue 
	 call rotoct(il  ,2,-th,-phi) ! Red
	 call rotoct(il+1,1,-th, phi) ! Yellow
	 call rotoct(il+1,2, th,-phi) ! green
	end do

	! get rotation matrices for Ylm for all four octahedral rotations
	! remember that if we have multiple orthorhombic cells along z; 
	! we cannot have independent tilt+rotations; 
	! so only octa rotations of first orthorhombic cell is required for Ylm.
	ths =  (/+1,-1,-1,+1/) * th;
	phis = (/+1,-1,+1,-1/) * phi;
	do i=1,4
	 call getRlmax(ths(i), phis(i), Rlmax(:,:,i))
	end do

	! unit cell rescales with the tilt/rotation:
	
	! determine lattice vectors
	! first determine position of oct2 with respect to oct 1:
	! notation: oct1: r's; oct2, s's
	! Oxygen atoms label: 1,2,3 along x,y,z; 4,5,6 along -x,-y,-z.
	il = 1;
	r2 = oct(il,1)%xor(:,2)
	r4 = oct(il,1)%xor(:,4)
	r5 = oct(il,1)%xor(:,5)

	s1 = oct(il,2)%xor(:,1)
	s2 = oct(il,2)%xor(:,2)
	s5 = oct(il,2)%xor(:,5)
	! pic of notebook calc saved in src dir for reference:	
	s = r4 - s1;
	avec(:,1) = r2 - (s5 + s); 
	avec(:,2) 	= s2 + s - r5;
		
	! for avec3
	r3 = oct(1,1)%xor(:,3)
	t6 = oct(2,1)%xor(:,6)
	t = r3 - t6;
	tt = t(3);
	avec(:,3) = (/0.0d0,  0.0d0, tt*nlayers/)

	!write(*,*)'main: TESTING: a1 <->  a2 avec:'
	!v = avec(:,1);
	!avec(:,1) = avec(:,2);
	!avec(:,2) 	= v;

!	write(*,*) avec(:,1)	
!	write(*,*) avec(:,2)	
!	write(*,*) avec(:,3)	


 	

	
	! calc ainv using updated avec:
	call r3minv(avec,ainv)



	

	if(1==0) then
	write(*,*) 'B-O distances: o1,o2,o3'
	il=1;io=1;
	v = oct(il,io)%xor(:,1);
	write(*,*) v(1)**2 + v(2)**2 + v(3)**2
	v = oct(il,io)%xor(:,2);
	write(*,*) v(1)**2 + v(2)**2 + v(3)**2
	v = oct(il,io)%xor(:,3);
	write(*,*) v(1)**2 + v(2)**2 + v(3)**2
	endif
	

	! atomic positions of TM atoms at the centre of octahedra
	! assuming a1=a2: if true, then its simple, 
	! otherwise rotation matrix has to be invoked
	do il=1,nlayers
	 oct(il,1)%rb = (/0.d0, 0d0 ,(il-1)*tt /);
	 oct(il,2)%rb = oct(il,1)%rb + s;
	enddo


	! using rescaled TM positions
	do il=1,nlayers
	 do io=1,2
	 
	 	! set abs value of oxygen position after tilt/rotation.
	  do i = 1,3
	   oct(il,io)%ro(:,i) =	oct(il,io)%rb(:) + oct(il,io)%xor(:,i);
	  end do

	 end do
	end do ! il
	
	do il=1,nlayers
	 do io=1,2

	  ! ro cartesian to fractional
	  do i=1,3
	   call r3mv(ainv,oct(il,io)%ro(:,i),v)
	   oct(il,io)%rof(:,i) = v
	  end do
	  
	  ! central B atom
	  call r3mv(ainv,oct(il,io)%rb,v)
	  oct(il,io)%rbf = v

	 end do
	end do ! il


	if(1==0) then
	il=1;io=1;
	write(*,*) 'cartesian:'
	write(*,*) oct(1,1)%ro(:,1)
	write(*,*) oct(1,1)%ro(:,2)
	write(*,*) oct(1,1)%ro(:,3)
	write(*,*) 'fractional:'
	write(*,*) oct(1,1)%rof(:,1)
	write(*,*) oct(1,1)%rof(:,2)
	write(*,*) oct(1,1)%rof(:,3)
	endif

	
	return
	end 	subroutine rotoctall






	subroutine getneighbours()
	implicit none
	
	call settmnn1()
	call setoxnn1()
	if(tmnn2) call settmnn2()
	if(oxnn2) call setoxnn2()

	return
	end 	subroutine getneighbours
!----------------------------------------------------------------------

	subroutine getmaps()
	implicit none
	
!	do il=1,nlayers
!	 do io=1,noctl ! noctl = 2 always
!	  	tm(il,io)%is = layersp(il)
!	 end do
!	end do
	! atom to orbitals map
	call mapatom2orbs()
	! atom to species (TM) map
	call mapatom2species()

	if(1==0) then
	write(*,'(a)')'TM:-------------------'
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always
	  ia = tm(il,io)%ia
	  write(*,'(a,i5,a,2i5)') 'ia=',ia, ' range: ',atom2orb(:,ia)
	 end do
	end do

	write(*,'(a)')'O:-------------------'
	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always
	  do i=1,3
	   ia = ox(il,io,i)%ia
	   write(*,'(a,i5,a,2i5)') 'ia=',ia, ' range: ',atom2orb(:,ia)
	  end do
	 end do
	end do

	endif

	return
	end 	subroutine getmaps
!----------------------------------------------------------------------
	! calculated on mathematica, first calc in complex Ylm, then basis changed to real Ylm, wiki's real Ylm.
	! gives eienvalues of H = L.S + Voct that are independent of orientation of the octahedra
	! verifying that L.S operator's matrix representation indeed corresponds to a scalar operator of 3D position space, i.e., the dot product L.S;

! siesta's Hsoc or its version with sign changed does not correspond to a 3D scalar operator... probably some mistakes in their matrix.

	subroutine mkhsocl2()
	use modmain, only: hsoc
	implicit none
	double complex, parameter :: i = (0.0d0, 1.0d0);
	double complex, parameter :: o = (1.0d0, 0.0d0);
	double complex, parameter :: t = (2.0d0, 0.0d0);
	double complex, parameter :: z = (0.0d0, 0.0d0);
	double complex, parameter :: d = zsqrt((3.d0,0.d0));

	double complex, dimension(10,10) :: hz, U
	
	! m= -2,-1,1,2,0 basis order
	hsoc = 0.0d0;
	hsoc(1,:) =(/z, z, z, -t*I, z, z, o, I, z, z/)
	hsoc(2,:) =(/z, z, -I, z, z, -o, z, z, I,   I*d/)
	hsoc(3,:) =(/z, I, z, z, z, -I, z, z, -o, d/)
	hsoc(4,:) =(/t*I, z, z, z, z, z, -I, o, z, z/)
	hsoc(5,:) =(/z, z, z, z, z, z, -I*d, -d, z,   z/)
	hsoc(6,:) =(/z, -o, I, z, z, z, z, z, t*I, z/)
	hsoc(7,:) =(/o, z, z, I, I*d, z, z, I, z, z/)
	hsoc(8,:) =(/-I, z, z, o, -d, z, -I, z, z, z/)
	hsoc(9,:) =(/z, -I, -o, z, z, -t*I, z, z, z, z/)
	hsoc(10,:) =(/z, -I*d, d, z, z, z, z, z, z, z/)


! m=-2:2 basis order
!	hsoc = 0.0d0;
!	hsoc(1,:) = (/z,z,z,z,-t*i,z,o,z,i,z/)
!	hsoc(2,:) = (/z,z,z,-i,z,-o,z,i*d,z,i/)
!	hsoc(3,:) = (/z,z,z,z,z,z,-i*d,z,-d,z/)
!	hsoc(4,:) = (/z,i,z,z,z,-i,z,d,z,-o/)
!	hsoc(5,:) = (/t*i,z,z,z,z,z,-i,z,o,z/)
!	hsoc(6,:) = (/z,-o,z,i,z,z,z,z,z,t*i/)
!	hsoc(7,:) = (/o,z,i*d,z,i,z,z,z,i,z/)
!	hsoc(8,:) = (/z,-i*d,z,d,z,z,z,z,z,z/)
!	hsoc(9,:) = (/-i,z,-d,z,o,z,-i,z,z,z/)
!	hsoc(10,:)=(/z,-i,z,-o,z,-t*i,z,z,z,z/)   


	! write Hsoc in complex spherical harmonics basis
	if(1==0) then
	 U = 0.0d0;
	 U(1:5,1:5) = Ur;
	 U(6:10,6:10) = Ur;
	 hz = matmul(U,hsoc)
	 hz = matmul(hz, conjg(transpose(U)))
	 write(*,'(10f8.4)') 	dble(hz)
	 write(*,*)'-----------------'
	 write(*,'(10f8.4)') 	dimag(hz)
	endif





	hsoc = 0.2d0 * hsoc ! to make the soc gap 1 electron volt at lambda=1


	
	return
	end 	subroutine mkhsocl2

!----------------------------------------------------------------------

!====================================================================
	subroutine writeB1ham4edrixs()
	implicit none
	integer ::i
	character*50 :: fname ='beast-edrixs-input-ham.in'
	double precision, dimension(5,5) :: h, R, RT
	double complex, dimension(5,5,3) :: hz
	
	R = Rlmax(5:9,5:9,1) ! Rlmax is rotation matrix for Ylm of atom 1
	RT = transpose(R);
	do i=1,3
	 ! transform to rotated frame of the octagedron 1 
	  h = B1hams(:,:,i)
	  ! R^T.h.R
	  h = matmul(RT,h)
	  h = matmul(h,R) ! 
	 ! convert dm to complex spherical harmonics
	 ! mind our order of real harmoncis: m2i list; i2m list?	 
	  call rtozflm(dcmplx(h),hz(:,:,i))
	end do
	
	open(10,file=trim(fname),action='write',position='append')
	do i=1,3
	 write(10,'(10f18.12)')  hz(:,:,i) ! write 5x5 matrices as an array of 25 elements
	end do
	close(10)
	
	return
	end 	subroutine writeB1ham4edrixs
!----------------------------------------------------------------





	end 	program perovskite
