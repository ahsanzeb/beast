
	program perovskite
	use modmain
	use skoster
	use hamiltonian
	use readinput
	use fermi
	use modgaunt, only: mkdgaunt
	use Hubbard, only: mkvee, mkvmat
	use scf, only: groundstate
	use pdos

	use estatic, only: setmadvar, initmadelung
	
	implicit none
	integer :: il, io, i, ia
	integer :: ik,ib
	double precision, dimension(3) :: kvec


	! read input file 'input.in'
	call input()

	! set nspin
	if(lspin) then
	 nspin=2
	elseif(lsoc .or. lhu) then
	 nspin = 2;
	 write(*,*)"main: setting nspin=2 as lsoc/lhu = T"
	else
	 nspin=1
	endif



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

	if(any(phi(:) > 45.0)) then
	 write(*,'("O atoms have to cross to make this big rotation!")')
	 stop 
	else
	 write(*,'(a)')
     .  'Octraherda size will be rescaled to keep "a" fixed.'
	endif

	!tmnn2 = .true.;
	!oxnn2 = .true.;


	!nlayers = 2;
	noctl = 2;
	noct = noctl*nlayers;
	natoms = noct*4;
	a = 7.0d0;


	!nsptm =1;
	nspecies = nsptm;
	
	norbtm = 5; norbtms = norbtm*nspin;
	norbo = 3; norbos = norbo*nspin;

	ntot = noct*norbtms + noct*3*norbos;
	ntottm = noct*norbtms;

	write(*,*)'ntot, ntottm, qtot = ', ntot, ntottm, qtot
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
	avec(:,1) = (/1.0d0, -1.0d0, 0.0d0/)*a
	avec(:,2) = (/1.0d0, +1.0d0, 0.0d0/)*a
	avec(:,3) = (/0.0d0,  0.0d0, 1.0d0/)*a*nlayers

	! calc ainv for coordinate transformations
	ainv = 0.0d0;
	call r3minv(avec,ainv)

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
	  	oct(il,io)%xo(1,:) = (/0.5d0,0.0d0,0.0d0/)*a
 !    .   +(/rand(0),rand(0),rand(0)/)*a
	  	oct(il,io)%xo(2,:) = (/0.0d0,0.5d0,0.0d0/)*a
!     .   +(/rand(0),rand(0),rand(0)/)*a
	  	oct(il,io)%xo(3,:) = (/0.0d0,0.0d0,0.5d0/)*a
!     .   +(/rand(0),rand(0),rand(0)/)*a
	  oct(il,io)%lo = 0.5d0 * a;
	 end do
	end do
	
	! rotate and set abs pos of oxygen
	do il=1,nlayers
	 do io=1, noctl
	  call rotate(il,io)
	 end do
	end do

	! write GEOMETRY.OUT for visualisation with VESTA
	call writegeom()

	!write(*,*)'-------------------- 1'

	! set nearest neighbours:
	call settmnn1()
	call setoxnn1()
	if(tmnn2) call settmnn2()
	if(oxnn2) call setoxnn2()

!	do il=1,nlayers
!	 do io=1,noctl ! noctl = 2 always
!	  	tm(il,io)%is = layersp(il)
!	 end do
!	end do
	! atom to orbitals map
	call mapatom2orbs()
	! atom to species (TM) map
	call mapatom2species()


	call reciplat(avec,bvec,omega,omegabz)

 !--------------------------------------------------------
 ! once before SCF cycle starts:
	call setmadvar()
	call initmadelung()
 !--------------------------------------------------------


	!write(*,*)'ia, i1, i2 '
	!do i=1,natoms
	!write(*,'(20i5)')i, atom2orb(:,i)
	!end do

	!write(*,*)'-------------------- 2'


	
	!write(*,*)'-------------------- 3'

	
	! dummy data set:
	! nsptm = number of species of TM atoms
	!allocate(skbo(nsptm,2)) 	! 2: sigma_pd, pi_pd
	!allocate(skbb(nsptm,nsptm,3))	! 3: sigma_dd, pi_dd, delta_dd
	!skbo = 0.0d0; skbb=0.0d0;
	!skbo(1,:) = (/1.0d0,0.5d0/);
		
	!skbb(1,1,:)	 = (/0.0d0,0.0d0,0.00d0/);
	!skoo = 0.0d0 ! some prob with o-o, calc... if finite skoo.
	
	!write(*,*)'-------------------- 1'

	! sets some basis transformation matrices:
	! some condition? inside an if block?
	call setUrUz() ! sets Uz and Ur matrices
	call setUlm2j() ! sets Ulm2j, uses Ur that is set in setUrUz() call above.


	! real hamiltonian matrix elements using SK method
	call realHij()

	! spin-orbit coupling universal hamiltonian in full d-orbital space of a TM atom
	call mkHsoc()

	!call realHii() ! sets onsite hamiltonian matrix elements

	!write(*,*)'------------real Hij done -------- '

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

	stop "main: stopping...."
	endif


	!call reciplat(avec,bvec,omega,omegabz)

	if(1==0) then
	write(*,'("Reciprocal lattive vectors:")')
	write(*,'(3f10.3)') bvec(:,1)
	write(*,'(3f10.3)') bvec(:,2)
	write(*,'(3f10.3)') bvec(:,3)
	endif
	


	if(lgs .or. 1==1) then
	 call groundstate()
	 call getpdos(nwplot)
	endif

	
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

	do ik= 1,np
	 !call getHk((/0.5d0,0.5d0,0.5d0/),hk)
	 !write(*,'(a,i10,3f10.4)')' ik, k = ',ik,vpl(:,ik)
	 kvec = bvec(:,1)*vpl(1,ik) + 
     .    bvec(:,2)*vpl(2,ik) +
     .    bvec(:,3)*vpl(3,ik)
	 call getHk(ik,kvec, hk, 0)
	 !write(*,*)' ham done... '


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
	!-------------------------------------------
	open(10,file='BAND.OUT',form='FORMATTED',action='write')
	do ib=1,ntot
	 do ik=1,np
	  write(10,'(2G20.8)') dp(ik), eval(ik,ib) -efermi
	 end do
	 write(10,'(2f15.8)')
	 write(10,'(2f15.8)') 
	end do
	close(10)
	!-------------------------------------------
	! band character:
	if(lbc) then
	 call getbc()
	endif
	!-------------------------------------------


	! output the vertex location lines
	open(50,file='BANDLINES.OUT',form='FORMATTED',action='write')
	do i=1,nv
	write(50,'(2G18.10)') dv(i),0.0 !emin
	write(50,'(2G18.10)') dv(i),0.0 !emax
	write(50,'("     ")')
	end do
	close(50)

	
	write(*,'(a)')'--------- complete -----------'

	end 	program perovskite
