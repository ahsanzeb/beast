
	program perovskite
	use modmain
	use skoster
	use hamiltonian
	use readinput
	
	implicit none
	integer :: il, io, i
	integer :: ik,ib
	double precision, dimension(3) :: kvec


	! read input file 'input.in'
	call input()
	
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
	
	norbtm = 5; norbtms = norbtm;
	norbo = 3; norbos = norbo;

	ntot = noct*norbtms + noct*3*norbos;
	write(*,*)'ntot = ', ntot
	i = noct*1 + noct*6*norbos
	write(*,*)'TM d^1: nelec = ', i
	write(*,*)'if only one spin, filled bands ~ ',0.5*i
	



	
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
	avec(1,:) = (/1.0d0, -1.0d0, 0.0d0/)*a
	avec(2,:) = (/1.0d0, +1.0d0, 0.0d0/)*a
	avec(3,:) = (/0.0d0,  0.0d0, 1.0d0/)*a*nlayers

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
	  	oct(il,io)%xo(1,:) = (/0.5d0,0.0d0,0.0d0/)*a;
	  	oct(il,io)%xo(2,:) = (/0.0d0,0.5d0,0.0d0/)*a;
	  	oct(il,io)%xo(3,:) = (/0.0d0,0.0d0,0.5d0/)*a;
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

	write(*,*)'-------------------- 1'

	! set nearest neighbours:
	call settmnn1()
	call setoxnn1()
	if(tmnn2) call settmnn2()
	if(oxnn2) call setoxnn2()

	write(*,*)'-------------------- 2'

	! atom to orbitals map
	call mapatom2orbs()
	! atom to species (TM) map
	call mapatom2species()

	
	write(*,*)'-------------------- 3'

	
	! dummy data set:
	! nsptm = number of species of TM atoms
	!allocate(skbo(nsptm,2)) 	! 2: sigma_pd, pi_pd
	!allocate(skbb(nsptm,nsptm,3))	! 3: sigma_dd, pi_dd, delta_dd
	!skbo = 0.0d0; skbb=0.0d0;
	!skbo(1,:) = (/1.0d0,0.5d0/);
		
	!skbb(1,1,:)	 = (/0.0d0,0.0d0,0.00d0/);
	!skoo = 0.0d0 ! some prob with o-o, calc... if finite skoo.
	
	write(*,*)'-------------------- 1'

	do il=1,nlayers
	 do io=1,noctl ! noctl = 2 always
	  	tm(il,io)%is = layersp(il)
	 end do
	end do




	! real hamiltonian matrix elements using SK method
	call realHij()


	write(*,*)'------------real Hij done -------- '



	call reciplat(transpose(avec),bvec,omega,omegabz)

	write(*,'("Reciprocal lattive vectors:")')
	write(*,'(3f10.3)') bvec(:,1)
	write(*,'(3f10.3)') bvec(:,2)
	write(*,'(3f10.3)') bvec(:,3)

	allocate(vpl(3,np))
	allocate(dv(nv))
	allocate(dp(np))

	call plotpt1d(bvec,nv,np,vvl,vpl,dv,dp)

	!write(*,*) vpl
	
	! vpl has kpoints along the bandlines....
	allocate(hk(ntot,ntot))
	allocate(eval(np,ntot))
	do ik=1,np
	 !call getHk((/0.5d0,0.5d0,0.5d0/),hk)
	 !write(*,'(a,i10,3f10.4)')' ik, k = ',ik,vpl(:,ik)
	 kvec = bvec(:,1)*vpl(1,ik) + 
     .    bvec(:,2)*vpl(2,ik) +
     .    bvec(:,3)*vpl(3,ik)
	 call getHk(kvec, hk)
	 !write(*,*)' ham done... '
	 call zdiag(ntot,hk,eval(ik,:),ik,ntot)
	end do
	!-------------------------------------------
	open(10,file='BAND.OUT',form='FORMATTED',action='write')
	do ib=1,ntot
	 do ik=1,np
	  write(10,'(2G20.8)') dp(ik), eval(ik,ib)
	 end do
	 write(10,'(2f15.8)')
	 write(10,'(2f15.8)') 
	end do
	close(10)
	!-------------------------------------------

	! output the vertex location lines
	open(50,file='BANDLINES.OUT',form='FORMATTED',action='write')
	do i=1,nv
	write(50,'(2G18.10)') dv(i),0.0 !emin
	write(50,'(2G18.10)') dv(i),0.0 !emax
	write(50,'("     ")')
	end do
	close(50)

	
	write(*,*)'-------------------- complete'

	end 	program perovskite
