

	module hamiltonian
	use modmain

	implicit none


	! O-TM  (1st neighbours) 
	! O-O   (2nd neighbours) 

	!...................................................................
	! TM-O  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	  ia = tm(il,io)%ia;
	  spd = ; ppd= ;
	  do k = 1,6 ! 1st nns O
	   allocate(tm(il,io)%nn1(k)%h(norbtm,norbo))
		 ja = tm(il,io)%nn1(k)%ia
		 lmn = getlmn(tm(il,io)%nn1(k)%dr)
	   do i = 1,norbtm
		  do j=1,norbo
	     tm(il,io)%nn1(k)%h(i,j) = slatkosdp(lmn,spd,ppd,i,j)
		  end do !j
		 end do ! k
	  end do ! i
	 end do ! io
	end do ! il
	!...................................................................
	! 	TM-TM (2nd neighbours) 
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	  ia = tm(il,io)%ia;
	  do k = 1,6 ! 2nd nns TM
	   allocate(tm(il,io)%nn2(k)%h(norbtm,norbtm))
		 ja = tm(il,io)%nn2(k)%ia
		 lmn = getlmn(tm(il,io)%nn2(k)%dr)
	   ddd = ; pdd = ; sdd = ; ! multiple TM species possible, in k loop
	   do i = 1,norbtm
		  do j=1,norbtm
	     tm(il,io)%nn2(k)%h(i,j) = slatkosdd(lmn,ddd,pdd,sdd,i,j)
		  end do !j
		 end do ! k
	  end do ! i
	 end do ! io
	end do ! il
	!...................................................................















	itm1 = 0;
	do il=1,nlayers
	 do io=1, noctl
	  ind1 = (il-1)*noctl+ io
	  itm1 = (ind1-1)*norbtm + 1;
	  	itm2 = itm1 + norbtm;
	  do jl = min(1,il-1), max(nlayers,il+1)
	   do jo = 1, noctl
	    ind2 = (jl-1)*noctl+ jo
	    itm1 = (ind1-1)*norbtm + 1;
	    	itm2 = itm1 + norbtm;

	    oct(il,io)%xo(1,:) = 

	    
	  	 end do
	  	end do
	 end do
	end do



























	end 	module hamiltonian
