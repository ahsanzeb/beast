

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
		 !ja = tm(il,io)%nn1(k)%ia
		 lmn = getlmn(tm(il,io)%nn1(k)%dr)
	   do i = 1,norbtm
		  do j=1,norbo
	     tm(il,io)%nn1(k)%h(i,j) = slatkosdp(lmn,spd,ppd,i,j)
		  end do !j
		 end do ! i
	  end do ! k
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
		 end do ! i
	  end do ! k
	 end do ! io
	end do ! il
	!...................................................................



	!...................................................................
	! O-TM  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  !ia = ox(il,io,ii)%ia;
	  do k = 1,2 ! 1st nns TM
	   allocate(ox(il,io,ii)%nn1(k)%h(norbo,norbtm))
		 ja = ox(il,io,ii)%nn1(k)%ia
		 lmn = getlmn(ox(il,io,ii)%nn1(k)%dr)
	   spd = ; ppd= ;
	   do i = 1,norbo
		  do j=1,norbtm
	     ox(il,io,ii)%nn1(k)%h(i,j) = slatkospd(lmn,spd,ppd,i,j)
		  end do !j
		 end do ! i
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	!...................................................................
	!...................................................................
	! O-O  (2nd neighbours)
	!...................................................................
	spp= ; ppp = ;
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  ia = ox(il,io,ii)%ia;
	  do k = 1,8 ! 2nd nns O
	   allocate(ox(il,io,ii)%nn2(k)%h(norbo,norbo))
		 !ja = ox(il,io,ii)%nn2(k)%ia
		 lmn = getlmn(ox(il,io,ii)%nn2(k)%dr)
	   do i = 1,norbo
		  do j=1,norbo
	     ox(il,io,ii)%nn2(k)%h(i,j) = slatkospp(lmn,spp,ppp,i,j)
		  end do !j
		 end do ! i
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	!...................................................................







	end 	module hamiltonian
