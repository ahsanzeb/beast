

	module slaterkoster
	use modmain
	implicit none

!------------------------------------------------------------------------
	subroutine realHij()
	implicit none
	
	integer :: il,io,is, js, ii
	double precision, dimension(3,norbtm) :: h	

	!...................................................................
	! TM-O  (1st neighbours)
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	  is = tm(il,io)%is; ! species index of TM atom
	  do k = 1,6 ! 1st nns O
	   allocate(tm(il,io)%nn1(k)%h(norbtm,norbo))
	   ! slatkospd computed (p,d): to get (d,p), r ===> -r and aux h as out
	   call slatkospd(-1.d0*tm(il,io)%nn1(k)%dr, skbo(is,:), h)
     tm(il,io)%nn1(k)%h = transpose(h);
	  end do ! k
	 end do ! io
	end do ! il
	!...................................................................
	! 	TM-TM (2nd neighbours) 
	!...................................................................
	do il=1,nlayers
	 do io=1, noctl
	  	is = tm(il,io)%is;
	  do k = 1,6 ! 2nd nns TM
	   allocate(tm(il,io)%nn2(k)%h(norbtm,norbtm))
		 js = tm(il,io)%nn2(k)%is;
		 lmn = getlmn(tm(il,io)%nn2(k)%dr)
	   call slatkosdd(tm(il,io)%nn2(k)%dr,
     .                   skbb(is,js,:),tm(il,io)%nn2(k)%h)
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
		 js = ox(il,io,ii)%nn1(k)%is
	   call slatkospd(ox(il,io,ii)%nn1(k)%dr,
     .                     skbo(js,:),ox(il,io,ii)%nn1(k)%h)
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	!...................................................................
	!...................................................................
	! O-O  (2nd neighbours)
	!...................................................................
	!skoo == (spp,ppp)
	do il=1,nlayers
	 do io=1, noctl
	 do ii=1,3
	  !ia = ox(il,io,ii)%ia;
	  do k = 1,8 ! 2nd nns O
	   allocate(ox(il,io,ii)%nn2(k)%h(norbo,norbo))
	   call slatkospp(ox(il,io,ii)%nn2(k)%dr,
     .                            skoo,ox(il,io,ii)%nn2(k)%h)
	  end do ! k
	 end do ! ii
	 end do ! io
	end do ! il
	!...................................................................



	return
	end 	subroutine realHij
!------------------------------------------------------------------------


	end 	module slaterkoster
