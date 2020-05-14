	module mtbshfl
	implicit none
	contains
	

      subroutine tbshfl(job,nm,n,m,dm)
C- Shuffle n x m matrix from MSM lm indices to TBE lm indices
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job: job = 0 MSM -> TBE
Ci        job = 1 TBE -> MSM
Ci   nm : leading dimension of dm
Ci   n,m: actual dimensions of dm
Ci   dm : a matrix to reshuffle
Cio Input/Outputs:
Cio   dm reordered
Cr Remarks:
Cr   Michael's structure constants are ordered according to the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     y     z     x    xy    yz 3z^2-r^2 zx   x^2-y^2
Cr   while the TBE programs use the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     x     y     z    xy    yz    zx  x^2-y^2  3z^2-r^2
Cr   This widget rearranges the matrix made by Michael's HBLSTR into
Cr   the TBE order.
Cr   The l>2 ordering is unchanged
Cr
Cr   the program is analogous to xxxind (see tb/tbesel.f) except that
Cr   it works with a rectangular matrix of arbitrary dimensions dm(n x m)
Cr   and in both directions
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer, intent(in)    :: job,nm,n,m
      real(8), intent(inout) :: dm(nm,m)
C Local Variables
      real(8) :: wk(nm,m)
!       real(8), allocatable :: wk(:,:)

      integer i, j, nn, mm
      integer, save, target :: ish(9) = (/1,4,2,3,5,6,8,9,7/)
      integer, save, target :: ishm(9) = (/1,3,4,2,5,6,9,7,8/)

      integer, pointer :: ind(:)

c... Checks
      if (n > nm) write(*,*)' tbshfl: incompatible dimensions'

!       allocate(wk(nm,m))

      if (job == 0) then
         ind => ish
      else if (job == 1) then
         ind => ishm
      endif

      nn = min(9,n)
      mm = min(9,m)
c... upper left corner
      forall ( i=1:nn, j=1:mm ) wk(i,j) = dm(ind(i),ind(j))
      dm(1:nn,1:mm) = wk(1:nn,1:mm)
c... remaining columns
      if (n > 9) then
        forall ( i=10:n, j=1:mm ) wk(i,j) = dm(i,ind(j))
        dm(10:n,1:mm) = wk(10:n,1:mm)
      endif
c... remaining rows
      if (m > 9) then
        forall ( i=1:nn, j=10:m) wk(i,j) = dm(ind(i),j)
        dm(1:nn,10:m) = wk(1:nn,10:m)
      endif

      end


      subroutine strfac(job,nm,n,m,dm)
      use esvar, only: ll
C- Include factor in n x m structure constant block
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job: job = 0 MSM -> TBE
Ci        job = 1 TBE -> MSM
Ci   nm : leading dimension of dm
Ci   n,m: actual dimensions of dm
Ci   dm : n x m matrix to rescale
Cio Input/Outputs:
Cio   dm * fac
Cr Remarks
Cr   job = 0:
Cr   To convert to Stone's convention the l_i,l_j element is multiplied
Cr   by fac = sqrt((2l_j+1)*(2l_i+1))/((2l_j+1)!!(2l_i+1)!!)
Cr   job = 1:
Cr   dm is multiplied by 1/fac
Cr
Cr   the program is analogous to xxxfac (see tb/tbesel.f) except that
Cr   it works with a rectangular matrix of arbitrary dimensions dm(n x m)
Cr   and in both directions
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer, intent(in)             :: job,nm,n,m
      double precision, intent(inout) :: dm(nm,m)
C Local Variables
      integer, parameter :: lmax = 5
      integer i,j,li,lj
      double precision arg,fac
      !external ll
C ... df(l) = (2l+1)!!
      integer, parameter :: df(0:lmax) = (/1,3,15,105,945,10395/)

c... Checks
      if (max(n,m) > (lmax+1)**2)
     .  write(*,*)' xsqfac: increase parameter lmax'
      if (n > nm) write(*,*)' strfac: incompatible dimensions'
      !call sanrg(.true.,job,0,1,' strfac:','job')

      do  i = 1, n
        li = ll(i)
        do  j = 1, m
          lj = ll(j)
          arg = dsqrt(dble((2*lj+1)*(2*li+1)))
          if (job == 0) then
            fac = arg/dble(df(lj)*df(li))
          else
            fac = dble(df(lj)*df(li))/arg
          endif
          dm(i,j) = dm(i,j)*fac
        enddo
      enddo

      end







C ----------------------------------------------------------------------
! A-site atoms with fixed monopoles: shuffl TM/O multipoles... 
      subroutine tbshflA(m,dm)
      implicit none
      integer, intent(in)    :: m
      real(8), intent(inout) :: dm(m)

      real(8) :: wk(min(9,m))
      integer j, mm
      integer, parameter:: ind(9) = (/1,4,2,3,5,6,8,9,7/)

	    wk = 0.0d0
      mm = min(9,m)
      forall (j=1:mm ) wk(j) = dm(ind(j))
      dm(1:mm) = wk(1:mm)

      end
C ----------------------------------------------------------------------
      subroutine strfacA(m,dm)
      use esvar, only: ll
C- Include factor in 1 x m structure constant block
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer, intent(in)             :: m
      double precision, intent(inout) :: dm(m)
C Local Variables
      integer, parameter :: lmax = 5
      integer j,lj
      double precision arg,fac
C ... df(l) = (2l+1)!!
      integer, parameter :: df(0:lmax) = (/1,3,15,105,945,10395/)

c... Checks
      if (m > (lmax+1)**2)
     .  write(*,*)' xsqfac: increase parameter lmax'

        ! li = 0; df(li)=1
        do  j = 1, m
          lj = ll(j)
          arg = dsqrt(dble(2*lj+1))
          fac = arg/dble(df(lj))
          dm(j) = dm(j)*fac
        enddo

      end
C ----------------------------------------------------------------------




      
	end 	module mtbshfl
