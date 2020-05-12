	module msylmnc
	implicit none
	contains
	

      subroutine sylmnc(c,lmx)
C- Normalization constants for the spherical harmonics
C ----------------------------------------------------------------
Ci Inputs
Ci   lmx
Co Outputs
Co   C
Cr Remarks
Cr   use together with sylm (from ASW package)
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lmx
      double precision c(*)
C ... Local parameters
      integer i,i1,i2,l,lav,lp1,m,n1,n2,n3
      double precision fn2,fpi,tlp1,tpi,y0

      tpi = 8d0*datan(1d0)
      fpi = 2d0*tpi
      y0 = 1d0/dsqrt(fpi)
      c(1) = y0
      do  l = 1, lmx
        lp1 = l+1
        tlp1 = l+lp1
        lav = l*lp1 + 1
        c(lav) = dsqrt(tlp1/fpi)
        do  m = 1, l
          n2 = lp1-m
          n1 = n2+1
          n3 = l+m
          fn2 = n2
          do  i = n1, n3
            fn2 = fn2*i
          enddo
          i1 = lav+m
          i2 = lav-m
          c(i1) = dsqrt(tlp1/(fn2*tpi))
          c(i2) = c(i1)
        enddo
      enddo
      end
